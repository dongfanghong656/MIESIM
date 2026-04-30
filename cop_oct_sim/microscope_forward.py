from __future__ import annotations
import numpy as np
from scipy.ndimage import convolve
from .config_schema import SimulationConfig
from .grids import centered_axis_um, lateral_sample_pitch_um, make_pupil_grid
from .measurement import add_camera_noise, integrate_camera_pixels
from .metrics import center_line_fwhm
from .optical_paths import apply_path_dichroic_to_pupil
from .propagation import fraunhofer_psf_from_pupil, simple_defocus_stack
from .pupil import build_shared_pupil
from .scatterers import (
    convolve_scatterer_stack,
    effective_scatterer_diameter_um,
    finite_target_fwhm_ratio,
    scattering_measurement_stack,
)
from .spectrometer import first_column, read_numeric_table

def _camera_qe(config: SimulationConfig, wavelength_nm: float) -> float:
    if not config.microscope.camera_qe_file:
        return 1.0
    table = read_numeric_table(config.microscope.camera_qe_file)
    wavelength = first_column(table, ("wavelength_nm", "lambda_nm", "wavelength", "lambda"))
    qe = first_column(table, ("qe", "quantum_efficiency", "response", "camera_qe"), fallback=1 if "col1" in table else 0)
    order = np.argsort(wavelength)
    return float(
        np.interp(
            float(wavelength_nm),
            wavelength[order],
            qe[order],
            left=qe[order][0],
            right=qe[order][-1],
        )
    )

def _resize_2d(arr: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float64)
    if arr.ndim == 0:
        return np.full(shape, float(arr), dtype=np.float32)
    if arr.ndim == 1:
        arr = arr[None, :]
    if arr.shape == shape:
        return arr.astype(np.float32)
    src_x = np.linspace(0.0, 1.0, arr.shape[1])
    dst_x = np.linspace(0.0, 1.0, shape[1])
    x_resized = np.vstack([np.interp(dst_x, src_x, row) for row in arr])
    src_y = np.linspace(0.0, 1.0, x_resized.shape[0])
    dst_y = np.linspace(0.0, 1.0, shape[0])
    out = np.vstack([np.interp(dst_y, src_y, x_resized[:, j]) for j in range(shape[1])]).T
    return out.astype(np.float32)

def _camera_flat_field(config: SimulationConfig, shape: tuple[int, int]) -> np.ndarray:
    if not config.microscope.camera_flat_field_file:
        return np.ones(shape, dtype=np.float32)
    flat = np.genfromtxt(config.microscope.camera_flat_field_file, delimiter=",", dtype=float)
    flat = np.nan_to_num(flat, nan=1.0, posinf=1.0, neginf=1.0)
    return np.clip(_resize_2d(flat, shape), 0.0, None).astype(np.float32)

def _illumination_pupil(config: SimulationConfig, grid) -> np.ndarray:
    illum_na = config.illumination.NA
    na_ratio = 1.0 if illum_na is None else float(illum_na) / max(float(config.objective.NA_nominal), 1e-12)
    aperture_radius = min(max(float(config.illumination.aperture_stop_radius), 0.0), 1.0)
    support_radius = min(max(na_ratio, 0.0), aperture_radius, 1.0)
    r2 = grid.rho_x**2 + grid.rho_y**2
    support = grid.mask & (r2 <= support_radius**2)
    amp = np.exp(-r2 / max(float(config.illumination.pupil_fill_ratio), 1e-6) ** 2)
    return np.where(support, amp, 0.0).astype(np.complex128)

def _illumination_envelope(config: SimulationConfig, n: int, dx_um: float) -> np.ndarray:
    field_stop = config.illumination.field_stop_diameter_um
    if field_stop is None or field_stop <= 0:
        return np.ones((n, n), dtype=np.float32)
    axis = centered_axis_um(n, dx_um)
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    radius = 0.5 * float(field_stop)
    edge = max(dx_um, radius * 0.05, 1e-12)
    rr = np.sqrt(xx**2 + yy**2)
    envelope = 0.5 * (1.0 - np.tanh((rr - radius) / edge))
    return envelope.astype(np.float32)

def _illumination_power_fraction(illumination_pupil: np.ndarray, grid) -> float:
    denominator = max(float(np.sum(grid.mask)), 1e-12)
    return float(np.sum(np.abs(illumination_pupil) ** 2) / denominator)

def _illumination_shape_active(config: SimulationConfig) -> bool:
    return bool(
        config.illumination.NA is not None
        or abs(float(config.illumination.pupil_fill_ratio) - 1.0) > 1e-12
        or float(config.illumination.aperture_stop_radius) < 1.0 - 1e-12
    )

def _delta_psf(shape: tuple[int, int]) -> np.ndarray:
    out = np.zeros(shape, dtype=np.float32)
    out[shape[0] // 2, shape[1] // 2] = 1.0
    return out

def _illumination_psf(config: SimulationConfig, illumination_pupil: np.ndarray, pad_factor: int) -> np.ndarray:
    shape = (illumination_pupil.shape[0] * pad_factor, illumination_pupil.shape[1] * pad_factor)
    if not _illumination_shape_active(config):
        return _delta_psf(shape)
    field = fraunhofer_psf_from_pupil(illumination_pupil, normalize=False, pad_factor=pad_factor)
    psf = np.abs(field) ** 2
    total = float(np.sum(psf))
    if total <= 0:
        return _delta_psf(psf.shape)
    return (psf / total).astype(np.float32)

def _apply_illumination_shape(stack: np.ndarray, illumination_psf: np.ndarray) -> np.ndarray:
    if illumination_psf.size == 1 or np.count_nonzero(illumination_psf) == 1:
        return stack
    shaped = np.empty_like(stack, dtype=np.float32)
    for i, frame in enumerate(np.asarray(stack, dtype=np.float32)):
        shaped[i] = convolve(frame, illumination_psf, mode="nearest")
    return shaped

def simulate_microscope_zstack(config: SimulationConfig, N: int = 256,
                               bead_diameter_um: float | None = None,
                               mode: str | None = None,
                               pad_factor: int = 2):
    rng = np.random.default_rng(config.random_seed)
    mode = mode or config.microscope.mode
    bead_diameter_um = effective_scatterer_diameter_um(config, bead_diameter_um)
    wavelength_um = config.microscope.wavelength_nm * 1e-3
    camera_qe = _camera_qe(config, config.microscope.wavelength_nm)
    grid = make_pupil_grid(N, config.objective.NA_nominal)
    illumination_pupil = _illumination_pupil(config, grid)
    illumination_psf = _illumination_psf(config, illumination_pupil, pad_factor)
    illumination_power_fraction = _illumination_power_fraction(illumination_pupil, grid)
    Pc = build_shared_pupil(grid, wavelength_um,
                            fill_ratio=config.objective.pupil_fill_ratio,
                            defocus_um=config.objective.defocus_um)
    Pm, dichroic_coeffs = apply_path_dichroic_to_pupil(
        Pc, config, "mic_transmission", config.microscope.wavelength_nm, input_pol="mean"
    )
    nz = config.microscope.z_count
    z_um = (np.arange(nz) - nz // 2) * config.microscope.z_step_um
    field_stack = simple_defocus_stack(
        Pm, z_um, wavelength_um, config.objective.NA_nominal, pad_factor=pad_factor
    )
    if mode == "complex":
        ideal = field_stack
    elif mode == "widefield":
        ideal = np.abs(field_stack) ** 2
    elif mode in {"reflectance_confocal", "confocal"}:
        ideal = np.abs(field_stack * np.conj(field_stack)) ** 2
    else:
        raise ValueError(f"Unknown microscope mode: {mode}")

    if np.iscomplexobj(ideal):
        measured = ideal.astype(np.complex64)
        background = np.zeros_like(measured, dtype=np.complex64)
        raw = measured
        corrected = measured
        bead_stack = measured
        camera_integrated = measured
        camera_flat_field = np.ones(measured.shape[-2:], dtype=np.float32)
        sample_um = lateral_sample_pitch_um(wavelength_um, config.objective.NA_nominal, pad_factor)
        camera_pixel_sample_um = sample_um
        illumination_envelope = _illumination_envelope(config, measured.shape[-1], sample_um)
        scatter_raw_stack = measured
        scatter_background_stack = background
        scatter_corrected_stack = corrected
        sample_scattering_amplitude = 1.0 + 0.0j
        target_fwhm_ratio = 0.0
    else:
        ideal = ideal / max(float(np.max(ideal)), 1e-12)
        sample_um = lateral_sample_pitch_um(wavelength_um, config.objective.NA_nominal, pad_factor)
        illumination_envelope = _illumination_envelope(config, ideal.shape[-1], sample_um)
        bead_stack = convolve_scatterer_stack(
            ideal,
            config.sample,
            bead_diameter_um,
            lateral_pixel_um=sample_um,
            axial_pixel_um=config.microscope.z_step_um,
        )
        scatter = scattering_measurement_stack(bead_stack, config, config.microscope.wavelength_nm)
        scatter_raw_stack = scatter["raw"]
        scatter_background_stack = scatter["background"]
        scatter_corrected_stack = (
            _apply_illumination_shape(scatter["corrected"], illumination_psf)
            * illumination_power_fraction
            * illumination_envelope[None, :, :]
        )
        sample_scattering_amplitude = scatter["amplitude"]
        target_fwhm_ratio = finite_target_fwhm_ratio(
            bead_diameter_um,
            center_line_fwhm(ideal[nz // 2], sample_um),
        )
        camera_pixel_sample_um = config.microscope.camera_pixel_um / config.microscope.magnification_total
        integrated_frames = []
        measured_frames = []
        camera_flat_field = None
        for i in range(nz):
            frame = scatter_corrected_stack[i] * camera_qe
            frame = integrate_camera_pixels(frame, camera_pixel_sample_um, sample_um).astype(np.float32)
            if camera_flat_field is None:
                camera_flat_field = _camera_flat_field(config, frame.shape)
            frame = frame * camera_flat_field
            integrated_frames.append(frame)
            frame = frame + config.microscope.background_adu + config.microscope.camera_dark_adu
            noisy = (
                add_camera_noise(
                    frame,
                    config.errors.camera_read_noise_e,
                    config.errors.photon_gain_e_per_adu,
                    rng,
                ).astype(np.float32)
            )
            if config.microscope.camera_saturation_adu is not None:
                noisy = np.clip(noisy, 0.0, float(config.microscope.camera_saturation_adu)).astype(np.float32)
            measured_frames.append(noisy)
        if camera_flat_field is None:
            camera_flat_field = np.ones_like(integrated_frames[0], dtype=np.float32)
        camera_integrated = np.stack(integrated_frames).astype(np.float32)
        measured = np.stack(measured_frames).astype(np.float32)
        background = np.full_like(
            measured,
            config.microscope.background_adu + config.microscope.camera_dark_adu,
            dtype=np.float32,
        )
        raw = measured
        corrected = raw - background

    nx = measured.shape[-1]
    ny = measured.shape[-2]
    dx_um = sample_um
    camera_dx_um = max(camera_pixel_sample_um, sample_um)
    return {
        "z_um": z_um.astype(np.float64),
        "x_um": centered_axis_um(ideal.shape[-1], dx_um),
        "y_um": centered_axis_um(ideal.shape[-2], dx_um),
        "camera_x_um": centered_axis_um(nx, camera_dx_um),
        "camera_y_um": centered_axis_um(ny, camera_dx_um),
        "dx_um": np.array(dx_um, dtype=np.float64),
        "camera_pixel_sample_um": np.array(camera_dx_um, dtype=np.float64),
        "ideal": ideal.astype(np.complex64 if np.iscomplexobj(ideal) else np.float32),
        "measured": measured,
        "raw": raw,
        "background": background,
        "corrected": corrected,
        "bead_convolved": bead_stack.astype(np.complex64 if np.iscomplexobj(bead_stack) else np.float32),
        "scattering_raw": scatter_raw_stack.astype(np.complex64 if np.iscomplexobj(scatter_raw_stack) else np.float32),
        "scattering_background": scatter_background_stack.astype(
            np.complex64 if np.iscomplexobj(scatter_background_stack) else np.float32
        ),
        "scattering_corrected": scatter_corrected_stack.astype(
            np.complex64 if np.iscomplexobj(scatter_corrected_stack) else np.float32
        ),
        "camera_integrated": camera_integrated,
        "camera_qe": np.array(camera_qe, dtype=np.float64),
        "camera_flat_field": camera_flat_field.astype(np.float32),
        "camera_dark_adu": np.array(config.microscope.camera_dark_adu, dtype=np.float64),
        "camera_saturation_adu": np.array(
            np.nan if config.microscope.camera_saturation_adu is None else config.microscope.camera_saturation_adu,
            dtype=np.float64,
        ),
        "illumination_pupil": illumination_pupil.astype(np.complex64),
        "illumination_psf": illumination_psf.astype(np.float32),
        "illumination_envelope": illumination_envelope.astype(np.float32),
        "illumination_power_fraction": np.array(illumination_power_fraction, dtype=np.float64),
        "illumination_na": np.array(
            config.objective.NA_nominal if config.illumination.NA is None else config.illumination.NA,
            dtype=np.float64,
        ),
        "illumination_pupil_fill_ratio": np.array(config.illumination.pupil_fill_ratio, dtype=np.float64),
        "illumination_aperture_stop_radius": np.array(config.illumination.aperture_stop_radius, dtype=np.float64),
        "pupil": Pm.astype(np.complex64),
        "mode": mode,
        "bead_diameter_um": np.array(bead_diameter_um, dtype=np.float64),
        "scatterer_model": config.sample.scatterer_model,
        "measurement_path": config.sample.measurement_path,
        "sample_scattering_amplitude": np.array(sample_scattering_amplitude, dtype=np.complex128),
        "finite_target_fwhm_ratio": np.array(target_fwhm_ratio, dtype=np.float64),
        "dichroic_path": "mic_transmission",
        "dichroic_s_amplitude": np.array(dichroic_coeffs.s_amplitude, dtype=np.float64),
        "dichroic_p_amplitude": np.array(dichroic_coeffs.p_amplitude, dtype=np.float64),
        "dichroic_s_phase_rad": np.array(dichroic_coeffs.s_phase_rad, dtype=np.float64),
        "dichroic_p_phase_rad": np.array(dichroic_coeffs.p_phase_rad, dtype=np.float64),
    }
