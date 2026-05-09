from __future__ import annotations
import numpy as np
from .config_schema import SimulationConfig
from .grids import (
    centered_axis_um,
    fft_depth_axis_um,
    lateral_sample_pitch_um,
    make_pupil_grid,
)
from .optical_paths import apply_path_dichroic_to_pupil, dichroic_path_coefficients
from .pupil import build_shared_pupil
from .propagation import fraunhofer_psf_from_pupil, physical_defocus_phase
from .reconstruction import window_vector
from .scatterers import effective_scatterer_diameter_um, sample_scattering_amplitude, scatterer_volume_offsets
from .spectrometer import apply_k_linearization_error, differential_dispersion_phase, effective_sensitivity_rolloff_per_um, make_oct_k_grid

def _chromatic_physical_defocus_um(config: SimulationConfig, wavelength_nm: float) -> float:
    span = max(config.oct.bandwidth_nm, 1e-12)
    relative = (wavelength_nm - config.oct.center_wavelength_nm) / span
    return config.errors.objective_focus_shift_um + relative * config.objective.focus_shift_um_850_vs_visible

def _tilt_waves_from_error(config: SimulationConfig) -> float:
    return np.deg2rad(config.errors.dichroic_tilt_deg) * 2.0

def _rolloff_amplitude(config: SimulationConfig, z_um: float) -> float:
    return float(np.exp(-abs(z_um) * effective_sensitivity_rolloff_per_um(config)))

def _make_oct_pupil(
    config: SimulationConfig,
    grid,
    wavelength_nm: float,
    physical_defocus_um: float = 0.0,
    zernike_defocus_opd_um: float | None = None,
):
    wl_um = wavelength_nm * 1e-3
    tilt = _tilt_waves_from_error(config)
    opd_um = config.objective.defocus_um if zernike_defocus_opd_um is None else float(zernike_defocus_opd_um)
    pc = build_shared_pupil(
        grid,
        wl_um,
        fill_ratio=config.objective.pupil_fill_ratio,
        defocus_um=opd_um,
        tilt_x_waves=tilt,
    )
    physical_um = _chromatic_physical_defocus_um(config, wavelength_nm) + float(physical_defocus_um)
    if physical_um != 0.0:
        pc = pc * physical_defocus_phase(grid, physical_um, wl_um, config.objective.NA_nominal)
    return apply_path_dichroic_to_pupil(pc, config, "oct_reflection", wavelength_nm, input_pol="mean")

def _linearize_stack(k: np.ndarray, stack: np.ndarray):
    order = np.argsort(k)
    k_sorted = k[order]
    stack_sorted = stack[order]
    k_lin = np.linspace(k_sorted.min(), k_sorted.max(), len(k_sorted))
    flat = stack_sorted.reshape(len(k_sorted), -1)
    out = np.empty_like(flat, dtype=np.complex128)
    for j in range(flat.shape[1]):
        out[:, j] = (
            np.interp(k_lin, k_sorted, flat[:, j].real)
            + 1j * np.interp(k_lin, k_sorted, flat[:, j].imag)
        )
    return k_lin, out.reshape(stack_sorted.shape)

def _sample_complex_bilinear(field: np.ndarray, x_um: float, y_um: float, dx_um: float) -> complex:
    ny, nx = field.shape
    x = nx // 2 + float(x_um) / max(float(dx_um), 1e-12)
    y = ny // 2 + float(y_um) / max(float(dx_um), 1e-12)
    if x < 0 or y < 0 or x > nx - 1 or y > ny - 1:
        return 0.0 + 0.0j

    x0 = int(np.floor(x))
    y0 = int(np.floor(y))
    x1 = min(x0 + 1, nx - 1)
    y1 = min(y0 + 1, ny - 1)
    wx = x - x0
    wy = y - y0
    return complex(
        (1.0 - wx) * (1.0 - wy) * field[y0, x0]
        + wx * (1.0 - wy) * field[y0, x1]
        + (1.0 - wx) * wy * field[y1, x0]
        + wx * wy * field[y1, x1]
    )

def _normalize_abs(arr: np.ndarray) -> np.ndarray:
    mag = np.abs(arr).astype(np.float64)
    peak = float(np.max(mag)) if mag.size else 0.0
    return mag / peak if peak > 0 else mag

def _normalize_profile(profile: np.ndarray) -> np.ndarray:
    arr = np.asarray(profile, dtype=np.float64)
    peak = float(np.max(np.abs(arr))) if arr.size else 0.0
    return arr / peak if peak > 0 else arr

def _axial_profile_products(volume: np.ndarray) -> dict:
    mag = np.abs(volume).astype(np.float64)
    cy = mag.shape[1] // 2
    cx = mag.shape[2] // 2
    return {
        "axial_profile_max": _normalize_profile(np.max(mag, axis=(1, 2))),
        "axial_profile_on_axis": _normalize_profile(mag[:, cy, cx]),
        "axial_profile_integrated": _normalize_profile(np.sum(mag, axis=(1, 2))),
    }

def _sample_scatterer_rci(
    config: SimulationConfig,
    rci_field: np.ndarray,
    k: float,
    x_um: float,
    y_um: float,
    dx_um: float,
    scatterer_diameter_um: float,
) -> complex:
    step_um = min(max(dx_um, 1e-6), max(float(scatterer_diameter_um) / 4.0, 1e-3))
    dz, dy, dx, weights = scatterer_volume_offsets(
        config.sample,
        scatterer_diameter_um,
        lateral_pixel_um=step_um,
        axial_pixel_um=step_um,
    )
    value = 0.0 + 0.0j
    for dzi, dyi, dxi, wi in zip(dz, dy, dx, weights):
        lateral = _sample_complex_bilinear(rci_field, x_um + float(dxi), y_um + float(dyi), dx_um)
        value += float(wi) * lateral * np.exp(1j * 2.0 * float(k) * float(dzi))
    return complex(value)

def _separability_nrmse(volume: np.ndarray) -> float:
    mag = _normalize_abs(volume)
    z_peak = int(np.argmax(np.max(mag, axis=(1, 2))))
    axial = np.max(mag, axis=(1, 2))
    axial = axial / max(float(np.max(axial)), 1e-12)
    spatial = mag[z_peak]
    spatial = spatial / max(float(np.max(spatial)), 1e-12)
    approx = axial[:, None, None] * spatial[None, :, :]
    return float(np.sqrt(np.mean((mag - approx) ** 2)))

def _through_focus_rci_stack(
    config: SimulationConfig,
    grid,
    depth_um: np.ndarray,
    wavelength_nm: float,
    pad_factor: int,
) -> np.ndarray:
    field_n = grid.mask.shape[0] * pad_factor
    out = np.empty((len(depth_um), field_n, field_n), dtype=np.complex64)
    for i, z in enumerate(depth_um):
        pupil_i, _ = _make_oct_pupil(config, grid, wavelength_nm, physical_defocus_um=float(z))
        pupil_d, _ = _make_oct_pupil(config, grid, wavelength_nm, physical_defocus_um=float(z))
        ui = fraunhofer_psf_from_pupil(pupil_i, pad_factor=pad_factor)
        ud = fraunhofer_psf_from_pupil(pupil_d, pad_factor=pad_factor)
        out[i] = (ui * np.conj(ud)).astype(np.complex64)
    peak = float(np.max(np.abs(out)))
    if peak > 0:
        out = (out / peak).astype(np.complex64)
    return out

def _depth_sample_indices(count: int, decimation: int) -> np.ndarray:
    if count <= 0:
        return np.array([], dtype=np.int64)
    step = max(int(decimation), 1)
    if step == 1:
        return np.arange(count, dtype=np.int64)
    indices = np.arange(0, count, step, dtype=np.int64)
    anchors = np.array([0, count // 2, count - 1], dtype=np.int64)
    return np.unique(np.concatenate([indices, anchors]))

def _interpolate_complex_stack(
    x_full: np.ndarray,
    x_sampled: np.ndarray,
    sampled: np.ndarray,
) -> np.ndarray:
    if len(x_sampled) == len(x_full):
        return sampled
    flat = sampled.reshape(len(x_sampled), -1)
    out = np.empty((len(x_full), flat.shape[1]), dtype=np.complex128)
    x_sampled = np.asarray(x_sampled, dtype=np.float64)
    x_full = np.asarray(x_full, dtype=np.float64)
    order = np.argsort(x_sampled)
    xs = x_sampled[order]
    vals = flat[order]
    for j in range(flat.shape[1]):
        out[:, j] = (
            np.interp(x_full, xs, vals[:, j].real)
            + 1j * np.interp(x_full, xs, vals[:, j].imag)
        )
    return out.reshape((len(x_full), *sampled.shape[1:]))

def _full_spectral_rci_direct_psf(
    config: SimulationConfig,
    grid,
    depth_um: np.ndarray,
    kgrid,
    k_measured: np.ndarray,
    dispersion: np.ndarray,
    rolloff: float,
    pad_factor: int,
    scatterer_z_um: float,
    depth_decimation: int,
) -> dict:
    """Low-resolution direct PSF using h_RCI(x,y,z;k) before k-space FFT.

    For each axial sample, this simulates the spectral interferogram produced by
    a point at that depth, including the wavelength-dependent RCI kernel, then
    reconstructs it on the same k-linearized grid and samples the matched depth.
    """
    field_n = grid.mask.shape[0] * pad_factor
    computed_indices = _depth_sample_indices(len(depth_um), depth_decimation)
    computed_depth_um = depth_um[computed_indices]
    sampled = np.empty((len(computed_indices), field_n, field_n), dtype=np.complex128)
    sampled_h_rci = np.empty((len(computed_indices), len(kgrid.k), field_n, field_n), dtype=np.complex64)
    window = window_vector(len(kgrid.k), config.oct.window).astype(np.float64)
    spectral_weight = (
        rolloff
        * kgrid.source_spectrum.astype(np.float64)
        * window
        * np.exp(1j * dispersion.astype(np.float64))
    )

    for zi, object_depth_um in enumerate(computed_depth_um):
        object_depth_um = float(object_depth_um)
        spectral = np.empty((len(kgrid.k), field_n, field_n), dtype=np.complex128)
        depth_phase_um = float(scatterer_z_um)
        for ki, (k, wl_nm) in enumerate(zip(kgrid.k, kgrid.wavelength_nm)):
            defocus_um = object_depth_um - float(scatterer_z_um)
            pupil_i, _ = _make_oct_pupil(config, grid, wl_nm, physical_defocus_um=defocus_um)
            pupil_d, _ = _make_oct_pupil(config, grid, wl_nm, physical_defocus_um=defocus_um)
            ui = fraunhofer_psf_from_pupil(pupil_i, pad_factor=pad_factor)
            ud = fraunhofer_psf_from_pupil(pupil_d, pad_factor=pad_factor)
            h_rci = ui * np.conj(ud)
            sampled_h_rci[zi, ki] = h_rci.astype(np.complex64)
            phase = np.exp(1j * 2.0 * float(k) * depth_phase_um)
            spectral[ki] = h_rci * spectral_weight[ki] * phase

        k_lin, spectral_lin = _linearize_stack(k_measured, spectral)
        reconstructed = np.fft.fftshift(np.fft.fft(spectral_lin, axis=0), axes=0)
        recon_depth_um = np.fft.fftshift(fft_depth_axis_um(k_lin)).astype(np.float64)
        match = int(np.argmin(np.abs(recon_depth_um - object_depth_um)))
        sampled[zi] = reconstructed[match]

    interpolated = len(computed_indices) != len(depth_um)
    out = _interpolate_complex_stack(depth_um, computed_depth_um, sampled) if interpolated else sampled

    max_abs = float(np.max(np.abs(out)))
    if max_abs > 0:
        out = out / max_abs
    dx_um_per_k = np.array(
        [
            lateral_sample_pitch_um(float(wl_nm) * 1e-3, config.objective.NA_nominal, pad_factor)
            for wl_nm in kgrid.wavelength_nm
        ],
        dtype=np.float64,
    )
    profiles = _axial_profile_products(out)
    return {
        "psf": out.astype(np.complex64),
        "sampled_h_rci": sampled_h_rci,
        "computed_depth_indices": computed_indices.astype(np.int64),
        "computed_depth_um": computed_depth_um.astype(np.float64),
        "depth_decimation": np.array(max(int(depth_decimation), 1), dtype=np.int64),
        "interpolated": np.array(interpolated, dtype=np.bool_),
        "k_true": kgrid.k.astype(np.float64),
        "k_measured": k_measured.astype(np.float64),
        "k_linear": k_lin.astype(np.float64),
        "wavelength_nm": kgrid.wavelength_nm.astype(np.float64),
        "dx_um_per_k": dx_um_per_k,
        "source": kgrid.source_spectrum.astype(np.float64),
        "window": window,
        "dispersion_phase_rad": dispersion.astype(np.float64),
        "phase_sign": "+exp(i 2 k scatterer_z)",
        "depth_phase_origin_um": np.array(scatterer_z_um, dtype=np.float64),
        "normalization": "h_RCI is stored before source/window/dispersion/depth phase; reconstructed PSF is globally peak normalized after optional depth interpolation",
        "model_contract": "low-resolution full spectral RCI direct PSF from h_RCI(x,y,z;k), source spectrum, fixed scatterer-depth phase, measured k linearization, dispersion phase, and k-space FFT reconstruction",
        **profiles,
    }

def simulate_oct_raw_direct(
    config: SimulationConfig,
    N: int = 128,
    scatterer_z_um: float = 0.0,
    scatterer_x_um: float = 0.0,
    scatterer_y_um: float = 0.0,
    k_samples: int | None = None,
):
    """Independent direct OCT forward model.

    This function intentionally does NOT import microscope_forward or any theory conversion module.
    It builds OCT fields from the shared pupil, source spectrum, reference phase and sample scatterer.
    """
    kgrid = make_oct_k_grid(config, k_samples or config.oct.spectrometer_pixels)
    grid = make_pupil_grid(N, config.objective.NA_nominal)
    Es = np.zeros_like(kgrid.k, dtype=np.complex128)
    h_rci_sample = np.zeros_like(kgrid.k, dtype=np.complex128)
    rolloff = _rolloff_amplitude(config, scatterer_z_um)
    dispersion = differential_dispersion_phase(config, kgrid.k)
    scatterer_diameter_um = effective_scatterer_diameter_um(config)
    scatter_amp = sample_scattering_amplitude(config, kgrid.wavelength_nm)

    for i, (k, wl_nm, s, phi) in enumerate(zip(kgrid.k, kgrid.wavelength_nm, kgrid.source_spectrum, dispersion)):
        Pi, _ = _make_oct_pupil(config, grid, wl_nm, physical_defocus_um=scatterer_z_um)
        Pd, _ = _make_oct_pupil(config, grid, wl_nm, physical_defocus_um=scatterer_z_um)
        ui = fraunhofer_psf_from_pupil(Pi, pad_factor=1)
        ud = fraunhofer_psf_from_pupil(Pd, pad_factor=1)
        dx_um = lateral_sample_pitch_um(wl_nm * 1e-3, config.objective.NA_nominal, pad_factor=1)
        h_rci_sample[i] = _sample_scatterer_rci(
            config,
            ui * np.conj(ud),
            float(k),
            scatterer_x_um,
            scatterer_y_um,
            dx_um,
            scatterer_diameter_um,
        )
        Es[i] = (
            config.oct.sample_amplitude
            * rolloff
            * s
            * scatter_amp[i]
            * h_rci_sample[i]
            * np.exp(1j * (2 * k * scatterer_z_um + phi))
        )

    Er = config.oct.reference_amplitude * np.ones_like(Es)
    background = np.abs(Er) ** 2 + np.abs(Es) ** 2
    interference = 2 * np.real(np.conj(Er) * Es)
    raw = background + interference
    k_measured = apply_k_linearization_error(config, kgrid.k)
    k_linear = np.linspace(float(np.min(k_measured)), float(np.max(k_measured)), len(k_measured))
    center_coeffs = dichroic_path_coefficients(config, "oct_reflection", config.oct.center_wavelength_nm)
    return {
        "raw": raw.astype(np.float64),
        "background": background.astype(np.float64),
        "interference": interference.astype(np.float64),
        "reference_field": Er.astype(np.complex128),
        "sample_field": Es.astype(np.complex128),
        "pixel_index": kgrid.pixel_index.astype(np.int64),
        "k": k_measured.astype(np.float64),
        "k_linear": k_linear.astype(np.float64),
        "k_true": kgrid.k.astype(np.float64),
        "wavelength_nm": kgrid.wavelength_nm,
        "source": kgrid.source_spectrum,
        "scatterer_z_um": np.array(scatterer_z_um, dtype=np.float64),
        "scatterer_x_um": np.array(scatterer_x_um, dtype=np.float64),
        "scatterer_y_um": np.array(scatterer_y_um, dtype=np.float64),
        "scatterer_diameter_um": np.array(scatterer_diameter_um, dtype=np.float64),
        "scatterer_model": config.sample.scatterer_model,
        "measurement_path": config.sample.measurement_path,
        "sample_scattering_amplitude": scatter_amp.astype(np.complex128),
        "h_rci_sample": h_rci_sample.astype(np.complex128),
        "dichroic_path": "oct_reflection",
        "dichroic_s_amplitude": np.array(center_coeffs.s_amplitude, dtype=np.float64),
        "dichroic_p_amplitude": np.array(center_coeffs.p_amplitude, dtype=np.float64),
        "dichroic_s_phase_rad": np.array(center_coeffs.s_phase_rad, dtype=np.float64),
        "dichroic_p_phase_rad": np.array(center_coeffs.p_phase_rad, dtype=np.float64),
        "rolloff_amplitude": np.array(rolloff, dtype=np.float64),
        "dispersion_phase_rad": dispersion.astype(np.float64),
    }

def simulate_oct_psf_direct(
    config: SimulationConfig,
    N: int = 64,
    pad_factor: int = 2,
    k_samples: int | None = None,
    scatterer_z_um: float = 0.0,
    full_spectral_rci: bool = False,
    full_spectral_rci_depth_decimation: int | None = None,
):
    """Generate an independent effective OCT PSF by stacking RCI kernels over k."""
    kgrid = make_oct_k_grid(config, k_samples or config.oct.spectrometer_pixels)
    grid = make_pupil_grid(N, config.objective.NA_nominal)
    field_n = N * pad_factor
    rci_stack = np.empty((len(kgrid.k), field_n, field_n), dtype=np.complex64)

    for i, (k, wl_nm) in enumerate(zip(kgrid.k, kgrid.wavelength_nm)):
        pupil_i, _ = _make_oct_pupil(config, grid, wl_nm)
        pupil_d, _ = _make_oct_pupil(config, grid, wl_nm)
        ui = fraunhofer_psf_from_pupil(pupil_i, pad_factor=pad_factor)
        ud = fraunhofer_psf_from_pupil(pupil_d, pad_factor=pad_factor)
        phase = np.exp(1j * 2 * k * scatterer_z_um)
        rci_stack[i] = (ui * np.conj(ud) * phase).astype(np.complex64)

    rolloff = _rolloff_amplitude(config, scatterer_z_um)
    dispersion = differential_dispersion_phase(config, kgrid.k)
    spectral = (
        rci_stack.astype(np.complex128)
        * rolloff
        * kgrid.source_spectrum[:, None, None]
        * window_vector(len(kgrid.k), config.oct.window)[:, None, None]
        * np.exp(1j * dispersion)[:, None, None]
    )
    k_measured = apply_k_linearization_error(config, kgrid.k)
    k_lin, spectral_lin = _linearize_stack(k_measured, spectral)
    spectral_psf = np.fft.fftshift(np.fft.fft(spectral_lin, axis=0), axes=0)
    depth_um = np.fft.fftshift(fft_depth_axis_um(k_lin)).astype(np.float64)
    rci_through_focus = _through_focus_rci_stack(
        config,
        grid,
        depth_um - float(scatterer_z_um),
        config.oct.center_wavelength_nm,
        pad_factor,
    ).astype(np.complex128)
    cy = spectral_psf.shape[1] // 2
    cx = spectral_psf.shape[2] // 2
    axial_gate = spectral_psf[:, cy, cx]
    psf = axial_gate[:, None, None] * rci_through_focus
    max_abs = float(np.max(np.abs(psf)))
    if max_abs > 0:
        psf = psf / max_abs
    spectral_rci_direct_psf = None
    spectral_rci_meta = None
    if full_spectral_rci:
        spectral_rci_meta = _full_spectral_rci_direct_psf(
            config,
            grid,
            depth_um,
            kgrid,
            k_measured,
            dispersion,
            rolloff,
            pad_factor,
            scatterer_z_um,
            full_spectral_rci_depth_decimation or config.oct.full_spectral_rci_depth_decimation,
        )
        spectral_rci_direct_psf = spectral_rci_meta["psf"]
        psf = spectral_rci_direct_psf.astype(np.complex128)

    wl_um = config.oct.center_wavelength_nm * 1e-3
    dx_um = lateral_sample_pitch_um(wl_um, config.objective.NA_nominal, pad_factor)
    result = {
        "psf": psf.astype(np.complex64),
        "magnitude": np.abs(psf).astype(np.float32),
        "rci_through_focus": rci_through_focus.astype(np.complex64),
        "separability_nrmse": np.array(_separability_nrmse(psf), dtype=np.float64),
        "full_spectral_rci": np.array(full_spectral_rci, dtype=np.bool_),
        "depth_um": depth_um,
        "x_um": centered_axis_um(field_n, dx_um),
        "y_um": centered_axis_um(field_n, dx_um),
        "dx_um": np.array(dx_um, dtype=np.float64),
        "k_linear": k_lin.astype(np.float64),
        "k_true": kgrid.k.astype(np.float64),
        "k_measured": k_measured.astype(np.float64),
        "source": np.interp(k_lin, np.sort(kgrid.k), kgrid.source_spectrum[np.argsort(kgrid.k)]).astype(np.float64),
        "rolloff_amplitude": np.array(rolloff, dtype=np.float64),
        "dispersion_phase_rad": dispersion.astype(np.float64),
    }
    if spectral_rci_direct_psf is not None:
        result["spectral_rci_direct_psf"] = spectral_rci_direct_psf.astype(np.complex64)
        result["spectral_rci_depth_um"] = depth_um.astype(np.float64)
        result["spectral_rci_sampled_h_rci"] = spectral_rci_meta["sampled_h_rci"]
        result["spectral_rci_computed_depth_indices"] = spectral_rci_meta["computed_depth_indices"]
        result["spectral_rci_computed_depth_um"] = spectral_rci_meta["computed_depth_um"]
        result["spectral_rci_k_true"] = spectral_rci_meta["k_true"]
        result["spectral_rci_k_measured"] = spectral_rci_meta["k_measured"]
        result["spectral_rci_k_linear"] = spectral_rci_meta["k_linear"]
        result["spectral_rci_wavelength_nm"] = spectral_rci_meta["wavelength_nm"]
        result["spectral_rci_dx_um_per_k"] = spectral_rci_meta["dx_um_per_k"]
        result["spectral_rci_source"] = spectral_rci_meta["source"]
        result["spectral_rci_window"] = spectral_rci_meta["window"]
        result["spectral_rci_dispersion_phase_rad"] = spectral_rci_meta["dispersion_phase_rad"]
        result["spectral_rci_axial_profile_max"] = spectral_rci_meta["axial_profile_max"]
        result["spectral_rci_axial_profile_on_axis"] = spectral_rci_meta["axial_profile_on_axis"]
        result["spectral_rci_axial_profile_integrated"] = spectral_rci_meta["axial_profile_integrated"]
        result["spectral_rci_phase_sign"] = spectral_rci_meta["phase_sign"]
        result["spectral_rci_depth_phase_origin_um"] = spectral_rci_meta["depth_phase_origin_um"]
        result["spectral_rci_normalization"] = spectral_rci_meta["normalization"]
        result["spectral_rci_model_contract"] = spectral_rci_meta["model_contract"]
        result["full_spectral_rci_depth_decimation"] = spectral_rci_meta["depth_decimation"]
        result["full_spectral_rci_interpolated"] = spectral_rci_meta["interpolated"]
    return result
