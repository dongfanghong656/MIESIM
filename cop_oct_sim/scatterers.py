from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import convolve

from .config_schema import SampleConfig, SimulationConfig
from .measurement import finite_sphere_kernel
from .spectrometer import first_column, read_numeric_table

SUPPORTED_SCATTERER_MODELS = {
    "delta",
    "gaussian",
    "uniform_sphere_projection",
    "rayleigh",
    "rayleigh_mie_lookup",
}
SUPPORTED_MEASUREMENT_PATHS = {"epi_backscatter", "oblique_darkfield"}

def _model_name(sample: SampleConfig) -> str:
    model = sample.scatterer_model.lower()
    if model not in SUPPORTED_SCATTERER_MODELS:
        raise ValueError(f"Unsupported scatterer_model: {sample.scatterer_model}")
    return model

def _measurement_path(sample: SampleConfig) -> str:
    path = sample.measurement_path.lower()
    if path not in SUPPORTED_MEASUREMENT_PATHS:
        raise ValueError(f"Unsupported measurement_path: {sample.measurement_path}")
    return path

def effective_scatterer_diameter_um(config: SimulationConfig, override_um: float | None = None) -> float:
    if override_um is not None:
        return float(override_um)
    if config.sample.scatterer_diameter_um is not None:
        return float(config.sample.scatterer_diameter_um)
    return float(config.microscope.bead_diameter_um)

def finite_target_fwhm_ratio(diameter_um: float, reference_fwhm_um: float) -> float:
    reference = max(float(reference_fwhm_um), 1e-12)
    return float(max(float(diameter_um), 0.0) / reference)

def _gaussian_kernel(
    diameter_um: float,
    lateral_pixel_um: float,
    axial_pixel_um: float,
) -> NDArray[np.float64]:
    sigma_um = max(float(diameter_um), 0.0) / 2.355
    if sigma_um <= 0:
        return np.ones((1, 1, 1), dtype=np.float64)
    rx = max(int(np.ceil(3.0 * sigma_um / max(lateral_pixel_um, 1e-12))), 1)
    rz = max(int(np.ceil(3.0 * sigma_um / max(axial_pixel_um, 1e-12))), 1)
    x = np.arange(-rx, rx + 1, dtype=np.float64) * lateral_pixel_um
    y = x.copy()
    z = np.arange(-rz, rz + 1, dtype=np.float64) * axial_pixel_um
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    kernel = np.exp(-0.5 * (xx**2 + yy**2 + zz**2) / max(sigma_um**2, 1e-24))
    return kernel / max(float(np.sum(kernel)), 1e-12)

def scatterer_kernel(
    sample: SampleConfig,
    diameter_um: float,
    lateral_pixel_um: float,
    axial_pixel_um: float,
) -> NDArray[np.float64]:
    model = _model_name(sample)
    if model == "delta" or diameter_um <= 0:
        return np.ones((1, 1, 1), dtype=np.float64)
    if model == "gaussian":
        return _gaussian_kernel(diameter_um, lateral_pixel_um, axial_pixel_um)
    return finite_sphere_kernel(diameter_um, lateral_pixel_um, axial_pixel_um)

def convolve_scatterer_stack(
    intensity_stack: NDArray,
    sample: SampleConfig,
    diameter_um: float,
    lateral_pixel_um: float,
    axial_pixel_um: float,
) -> NDArray:
    kernel = scatterer_kernel(sample, diameter_um, lateral_pixel_um, axial_pixel_um)
    if kernel.size == 1:
        return intensity_stack
    return convolve(intensity_stack, kernel, mode="nearest")

def scatterer_volume_offsets(
    sample: SampleConfig,
    diameter_um: float,
    lateral_pixel_um: float,
    axial_pixel_um: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    kernel = scatterer_kernel(sample, diameter_um, lateral_pixel_um, axial_pixel_um)
    if kernel.size == 1:
        return (
            np.array([0.0], dtype=np.float64),
            np.array([0.0], dtype=np.float64),
            np.array([0.0], dtype=np.float64),
            np.array([1.0], dtype=np.float64),
        )
    threshold = float(np.max(kernel)) * 1e-4
    zi, yi, xi = np.nonzero(kernel >= threshold)
    weights = kernel[zi, yi, xi].astype(np.float64)
    weights /= max(float(np.sum(weights)), 1e-12)
    cz = kernel.shape[0] // 2
    cy = kernel.shape[1] // 2
    cx = kernel.shape[2] // 2
    dz = (zi.astype(np.float64) - cz) * axial_pixel_um
    dy = (yi.astype(np.float64) - cy) * lateral_pixel_um
    dx = (xi.astype(np.float64) - cx) * lateral_pixel_um
    return dz, dy, dx, weights

def _rayleigh_relative_amplitude(config: SimulationConfig, wavelength_nm: NDArray[np.float64]) -> NDArray[np.float64]:
    diameter_um = effective_scatterer_diameter_um(config)
    radius_um = max(0.5 * diameter_um, 1e-12)
    wavelength_um = np.maximum(wavelength_nm * 1e-3, 1e-12)
    medium_n = max(float(config.sample.medium_refractive_index), 1e-12)
    rel_n = float(config.sample.particle_refractive_index) / medium_n
    contrast = abs((rel_n**2 - 1.0) / (rel_n**2 + 2.0))
    k_medium = 2.0 * np.pi * medium_n / wavelength_um
    amp = contrast * (radius_um**3) * (k_medium**2)
    peak = float(np.max(np.abs(amp))) if amp.size else 0.0
    return amp / peak if peak > 0 else np.ones_like(wavelength_nm, dtype=np.float64)

def _lookup_scattering_amplitude(path: str, wavelength_nm: NDArray[np.float64]) -> NDArray[np.complex128]:
    table = read_numeric_table(path)
    source_wavelength = first_column(table, ("wavelength_nm", "lambda_nm", "wavelength", "lambda"))
    amplitude = first_column(table, ("amplitude", "amp", "scattering_amplitude"), fallback=1 if "col1" in table else 0)
    if any(name in table for name in ("phase_rad", "phase", "phase_radian")):
        phase = first_column(table, ("phase_rad", "phase", "phase_radian"))
    else:
        phase = np.zeros_like(amplitude, dtype=np.float64)
    order = np.argsort(source_wavelength)
    amp_i = np.interp(
        wavelength_nm,
        source_wavelength[order],
        amplitude[order],
        left=amplitude[order][0],
        right=amplitude[order][-1],
    )
    phase_i = np.interp(
        wavelength_nm,
        source_wavelength[order],
        phase[order],
        left=phase[order][0],
        right=phase[order][-1],
    )
    return amp_i.astype(np.float64) * np.exp(1j * phase_i.astype(np.float64))

def sample_scattering_amplitude(config: SimulationConfig, wavelength_nm) -> NDArray[np.complex128]:
    wavelength = np.asarray(wavelength_nm, dtype=np.float64)
    model = _model_name(config.sample)
    if model == "rayleigh_mie_lookup":
        if not config.sample.scattering_lookup_file:
            raise NotImplementedError("rayleigh_mie_lookup requires sample.scattering_lookup_file.")
        return _lookup_scattering_amplitude(config.sample.scattering_lookup_file, wavelength)
    if model == "rayleigh":
        return _rayleigh_relative_amplitude(config, wavelength).astype(np.complex128)
    return np.ones_like(wavelength, dtype=np.complex128)

def scattering_measurement_stack(
    intensity_stack: NDArray,
    config: SimulationConfig,
    wavelength_nm: float,
) -> dict[str, NDArray | complex | str]:
    path = _measurement_path(config.sample)
    amp = complex(sample_scattering_amplitude(config, np.array([wavelength_nm], dtype=np.float64))[0])
    field = np.sqrt(np.clip(np.asarray(intensity_stack, dtype=np.float64), 0.0, None)) * amp
    background = complex(float(config.sample.background_amplitude))

    if path == "epi_backscatter":
        bg_field = np.full_like(field, background, dtype=np.complex128)
        raw = np.abs(bg_field + field) ** 2
        bg = np.full_like(raw, abs(background) ** 2, dtype=np.float64)
    else:
        leak = background * float(config.sample.darkfield_background_fraction)
        raw = np.abs(field) ** 2 + abs(leak) ** 2
        bg = np.full_like(raw, abs(leak) ** 2, dtype=np.float64)

    corrected = raw - bg if config.sample.background_subtract else raw
    return {
        "raw": raw.astype(np.float32),
        "background": bg.astype(np.float32),
        "corrected": corrected.astype(np.float32),
        "amplitude": amp,
        "measurement_path": path,
    }
