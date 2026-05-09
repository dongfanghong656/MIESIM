from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from .config_schema import SimulationConfig
from .grids import fft_depth_axis_um
from .metrics import normalize
from .reconstruction import window_vector
from .spectrometer import apply_k_linearization_error, differential_dispersion_phase, make_oct_k_grid

def focal_plane_from_zstack(stack: NDArray) -> NDArray:
    arr = np.asarray(stack)
    if arr.ndim != 3:
        raise ValueError("Expected a z-stack with shape [z, y, x].")
    z = int(np.argmax(np.max(np.abs(arr), axis=(1, 2))))
    return arr[z]

def path_e_separable_psf(
    microscope_spatial: NDArray,
    axial_profile: NDArray,
    microscope_mode: str = "widefield",
) -> NDArray[np.complex64]:
    """Path E baseline: microscope spatial term times OCT axial gate.

    This module is intentionally outside `oct_forward.py`; it is a model under
    test, not a direct-truth generator.
    """
    spatial = np.asarray(microscope_spatial)
    if spatial.ndim == 3:
        spatial = focal_plane_from_zstack(spatial)
    spatial = np.abs(spatial).astype(np.float64)
    if microscope_mode in {"reflectance_confocal", "confocal"}:
        spatial = np.sqrt(np.maximum(spatial, 0.0))
    spatial = normalize(spatial)
    axial = normalize(np.asarray(axial_profile).astype(np.complex128))
    converted = axial[:, None, None] * spatial[None, :, :]
    max_abs = float(np.max(np.abs(converted)))
    if max_abs > 0:
        converted = converted / max_abs
    return converted.astype(np.complex64)

def axial_profile_from_direct_psf(psf: NDArray) -> NDArray[np.complex128]:
    arr = np.asarray(psf)
    if arr.ndim != 3:
        raise ValueError("Expected direct PSF volume [z, y, x].")
    mag = np.abs(arr)
    y, x = np.unravel_index(int(np.argmax(np.max(mag, axis=0))), mag.shape[1:])
    return arr[:, y, x].astype(np.complex128)

def predict_axial_gate_from_source(
    config: SimulationConfig,
    k_samples: int | None = None,
    measured_k: NDArray | None = None,
    window: str | None = None,
) -> dict:
    """Predict the OCT axial coherence gate without reading direct OCT PSF truth."""
    kgrid = make_oct_k_grid(config, k_samples or config.oct.spectrometer_pixels)
    k_measured = np.asarray(measured_k, dtype=np.float64) if measured_k is not None else apply_k_linearization_error(config, kgrid.k)
    order = np.argsort(k_measured)
    k_sorted = k_measured[order]
    k_linear = np.linspace(float(k_sorted.min()), float(k_sorted.max()), len(k_sorted))
    source_sorted = kgrid.source_spectrum[order]
    dispersion_sorted = differential_dispersion_phase(config, kgrid.k)[order]
    source_linear = np.interp(k_linear, k_sorted, source_sorted)
    dispersion_linear = np.interp(k_linear, k_sorted, dispersion_sorted)
    spectral = (
        source_linear.astype(np.complex128)
        * window_vector(len(k_linear), window or config.oct.window)
        * np.exp(1j * dispersion_linear)
    )
    axial = np.fft.fftshift(np.fft.fft(spectral))
    return {
        "axial_profile": axial.astype(np.complex128),
        "magnitude": np.abs(axial).astype(np.float32),
        "depth_um": np.fft.fftshift(fft_depth_axis_um(k_linear)).astype(np.float64),
        "k_linear": k_linear.astype(np.float64),
        "source_linear": source_linear.astype(np.float64),
        "dispersion_phase_rad": dispersion_linear.astype(np.float64),
    }
