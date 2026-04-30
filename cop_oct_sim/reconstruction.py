from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from scipy.signal.windows import tukey
from .grids import fft_depth_axis_um

def window_vector(n: int, name: str = "hann") -> NDArray[np.float64]:
    name = name.lower()
    if name in {"rect", "rectangle", "none"}:
        return np.ones(n, dtype=np.float64)
    if name == "hann":
        return np.hanning(n)
    if name == "tukey":
        return tukey(n, alpha=0.25)
    raise ValueError(f"Unknown window: {name}")

def _interp_complex(x_new: NDArray, x: NDArray, y: NDArray) -> NDArray:
    if np.iscomplexobj(y):
        return np.interp(x_new, x, np.real(y)) + 1j * np.interp(x_new, x, np.imag(y))
    return np.interp(x_new, x, y)

def reconstruct_sd_oct(raw: NDArray, k: NDArray, dispersion_phase: NDArray | None = None,
                       window: str = "hann",
                       source_spectrum: NDArray | None = None,
                       return_depth: bool = False) -> NDArray[np.complex128] | dict:
    order = np.argsort(k)
    k_sorted = k[order]
    raw_sorted = raw[order]
    k_lin = np.linspace(k_sorted.min(), k_sorted.max(), len(k_sorted))
    raw_lin = _interp_complex(k_lin, k_sorted, raw_sorted).astype(np.complex128)
    if source_spectrum is not None:
        source_sorted = np.asarray(source_spectrum)[order]
        source_lin = np.interp(k_lin, k_sorted, source_sorted)
        raw_lin = raw_lin / np.maximum(source_lin, 1e-9)
    raw_lin = raw_lin - np.mean(raw_lin)
    if dispersion_phase is not None:
        phase = np.interp(k_lin, k_sorted, dispersion_phase[order])
        raw_lin = raw_lin * np.exp(1j * phase)
    raw_lin = raw_lin * window_vector(len(raw_lin), window)
    complex_depth = np.fft.fftshift(np.fft.fft(raw_lin))
    if not return_depth:
        return complex_depth
    depth_um = np.fft.fftshift(fft_depth_axis_um(k_lin))
    return {"complex": complex_depth, "depth_um": depth_um, "k_linear": k_lin}
