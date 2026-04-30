from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray

@dataclass(frozen=True)
class PupilGrid:
    rho_x: NDArray[np.float64]
    rho_y: NDArray[np.float64]
    mask: NDArray[np.bool_]
    dx: float
    NA: float

@dataclass(frozen=True)
class ImageGrid:
    x_um: NDArray[np.float64]
    y_um: NDArray[np.float64]
    z_um: NDArray[np.float64]

@dataclass(frozen=True)
class KGrid:
    k: NDArray[np.float64]
    wavelength_nm: NDArray[np.float64]
    source_spectrum: NDArray[np.float64]
    pixel_index: NDArray[np.int64]

def make_pupil_grid(N: int, NA: float) -> PupilGrid:
    if N < 8:
        raise ValueError("Pupil grid N must be at least 8.")
    if NA <= 0:
        raise ValueError("NA must be positive.")
    axis = np.linspace(-1.0, 1.0, N, endpoint=False)
    rx, ry = np.meshgrid(axis, axis, indexing="xy")
    mask = (rx**2 + ry**2) <= 1.0
    return PupilGrid(rx, ry, mask, dx=float(axis[1] - axis[0]), NA=NA)

def make_image_grid(nx=128, ny=128, nz=81, dx_um=0.25, dz_um=0.25) -> ImageGrid:
    x = (np.arange(nx) - nx // 2) * dx_um
    y = (np.arange(ny) - ny // 2) * dx_um
    z = (np.arange(nz) - nz // 2) * dz_um
    return ImageGrid(x, y, z)

def make_gaussian_k_grid(center_nm=850.0, bandwidth_nm=55.0, n=2048) -> KGrid:
    if center_nm <= 0 or bandwidth_nm <= 0:
        raise ValueError("Center wavelength and bandwidth must be positive.")
    if n < 16:
        raise ValueError("Spectrometer grid must contain at least 16 samples.")
    sigma_nm = bandwidth_nm / 2.355
    wl = np.linspace(center_nm - 3.0 * sigma_nm, center_nm + 3.0 * sigma_nm, n)
    k = 2 * np.pi / (wl * 1e-3)  # rad/um
    s = np.exp(-0.5 * ((wl - center_nm) / sigma_nm) ** 2)
    return KGrid(k=k, wavelength_nm=wl, source_spectrum=s / np.max(s), pixel_index=np.arange(n))

def lateral_sample_pitch_um(wavelength_um: float, NA: float, pad_factor: int = 1) -> float:
    if wavelength_um <= 0 or NA <= 0 or pad_factor <= 0:
        raise ValueError("wavelength, NA, and pad_factor must be positive.")
    return wavelength_um / (2.0 * NA * pad_factor)

def centered_axis_um(n: int, dx_um: float) -> NDArray[np.float64]:
    return (np.arange(n) - n // 2) * dx_um

def fft_depth_axis_um(k_linear: NDArray[np.float64]) -> NDArray[np.float64]:
    if k_linear.ndim != 1 or len(k_linear) < 2:
        raise ValueError("k_linear must be a one-dimensional grid.")
    dk = float(np.mean(np.diff(k_linear)))
    return np.fft.fftfreq(len(k_linear), d=abs(dk)) * np.pi
