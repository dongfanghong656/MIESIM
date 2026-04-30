from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from .grids import PupilGrid

def clear_circular_pupil(grid: PupilGrid) -> NDArray[np.complex128]:
    return np.where(grid.mask, 1.0 + 0.0j, 0.0 + 0.0j)

def zernike_defocus(grid: PupilGrid, coeff_um: float, wavelength_um: float) -> NDArray[np.complex128]:
    r2 = grid.rho_x**2 + grid.rho_y**2
    phase = 2 * np.pi * coeff_um * (2 * r2 - 1) / wavelength_um
    out = np.exp(1j * phase)
    return np.where(grid.mask, out, 0.0)

def zernike_astigmatism_45(grid: PupilGrid, coeff_um: float, wavelength_um: float) -> NDArray[np.complex128]:
    phase = 2 * np.pi * coeff_um * (2.0 * grid.rho_x * grid.rho_y) / wavelength_um
    return np.where(grid.mask, np.exp(1j * phase), 0.0)

def gaussian_underfill(grid: PupilGrid, fill_ratio: float) -> NDArray[np.float64]:
    r2 = grid.rho_x**2 + grid.rho_y**2
    amp = np.exp(-r2 / max(fill_ratio, 1e-6) ** 2)
    return np.where(grid.mask, amp, 0.0)

def clipped_gaussian_underfill(
    grid: PupilGrid, fill_ratio: float, clip_radius: float = 1.0
) -> NDArray[np.float64]:
    r2 = grid.rho_x**2 + grid.rho_y**2
    support = grid.mask & (r2 <= clip_radius**2)
    amp = np.exp(-r2 / max(fill_ratio, 1e-6) ** 2)
    return np.where(support, amp, 0.0)

def tilt_phase(grid: PupilGrid, tilt_x_waves: float = 0.0, tilt_y_waves: float = 0.0) -> NDArray[np.complex128]:
    phase = 2 * np.pi * (tilt_x_waves * grid.rho_x + tilt_y_waves * grid.rho_y)
    return np.where(grid.mask, np.exp(1j * phase), 0.0)

def build_shared_pupil(grid: PupilGrid, wavelength_um: float, fill_ratio: float = 0.82,
                       defocus_um: float = 0.0,
                       astig45_um: float = 0.0,
                       tilt_x_waves: float = 0.0,
                       tilt_y_waves: float = 0.0) -> NDArray[np.complex128]:
    amp = gaussian_underfill(grid, fill_ratio)
    ph = zernike_defocus(grid, defocus_um, wavelength_um)
    ph *= zernike_astigmatism_45(grid, astig45_um, wavelength_um)
    ph *= tilt_phase(grid, tilt_x_waves, tilt_y_waves)
    return (amp * ph).astype(np.complex128)
