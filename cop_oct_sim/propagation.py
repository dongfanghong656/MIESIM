from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

def _center_pad(pupil: NDArray[np.complex128], pad_factor: int) -> NDArray[np.complex128]:
    if pad_factor < 1:
        raise ValueError("pad_factor must be >= 1.")
    if pad_factor == 1:
        return pupil
    ny, nx = pupil.shape
    out = np.zeros((ny * pad_factor, nx * pad_factor), dtype=pupil.dtype)
    y0 = (out.shape[0] - ny) // 2
    x0 = (out.shape[1] - nx) // 2
    out[y0:y0 + ny, x0:x0 + nx] = pupil
    return out

def fraunhofer_psf_from_pupil(
    pupil: NDArray[np.complex128],
    normalize: bool = False,
    pad_factor: int = 1,
) -> NDArray[np.complex128]:
    padded = _center_pad(pupil, pad_factor)
    field = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(padded), norm="ortho"))
    if normalize:
        m = np.max(np.abs(field))
        if m > 0:
            field = field / m
    return field

def physical_defocus_phase(
    grid,
    z_um: float,
    wavelength_um: float,
    NA: float,
) -> NDArray[np.complex128]:
    r2 = grid.rho_x**2 + grid.rho_y**2
    phase = np.pi * float(z_um) * (float(NA) ** 2) * r2 / max(float(wavelength_um), 1e-12)
    return np.where(grid.mask, np.exp(1j * phase), 0.0).astype(np.complex128)

def simple_defocus_stack(pupil: NDArray[np.complex128], z_um, wavelength_um: float,
                         NA: float, pad_factor: int = 1) -> NDArray[np.complex128]:
    ny, nx = pupil.shape
    from .grids import make_pupil_grid

    grid = make_pupil_grid(ny, NA)
    z_um = np.asarray(z_um, dtype=float)
    out_shape = (len(z_um), ny * pad_factor, nx * pad_factor)
    out = np.empty(out_shape, dtype=np.complex128)
    for i, z in enumerate(z_um):
        phase = physical_defocus_phase(grid, float(z), wavelength_um, NA)
        out[i] = fraunhofer_psf_from_pupil(pupil * phase, normalize=False, pad_factor=pad_factor)
    return out

def propagation_energy_ratio(
    pupil: NDArray[np.complex128], field: NDArray[np.complex128]
) -> float:
    denom = float(np.sum(np.abs(pupil) ** 2))
    if denom == 0:
        return float("nan")
    return float(np.sum(np.abs(field) ** 2) / denom)
