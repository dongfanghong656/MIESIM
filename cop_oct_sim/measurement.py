from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import convolve

def finite_sphere_kernel(
    bead_diameter_um: float,
    lateral_pixel_um: float,
    axial_pixel_um: float | None = None,
) -> NDArray[np.float64]:
    """Return a normalized finite sphere kernel on the simulation grid.

    With ``axial_pixel_um`` this is a 3D uniform-volume bead. Without it, the
    kernel is the lateral projected path length through the bead.
    """
    radius_um = 0.5 * float(bead_diameter_um)
    lateral_pixel_um = max(float(lateral_pixel_um), 1e-12)
    if radius_um <= 0:
        shape = (1, 1, 1) if axial_pixel_um is not None else (1, 1)
        return np.ones(shape, dtype=np.float64)

    rx = max(int(np.ceil(radius_um / lateral_pixel_um)), 1)
    x = np.arange(-rx, rx + 1, dtype=np.float64) * lateral_pixel_um
    y = x.copy()

    if axial_pixel_um is None:
        yy, xx = np.meshgrid(y, x, indexing="ij")
        r2 = xx**2 + yy**2
        path_length = 2.0 * np.sqrt(np.maximum(radius_um**2 - r2, 0.0))
        kernel = np.where(r2 <= radius_um**2, path_length, 0.0)
    else:
        axial_pixel_um = max(float(axial_pixel_um), 1e-12)
        rz = max(int(np.ceil(radius_um / axial_pixel_um)), 1)
        z = np.arange(-rz, rz + 1, dtype=np.float64) * axial_pixel_um
        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        kernel = (xx**2 + yy**2 + zz**2 <= radius_um**2).astype(np.float64)

    total = float(np.sum(kernel))
    if total <= 0:
        center = tuple(s // 2 for s in kernel.shape)
        kernel[center] = 1.0
        total = 1.0
    return kernel / total

def convolve_finite_bead(intensity: NDArray, bead_diameter_um: float, pixel_um: float) -> NDArray:
    kernel = finite_sphere_kernel(bead_diameter_um, pixel_um)
    if kernel.size == 1:
        return intensity
    return convolve(intensity, kernel, mode="nearest")

def convolve_finite_sphere_stack(
    intensity_stack: NDArray,
    bead_diameter_um: float,
    lateral_pixel_um: float,
    axial_pixel_um: float,
) -> NDArray:
    kernel = finite_sphere_kernel(bead_diameter_um, lateral_pixel_um, axial_pixel_um)
    if kernel.size == 1:
        return intensity_stack
    return convolve(intensity_stack, kernel, mode="nearest")

def integrate_camera_pixels(image: NDArray, pixel_um: float, sample_um: float) -> NDArray:
    width_px = max(int(round(pixel_um / max(sample_um, 1e-12))), 1)
    if width_px <= 1:
        return image
    arr = np.asarray(image)
    ny, nx = arr.shape
    ny_crop = max((ny // width_px) * width_px, width_px)
    nx_crop = max((nx // width_px) * width_px, width_px)
    ny_crop = min(ny_crop, ny)
    nx_crop = min(nx_crop, nx)
    cropped = arr[:ny_crop, :nx_crop]
    return cropped.reshape(ny_crop // width_px, width_px, nx_crop // width_px, width_px).mean(axis=(1, 3))

def add_camera_noise(image: NDArray, read_noise: float = 2.0, gain: float = 1.0,
                     rng: np.random.Generator | None = None) -> NDArray:
    rng = rng or np.random.default_rng()
    gain = max(float(gain), 1e-12)
    electrons = np.clip(image * gain, 0, None)
    noisy_e = rng.poisson(electrons).astype(np.float64)
    noisy_e += rng.normal(0.0, read_noise, size=image.shape)
    return np.clip(noisy_e / gain, 0.0, None)
