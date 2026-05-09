from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

def normalize(x: NDArray) -> NDArray:
    m = np.max(np.abs(x))
    return x / m if m > 0 else x

def nrmse(a: NDArray, b: NDArray) -> float:
    a = normalize(a); b = normalize(b)
    return float(np.sqrt(np.mean(np.abs(a - b) ** 2)))

def strehl_ratio(psf_actual: NDArray, psf_ideal: NDArray) -> float:
    """Return the energy-normalized on-axis Strehl ratio.

    The on-axis voxel is the center sample of each PSF volume. Each PSF is
    first normalized by integrated intensity so the metric is insensitive to
    arbitrary complex-amplitude scale, but still detects focus-quality loss at
    the common-path on-axis reference point. Independent peak normalization
    would make the already peak-normalized direct PSF insensitive to defocus.
    """
    actual_i = np.abs(np.asarray(psf_actual)) ** 2
    ideal_i = np.abs(np.asarray(psf_ideal)) ** 2
    if actual_i.shape != ideal_i.shape:
        raise ValueError(f"PSF shapes differ: actual={actual_i.shape}, ideal={ideal_i.shape}")
    actual_energy = float(np.sum(actual_i)) if actual_i.size else 0.0
    ideal_energy = float(np.sum(ideal_i)) if ideal_i.size else 0.0
    if actual_energy <= 0.0 or ideal_energy <= 0.0:
        return float("nan")
    center = tuple(dim // 2 for dim in actual_i.shape)
    actual_center = float(actual_i[center] / actual_energy)
    ideal_center = float(ideal_i[center] / ideal_energy)
    if ideal_center <= 0.0:
        return float("nan")
    return float(actual_center / ideal_center)

def fwhm_1d(profile: NDArray, dx: float = 1.0) -> float:
    p = np.abs(profile)
    if p.max() <= 0:
        return float("nan")
    half = p.max() / 2
    idx = np.where(p >= half)[0]
    if len(idx) == 0:
        return float("nan")
    left_i = int(idx[0])
    right_i = int(idx[-1])

    left = float(left_i)
    if left_i > 0 and p[left_i] != p[left_i - 1]:
        left = left_i - 1 + (half - p[left_i - 1]) / (p[left_i] - p[left_i - 1])

    right = float(right_i)
    if right_i < len(p) - 1 and p[right_i] != p[right_i + 1]:
        right = right_i + (half - p[right_i]) / (p[right_i + 1] - p[right_i])

    return float(max(right - left, 0.0) * dx)

def central_lobe_fwhm_1d(profile: NDArray, dx: float = 1.0, peak_index: int | None = None) -> float:
    p = np.abs(profile)
    if p.max() <= 0:
        return float("nan")
    peak = int(np.argmax(p)) if peak_index is None else int(peak_index)
    if peak < 0 or peak >= len(p) or p[peak] <= 0:
        return float("nan")
    half = p[peak] / 2.0

    left_i = peak
    while left_i > 0 and p[left_i - 1] >= half:
        left_i -= 1
    right_i = peak
    while right_i < len(p) - 1 and p[right_i + 1] >= half:
        right_i += 1

    left = float(left_i)
    if left_i > 0 and p[left_i] != p[left_i - 1]:
        left = left_i - 1 + (half - p[left_i - 1]) / (p[left_i] - p[left_i - 1])

    right = float(right_i)
    if right_i < len(p) - 1 and p[right_i] != p[right_i + 1]:
        right = right_i + (half - p[right_i]) / (p[right_i + 1] - p[right_i])

    return float(max(right - left, 0.0) * dx)

def center_line_fwhm(image: NDArray, dx: float = 1.0, axis: int = 1) -> float:
    arr = np.abs(image)
    cy, cx = np.unravel_index(int(np.argmax(arr)), arr.shape)
    line = arr[cy, :] if axis == 1 else arr[:, cx]
    return central_lobe_fwhm_1d(line, dx)

def axial_fwhm(psf: NDArray, depth_um: NDArray) -> float:
    mag = np.abs(psf)
    profile = np.max(mag, axis=tuple(range(1, mag.ndim)))
    dx = float(np.mean(np.diff(depth_um))) if len(depth_um) > 1 else 1.0
    return fwhm_1d(profile, abs(dx))

def on_axis_axial_fwhm(psf: NDArray, depth_um: NDArray) -> float:
    mag = np.abs(psf)
    if mag.ndim < 3:
        return axial_fwhm(mag, depth_um)
    cy = mag.shape[-2] // 2
    cx = mag.shape[-1] // 2
    profile = mag[:, cy, cx]
    dx = float(np.mean(np.diff(depth_um))) if len(depth_um) > 1 else 1.0
    return central_lobe_fwhm_1d(profile, abs(dx))

def integrated_axial_fwhm(psf: NDArray, depth_um: NDArray) -> float:
    mag = np.abs(psf)
    if mag.ndim < 3:
        return axial_fwhm(mag, depth_um)
    profile = np.sum(mag, axis=tuple(range(1, mag.ndim)))
    dx = float(np.mean(np.diff(depth_um))) if len(depth_um) > 1 else 1.0
    return central_lobe_fwhm_1d(profile, abs(dx))

def lateral_fwhm_from_volume(psf: NDArray, dx_um: float) -> float:
    mag = np.abs(psf)
    if mag.ndim == 3:
        z = int(np.argmax(np.max(mag, axis=(1, 2))))
        plane = mag[z]
    else:
        plane = mag
    return center_line_fwhm(plane, dx_um, axis=1)

def peak_shift_um(a: NDArray, b: NDArray, spacings: tuple[float, ...]) -> float:
    ia = np.array(np.unravel_index(int(np.argmax(np.abs(a))), np.shape(a)))
    ib = np.array(np.unravel_index(int(np.argmax(np.abs(b))), np.shape(b)))
    delta = (ia - ib).astype(float) * np.array(spacings)
    return float(np.sqrt(np.sum(delta**2)))

def sidelobe_ratio_1d(profile: NDArray, guard_px: int = 3) -> float:
    p = np.abs(profile).astype(float)
    if p.max() <= 0:
        return float("nan")
    peak = int(np.argmax(p))
    mask = np.ones_like(p, dtype=bool)
    mask[max(0, peak - guard_px):min(len(p), peak + guard_px + 1)] = False
    if not np.any(mask):
        return 0.0
    return float(np.max(p[mask]) / np.max(p))

def compare_psf_volumes(direct: NDArray, converted: NDArray, depth_um: NDArray, dx_um: float) -> dict:
    direct_mag = normalize(np.abs(direct))
    converted_mag = normalize(np.abs(converted))
    return {
        "nrmse_3d": nrmse(direct_mag, converted_mag),
        "fwhm_ax_direct_um": axial_fwhm(direct_mag, depth_um),
        "fwhm_ax_converted_um": axial_fwhm(converted_mag, depth_um),
        "fwhm_ax_on_axis_direct_um": on_axis_axial_fwhm(direct_mag, depth_um),
        "fwhm_ax_on_axis_converted_um": on_axis_axial_fwhm(converted_mag, depth_um),
        "fwhm_ax_integrated_direct_um": integrated_axial_fwhm(direct_mag, depth_um),
        "fwhm_ax_integrated_converted_um": integrated_axial_fwhm(converted_mag, depth_um),
        "fwhm_lat_direct_um": lateral_fwhm_from_volume(direct_mag, dx_um),
        "fwhm_lat_converted_um": lateral_fwhm_from_volume(converted_mag, dx_um),
        "peak_shift_um": peak_shift_um(direct_mag, converted_mag, (abs(float(np.mean(np.diff(depth_um)))), dx_um, dx_um)),
    }
