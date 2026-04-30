from __future__ import annotations
from pathlib import Path

import numpy as np

from .config_schema import SimulationConfig
from .grids import KGrid, make_gaussian_k_grid

def normalize_source(source: np.ndarray) -> np.ndarray:
    source = np.clip(np.asarray(source, dtype=np.float64), 0.0, None)
    peak = float(np.max(source)) if source.size else 0.0
    if peak <= 0:
        return np.ones_like(source, dtype=np.float64)
    return source / peak

def read_numeric_table(path: str | Path) -> dict[str, np.ndarray]:
    path = Path(path)
    named = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    if named.dtype.names:
        return {name.lower(): np.asarray(named[name], dtype=np.float64) for name in named.dtype.names}
    raw = np.genfromtxt(path, delimiter=",", dtype=float)
    raw = np.atleast_2d(raw)
    return {f"col{i}": raw[:, i].astype(np.float64) for i in range(raw.shape[1])}

def first_column(table: dict[str, np.ndarray], names: tuple[str, ...], fallback: int = 0) -> np.ndarray:
    for name in names:
        if name in table:
            return table[name]
    return table[f"col{fallback}"]

def resample_by_index(values: np.ndarray, n: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if len(values) == n:
        return values.copy()
    src = np.arange(len(values), dtype=np.float64)
    dst = np.linspace(0.0, len(values) - 1.0, n)
    return np.interp(dst, src, values)

def load_wavelength_nm_or_k(path: str | Path, n: int) -> tuple[np.ndarray, np.ndarray]:
    table = read_numeric_table(path)
    if any(name in table for name in ("k_rad_per_um", "k")):
        k = first_column(table, ("k_rad_per_um", "k"))
        k = resample_by_index(k, n)
        wavelength_nm = 2.0 * np.pi / k * 1e3
    elif any(name in table for name in ("wavelength_nm", "lambda_nm", "wavelength", "lambda")):
        wavelength_nm = first_column(table, ("wavelength_nm", "lambda_nm", "wavelength", "lambda"))
        wavelength_nm = resample_by_index(wavelength_nm, n)
        k = 2.0 * np.pi / (wavelength_nm * 1e-3)
    elif "col1" in table:
        wavelength_nm = resample_by_index(table["col1"], n)
        k = 2.0 * np.pi / (wavelength_nm * 1e-3)
    else:
        wavelength_nm = resample_by_index(table["col0"], n)
        k = 2.0 * np.pi / (wavelength_nm * 1e-3)
    if len(k) < 16 or np.any(k <= 0) or np.any(wavelength_nm <= 0):
        raise ValueError("Spectrometer calibration must contain at least 16 positive k/wavelength samples.")
    return wavelength_nm.astype(np.float64), k.astype(np.float64)

def load_source_spectrum(path: str | Path, wavelength_nm: np.ndarray, k: np.ndarray | None = None) -> np.ndarray:
    table = read_numeric_table(path)
    y = first_column(table, ("source", "intensity", "power", "spectrum"), fallback=1 if "col1" in table else 0)

    if any(name in table for name in ("wavelength_nm", "lambda_nm", "wavelength", "lambda")):
        x = first_column(table, ("wavelength_nm", "lambda_nm", "wavelength", "lambda"))
        order = np.argsort(x)
        source = np.interp(wavelength_nm, x[order], y[order], left=y[order][0], right=y[order][-1])
    elif k is not None and any(name in table for name in ("k_rad_per_um", "k")):
        x = first_column(table, ("k_rad_per_um", "k"))
        order = np.argsort(x)
        source = np.interp(k, x[order], y[order], left=y[order][0], right=y[order][-1])
    elif "pixel" in table:
        x = table["pixel"]
        source = np.interp(np.arange(len(wavelength_nm), dtype=np.float64), x, y, left=y[0], right=y[-1])
    elif len(y) == len(wavelength_nm):
        source = y
    else:
        source = resample_by_index(y, len(wavelength_nm))
    return normalize_source(source)

def make_oct_k_grid(config: SimulationConfig, n: int) -> KGrid:
    if config.oct.pixel_to_lambda_file:
        wavelength_nm, k = load_wavelength_nm_or_k(config.oct.pixel_to_lambda_file, n)
        source = (
            load_source_spectrum(config.oct.source_spectrum_file, wavelength_nm, k)
            if config.oct.source_spectrum_file
            else normalize_source(
                np.exp(
                    -0.5
                    * ((wavelength_nm - config.oct.center_wavelength_nm) / max(config.oct.bandwidth_nm / 2.355, 1e-12)) ** 2
                )
            )
        )
        return KGrid(k=k, wavelength_nm=wavelength_nm, source_spectrum=source, pixel_index=np.arange(len(k), dtype=np.int64))

    kgrid = make_gaussian_k_grid(config.oct.center_wavelength_nm, config.oct.bandwidth_nm, n)
    if not config.oct.source_spectrum_file:
        return kgrid
    return KGrid(
        k=kgrid.k,
        wavelength_nm=kgrid.wavelength_nm,
        source_spectrum=load_source_spectrum(config.oct.source_spectrum_file, kgrid.wavelength_nm, kgrid.k),
        pixel_index=kgrid.pixel_index,
    )

def spectral_coordinate(k: np.ndarray) -> np.ndarray:
    span = float(np.max(k) - np.min(k))
    if span <= 0:
        return np.zeros_like(k)
    return 2.0 * (k - np.mean(k)) / span

def dispersion_phase(config: SimulationConfig, k: np.ndarray) -> np.ndarray:
    xi = spectral_coordinate(k)
    return config.errors.dispersion_quadratic_rad * xi**2

def apply_k_linearization_error(config: SimulationConfig, k: np.ndarray) -> np.ndarray:
    rms_px = float(config.errors.k_linearization_rms_pixel)
    if rms_px == 0:
        return k.copy()
    n = len(k)
    t = np.linspace(-1.0, 1.0, n)
    residual = np.sin(np.pi * t) + 0.5 * np.sin(2.0 * np.pi * t + 0.35)
    residual = residual - np.mean(residual)
    residual = residual / max(float(np.sqrt(np.mean(residual**2))), 1e-12)
    dk_px = float(np.mean(np.abs(np.diff(k))))
    return k + rms_px * dk_px * residual
