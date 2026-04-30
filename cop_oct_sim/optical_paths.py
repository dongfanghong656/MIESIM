from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config_schema import SimulationConfig
from .jones import apply_jones_to_scalar_pupil
from .spectrometer import first_column, read_numeric_table

@dataclass(frozen=True)
class PathDichroicCoefficients:
    path: str
    s_amplitude: float
    p_amplitude: float
    s_phase_rad: float
    p_phase_rad: float

def _nearest_aoi_rows(table: dict[str, np.ndarray], aoi_deg: float) -> np.ndarray:
    if "aoi" not in table:
        return np.ones_like(next(iter(table.values())), dtype=bool)
    aoi = table["aoi"]
    nearest = float(aoi[np.argmin(np.abs(aoi - aoi_deg))])
    return np.isclose(aoi, nearest)

def _interp_table_column(
    table: dict[str, np.ndarray],
    wavelength_nm: float,
    names: tuple[str, ...],
    default: float,
    mask: np.ndarray,
) -> float:
    for name in names:
        if name in table:
            y = table[name][mask]
            if any(col in table for col in ("wavelength_nm", "lambda_nm", "wavelength", "lambda")):
                x = first_column(table, ("wavelength_nm", "lambda_nm", "wavelength", "lambda"))[mask]
                order = np.argsort(x)
                return float(np.interp(wavelength_nm, x[order], y[order], left=y[order][0], right=y[order][-1]))
            return float(y[0])
    return default

def _coefficients_from_table(config: SimulationConfig, path: str, wavelength_nm: float) -> PathDichroicCoefficients:
    table = read_numeric_table(Path(config.dichroic.table_file))
    mask = _nearest_aoi_rows(table, config.dichroic.angle_deg)
    if path == "mic_transmission":
        s_power = _interp_table_column(table, wavelength_nm, ("ts", "t_s", "transmission_s"), 1.0, mask)
        p_power = _interp_table_column(table, wavelength_nm, ("tp", "t_p", "transmission_p"), 1.0, mask)
    elif path == "oct_reflection":
        s_power = _interp_table_column(table, wavelength_nm, ("rs", "r_s", "reflection_s"), 1.0, mask)
        p_power = _interp_table_column(table, wavelength_nm, ("rp", "r_p", "reflection_p"), 1.0, mask)
    else:
        raise ValueError(f"Unknown optical path: {path}")
    s_phase = _interp_table_column(table, wavelength_nm, ("phase_s", "s_phase_rad"), 0.0, mask)
    p_phase = _interp_table_column(table, wavelength_nm, ("phase_p", "p_phase_rad"), 0.0, mask)
    return PathDichroicCoefficients(
        path=path,
        s_amplitude=float(np.sqrt(max(s_power, 0.0))),
        p_amplitude=float(np.sqrt(max(p_power, 0.0))),
        s_phase_rad=s_phase,
        p_phase_rad=p_phase,
    )

def dichroic_path_coefficients(
    config: SimulationConfig,
    path: str,
    wavelength_nm: float,
) -> PathDichroicCoefficients:
    if config.dichroic.table_file:
        return _coefficients_from_table(config, path, wavelength_nm)
    if path == "mic_transmission":
        return PathDichroicCoefficients(
            path=path,
            s_amplitude=config.dichroic.amplitude_s * config.dichroic.mic_transmission_amplitude_s,
            p_amplitude=config.dichroic.amplitude_p * config.dichroic.mic_transmission_amplitude_p,
            s_phase_rad=config.dichroic.phase_s_rad + config.dichroic.mic_transmission_phase_s_rad,
            p_phase_rad=config.dichroic.phase_p_rad + config.dichroic.mic_transmission_phase_p_rad,
        )
    if path == "oct_reflection":
        return PathDichroicCoefficients(
            path=path,
            s_amplitude=config.dichroic.amplitude_s * config.dichroic.oct_reflection_amplitude_s,
            p_amplitude=config.dichroic.amplitude_p * config.dichroic.oct_reflection_amplitude_p,
            s_phase_rad=config.dichroic.phase_s_rad + config.dichroic.oct_reflection_phase_s_rad,
            p_phase_rad=config.dichroic.phase_p_rad + config.dichroic.oct_reflection_phase_p_rad,
        )
    raise ValueError(f"Unknown optical path: {path}")

def apply_path_dichroic_to_pupil(
    pupil: np.ndarray,
    config: SimulationConfig,
    path: str,
    wavelength_nm: float,
    input_pol: str = "mean",
) -> tuple[np.ndarray, PathDichroicCoefficients]:
    coeffs = dichroic_path_coefficients(config, path, wavelength_nm)
    return (
        apply_jones_to_scalar_pupil(
            pupil,
            s_amplitude=coeffs.s_amplitude,
            p_amplitude=coeffs.p_amplitude,
            s_phase_rad=coeffs.s_phase_rad,
            p_phase_rad=coeffs.p_phase_rad,
            input_pol=input_pol,
        ),
        coeffs,
    )
