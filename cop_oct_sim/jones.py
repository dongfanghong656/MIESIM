from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray

@dataclass(frozen=True)
class DichroicPowerBalance:
    s_power_sum: float
    p_power_sum: float
    s_absorption: float
    p_absorption: float

def scalar_dichroic_operator(shape, amplitude=1.0, phase_rad=0.0) -> NDArray[np.complex128]:
    return amplitude * np.exp(1j * phase_rad) * np.ones(shape, dtype=np.complex128)

def apply_dichroic_scalar(pupil: NDArray[np.complex128], amplitude: float = 1.0,
                          phase_rad: float = 0.0) -> NDArray[np.complex128]:
    return pupil * scalar_dichroic_operator(pupil.shape, amplitude, phase_rad)

def jones_dichroic_operator(
    shape,
    s_amplitude: float = 1.0,
    p_amplitude: float = 1.0,
    s_phase_rad: float = 0.0,
    p_phase_rad: float = 0.0,
) -> NDArray[np.complex128]:
    op = np.zeros((2, 2, *shape), dtype=np.complex128)
    op[0, 0] = scalar_dichroic_operator(shape, s_amplitude, s_phase_rad)
    op[1, 1] = scalar_dichroic_operator(shape, p_amplitude, p_phase_rad)
    return op

def apply_jones_to_scalar_pupil(
    pupil: NDArray[np.complex128],
    s_amplitude: float = 1.0,
    p_amplitude: float = 1.0,
    s_phase_rad: float = 0.0,
    p_phase_rad: float = 0.0,
    input_pol: str = "s",
) -> NDArray[np.complex128]:
    if input_pol not in {"s", "p", "mean"}:
        raise ValueError("input_pol must be 's', 'p', or 'mean'.")
    if input_pol == "s":
        return apply_dichroic_scalar(pupil, s_amplitude, s_phase_rad)
    if input_pol == "p":
        return apply_dichroic_scalar(pupil, p_amplitude, p_phase_rad)
    s_coeff = s_amplitude * np.exp(1j * s_phase_rad)
    p_coeff = p_amplitude * np.exp(1j * p_phase_rad)
    intensity_mean = (abs(s_coeff) ** 2 + abs(p_coeff) ** 2) / 2.0
    coherent_hint = s_coeff + p_coeff
    phase = np.angle(coherent_hint) if abs(coherent_hint) > 1e-12 else 0.0
    return pupil * np.sqrt(intensity_mean) * np.exp(1j * phase)

def dichroic_energy_balance(
    r_s_amplitude: float,
    t_s_amplitude: float,
    r_p_amplitude: float | None = None,
    t_p_amplitude: float | None = None,
) -> DichroicPowerBalance:
    r_p_amplitude = r_s_amplitude if r_p_amplitude is None else r_p_amplitude
    t_p_amplitude = t_s_amplitude if t_p_amplitude is None else t_p_amplitude
    s_sum = float(abs(r_s_amplitude) ** 2 + abs(t_s_amplitude) ** 2)
    p_sum = float(abs(r_p_amplitude) ** 2 + abs(t_p_amplitude) ** 2)
    return DichroicPowerBalance(
        s_power_sum=s_sum,
        p_power_sum=p_sum,
        s_absorption=float(max(0.0, 1.0 - s_sum)),
        p_absorption=float(max(0.0, 1.0 - p_sum)),
    )

def jones_table_placeholder(k: float, aoi_deg: float, pol: str = "s") -> complex:
    """Replace with interpolation over Rs/Rp/Ts/Tp and phase table."""
    return 1.0 + 0j
