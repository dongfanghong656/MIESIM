from __future__ import annotations
from .config_schema import SimulationConfig
from .microscope_forward import simulate_microscope_zstack
from .oct_forward import simulate_oct_psf_direct, simulate_oct_raw_direct
from .reconstruction import reconstruct_sd_oct
from .metrics import compare_psf_volumes, fwhm_1d
from .theory_conversion import axial_profile_from_direct_psf, path_e_separable_psf, predict_axial_gate_from_source

SUPPORTED_DIRECT_PSF_MODELS = {"hybrid_rci", "full_spectral_rci"}

def _direct_psf_model(config: SimulationConfig) -> str:
    model = config.oct.direct_psf_model.lower()
    if model not in SUPPORTED_DIRECT_PSF_MODELS:
        raise ValueError(f"Unsupported oct.direct_psf_model: {config.oct.direct_psf_model}")
    return model

def run_minimal(config: SimulationConfig, N: int = 64, k_samples: int = 256, pad_factor: int = 2):
    mic = simulate_microscope_zstack(config, N=N, pad_factor=pad_factor)
    oct_raw = simulate_oct_raw_direct(config, N=N, k_samples=k_samples)
    oct_complex = reconstruct_sd_oct(oct_raw["interference"], oct_raw["k"], window=config.oct.window)
    direct_model = _direct_psf_model(config)
    oct_direct = simulate_oct_psf_direct(
        config,
        N=N,
        pad_factor=pad_factor,
        k_samples=k_samples,
        full_spectral_rci=direct_model == "full_spectral_rci",
    )
    predictive_gate = predict_axial_gate_from_source(config, k_samples=k_samples, measured_k=oct_raw["k"])
    measured_spatial = mic.get("corrected", mic["measured"])
    measured_predictive = path_e_separable_psf(measured_spatial, predictive_gate["axial_profile"], mic["mode"])
    oracle_axial = axial_profile_from_direct_psf(oct_direct["psf"])
    ideal_upper_bound = path_e_separable_psf(mic["ideal"], oracle_axial, mic["mode"])
    metrics_measured_predictive = compare_psf_volumes(
        oct_direct["psf"],
        measured_predictive,
        oct_direct["depth_um"],
        float(oct_direct["dx_um"]),
    )
    metrics_ideal_upper_bound = compare_psf_volumes(
        oct_direct["psf"],
        ideal_upper_bound,
        oct_direct["depth_um"],
        float(oct_direct["dx_um"]),
    )
    return {
        "mic": mic,
        "oct_raw": oct_raw,
        "oct_complex": oct_complex,
        "oct_direct": oct_direct,
        "direct_psf_model": direct_model,
        "predictive_axial_gate": predictive_gate,
        "path_e_measured_predictive": measured_predictive,
        "path_e_ideal_upper_bound": ideal_upper_bound,
        "path_e_predictive": measured_predictive,
        "path_e_oracle": ideal_upper_bound,
        "path_e_converted": measured_predictive,
        "metrics_measured_predictive": metrics_measured_predictive,
        "metrics_ideal_upper_bound": metrics_ideal_upper_bound,
        "metrics_predictive": metrics_measured_predictive,
        "metrics_oracle": metrics_ideal_upper_bound,
        "metrics": metrics_measured_predictive,
        "oct_axial_fwhm_px": fwhm_1d(abs(oct_complex)),
    }
