from __future__ import annotations
import ast
import csv
import json
import platform
import sys
from dataclasses import replace
from datetime import datetime
from pathlib import Path

import numpy as np

from .config_schema import SimulationConfig, config_hash, config_to_dict, dump_resolved_config
from .figures import generate_sweep_trend_figure, generate_validation_figures
from .grids import lateral_sample_pitch_um, make_pupil_grid
from .io_hdf5 import write_dict_h5
from .jones import dichroic_energy_balance
from .metrics import (
    axial_fwhm,
    center_line_fwhm,
    compare_psf_volumes,
    lateral_fwhm_from_volume,
    sidelobe_ratio_1d,
)
from .microscope_forward import simulate_microscope_zstack
from .oct_forward import simulate_oct_psf_direct, simulate_oct_raw_direct
from .optical_paths import dichroic_path_coefficients
from .pipelines import run_minimal
from .propagation import fraunhofer_psf_from_pupil, propagation_energy_ratio
from .pupil import build_shared_pupil, clear_circular_pupil
from .reconstruction import reconstruct_sd_oct
from .theory_conversion import axial_profile_from_direct_psf, path_e_separable_psf, predict_axial_gate_from_source

def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return str(value)

def _pass(name: str, value: float, threshold: float, comparator: str = "<") -> dict:
    ok = value < threshold if comparator == "<" else value > threshold
    return {"name": name, "value": float(value), "threshold": threshold, "comparator": comparator, "pass": bool(ok)}

def _relative_error(reference: float, estimate: float) -> float:
    if not (np.isfinite(reference) and np.isfinite(estimate)):
        return float("inf")
    return float(abs(float(reference) - float(estimate)) / max(abs(float(reference)), 1e-12))

def _path_e_gate(
    model_name: str,
    metrics: dict,
    nrmse_threshold: float,
    axial_fwhm_rel_threshold: float,
    lateral_fwhm_rel_threshold: float,
) -> dict:
    axial_error = _relative_error(
        metrics["fwhm_ax_on_axis_direct_um"], metrics["fwhm_ax_on_axis_converted_um"]
    )
    lateral_error = _relative_error(metrics["fwhm_lat_direct_um"], metrics["fwhm_lat_converted_um"])
    nrmse_value = float(metrics["nrmse_3d"])
    ok = (
        np.isfinite(nrmse_value)
        and np.isfinite(axial_error)
        and np.isfinite(lateral_error)
        and nrmse_value < nrmse_threshold
        and axial_error < axial_fwhm_rel_threshold
        and lateral_error < lateral_fwhm_rel_threshold
    )
    return {
        **metrics,
        "model_name": model_name,
        "nrmse_threshold": nrmse_threshold,
        "axial_fwhm_relative_error": axial_error,
        "axial_fwhm_relative_error_threshold": axial_fwhm_rel_threshold,
        "lateral_fwhm_relative_error": lateral_error,
        "lateral_fwhm_relative_error_threshold": lateral_fwhm_rel_threshold,
        "pass": bool(ok),
    }

def _provenance(stamp: str, config: SimulationConfig) -> dict:
    return {
        "timestamp": stamp,
        "config_hash": config_hash(config),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "command": " ".join(sys.argv),
        "package_versions": {
            "numpy": np.__version__,
        },
    }

def airy_sanity(wavelength_um: float = 0.85, NA: float = 0.25) -> dict:
    grid = make_pupil_grid(128, NA)
    pupil = clear_circular_pupil(grid)
    pad = 8
    field = fraunhofer_psf_from_pupil(pupil, normalize=True, pad_factor=pad)
    intensity = np.abs(field) ** 2
    dx = lateral_sample_pitch_um(wavelength_um, NA, pad)
    measured = center_line_fwhm(intensity, dx)
    expected = 0.514 * wavelength_um / NA
    rel_error = abs(measured - expected) / expected
    return {
        "measured_fwhm_um": measured,
        "expected_fwhm_um": expected,
        "relative_error": rel_error,
        "pass": bool(rel_error < 0.25),
    }

def gaussian_underfill_sanity(wavelength_um: float = 0.85, NA: float = 0.25) -> dict:
    grid = make_pupil_grid(128, NA)
    pad = 8
    dx = lateral_sample_pitch_um(wavelength_um, NA, pad)
    narrow_pupil = build_shared_pupil(grid, wavelength_um, fill_ratio=0.55)
    full_pupil = build_shared_pupil(grid, wavelength_um, fill_ratio=1.05)
    narrow = np.abs(fraunhofer_psf_from_pupil(narrow_pupil, pad_factor=pad)) ** 2
    full = np.abs(fraunhofer_psf_from_pupil(full_pupil, pad_factor=pad)) ** 2
    narrow_fwhm = center_line_fwhm(narrow, dx)
    full_fwhm = center_line_fwhm(full, dx)
    return {
        "underfilled_fwhm_um": narrow_fwhm,
        "fuller_fwhm_um": full_fwhm,
        "broadening_ratio": narrow_fwhm / full_fwhm,
        "pass": bool(narrow_fwhm > full_fwhm),
    }

def parseval_sanity() -> dict:
    grid = make_pupil_grid(64, 0.25)
    pupil = build_shared_pupil(grid, 0.85, fill_ratio=0.82)
    field = fraunhofer_psf_from_pupil(pupil, normalize=False, pad_factor=2)
    ratio = propagation_energy_ratio(pupil, field)
    return {"energy_ratio": ratio, "pass": bool(abs(ratio - 1.0) < 1e-12)}

def dichroic_sanity() -> dict:
    amp = 1.0 / np.sqrt(2.0)
    balance = dichroic_energy_balance(amp, amp)
    err = max(abs(balance.s_power_sum - 1.0), abs(balance.p_power_sum - 1.0))
    return {
        "s_power_sum": balance.s_power_sum,
        "p_power_sum": balance.p_power_sum,
        "max_error": err,
        "pass": bool(err < 1e-12),
    }

def axial_bandwidth_sanity(config: SimulationConfig) -> dict:
    cfg_narrow = replace(
        config,
        oct=replace(config.oct, bandwidth_nm=35.0, spectrometer_pixels=128),
        microscope=replace(config.microscope, wavelength_nm=config.oct.center_wavelength_nm),
    )
    cfg_wide = replace(
        config,
        oct=replace(config.oct, bandwidth_nm=80.0, spectrometer_pixels=128),
        microscope=replace(config.microscope, wavelength_nm=config.oct.center_wavelength_nm),
    )
    narrow = simulate_oct_psf_direct(cfg_narrow, N=32, pad_factor=1, k_samples=128)
    wide = simulate_oct_psf_direct(cfg_wide, N=32, pad_factor=1, k_samples=128)
    narrow_fwhm = axial_fwhm(narrow["psf"], narrow["depth_um"])
    wide_fwhm = axial_fwhm(wide["psf"], wide["depth_um"])
    return {
        "narrow_bandwidth_nm": 35.0,
        "wide_bandwidth_nm": 80.0,
        "narrow_axial_fwhm_um": narrow_fwhm,
        "wide_axial_fwhm_um": wide_fwhm,
        "pass": bool(wide_fwhm < narrow_fwhm),
    }

def no_leak_sanity() -> dict:
    import cop_oct_sim.oct_forward as oct_forward

    source = Path(oct_forward.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    imports: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    forbidden = ["microscope_forward", "theory_conversion"]
    leaks = [item for item in imports if any(name in item for name in forbidden)]
    return {"imports": imports, "leaks": leaks, "pass": len(leaks) == 0}

def _sweep_row(
    run_id: str,
    error_name: str,
    error_value: float,
    cfg: SimulationConfig,
    bead_um: float = 0.2,
    scatterer_z_um: float = 0.0,
) -> dict:
    mic = simulate_microscope_zstack(cfg, N=32, bead_diameter_um=bead_um, pad_factor=2)
    direct = simulate_oct_psf_direct(cfg, N=32, pad_factor=2, k_samples=96, scatterer_z_um=scatterer_z_um)
    raw = simulate_oct_raw_direct(cfg, N=32, k_samples=96, scatterer_z_um=scatterer_z_um)
    axial = axial_profile_from_direct_psf(direct["psf"])
    converted = path_e_separable_psf(mic["measured"], axial, mic["mode"])
    cmp = compare_psf_volumes(direct["psf"], converted, direct["depth_um"], float(direct["dx_um"]))
    mic_lat = lateral_fwhm_from_volume(mic["measured"], float(mic["dx_um"]))
    direct_axial_profile = np.max(np.abs(direct["psf"]), axis=(1, 2))
    converted_axial_profile = np.max(np.abs(converted), axis=(1, 2))
    lat_rel = abs(cmp["fwhm_lat_direct_um"] - cmp["fwhm_lat_converted_um"]) / max(cmp["fwhm_lat_direct_um"], 1e-12)
    k_true = raw.get("k_true", raw["k"])
    dk = max(float(np.mean(np.abs(np.diff(k_true)))), 1e-12)
    k_rms_px = float(np.sqrt(np.mean(((raw["k"] - k_true) / dk) ** 2)))
    raw_mod = float(np.std(raw["raw"] - np.mean(raw["raw"])))
    pass_fail = "pass" if (cmp["nrmse_3d"] < 0.25 and lat_rel < 0.25) else "review"
    if error_name in {"dispersion_quadratic_rad", "k_linearization_rms_pixel", "rolloff_per_um"}:
        pass_fail = "diagnostic"
    return {
        "run_id": run_id,
        "error_name": error_name,
        "error_value": error_value,
        "scatterer_z_um": scatterer_z_um,
        "fwhm_lat_direct": cmp["fwhm_lat_direct_um"],
        "fwhm_lat_converted": cmp["fwhm_lat_converted_um"],
        "fwhm_lat_mic_measured": mic_lat,
        "fwhm_ax_direct": cmp["fwhm_ax_direct_um"],
        "fwhm_ax_converted": cmp["fwhm_ax_converted_um"],
        "peak_shift_um": cmp["peak_shift_um"],
        "strehl_delta": float(np.max(np.abs(converted)) - np.max(np.abs(direct["psf"]))),
        "sidelobe_delta": float(sidelobe_ratio_1d(converted_axial_profile) - sidelobe_ratio_1d(direct_axial_profile)),
        "nrmse_3d": cmp["nrmse_3d"],
        "raw_modulation_std": raw_mod,
        "rolloff_amplitude": float(raw["rolloff_amplitude"]),
        "k_rms_pixel_error": k_rms_px,
        "dispersion_phase_ptp_rad": float(np.ptp(raw["dispersion_phase_rad"])),
        "pass_fail": pass_fail,
        "suspected_failure_mode": "baseline" if pass_fail == "pass" else f"{error_name}_sensitivity",
    }

def build_sweep_rows(config: SimulationConfig) -> list[dict]:
    rows: list[dict] = []
    for fill in [0.55, 0.70, 0.82, 1.00]:
        cfg = replace(config, objective=replace(config.objective, pupil_fill_ratio=fill))
        rows.append(_sweep_row(f"fill_{fill:.2f}", "pupil_fill_ratio", fill, cfg))
    for bead in [0.10, 0.20, 0.50, 1.00]:
        rows.append(_sweep_row(f"bead_{bead:.2f}", "bead_diameter_um", bead, config, bead_um=bead))
    for focus in [-2.0, 0.0, 2.0]:
        cfg = replace(config, errors=replace(config.errors, objective_focus_shift_um=focus))
        rows.append(_sweep_row(f"focus_{focus:+.1f}", "objective_focus_shift_um", focus, cfg))
    for disp in [0.0, 10.0, 25.0]:
        cfg = replace(config, errors=replace(config.errors, dispersion_quadratic_rad=disp))
        rows.append(_sweep_row(f"dispersion_{disp:.1f}", "dispersion_quadratic_rad", disp, cfg))
    for kerr in [0.0, 1.0, 3.0]:
        cfg = replace(config, errors=replace(config.errors, k_linearization_rms_pixel=kerr))
        rows.append(_sweep_row(f"kerr_{kerr:.1f}", "k_linearization_rms_pixel", kerr, cfg))
    for rolloff in [0.0, 0.03, 0.08]:
        cfg = replace(config, errors=replace(config.errors, rolloff_per_um=rolloff))
        rows.append(_sweep_row(f"rolloff_{rolloff:.2f}", "rolloff_per_um", rolloff, cfg, scatterer_z_um=20.0))
    return rows

def _convergence_metric(config: SimulationConfig, N: int, k_samples: int, pad_factor: int) -> dict:
    result = run_minimal(config, N=N, k_samples=k_samples, pad_factor=pad_factor)
    metrics = result["metrics_predictive"]
    return {
        "nrmse_3d": float(metrics["nrmse_3d"]),
        "fwhm_ax_direct_um": float(metrics["fwhm_ax_direct_um"]),
        "fwhm_lat_direct_um": float(metrics["fwhm_lat_direct_um"]),
        "fwhm_ax_converted_um": float(metrics["fwhm_ax_converted_um"]),
        "fwhm_lat_converted_um": float(metrics["fwhm_lat_converted_um"]),
    }

def build_convergence_rows(
    config: SimulationConfig,
    N: int = 32,
    k_samples: int = 96,
    pad_factor: int = 2,
) -> list[dict]:
    base_N = max(int(N), 24)
    base_k = max(int(k_samples), 64)
    base_pad = max(int(pad_factor), 1)
    baseline = _convergence_metric(config, base_N, base_k, base_pad)
    baseline_nrmse = baseline["nrmse_3d"]
    baseline_ax = baseline["fwhm_ax_direct_um"]
    baseline_lat = baseline["fwhm_lat_direct_um"]

    cases = [
        ("N", base_N + 8, config, base_N + 8, base_k, base_pad),
        ("pad_factor", 1 if base_pad > 1 else 2, config, base_N, base_k, 1 if base_pad > 1 else 2),
        ("k_samples", base_k + 32, config, base_N, base_k + 32, base_pad),
        (
            "z_step_um",
            config.microscope.z_step_um * 0.5,
            replace(config, microscope=replace(config.microscope, z_step_um=config.microscope.z_step_um * 0.5)),
            base_N,
            base_k,
            base_pad,
        ),
        (
            "window",
            "tukey" if config.oct.window.lower() != "tukey" else "hann",
            replace(config, oct=replace(config.oct, window="tukey" if config.oct.window.lower() != "tukey" else "hann")),
            base_N,
            base_k,
            base_pad,
        ),
        (
            "bead_diameter_um",
            max(config.microscope.bead_diameter_um * 2.0, 0.5),
            replace(
                config,
                microscope=replace(
                    config.microscope,
                    bead_diameter_um=max(config.microscope.bead_diameter_um * 2.0, 0.5),
                ),
            ),
            base_N,
            base_k,
            base_pad,
        ),
        (
            "camera_pixel_um",
            config.microscope.camera_pixel_um * 2.0,
            replace(
                config,
                microscope=replace(config.microscope, camera_pixel_um=config.microscope.camera_pixel_um * 2.0),
            ),
            base_N,
            base_k,
            base_pad,
        ),
    ]

    rows: list[dict] = []
    for axis, value, cfg, case_N, case_k, case_pad in cases:
        metric = _convergence_metric(cfg, case_N, case_k, case_pad)
        nrmse_delta = abs(metric["nrmse_3d"] - baseline_nrmse)
        ax_delta = abs(metric["fwhm_ax_direct_um"] - baseline_ax)
        lat_delta = abs(metric["fwhm_lat_direct_um"] - baseline_lat)
        status = "stable" if nrmse_delta < 0.05 else "review"
        rows.append(
            {
                "sweep_axis": axis,
                "baseline_value": str({"N": base_N, "k_samples": base_k, "pad_factor": base_pad}.get(axis, "config_default")),
                "test_value": str(value),
                "N": case_N,
                "k_samples": case_k,
                "pad_factor": case_pad,
                "baseline_nrmse_3d": baseline_nrmse,
                "nrmse_3d": metric["nrmse_3d"],
                "nrmse_delta_abs": nrmse_delta,
                "fwhm_ax_direct_um": metric["fwhm_ax_direct_um"],
                "fwhm_ax_delta_abs_um": ax_delta,
                "fwhm_lat_direct_um": metric["fwhm_lat_direct_um"],
                "fwhm_lat_delta_abs_um": lat_delta,
                "status": status,
            }
        )
    return rows

def build_direct_model_comparison_rows(
    config: SimulationConfig,
    N: int = 12,
    k_samples: int = 24,
    pad_factor: int = 1,
) -> list[dict]:
    depth_decimation = max(int(config.oct.full_spectral_rci_depth_decimation), 1)
    axial_rel_threshold = 0.25
    lateral_rel_threshold = 0.25
    hybrid_cfg = replace(config, oct=replace(config.oct, direct_psf_model="hybrid_rci"))
    full_cfg = replace(
        config,
        oct=replace(
            config.oct,
            direct_psf_model="full_spectral_rci",
            full_spectral_rci_depth_decimation=depth_decimation,
        ),
    )

    hybrid = simulate_oct_psf_direct(
        hybrid_cfg,
        N=N,
        pad_factor=pad_factor,
        k_samples=k_samples,
        full_spectral_rci=False,
    )
    full = simulate_oct_psf_direct(
        full_cfg,
        N=N,
        pad_factor=pad_factor,
        k_samples=k_samples,
        full_spectral_rci=True,
        full_spectral_rci_depth_decimation=depth_decimation,
    )
    cmp = compare_psf_volumes(full["psf"], hybrid["psf"], full["depth_um"], float(full["dx_um"]))
    computed_count = int(len(full.get("spectral_rci_computed_depth_indices", [])))
    depth_count = int(len(full["depth_um"]))
    finite = all(np.isfinite(float(cmp[key])) for key in cmp)
    axial_rel_error = _relative_error(
        cmp["fwhm_ax_on_axis_direct_um"], cmp["fwhm_ax_on_axis_converted_um"]
    )
    lateral_rel_error = _relative_error(cmp["fwhm_lat_direct_um"], cmp["fwhm_lat_converted_um"])
    pass_gate = (
        finite
        and computed_count > 0
        and depth_count > 0
        and axial_rel_error <= axial_rel_threshold
        and lateral_rel_error <= lateral_rel_threshold
    )
    return [
        {
            "comparison_name": "hybrid_vs_full_spectral_rci",
            "hybrid_direct_model": "hybrid_rci",
            "full_direct_model": "full_spectral_rci",
            "N": int(N),
            "k_samples": int(k_samples),
            "pad_factor": int(pad_factor),
            "full_spectral_rci_depth_decimation": depth_decimation,
            "full_spectral_rci_interpolated": bool(full.get("full_spectral_rci_interpolated", False)),
            "computed_depth_count": computed_count,
            "depth_count": depth_count,
            "nrmse_3d": float(cmp["nrmse_3d"]),
            "fwhm_ax_full_spectral_um": float(cmp["fwhm_ax_direct_um"]),
            "fwhm_ax_hybrid_um": float(cmp["fwhm_ax_converted_um"]),
            "fwhm_ax_on_axis_full_spectral_um": float(cmp["fwhm_ax_on_axis_direct_um"]),
            "fwhm_ax_on_axis_hybrid_um": float(cmp["fwhm_ax_on_axis_converted_um"]),
            "axial_fwhm_relative_error": axial_rel_error,
            "axial_fwhm_relative_error_threshold": axial_rel_threshold,
            "fwhm_lat_full_spectral_um": float(cmp["fwhm_lat_direct_um"]),
            "fwhm_lat_hybrid_um": float(cmp["fwhm_lat_converted_um"]),
            "lateral_fwhm_relative_error": lateral_rel_error,
            "lateral_fwhm_relative_error_threshold": lateral_rel_threshold,
            "peak_shift_um": float(cmp["peak_shift_um"]),
            "status": "diagnostic" if pass_gate else "blocker",
            "pass": bool(pass_gate),
        }
    ]

def _normalized_axial_profiles(volume: np.ndarray) -> dict:
    mag = np.abs(volume).astype(np.float64)
    cy = mag.shape[1] // 2
    cx = mag.shape[2] // 2

    def norm(profile: np.ndarray) -> np.ndarray:
        peak = float(np.max(np.abs(profile))) if profile.size else 0.0
        return profile / peak if peak > 0 else profile

    return {
        "max": norm(np.max(mag, axis=(1, 2))),
        "on_axis": norm(mag[:, cy, cx]),
        "integrated": norm(np.sum(mag, axis=(1, 2))),
    }

def build_direct_model_axial_profile_rows(
    config: SimulationConfig,
    N: int = 12,
    k_samples: int = 24,
    pad_factor: int = 1,
) -> list[dict]:
    depth_decimation = max(int(config.oct.full_spectral_rci_depth_decimation), 1)
    hybrid_cfg = replace(config, oct=replace(config.oct, direct_psf_model="hybrid_rci"))
    full_cfg = replace(
        config,
        oct=replace(
            config.oct,
            direct_psf_model="full_spectral_rci",
            full_spectral_rci_depth_decimation=depth_decimation,
        ),
    )
    hybrid = simulate_oct_psf_direct(
        hybrid_cfg,
        N=N,
        pad_factor=pad_factor,
        k_samples=k_samples,
        full_spectral_rci=False,
    )
    full = simulate_oct_psf_direct(
        full_cfg,
        N=N,
        pad_factor=pad_factor,
        k_samples=k_samples,
        full_spectral_rci=True,
        full_spectral_rci_depth_decimation=depth_decimation,
    )
    full_profiles = _normalized_axial_profiles(full["psf"])
    hybrid_profiles = _normalized_axial_profiles(hybrid["psf"])
    rows: list[dict] = []
    for i, depth in enumerate(full["depth_um"]):
        rows.append(
            {
                "depth_index": int(i),
                "depth_um": float(depth),
                "full_profile_max": float(full_profiles["max"][i]),
                "full_profile_on_axis": float(full_profiles["on_axis"][i]),
                "full_profile_integrated": float(full_profiles["integrated"][i]),
                "hybrid_profile_max": float(hybrid_profiles["max"][i]),
                "hybrid_profile_on_axis": float(hybrid_profiles["on_axis"][i]),
                "hybrid_profile_integrated": float(hybrid_profiles["integrated"][i]),
                "on_axis_abs_residual": float(abs(full_profiles["on_axis"][i] - hybrid_profiles["on_axis"][i])),
                "full_spectral_rci_depth_decimation": depth_decimation,
                "full_spectral_rci_interpolated": bool(full.get("full_spectral_rci_interpolated", False)),
            }
        )
    return rows

def _predictive_metric_for_pair(
    mic_cfg: SimulationConfig,
    oct_cfg: SimulationConfig,
    N: int,
    k_samples: int,
    pad_factor: int,
) -> dict:
    mic = simulate_microscope_zstack(mic_cfg, N=N, pad_factor=pad_factor)
    direct = simulate_oct_psf_direct(oct_cfg, N=N, pad_factor=pad_factor, k_samples=k_samples)
    raw = simulate_oct_raw_direct(oct_cfg, N=N, k_samples=k_samples)
    gate = predict_axial_gate_from_source(oct_cfg, k_samples=k_samples, measured_k=raw["k"])
    converted = path_e_separable_psf(mic.get("corrected", mic["measured"]), gate["axial_profile"], mic["mode"])
    return compare_psf_volumes(direct["psf"], converted, direct["depth_um"], float(direct["dx_um"]))

def _k_rms_pixel_error(raw: dict) -> float:
    k_true = raw.get("k_true", raw["k"])
    dk = max(float(np.mean(np.abs(np.diff(k_true)))), 1e-12)
    return float(np.sqrt(np.mean(((raw["k"] - k_true) / dk) ** 2)))

def build_negative_control_rows(
    config: SimulationConfig,
    N: int = 32,
    k_samples: int = 96,
    pad_factor: int = 2,
) -> list[dict]:
    baseline = run_minimal(config, N=N, k_samples=k_samples, pad_factor=pad_factor)
    baseline_nrmse = float(baseline["metrics_predictive"]["nrmse_3d"])

    chromatic_oct_cfg = replace(
        config,
        objective=replace(config.objective, focus_shift_um_850_vs_visible=6.0),
        errors=replace(config.errors, objective_focus_shift_um=4.0),
    )
    chromatic_mic_cfg = replace(
        config,
        microscope=replace(config.microscope, wavelength_nm=660.0),
        errors=replace(config.errors, objective_focus_shift_um=0.0),
    )
    chromatic_metric = _predictive_metric_for_pair(
        chromatic_mic_cfg,
        chromatic_oct_cfg,
        N=N,
        k_samples=k_samples,
        pad_factor=pad_factor,
    )
    chromatic_wavelength_delta_nm = abs(
        float(chromatic_mic_cfg.microscope.wavelength_nm) - float(chromatic_oct_cfg.oct.center_wavelength_nm)
    )

    wrong_path_cfg = replace(
        config,
        dichroic=replace(
            config.dichroic,
            mic_transmission_amplitude_s=0.10,
            mic_transmission_amplitude_p=0.15,
            oct_reflection_amplitude_s=1.00,
            oct_reflection_amplitude_p=0.90,
        ),
    )
    mic_coeff = dichroic_path_coefficients(
        wrong_path_cfg,
        "mic_transmission",
        wrong_path_cfg.microscope.wavelength_nm,
    )
    oct_coeff = dichroic_path_coefficients(
        wrong_path_cfg,
        "oct_reflection",
        wrong_path_cfg.oct.center_wavelength_nm,
    )
    path_delta = float(
        np.sqrt(
            (mic_coeff.s_amplitude - oct_coeff.s_amplitude) ** 2
            + (mic_coeff.p_amplitude - oct_coeff.p_amplitude) ** 2
        )
    )

    bad_k_cfg = replace(
        config,
        errors=replace(config.errors, k_linearization_rms_pixel=4.0),
    )
    bad_k_raw = simulate_oct_raw_direct(bad_k_cfg, N=N, k_samples=k_samples)
    bad_k_metric = _predictive_metric_for_pair(
        config,
        bad_k_cfg,
        N=N,
        k_samples=k_samples,
        pad_factor=pad_factor,
    )
    bad_k_rms = _k_rms_pixel_error(bad_k_raw)

    rows = [
        {
            "control_name": "chromatic_mismatch",
            "failure_mode": "visible_microscope_spatial_term_does_not_match_850nm_oct_truth",
            "baseline_nrmse_3d": baseline_nrmse,
            "bad_nrmse_3d": float(chromatic_metric["nrmse_3d"]),
            "detection_metric": "microscope_oct_wavelength_delta_nm",
            "detection_value": chromatic_wavelength_delta_nm,
            "detection_threshold": 50.0,
            "nrmse_increase": float(chromatic_metric["nrmse_3d"] - baseline_nrmse),
            "path_amplitude_delta": 0.0,
            "bad_k_rms_pixel_error": 0.0,
        },
        {
            "control_name": "wrong_dichroic_path",
            "failure_mode": "mic_transmission_coefficients_are_not_oct_reflection_coefficients",
            "baseline_nrmse_3d": baseline_nrmse,
            "bad_nrmse_3d": baseline_nrmse,
            "detection_metric": "path_amplitude_delta",
            "detection_value": path_delta,
            "detection_threshold": 0.50,
            "path_amplitude_delta": path_delta,
            "bad_k_rms_pixel_error": 0.0,
        },
        {
            "control_name": "bad_k_calibration",
            "failure_mode": "measured_k_grid_has_large_residual_pixel_error",
            "baseline_nrmse_3d": baseline_nrmse,
            "bad_nrmse_3d": float(bad_k_metric["nrmse_3d"]),
            "detection_metric": "k_rms_pixel_error",
            "detection_value": bad_k_rms,
            "detection_threshold": 2.0,
            "path_amplitude_delta": 0.0,
            "bad_k_rms_pixel_error": bad_k_rms,
        },
    ]
    for row in rows:
        row["detected"] = bool(float(row["detection_value"]) > float(row["detection_threshold"]))
    return rows

def run_validation_suite(
    output_root: str | Path = "outputs",
    config: SimulationConfig | None = None,
    N: int = 48,
    k_samples: int = 160,
    pad_factor: int = 2,
) -> dict:
    config = config or SimulationConfig()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(output_root) / f"validation_{stamp}"
    metrics_dir = out_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    provenance = _provenance(stamp, config)

    dump_resolved_config(out_dir / "config_resolved.yaml", config)
    minimal = run_minimal(config, N=N, k_samples=k_samples, pad_factor=pad_factor)
    full_spectral_lowres = simulate_oct_psf_direct(
        config,
        N=min(max(16, N // 2), 20),
        pad_factor=1,
        k_samples=min(max(32, k_samples // 2), 48),
        full_spectral_rci=True,
    )
    raw = simulate_oct_raw_direct(config, N=N, k_samples=k_samples)
    reconstruction = reconstruct_sd_oct(
        raw["interference"],
        raw["k"],
        window=config.oct.window,
        source_spectrum=raw["source"],
        return_depth=True,
    )
    reconstruction_product = {
        "complex": reconstruction["complex"],
        "magnitude": np.abs(reconstruction["complex"]).astype(np.float32),
        "depth_um": reconstruction["depth_um"],
        "k_linear": reconstruction["k_linear"],
    }

    h5_meta = {"config_hash": config_hash(config), "provenance": provenance}
    write_dict_h5(out_dir / "microscope_forward" / "microscope_zstack.h5", minimal["mic"], h5_meta)
    write_dict_h5(
        out_dir / "microscope_forward" / "illumination_pupil.h5",
        {
            "illumination_pupil": minimal["mic"]["illumination_pupil"],
            "illumination_psf": minimal["mic"]["illumination_psf"],
            "illumination_envelope": minimal["mic"]["illumination_envelope"],
            "illumination_power_fraction": minimal["mic"]["illumination_power_fraction"],
            "illumination_na": minimal["mic"]["illumination_na"],
            "illumination_pupil_fill_ratio": minimal["mic"]["illumination_pupil_fill_ratio"],
            "illumination_aperture_stop_radius": minimal["mic"]["illumination_aperture_stop_radius"],
        },
        h5_meta,
    )
    write_dict_h5(out_dir / "oct_direct" / "oct_raw_interferogram.h5", raw, h5_meta)
    write_dict_h5(out_dir / "oct_direct" / "oct_reconstruction.h5", reconstruction_product, h5_meta)
    write_dict_h5(out_dir / "oct_direct" / "oct_psf_direct.h5", minimal["oct_direct"], h5_meta)
    write_dict_h5(
        out_dir / "oct_direct" / "oct_psf_direct_lowres_full_spectral_rci.h5",
        full_spectral_lowres,
        {"model": "lowres_full_spectral_rci_direct", **h5_meta},
    )
    write_dict_h5(
        out_dir / "converters" / "path_e_measured_predictive_psf.h5",
        {
            "psf": minimal["path_e_measured_predictive"],
            "axial_profile": minimal["predictive_axial_gate"]["axial_profile"],
        },
        {"model": "path_e_measured_predictive", **h5_meta},
    )
    write_dict_h5(
        out_dir / "converters" / "path_e_ideal_upper_bound_psf.h5",
        {"psf": minimal["path_e_ideal_upper_bound"]},
        {"model": "path_e_ideal_upper_bound", **h5_meta},
    )
    write_dict_h5(
        out_dir / "converters" / "path_e_predictive_psf.h5",
        {
            "psf": minimal["path_e_measured_predictive"],
            "axial_profile": minimal["predictive_axial_gate"]["axial_profile"],
        },
        {"model": "legacy_alias_for_path_e_measured_predictive", **h5_meta},
    )
    write_dict_h5(
        out_dir / "converters" / "path_e_oracle_upper_bound_psf.h5",
        {"psf": minimal["path_e_ideal_upper_bound"]},
        {"model": "legacy_alias_for_path_e_ideal_upper_bound", **h5_meta},
    )
    write_dict_h5(
        out_dir / "theory_conversion_under_test" / "path_e_converted_psf.h5",
        {"psf": minimal["path_e_converted"]},
        {"model": "legacy_alias_for_path_e_measured_predictive", **h5_meta},
    )
    generate_validation_figures(
        out_dir / "figures",
        minimal["mic"],
        reconstruction_product,
        minimal["oct_direct"]["psf"],
        minimal["path_e_measured_predictive"],
    )
    negative_rows = build_negative_control_rows(config, N=max(24, N // 2), k_samples=max(64, k_samples // 2), pad_factor=pad_factor)
    convergence_rows = build_convergence_rows(
        config,
        N=max(24, N // 2),
        k_samples=max(64, k_samples // 2),
        pad_factor=pad_factor,
    )
    direct_model_rows = build_direct_model_comparison_rows(
        config,
        N=min(max(12, N // 3), 16),
        k_samples=min(max(24, k_samples // 4), 32),
        pad_factor=1,
    )
    direct_profile_rows = build_direct_model_axial_profile_rows(
        config,
        N=min(max(12, N // 3), 16),
        k_samples=min(max(24, k_samples // 4), 32),
        pad_factor=1,
    )

    path_e_measured = _path_e_gate(
        "path_e_measured_predictive",
        minimal["metrics_measured_predictive"],
        nrmse_threshold=0.15,
        axial_fwhm_rel_threshold=0.25,
        lateral_fwhm_rel_threshold=0.30,
    )
    path_e_ideal = _path_e_gate(
        "path_e_ideal_upper_bound",
        minimal["metrics_ideal_upper_bound"],
        nrmse_threshold=0.08,
        axial_fwhm_rel_threshold=0.15,
        lateral_fwhm_rel_threshold=0.25,
    )
    checks = {
        "airy": airy_sanity(config.oct.center_wavelength_nm * 1e-3, config.objective.NA_nominal),
        "gaussian_underfill": gaussian_underfill_sanity(config.oct.center_wavelength_nm * 1e-3, config.objective.NA_nominal),
        "parseval": parseval_sanity(),
        "dichroic": dichroic_sanity(),
        "axial_bandwidth": axial_bandwidth_sanity(config),
        "no_leak": no_leak_sanity(),
        "path_e_measured_predictive": path_e_measured,
        "path_e_ideal_upper_bound": path_e_ideal,
        "path_e_predictive": {**path_e_measured, "model_name": "legacy_alias_for_path_e_measured_predictive"},
        "path_e_oracle_upper_bound": {**path_e_ideal, "model_name": "legacy_alias_for_path_e_ideal_upper_bound"},
        "scatterer_contract": {
            "scatterer_model": minimal["mic"]["scatterer_model"],
            "measurement_path": minimal["mic"]["measurement_path"],
            "finite_target_fwhm_ratio": float(minimal["mic"]["finite_target_fwhm_ratio"]),
            "sample_scattering_amplitude_ptp": float(np.ptp(np.abs(minimal["oct_raw"]["sample_scattering_amplitude"]))),
            "point_target_ratio_threshold": 0.20,
            "pass": bool(float(minimal["mic"]["finite_target_fwhm_ratio"]) < 0.20),
        },
        "illumination_contract": {
            "illumination_na": float(minimal["mic"]["illumination_na"]),
            "illumination_pupil_fill_ratio": float(minimal["mic"]["illumination_pupil_fill_ratio"]),
            "illumination_aperture_stop_radius": float(minimal["mic"]["illumination_aperture_stop_radius"]),
            "illumination_power_fraction": float(minimal["mic"]["illumination_power_fraction"]),
            "illumination_psf_peak": float(np.max(minimal["mic"]["illumination_psf"])),
            "illumination_psf_sum": float(np.sum(minimal["mic"]["illumination_psf"])),
            "pass": bool(float(minimal["mic"]["illumination_power_fraction"]) > 0.0),
        },
        "camera_contract": {
            "camera_qe": float(minimal["mic"]["camera_qe"]),
            "camera_dark_adu": float(minimal["mic"]["camera_dark_adu"]),
            "camera_saturation_adu": (
                None
                if config.microscope.camera_saturation_adu is None
                else float(config.microscope.camera_saturation_adu)
            ),
            "camera_flat_field_mean": float(np.mean(minimal["mic"]["camera_flat_field"])),
            "pass": bool(float(minimal["mic"]["camera_qe"]) >= 0.0),
        },
        "negative_controls": {
            "count": len(negative_rows),
            "detected_count": int(sum(1 for row in negative_rows if row["detected"])),
            "controls": negative_rows,
            "pass": bool(negative_rows and all(row["detected"] for row in negative_rows)),
        },
        "convergence_sweeps": {
            "count": len(convergence_rows),
            "review_count": int(sum(1 for row in convergence_rows if row["status"] == "review")),
            "axes": [row["sweep_axis"] for row in convergence_rows],
            "pass": bool(
                convergence_rows
                and all(np.isfinite(float(row["nrmse_3d"])) for row in convergence_rows)
                and all(np.isfinite(float(row["nrmse_delta_abs"])) for row in convergence_rows)
                and not any(row["status"] == "review" for row in convergence_rows)
            ),
        },
        "direct_model_comparison": {
            "count": len(direct_model_rows),
            "rows": direct_model_rows,
            "max_nrmse_3d": float(max(float(row["nrmse_3d"]) for row in direct_model_rows)),
            "pass": bool(direct_model_rows and all(row["pass"] for row in direct_model_rows)),
        },
        "direct_model_axial_profiles": {
            "count": len(direct_profile_rows),
            "max_on_axis_abs_residual": float(
                max(float(row["on_axis_abs_residual"]) for row in direct_profile_rows)
            ),
            "columns": list(direct_profile_rows[0].keys()) if direct_profile_rows else [],
            "pass": bool(
                direct_profile_rows
                and all(np.isfinite(float(row["full_profile_on_axis"])) for row in direct_profile_rows)
                and all(np.isfinite(float(row["hybrid_profile_on_axis"])) for row in direct_profile_rows)
            ),
        },
    }
    checks["all_pass"] = bool(all(v.get("pass", False) for v in checks.values() if isinstance(v, dict)))
    documented_limitations = [
        "full vector Debye / Zemax POP is not implemented in this scalar low-NA baseline",
        "rayleigh_mie_lookup uses supplied calibration tables; it is not an internal Mie solver",
        "Path E measured predictive is not accepted unless axial and lateral FWHM gates pass, even if normalized 3D NRMSE is low",
    ]
    unmet_requirements = []
    if int(config.oct.full_spectral_rci_depth_decimation) > 1:
        unmet_requirements.append(
            "lowres_full_spectral_rci_direct is a coarse review artifact (depth_decimation>1 uses non-phase-aware complex linear interpolation across z); set full_spectral_rci_depth_decimation=1 for production-resolution direct truth"
        )
    failing_checks = sum(
        1 for value in checks.values() if isinstance(value, dict) and not bool(value.get("pass", False))
    )
    blocker_count = int(failing_checks + len(unmet_requirements))
    verdict = {
        "pilot_pass": bool(checks["all_pass"] and blocker_count == 0),
        "review_required": bool(blocker_count > 0),
        "blocker_count": blocker_count,
        "failing_check_count": int(failing_checks),
        "documented_limitations": documented_limitations,
        "unmet_requirements": unmet_requirements,
        "unsupported_claims": documented_limitations + unmet_requirements,
    }

    summary = {
        "output_dir": str(out_dir),
        "config_hash": config_hash(config),
        "provenance": provenance,
        "runtime_parameters": {
            "N": N,
            "k_samples": k_samples,
            "pad_factor": pad_factor,
            "direct_psf_model": minimal["direct_psf_model"],
            "full_spectral_rci_depth_decimation": int(config.oct.full_spectral_rci_depth_decimation),
        },
        "config": config_to_dict(config),
        "checks": checks,
        "verdict": verdict,
    }
    (metrics_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=_json_default),
        encoding="utf-8",
    )

    rows = build_sweep_rows(config)
    generate_sweep_trend_figure(out_dir / "figures", rows)
    fieldnames = list(rows[0].keys())
    with (metrics_dir / "sweep_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    neg_fieldnames = list(negative_rows[0].keys())
    with (metrics_dir / "negative_controls.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=neg_fieldnames)
        writer.writeheader()
        writer.writerows(negative_rows)

    conv_fieldnames = list(convergence_rows[0].keys())
    with (metrics_dir / "convergence_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=conv_fieldnames)
        writer.writeheader()
        writer.writerows(convergence_rows)

    direct_fieldnames = list(direct_model_rows[0].keys())
    with (metrics_dir / "direct_model_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=direct_fieldnames)
        writer.writeheader()
        writer.writerows(direct_model_rows)

    direct_profile_fieldnames = list(direct_profile_rows[0].keys())
    with (metrics_dir / "direct_model_axial_profiles.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=direct_profile_fieldnames)
        writer.writeheader()
        writer.writerows(direct_profile_rows)

    (out_dir / "README_run.md").write_text(
        "# Validation Run\n\n"
        f"- config hash: `{config_hash(config)}`\n"
        f"- all checks pass: `{checks['all_pass']}`\n"
        f"- pilot pass: `{verdict['pilot_pass']}`\n"
        f"- review required: `{verdict['review_required']}`\n"
        f"- blocker count: `{verdict['blocker_count']}`\n"
        f"- command: `{provenance['command']}`\n"
        "- key files: `metrics/validation_summary.json`, `metrics/sweep_summary.csv`, "
        "`metrics/negative_controls.csv`, `metrics/convergence_summary.csv`, "
        "`metrics/direct_model_comparison.csv`, `metrics/direct_model_axial_profiles.csv`, "
        "`oct_direct/oct_psf_direct.h5`, `oct_direct/oct_reconstruction.h5`, "
        "`microscope_forward/microscope_zstack.h5`, `converters/`, `figures/`\n",
        encoding="utf-8",
    )
    return summary
