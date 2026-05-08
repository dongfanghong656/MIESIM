from dataclasses import replace
from pathlib import Path
import shutil
import uuid

import h5py
import numpy as np

from cop_oct_sim.config_schema import IlluminationConfig, SampleConfig, SimulationConfig, load_config
from cop_oct_sim.grids import make_pupil_grid
from cop_oct_sim.jones import apply_jones_to_scalar_pupil
from cop_oct_sim.microscope_forward import simulate_microscope_zstack
from cop_oct_sim.oct_forward import simulate_oct_psf_direct, simulate_oct_raw_direct
from cop_oct_sim.pipelines import run_minimal
from cop_oct_sim.propagation import fraunhofer_psf_from_pupil, physical_defocus_phase, propagation_energy_ratio, simple_defocus_stack
from cop_oct_sim.pupil import build_shared_pupil
from cop_oct_sim.reconstruction import reconstruct_sd_oct
from cop_oct_sim.scatterers import sample_scattering_amplitude
from cop_oct_sim.theory_conversion import predict_axial_gate_from_source
from cop_oct_sim.validation import (
    airy_sanity,
    axial_bandwidth_sanity,
    build_convergence_rows,
    build_direct_model_axial_profile_rows,
    build_direct_model_comparison_rows,
    build_negative_control_rows,
    dichroic_sanity,
    gaussian_underfill_sanity,
    no_leak_sanity,
    parseval_sanity,
    run_validation_suite,
)
from cop_oct_sim.metrics import (
    axial_fwhm,
    central_lobe_fwhm_1d,
    lateral_fwhm_from_volume,
    nrmse,
    normalize,
    on_axis_axial_fwhm,
)

def small_config() -> SimulationConfig:
    cfg = SimulationConfig()
    return replace(
        cfg,
        microscope=replace(cfg.microscope, wavelength_nm=cfg.oct.center_wavelength_nm, z_count=15),
        oct=replace(cfg.oct, spectrometer_pixels=96),
        errors=replace(cfg.errors, photon_gain_e_per_adu=5000.0),
    )

def test_microscope_shapes_and_reproducibility():
    cfg = small_config()
    mic1 = simulate_microscope_zstack(cfg, N=32, pad_factor=2)
    mic2 = simulate_microscope_zstack(cfg, N=32, pad_factor=2)
    assert mic1["measured"].ndim == 3
    assert mic1["ideal"].shape == mic1["measured"].shape
    np.testing.assert_allclose(mic1["measured"], mic2["measured"])
    assert mic1["raw"].shape == mic1["measured"].shape
    assert mic1["background"].shape == mic1["measured"].shape
    np.testing.assert_allclose(mic1["raw"], mic1["measured"])
    np.testing.assert_allclose(mic1["corrected"], mic1["raw"] - mic1["background"])
    assert mic1["camera_x_um"].shape[0] == mic1["measured"].shape[-1]
    assert mic1["camera_y_um"].shape[0] == mic1["measured"].shape[-2]

def test_camera_pixel_integration_can_downsample_to_camera_grid():
    base = small_config()
    cfg = replace(
        base,
        microscope=replace(
            base.microscope,
            camera_pixel_um=10.0,
            magnification_total=1.0,
            z_count=5,
        ),
        errors=replace(base.errors, photon_gain_e_per_adu=1.0e8, camera_read_noise_e=0.0),
    )
    mic = simulate_microscope_zstack(cfg, N=32, pad_factor=2)

    assert mic["camera_integrated"].shape == mic["measured"].shape
    assert mic["measured"].shape[-1] < mic["ideal"].shape[-1]
    assert float(mic["camera_pixel_sample_um"]) > float(mic["dx_um"])
    assert mic["camera_x_um"].shape[0] == mic["measured"].shape[-1]

def test_camera_qe_file_scales_microscope_signal():
    work = Path.cwd() / ".test-output" / f"camera-qe-{uuid.uuid4().hex}"
    work.mkdir(parents=True, exist_ok=True)
    try:
        qe_path = work / "camera_qe.csv"
        qe_path.write_text(
            "wavelength_nm,qe\n"
            "540,0.50\n"
            "550,0.50\n"
            "560,0.50\n",
            encoding="utf-8",
        )
        base = small_config()
        cfg_ref = replace(
            base,
            microscope=replace(base.microscope, wavelength_nm=550.0),
            errors=replace(base.errors, photon_gain_e_per_adu=1.0e8, camera_read_noise_e=0.0),
        )
        cfg_qe = replace(
            cfg_ref,
            microscope=replace(cfg_ref.microscope, camera_qe_file=str(qe_path)),
        )

        ref = simulate_microscope_zstack(cfg_ref, N=24, pad_factor=1)
        qe = simulate_microscope_zstack(cfg_qe, N=24, pad_factor=1)

        assert float(qe["camera_qe"]) == 0.5
        np.testing.assert_allclose(np.max(qe["camera_integrated"]), 0.5 * np.max(ref["camera_integrated"]), rtol=0.02)
    finally:
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)

def test_camera_dark_flat_field_and_saturation_are_applied():
    work = Path.cwd() / ".test-output" / f"camera-chain-{uuid.uuid4().hex}"
    work.mkdir(parents=True, exist_ok=True)
    try:
        flat_path = work / "flat_field.csv"
        np.savetxt(flat_path, np.full((4, 4), 0.5), delimiter=",")
        base = small_config()
        cfg = replace(
            base,
            microscope=replace(
                base.microscope,
                camera_flat_field_file=str(flat_path),
                camera_dark_adu=0.10,
                camera_saturation_adu=0.40,
            ),
            errors=replace(base.errors, photon_gain_e_per_adu=1.0e8, camera_read_noise_e=0.0),
        )

        mic = simulate_microscope_zstack(cfg, N=24, pad_factor=1)

        assert "camera_flat_field" in mic
        assert float(mic["camera_dark_adu"]) == 0.10
        assert float(mic["camera_saturation_adu"]) == 0.40
        np.testing.assert_allclose(np.mean(mic["camera_flat_field"]), 0.5)
        assert np.max(mic["raw"]) <= 0.400001
        assert np.min(mic["background"]) >= 0.10
    finally:
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)

def test_illumination_pupil_limits_signal_and_is_reported():
    base = small_config()
    cfg_full = replace(
        base,
        errors=replace(base.errors, photon_gain_e_per_adu=1.0e8, camera_read_noise_e=0.0),
    )
    cfg_limited = replace(
        cfg_full,
        illumination=IlluminationConfig(
            NA=0.08,
            pupil_fill_ratio=0.45,
            aperture_stop_radius=0.55,
            field_stop_diameter_um=10.0,
        ),
    )

    full = simulate_microscope_zstack(cfg_full, N=32, pad_factor=2)
    limited = simulate_microscope_zstack(cfg_limited, N=32, pad_factor=2)

    assert "illumination_pupil" in limited
    assert "illumination_envelope" in limited
    assert limited["illumination_envelope"].shape == limited["ideal"].shape[-2:]
    assert float(limited["illumination_power_fraction"]) < float(full["illumination_power_fraction"])
    assert np.sum(np.abs(limited["illumination_pupil"]) ** 2) < np.sum(np.abs(full["illumination_pupil"]) ** 2)
    assert np.max(limited["camera_integrated"]) < np.max(full["camera_integrated"])

def test_illumination_na_changes_epi_response_shape():
    base = small_config()
    quiet = replace(base.errors, photon_gain_e_per_adu=1.0e8, camera_read_noise_e=0.0)
    cfg_full = replace(base, errors=quiet)
    cfg_limited = replace(
        cfg_full,
        illumination=IlluminationConfig(NA=0.07, pupil_fill_ratio=0.45, aperture_stop_radius=0.5),
    )

    full = simulate_microscope_zstack(cfg_full, N=32, pad_factor=2)
    limited = simulate_microscope_zstack(cfg_limited, N=32, pad_factor=2)

    assert "illumination_psf" in limited
    full_norm = normalize(full["scattering_corrected"])
    limited_norm = normalize(limited["scattering_corrected"])
    assert not np.allclose(full_norm, limited_norm)
    assert lateral_fwhm_from_volume(limited_norm, float(limited["dx_um"])) > lateral_fwhm_from_volume(
        full_norm,
        float(full["dx_um"]),
    )

def test_fraunhofer_default_preserves_energy():
    grid = make_pupil_grid(48, 0.25)
    pupil = build_shared_pupil(grid, 0.85, fill_ratio=0.82)
    field = fraunhofer_psf_from_pupil(pupil, pad_factor=2)

    assert abs(propagation_energy_ratio(pupil, field) - 1.0) < 1e-12

def test_jones_mean_does_not_double_power():
    pupil = np.ones((8, 8), dtype=np.complex128)
    mean = apply_jones_to_scalar_pupil(pupil, s_amplitude=1.0, p_amplitude=1.0, input_pol="mean")
    attenuated = apply_jones_to_scalar_pupil(pupil, s_amplitude=0.5, p_amplitude=1.0, input_pol="mean")

    np.testing.assert_allclose(np.sum(np.abs(mean) ** 2), np.sum(np.abs(pupil) ** 2))
    expected_attenuated_power = np.sum(np.abs(pupil) ** 2) * (0.5**2 + 1.0**2) / 2.0
    np.testing.assert_allclose(np.sum(np.abs(attenuated) ** 2), expected_attenuated_power)

def test_defocus_stack_preserves_energy_without_plane_peak_normalization():
    grid = make_pupil_grid(48, 0.25)
    pupil = build_shared_pupil(grid, 0.85, fill_ratio=0.82)
    z_um = np.linspace(-8.0, 8.0, 9)

    stack = simple_defocus_stack(pupil, z_um, 0.85, 0.25, pad_factor=2)
    peaks = np.max(np.abs(stack), axis=(1, 2))
    energy = np.sum(np.abs(stack) ** 2, axis=(1, 2))

    assert np.ptp(peaks) > 0.05 * np.max(peaks)
    np.testing.assert_allclose(energy / energy[len(energy) // 2], np.ones_like(energy), rtol=1e-12)

def test_direct_oct_no_import_leak():
    assert no_leak_sanity()["pass"]

def test_oct_reconstruction_runs():
    cfg = small_config()
    raw = simulate_oct_raw_direct(cfg, N=24, k_samples=96)
    z = reconstruct_sd_oct(raw["raw"], raw["k"])
    assert z.ndim == 1
    assert len(z) == 96
    for key in ["reference_field", "sample_field", "background", "interference", "pixel_index", "k_linear"]:
        assert key in raw
    np.testing.assert_allclose(raw["raw"], raw["background"] + raw["interference"])
    assert np.all(np.diff(raw["k_linear"]) > 0)

def test_oct_raw_direct_samples_spatial_rci_not_scalar_peak():
    cfg = small_config()
    centered = simulate_oct_raw_direct(cfg, N=32, k_samples=96, scatterer_x_um=0.0)
    off_axis = simulate_oct_raw_direct(cfg, N=32, k_samples=96, scatterer_x_um=5.0)

    centered_mod = np.std(centered["raw"] - np.mean(centered["raw"]))
    off_axis_mod = np.std(off_axis["raw"] - np.mean(off_axis["raw"]))

    assert "h_rci_sample" in centered
    assert centered_mod > 1.2 * off_axis_mod
    assert not np.allclose(centered["h_rci_sample"], off_axis["h_rci_sample"])

def test_path_e_predictive_axial_gate_is_separate_from_oracle():
    cfg = small_config()
    gate = predict_axial_gate_from_source(cfg, k_samples=96)
    result = run_minimal(cfg, N=24, k_samples=96, pad_factor=2)

    assert gate["axial_profile"].shape == (96,)
    assert gate["depth_um"].shape == (96,)
    assert "path_e_predictive" in result
    assert "path_e_oracle" in result
    assert "path_e_measured_predictive" in result
    assert "path_e_ideal_upper_bound" in result
    assert "metrics_predictive" in result
    assert "metrics_oracle" in result
    assert "metrics_measured_predictive" in result
    assert "metrics_ideal_upper_bound" in result
    np.testing.assert_allclose(result["path_e_converted"], result["path_e_measured_predictive"])
    np.testing.assert_allclose(result["path_e_oracle"], result["path_e_ideal_upper_bound"])

def test_path_e_predictive_uses_measured_microscope_stack():
    cfg = replace(
        small_config(),
        microscope=replace(small_config().microscope, background_adu=0.2),
        errors=replace(small_config().errors, photon_gain_e_per_adu=10.0, camera_read_noise_e=0.0),
    )
    result = run_minimal(cfg, N=24, k_samples=96, pad_factor=2)

    assert not np.allclose(result["mic"]["ideal"], result["mic"]["corrected"])
    np.testing.assert_allclose(result["path_e_predictive"], result["path_e_measured_predictive"])
    assert not np.allclose(result["path_e_measured_predictive"], result["path_e_ideal_upper_bound"])

def test_oct_direct_uses_pixel_lambda_and_source_files():
    work = Path.cwd() / ".test-output" / f"spectral-{uuid.uuid4().hex}"
    work.mkdir(parents=True, exist_ok=True)
    try:
        n = 32
        pixel = np.arange(n)
        wavelength_nm = 840.0 + 0.15 * pixel + 0.002 * pixel**2
        source = 0.2 + 0.8 * np.exp(-0.5 * ((wavelength_nm - wavelength_nm.mean()) / 1.2) ** 2)
        pixel_lambda_path = work / "pixel_lambda.csv"
        source_path = work / "source.csv"
        np.savetxt(
            pixel_lambda_path,
            np.column_stack([pixel, wavelength_nm]),
            delimiter=",",
            header="pixel,wavelength_nm",
            comments="",
        )
        np.savetxt(
            source_path,
            np.column_stack([wavelength_nm, source]),
            delimiter=",",
            header="wavelength_nm,source",
            comments="",
        )

        base = small_config()
        cfg = replace(
            base,
            oct=replace(
                base.oct,
                spectrometer_pixels=n,
                pixel_to_lambda_file=str(pixel_lambda_path),
                source_spectrum_file=str(source_path),
            ),
        )

        raw = simulate_oct_raw_direct(cfg, N=24)

        np.testing.assert_allclose(raw["wavelength_nm"], wavelength_nm)
        np.testing.assert_allclose(raw["source"], source / np.max(source))
        assert not np.allclose(np.diff(raw["k_true"]), np.diff(raw["k_true"])[0])
    finally:
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)

def test_dichroic_table_separates_mic_transmission_and_oct_reflection():
    work = Path.cwd() / ".test-output" / f"dichroic-{uuid.uuid4().hex}"
    work.mkdir(parents=True, exist_ok=True)
    try:
        table = work / "dichroic.csv"
        table.write_text(
            "wavelength_nm,AOI,Rs,Rp,Ts,Tp,phase_s,phase_p\n"
            "660,45,0.04,0.09,0.25,0.36,0.10,0.20\n"
            "850,45,0.81,0.64,0.01,0.04,0.30,0.40\n",
            encoding="utf-8",
        )
        base = small_config()
        cfg = replace(
            base,
            microscope=replace(base.microscope, wavelength_nm=660.0),
            dichroic=replace(base.dichroic, table_file=str(table)),
        )

        mic = simulate_microscope_zstack(cfg, N=24, pad_factor=1)
        raw = simulate_oct_raw_direct(cfg, N=24, k_samples=96)

        assert mic["dichroic_path"] == "mic_transmission"
        assert raw["dichroic_path"] == "oct_reflection"
        np.testing.assert_allclose(float(mic["dichroic_s_amplitude"]), 0.5)
        np.testing.assert_allclose(float(mic["dichroic_p_amplitude"]), 0.6)
        np.testing.assert_allclose(float(raw["dichroic_s_amplitude"]), 0.9)
        np.testing.assert_allclose(float(raw["dichroic_p_amplitude"]), 0.8)
    finally:
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)

def test_legacy_dichroic_coefficients_are_common_path_factors():
    base = small_config()
    cfg = replace(
        base,
        dichroic=replace(
            base.dichroic,
            amplitude_s=0.5,
            amplitude_p=0.25,
            mic_transmission_amplitude_s=0.8,
            mic_transmission_amplitude_p=0.6,
            oct_reflection_amplitude_s=0.7,
            oct_reflection_amplitude_p=0.4,
        ),
    )

    mic = simulate_microscope_zstack(cfg, N=24, pad_factor=1)
    raw = simulate_oct_raw_direct(cfg, N=24, k_samples=96)

    np.testing.assert_allclose(float(mic["dichroic_s_amplitude"]), 0.4)
    np.testing.assert_allclose(float(mic["dichroic_p_amplitude"]), 0.15)
    np.testing.assert_allclose(float(raw["dichroic_s_amplitude"]), 0.35)
    np.testing.assert_allclose(float(raw["dichroic_p_amplitude"]), 0.10)

def test_direct_oct_psf_shapes():
    cfg = small_config()
    direct = simulate_oct_psf_direct(cfg, N=24, pad_factor=2, k_samples=96)
    assert direct["psf"].shape == (96, 48, 48)
    assert direct["magnitude"].shape == direct["psf"].shape
    assert direct["rci_through_focus"].shape == direct["psf"].shape
    assert np.isfinite(direct["magnitude"]).all()

def test_direct_oct_psf_nonseparable_with_through_focus_defocus():
    base = small_config()
    cfg = replace(
        base,
        objective=replace(base.objective, defocus_um=2.5),
        errors=replace(base.errors, objective_focus_shift_um=2.0),
    )
    direct = simulate_oct_psf_direct(cfg, N=24, pad_factor=2, k_samples=96)
    vol = normalize(np.abs(direct["psf"]))
    z = int(np.argmax(np.max(vol, axis=(1, 2))))
    axial = normalize(np.max(vol, axis=(1, 2)))
    spatial = normalize(vol[z])
    separable = axial[:, None, None] * spatial[None, :, :]

    assert float(direct["separability_nrmse"]) > 1e-3
    assert nrmse(vol, separable) > 1e-3

def test_lowres_full_spectral_rci_direct_psf_is_reported():
    cfg = small_config()
    direct = simulate_oct_psf_direct(cfg, N=20, pad_factor=1, k_samples=32, full_spectral_rci=True)

    assert direct["psf"].shape == (32, 20, 20)
    assert "spectral_rci_direct_psf" in direct
    assert "spectral_rci_depth_um" in direct
    assert bool(direct["full_spectral_rci"])
    assert direct["spectral_rci_direct_psf"].shape == direct["psf"].shape
    assert np.isfinite(direct["spectral_rci_direct_psf"]).all()
    assert direct["spectral_rci_sampled_h_rci"].shape == (
        len(direct["spectral_rci_computed_depth_indices"]),
        32,
        20,
        20,
    )
    for key in [
        "spectral_rci_k_true",
        "spectral_rci_k_measured",
        "spectral_rci_wavelength_nm",
        "spectral_rci_dx_um_per_k",
        "spectral_rci_source",
        "spectral_rci_window",
        "spectral_rci_dispersion_phase_rad",
        "spectral_rci_axial_profile_max",
        "spectral_rci_axial_profile_on_axis",
        "spectral_rci_axial_profile_integrated",
    ]:
        assert key in direct
        assert np.isfinite(direct[key]).all()
    for key in [
        "spectral_rci_phase_sign",
        "spectral_rci_normalization",
        "spectral_rci_model_contract",
    ]:
        assert key in direct
        assert isinstance(direct[key], str)
    assert "scatterer_z" in direct["spectral_rci_phase_sign"]
    assert float(direct["spectral_rci_depth_phase_origin_um"]) == 0.0

def test_pipeline_can_promote_full_spectral_rci_to_direct_truth():
    base = small_config()
    cfg = replace(base, oct=replace(base.oct, direct_psf_model="full_spectral_rci"))
    result = run_minimal(cfg, N=12, k_samples=24, pad_factor=1)

    assert result["direct_psf_model"] == "full_spectral_rci"
    assert bool(result["oct_direct"]["full_spectral_rci"])
    assert "spectral_rci_direct_psf" in result["oct_direct"]
    np.testing.assert_allclose(result["oct_direct"]["psf"], result["oct_direct"]["spectral_rci_direct_psf"])

def test_full_spectral_rci_depth_decimation_is_reported():
    base = small_config()
    cfg = replace(
        base,
        oct=replace(
            base.oct,
            direct_psf_model="full_spectral_rci",
            full_spectral_rci_depth_decimation=4,
        ),
    )
    direct = simulate_oct_psf_direct(cfg, N=12, pad_factor=1, k_samples=24, full_spectral_rci=True)

    indices = direct["spectral_rci_computed_depth_indices"]
    assert direct["psf"].shape == (24, 12, 12)
    assert int(direct["full_spectral_rci_depth_decimation"]) == 4
    assert bool(direct["full_spectral_rci_interpolated"])
    assert len(indices) < 24
    assert 0 in set(indices.tolist())
    assert 12 in set(indices.tolist())
    assert 23 in set(indices.tolist())

def test_physical_defocus_small_z_keeps_on_axis_rci_reasonable():
    grid = make_pupil_grid(64, 0.25)
    pupil = build_shared_pupil(grid, 0.85, fill_ratio=1.0, defocus_um=0.0)
    focused = fraunhofer_psf_from_pupil(pupil, pad_factor=2)
    small_z = fraunhofer_psf_from_pupil(
        pupil * physical_defocus_phase(grid, z_um=2.5, wavelength_um=0.85, NA=0.25),
        pad_factor=2,
    )
    cy = focused.shape[0] // 2
    cx = focused.shape[1] // 2
    focused_rci = focused * np.conj(focused)
    small_z_rci = small_z * np.conj(small_z)
    ratio = abs(small_z_rci[cy, cx]) / max(abs(focused_rci[cy, cx]), 1e-12)

    assert ratio > 0.5

def test_central_lobe_fwhm_ignores_mirror_peak():
    profile = np.zeros(101, dtype=float)
    profile[20] = 1.0
    profile[19] = 0.75
    profile[21] = 0.75
    profile[80] = 0.95
    profile[79] = 0.70
    profile[81] = 0.70

    assert central_lobe_fwhm_1d(profile, dx=1.0, peak_index=20) < 4.0
    assert axial_fwhm(profile[:, None, None], np.arange(101, dtype=float)) > 50.0

def test_on_axis_axial_fwhm_uses_center_pixel_not_off_axis_max():
    depth = np.arange(101, dtype=float)
    volume = np.zeros((101, 5, 5), dtype=float)
    volume[50, 2, 2] = 1.0
    volume[49, 2, 2] = 0.75
    volume[51, 2, 2] = 0.75
    volume[15, 0, 0] = 0.99
    volume[85, 4, 4] = 0.99

    assert on_axis_axial_fwhm(volume, depth) < 4.0
    assert axial_fwhm(volume, depth) > 60.0

def test_dispersion_residual_broadens_axial_psf():
    cfg = small_config()
    clean = simulate_oct_psf_direct(cfg, N=24, pad_factor=1, k_samples=96)
    dispersed_cfg = replace(
        cfg,
        errors=replace(cfg.errors, dispersion_quadratic_rad=25.0),
    )
    dispersed = simulate_oct_psf_direct(dispersed_cfg, N=24, pad_factor=1, k_samples=96)
    assert axial_fwhm(dispersed["psf"], dispersed["depth_um"]) > axial_fwhm(clean["psf"], clean["depth_um"])

def test_k_linearization_error_returns_perturbed_measured_grid():
    cfg = small_config()
    errored_cfg = replace(
        cfg,
        errors=replace(cfg.errors, k_linearization_rms_pixel=2.0),
    )
    raw = simulate_oct_raw_direct(errored_cfg, N=24, k_samples=96)
    assert "k_true" in raw
    assert not np.allclose(raw["k"], raw["k_true"])
    assert np.all(np.diff(np.sort(raw["k"])) > 0)

def test_rolloff_reduces_deep_raw_modulation():
    cfg = replace(
        small_config(),
        errors=replace(small_config().errors, rolloff_per_um=0.08),
    )
    shallow = simulate_oct_raw_direct(cfg, N=24, k_samples=96, scatterer_z_um=0.0)
    deep = simulate_oct_raw_direct(cfg, N=24, k_samples=96, scatterer_z_um=20.0)
    shallow_mod = np.std(shallow["raw"] - np.mean(shallow["raw"]))
    deep_mod = np.std(deep["raw"] - np.mean(deep["raw"]))
    assert deep_mod < shallow_mod

def test_configured_finite_bead_broadens_microscope_stack_axially():
    base = small_config()
    cfg_point = replace(
        base,
        microscope=replace(base.microscope, z_count=121, z_step_um=0.5, bead_diameter_um=0.0),
        errors=replace(base.errors, photon_gain_e_per_adu=1.0e10, camera_read_noise_e=0.0),
    )
    cfg_bead = replace(
        cfg_point,
        microscope=replace(cfg_point.microscope, bead_diameter_um=4.0),
    )

    point = simulate_microscope_zstack(cfg_point, N=32, pad_factor=2)
    bead = simulate_microscope_zstack(cfg_bead, N=32, pad_factor=2)

    assert float(bead["bead_diameter_um"]) == 4.0
    assert axial_fwhm(bead["measured"], bead["z_um"]) > 1.1 * axial_fwhm(point["measured"], point["z_um"])

def test_scattering_measurement_paths_are_reported_and_different():
    base = small_config()
    epi_cfg = replace(
        base,
        sample=SampleConfig(
            scatterer_model="uniform_sphere_projection",
            scatterer_diameter_um=0.6,
            measurement_path="epi_backscatter",
            background_amplitude=0.4,
            background_subtract=True,
        ),
        errors=replace(base.errors, photon_gain_e_per_adu=1.0e8, camera_read_noise_e=0.0),
    )
    dark_cfg = replace(
        epi_cfg,
        sample=replace(epi_cfg.sample, measurement_path="oblique_darkfield", darkfield_background_fraction=0.02),
    )

    epi = simulate_microscope_zstack(epi_cfg, N=32, pad_factor=2)
    dark = simulate_microscope_zstack(dark_cfg, N=32, pad_factor=2)

    assert epi["measurement_path"] == "epi_backscatter"
    assert dark["measurement_path"] == "oblique_darkfield"
    assert epi["scatterer_model"] == "uniform_sphere_projection"
    assert float(epi["finite_target_fwhm_ratio"]) > 0.0
    assert "scattering_raw" in epi
    assert "scattering_background" in epi
    assert not np.allclose(epi["scattering_raw"], dark["scattering_raw"])

def test_rayleigh_mie_lookup_modulates_oct_raw_and_reports_amplitude():
    work = Path.cwd() / ".test-output" / f"scatter-{uuid.uuid4().hex}"
    work.mkdir(parents=True, exist_ok=True)
    try:
        lookup = work / "scattering_lookup.csv"
        lookup.write_text(
            "wavelength_nm,amplitude,phase_rad\n"
            "830,0.25,0.0\n"
            "850,1.00,0.2\n"
            "870,0.50,0.4\n",
            encoding="utf-8",
        )
        base = small_config()
        cfg = replace(
            base,
            sample=SampleConfig(
                scatterer_model="rayleigh_mie_lookup",
                scatterer_diameter_um=0.8,
                measurement_path="epi_backscatter",
                scattering_lookup_file=str(lookup),
            ),
        )

        amp = sample_scattering_amplitude(cfg, np.array([830.0, 850.0, 870.0]))
        raw = simulate_oct_raw_direct(cfg, N=24, k_samples=96)

        np.testing.assert_allclose(np.abs(amp), [0.25, 1.0, 0.5])
        assert raw["scatterer_model"] == "rayleigh_mie_lookup"
        assert "sample_scattering_amplitude" in raw
        assert np.ptp(np.abs(raw["sample_scattering_amplitude"])) > 0.1
        assert float(raw["scatterer_diameter_um"]) == 0.8
    finally:
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)

def test_physics_sanity_gates():
    cfg = small_config()
    assert airy_sanity()["pass"]
    assert gaussian_underfill_sanity()["pass"]
    assert parseval_sanity()["pass"]
    assert dichroic_sanity()["pass"]
    assert axial_bandwidth_sanity(cfg)["pass"]

def test_negative_controls_are_detected():
    rows = build_negative_control_rows(small_config(), N=24, k_samples=64, pad_factor=2)
    by_name = {row["control_name"]: row for row in rows}

    assert {"chromatic_mismatch", "wrong_dichroic_path", "bad_k_calibration"} <= set(by_name)
    assert all(row["detected"] for row in rows)
    assert by_name["chromatic_mismatch"]["bad_nrmse_3d"] > by_name["chromatic_mismatch"]["baseline_nrmse_3d"]
    assert by_name["wrong_dichroic_path"]["path_amplitude_delta"] > 0.5
    assert by_name["bad_k_calibration"]["bad_k_rms_pixel_error"] > 2.0

def test_convergence_sweeps_cover_required_axes():
    rows = build_convergence_rows(small_config(), N=24, k_samples=64, pad_factor=2)
    by_axis = {row["sweep_axis"]: row for row in rows}

    assert {"N", "pad_factor", "k_samples", "z_step_um", "window", "bead_diameter_um", "camera_pixel_um"} <= set(by_axis)
    assert all(np.isfinite(float(row["nrmse_3d"])) for row in rows)
    assert all(np.isfinite(float(row["nrmse_delta_abs"])) for row in rows)
    assert all(row["status"] in {"stable", "review"} for row in rows)

def test_direct_model_comparison_rows_compare_hybrid_and_full_spectral():
    base = small_config()
    cfg = replace(base, oct=replace(base.oct, full_spectral_rci_depth_decimation=3))
    rows = build_direct_model_comparison_rows(cfg, N=12, k_samples=24, pad_factor=1)

    assert len(rows) == 1
    row = rows[0]
    assert row["comparison_name"] == "hybrid_vs_full_spectral_rci_lowres"
    assert row["hybrid_direct_model"] == "hybrid_rci"
    assert row["full_direct_model"] == "full_spectral_rci"
    assert int(row["full_spectral_rci_depth_decimation"]) == 3
    assert int(row["computed_depth_count"]) < int(row["depth_count"])
    assert row["full_spectral_rci_interpolated"] is True
    assert np.isfinite(float(row["nrmse_3d"]))
    assert float(row["axial_fwhm_relative_error"]) <= 0.25
    assert row["status"] == "diagnostic"
    assert row["pass"] is True

def test_direct_model_axial_profile_rows_report_hybrid_and_full_spectral_profiles():
    base = small_config()
    cfg = replace(base, oct=replace(base.oct, full_spectral_rci_depth_decimation=3))
    rows = build_direct_model_axial_profile_rows(cfg, N=12, k_samples=24, pad_factor=1)

    assert len(rows) == 24
    first = rows[0]
    for key in [
        "depth_um",
        "full_profile_max",
        "full_profile_on_axis",
        "full_profile_integrated",
        "hybrid_profile_max",
        "hybrid_profile_on_axis",
        "hybrid_profile_integrated",
        "on_axis_abs_residual",
    ]:
        assert key in first
        assert np.isfinite(float(first[key]))
    assert max(float(row["full_profile_max"]) for row in rows) == 1.0
    assert max(float(row["hybrid_profile_max"]) for row in rows) == 1.0

def test_config_paths_resolve_relative_to_project_root_from_subdir():
    original_cwd = Path.cwd()
    try:
        subdir = original_cwd / "cop_oct_sim"
        import os

        os.chdir(subdir)
        cfg = load_config("../configs/config_visible_660.yaml")
    finally:
        os.chdir(original_cwd)

    assert Path(cfg.microscope.camera_qe_file).is_absolute()
    assert Path(cfg.microscope.camera_qe_file).exists()
    assert Path(cfg.oct.source_spectrum_file).is_absolute()
    assert Path(cfg.oct.source_spectrum_file).exists()
    assert Path(cfg.dichroic.table_file).is_absolute()
    assert Path(cfg.dichroic.table_file).exists()

def test_pilot_pass_false_when_review_required():
    base = small_config()
    cfg = replace(
        base,
        oct=replace(base.oct, direct_psf_model="full_spectral_rci", full_spectral_rci_depth_decimation=3),
    )
    output_root = Path.cwd() / ".test-output" / f"pytest-review-required-{uuid.uuid4().hex}"
    try:
        summary = run_validation_suite(output_root, cfg, N=12, k_samples=24, pad_factor=1)
        assert summary["verdict"]["review_required"] is True
        assert summary["verdict"]["pilot_pass"] is False
        assert summary["verdict"]["blocker_count"] > 0
    finally:
        if output_root.exists():
            shutil.rmtree(output_root, ignore_errors=True)

def test_validation_runner_writes_data():
    output_root = Path.cwd() / ".test-output" / f"pytest-{uuid.uuid4().hex}"
    try:
        summary = run_validation_suite(output_root, small_config(), N=24, k_samples=64, pad_factor=2)
        out = Path(summary["output_dir"])
        assert summary["config"]["microscope"]["wavelength_nm"] == small_config().microscope.wavelength_nm
        assert summary["runtime_parameters"]["direct_psf_model"] == small_config().oct.direct_psf_model
        assert summary["runtime_parameters"]["full_spectral_rci_depth_decimation"] == 1
        assert summary["checks"]["all_pass"] is True
        assert summary["verdict"]["pilot_pass"] is True
        assert summary["verdict"]["review_required"] is False
        assert summary["verdict"]["blocker_count"] == 0
        assert summary["verdict"]["failing_check_count"] == 0
        assert summary["verdict"]["unmet_requirements"] == []
        assert isinstance(summary["verdict"]["documented_limitations"], list)
        assert isinstance(summary["verdict"]["unsupported_claims"], list)
        assert summary["checks"]["path_e_measured_predictive"]["pass"] is True
        assert summary["checks"]["direct_model_comparison"]["pass"] is True
        assert "axial_fwhm_relative_error" in summary["checks"]["path_e_measured_predictive"]
        assert "path_e_ideal_upper_bound" in summary["checks"]
        assert "scatterer_contract" in summary["checks"]
        assert "illumination_contract" in summary["checks"]
        assert "camera_contract" in summary["checks"]
        assert "negative_controls" in summary["checks"]
        assert summary["checks"]["negative_controls"]["pass"] is True
        assert "convergence_sweeps" in summary["checks"]
        assert summary["checks"]["convergence_sweeps"]["pass"] is True
        assert (out / "metrics" / "validation_summary.json").exists()
        assert (out / "metrics" / "sweep_summary.csv").exists()
        assert (out / "metrics" / "negative_controls.csv").exists()
        assert (out / "metrics" / "convergence_summary.csv").exists()
        assert (out / "metrics" / "direct_model_comparison.csv").exists()
        assert (out / "metrics" / "direct_model_axial_profiles.csv").exists()
        assert "direct_model_comparison" in summary["checks"]
        assert summary["checks"]["direct_model_comparison"]["pass"] is True
        assert "direct_model_axial_profiles" in summary["checks"]
        assert summary["checks"]["direct_model_axial_profiles"]["pass"] is True
        assert (out / "oct_direct" / "oct_psf_direct.h5").exists()
        assert (out / "oct_direct" / "oct_psf_direct_lowres_full_spectral_rci.h5").exists()
        assert (out / "oct_direct" / "oct_reconstruction.h5").exists()
        assert (out / "microscope_forward" / "illumination_pupil.h5").exists()
        assert (out / "converters" / "path_e_measured_predictive_psf.h5").exists()
        assert (out / "converters" / "path_e_ideal_upper_bound_psf.h5").exists()
        for name in [
            "mic_central_xy.png",
            "mic_xz.png",
            "mic_yz.png",
            "oct_a_scan.png",
            "path_e_residual_xy.png",
            "sweep_trend.png",
        ]:
            assert (out / "figures" / name).exists()
        with h5py.File(out / "oct_direct" / "oct_raw_interferogram.h5", "r") as f:
            for key in ["raw", "background", "interference", "reference_field", "sample_field", "k_linear"]:
                assert key in f
        with h5py.File(out / "microscope_forward" / "microscope_zstack.h5", "r") as f:
            for key in ["raw", "background", "corrected", "camera_integrated", "bead_convolved", "camera_flat_field"]:
                assert key in f
        with h5py.File(out / "microscope_forward" / "illumination_pupil.h5", "r") as f:
            for key in ["illumination_pupil", "illumination_psf", "illumination_envelope", "illumination_power_fraction"]:
                assert key in f
        with h5py.File(out / "oct_direct" / "oct_psf_direct_lowres_full_spectral_rci.h5", "r") as f:
            for key in [
                "spectral_rci_sampled_h_rci",
                "spectral_rci_k_true",
                "spectral_rci_k_measured",
                "spectral_rci_wavelength_nm",
                "spectral_rci_axial_profile_on_axis",
                "spectral_rci_depth_phase_origin_um",
            ]:
                assert key in f
    finally:
        if output_root.exists():
            shutil.rmtree(output_root, ignore_errors=True)
