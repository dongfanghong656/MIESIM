# Round 4 Plan — Verify 222556 Claims, Then Promote full_spectral_rci to Production Truth

Owner: Claude (planner+executor) · Reviewer: PRO
Date: 2026-04-30
Replaces: prior round4_plan.md (was based on hallucinated paths; deleted)

## Context

PRO audited package `20260428-095635` and flagged 6 P0s (defocus unit, hybrid-as-truth, full-spectral-not-truth, pass-gate, config paths, FWHM metrics).

Codex's later package `20260428-222556` CHANGELOG **claims all P0s addressed**. PRO has NOT re-audited 222556. The claims must be independently verified before any new work is started.

Real package layout (flat): `cop_oct_sim/{config_schema, errors, figures, grids, io_hdf5, jones, measurement, metrics, microscope_forward, oct_forward, optical_paths, pipelines, propagation, pupil, reconstruction, scatterers, spectrometer, theory_conversion, validation}.py`.

## Phase 0 — Independent verification of 222556 claims (BLOCKING)

Do not start new physics work until each below is checked. For every claim, write a short note in `status/round4_verification_log.md` with file:line evidence + a passing/failing test command.

### V0.1 — Defocus split is correct
- Confirm `_make_oct_pupil` (cop_oct_sim/oct_forward.py:28-48) sends physical z through `physical_defocus_phase` and Zernike OPD through `build_shared_pupil(defocus_um=...)`.
- Confirm `physical_defocus_phase` formula in `cop_oct_sim/propagation.py` matches PHYSICS_CONTRACT.md line 18-20: `phase(rho) = pi * z_um * NA^2 * rho^2 / wavelength_um`.
- Run any new regression test Codex shipped (search `tests/` for `defocus`).
- Acceptance: at NA=0.05, λ=850 nm, z=2.5 µm physical defocus produces an analytic Strehl drop matching the paraxial formula within 5 %.

### V0.2 — Gate semantics fixed
- Read `cop_oct_sim/validation.py` near the verdict assembly. Confirm `pilot_pass = checks.all_pass AND blocker_count == 0 AND not unsupported_review_required`.
- Run `python scripts/run_validation_suite.py --config configs/config_full_spectral_rci_production_smoke.yaml --strict-pass`. Verify exit code 1 when any blocker is present.
- Acceptance: the previously contradictory `all_pass=True ∧ review_required=True ∧ pilot_pass=True` state is impossible.

### V0.3 — direct_model_comparison gates on FWHM, not just NRMSE
- Locate the comparison logic. Confirm both 3D NRMSE AND axial-FWHM relative error (≤10 % per spec/PRO) gate it.
- Acceptance: hand-craft a synthetic case where NRMSE<0.05 but axial FWHM differs 3× → must FAIL.

### V0.4 — Config relative paths
- Confirm load_config in `cop_oct_sim/config_schema.py` resolves CSV/YAML relatives against `Path(config_path).parent`.
- Run `cd review_packages/.../20260428-222556/project && python -m scripts.run_validation_suite --config ../configs/config_visible_660.yaml`. Should NOT FileNotFoundError on `calibration_data/camera_qe.csv`.
- Acceptance: same config invoked from two different CWDs returns the same resolved absolute paths.

### V0.5 — Central-lobe / on-axis FWHM metrics exist and are wired
- Find `central_lobe_fwhm_1d`, `on_axis_axial_fwhm`, `integrated_axial_width` in `cop_oct_sim/metrics.py`.
- Confirm they are exercised in `validation.py` axial-FWHM gate, not just `axial_fwhm` (max-over-XY).
- Acceptance: a synthetic two-peak axial profile (mirror peak) is correctly handled by central-lobe extractor.

### V0.6 — Run pytest
- `python -m pytest -q` over the full repo. Capture pass count + elapsed.
- Run `python scripts/run_validation_suite.py --config configs/config_minimal.yaml` and the full-spectral smoke. Capture verdict JSON.

## Phase 1 — Promote full_spectral_rci to production direct truth

Only after Phase 0 verification log is clean.

### T1.1 — Remove low-resolution caveats
- The full_spectral_rci output currently uses depth_decimation>1 + complex interpolation across z. PRO flagged this as not physical truth (mixes phase and amplitude). Make production path require decimation=1 OR add an explicit phase-aware reinterpolation justified by a unit test.
- Files: `cop_oct_sim/oct_forward.py:_full_spectral_rci_direct_psf` (line 187-273), `_interpolate_complex_stack` (line 166-185).

### T1.2 — One direct truth API
- `simulate_oct_psf_direct(..., full_spectral_rci=True)` becomes the documented truth. `hybrid_rci` (the post-FFT axial_gate × through_focus product, line 400-406) is renamed/marked diagnostic-only and excluded from `direct_model_comparison` pass criterion (it can still be reported as a diagnostic).

### T1.3 — Negative test against hybrid
- Add a test that intentionally injects a low-NA aberration and asserts `hybrid_rci` and `full_spectral_rci` axial FWHM disagree by >25 % → blocker fires.

## Phase 2 — Remaining V-gate completeness (only after Phase 1 lands)

- Add Strehl ratio gate (PRO 04 mentions it as missing).
- Tighten dispersion compensation: ensure differential, not absolute, dispersion in common-path geometry.
- Replace synthetic-noise input in phase stability test with propagated reflector signal.
- Sensitivity roll-off vs spec depth-of-focus formula.

## Phase 3 — Closed-loop measurement (deferred to round 5)

- Finite-size scattering ball with dispersion + spherical aberration; round-trip a measured PSF.
- Köhler illumination separation: aperture stop vs field stop as independent params.
- Pupil-angle-resolved dichroic with AOI table.

## Out of scope this round

- Polarization, full Jones, Mie internal solver, Zemax POP import, FDTD spot checks. CHANGELOG says these are explicit boundaries; do not start.

## Exit criteria for round-4 review package

1. Phase 0 verification log committed under `status/round4_verification_log.md` with all 6 items signed off.
2. Phase 1 promoted full_spectral_rci is the default in `config_full_spectral_rci_production_smoke.yaml`; hybrid is diagnostic only.
3. Round-4 response DOCX (≤6 pages) showing axial-line plots at z ∈ {0, 5, 10, 20} µm for the corrected defocus path, plus the negative-test result.
4. PRO sign-off on each of their 6 original P0s.
