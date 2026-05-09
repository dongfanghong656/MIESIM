## T5.1 — Strehl-ratio V-gate (2026-05-09, SHA 60565b1)

### What changed
- `cop_oct_sim/metrics.py:13` adds `strehl_ratio(psf_actual, psf_ideal) -> float`.
- `cop_oct_sim/validation.py:36` adds `ValidationConfig.strehl_threshold = 0.80`.
- `cop_oct_sim/validation.py:104` adds `_strehl_v_gate(...)`, and `cop_oct_sim/validation.py:854` registers `checks["strehl_v_gate"]`.
- `tests/test_sanity.py:222`, `tests/test_sanity.py:735`, and `tests/test_sanity.py:747` add the three T5.1 regression tests.
- The gate records `strehl_normalization = "integrated_intensity"` because the current direct-PSF output is already peak-normalized; applying the task-book's literal peak-normalized formula to this output gives heavy-defocus Strehl = 1.0 and cannot satisfy the specified fail test.

### Verification
- Local AST/in-memory compile passed for `cop_oct_sim/metrics.py`, `cop_oct_sim/validation.py`, and `tests/test_sanity.py`.
- Local direct execution passed:
  - `test_strehl_ratio_is_one_for_identical_psfs`
  - `test_strehl_v_gate_passes_for_unaberrated_minimal_config`
  - `test_strehl_v_gate_fails_when_defocus_is_heavy`
- CI run: https://github.com/dongfanghong656/MIESIM/actions/runs/25597551164 — success.
- CI artifact `validation-outputs-py3.11` includes two `validation_summary.json` files with `checks.strehl_v_gate.strehl_ratio = 1`, `checks.strehl_v_gate.pass = true`, `checks.strehl_v_gate.strehl_threshold = 0.8`, and `verdict.pilot_pass = true`.

### Open items carried forward
- None for T5.1. If PRO requires the literal peak-normalized Strehl formula, the direct PSF pipeline must expose an unnormalized PSF or center-intensity metadata; using the already peak-normalized returned PSF is mathematically insensitive to defocus.

## T5.2 — Differential dispersion in common-path (2026-05-09, SHA 5544a37)

### What changed
- `cop_oct_sim/spectrometer.py:111` renames the absolute sample-arm quadratic phase helper to `absolute_dispersion_phase(...)`.
- `cop_oct_sim/spectrometer.py:115` adds `differential_dispersion_phase(...)`, using `(dispersion_quadratic_rad - reference_dispersion_quadratic_rad) * xi(k)^2`.
- `cop_oct_sim/config_schema.py:89` adds `ErrorConfig.reference_dispersion_quadratic_rad = 0.0`.
- `cop_oct_sim/oct_forward.py:293` and `cop_oct_sim/oct_forward.py:381` switch direct OCT field assembly to `differential_dispersion_phase(...)`.
- `cop_oct_sim/theory_conversion.py:62` switches the predicted axial gate to the same differential phase convention.
- `PHYSICS_CONTRACT.md:31` documents the common-path dispersion identity and the sample/reference mismatch definition.
- `tests/test_sanity.py:521` and `tests/test_sanity.py:535` add the cancellation and mismatch-broadening tests.

### Verification
- Local full pytest: `python -m pytest tests/test_sanity.py -q` -> 44 passed.
- Local mismatch check: clean axial FWHM = 7.332735312545693 um; mismatch axial FWHM = 8.453789106279611 um; broadening ratio = 1.152883439255729.
- CI run: https://github.com/dongfanghong656/MIESIM/actions/runs/25603263330 — success.
- CI artifact `validation-outputs-py3.11` includes two `validation_summary.json` files with `verdict.pilot_pass = true` and `checks.strehl_v_gate.pass = true`.
- CI `direct_model_comparison.csv` retains the Round-4 faithfulness floor: `axial_fwhm_relative_error = 7.232718201331539e-05`, `nrmse_3d = 3.851839665003354e-06`, `full_spectral_rci_interpolated = False`.

### Open items carried forward
- None for T5.2.

## T5.3 — Sensitivity-rolloff V-gate (2026-05-09, SHA 3a48597)

### What changed
- `cop_oct_sim/spectrometer.py:125` adds amplitude roll-off unit conversion helpers.
- `cop_oct_sim/spectrometer.py:131` adds `theoretical_sensitivity_rolloff_db_per_mm(config)` using the Round 5 6.7 dB-at-half-range convention.
- `cop_oct_sim/spectrometer.py:140` adds `effective_sensitivity_rolloff_per_um(config)`.
- `cop_oct_sim/oct_forward.py:26` applies the theoretical spectrometer roll-off plus `errors.rolloff_per_um` as an additional empirical error term.
- `cop_oct_sim/validation.py:138` adds measured roll-off fitting from forward-model `oct_raw` roll-off amplitudes.
- `cop_oct_sim/validation.py:167` adds `checks["sensitivity_rolloff_v_gate"]`.
- `tests/test_sanity.py:574`, `tests/test_sanity.py:584`, and `tests/test_sanity.py:589` add the three T5.3 tests.
- `PHYSICS_CONTRACT.md:46` documents the sensitivity roll-off convention.

### Verification
- Local full pytest: `python -m pytest tests/test_sanity.py -q` -> 47 passed.
- CI run: https://github.com/dongfanghong656/MIESIM/actions/runs/25603613595 — success.
- CI artifact `validation-outputs-py3.11` includes two `validation_summary.json` files with `checks.sensitivity_rolloff_v_gate.pass = true`.
- CI roll-off slopes:
  - minimal config: theoretical = 5.11320193525492 dB/mm, measured = 5.11320193525492 dB/mm, relative error = 1.38962384970754e-15.
  - production smoke: theoretical = 110.176045910177 dB/mm, measured = 110.176045910177 dB/mm, relative error = 1.28983161428643e-16.
- CI `direct_model_comparison.csv` retains the Round-4 faithfulness floor: `axial_fwhm_relative_error = 7.232718201331539e-05`, `nrmse_3d = 3.851839665003354e-06`, `full_spectral_rci_interpolated = False`.

### Open items carried forward
- None for T5.3.

## T5.4 — Phase-stability via propagated reflector (2026-05-09, SHA ca8f262)

### What changed
- `cop_oct_sim/validation.py` adds `ValidationConfig.phase_stability_max_std_rad`, `phase_stability_repeats`, and `phase_stability_reflector_z_um`.
- `cop_oct_sim/validation.py` adds `_phase_stability_v_gate(...)`, which simulates a single reflector with `simulate_oct_raw_direct`, adds camera read noise through the existing `ErrorConfig.camera_read_noise_e` knob, reconstructs with `reconstruct_sd_oct`, and measures phase standard deviation across repeats.
- `cop_oct_sim/validation.py` registers `checks["phase_stability_v_gate"]`.
- `tests/test_sanity.py` adds `test_phase_stability_uses_propagated_reflector` and `test_phase_stability_passes_on_minimal`.
- `PHYSICS_CONTRACT.md` documents the propagated-reflector phase-stability convention.

### Verification
- Local AST parse passed for `cop_oct_sim/validation.py` and `tests/test_sanity.py`.
- Local targeted pytest: `python -m pytest tests/test_sanity.py -q -k "phase_stability"` -> 2 passed.
- Local full pytest: `python -m pytest tests/test_sanity.py -q` -> 49 passed.
- CI run: https://github.com/dongfanghong656/MIESIM/actions/runs/25604099700 — success.
- CI artifact `validation-outputs-py3.11` includes two `validation_summary.json` files with `checks.phase_stability_v_gate.pass = true`.
- CI phase-stability values:
  - minimal config: `phase_std_rad = 2.5616344142886695e-06`, `reflector_z_um = 20.0`, `n_repeats = 8`, `camera_read_noise_e = 2.0`.
  - production smoke: `phase_std_rad = 3.262574608972163e-06`, `reflector_z_um = 20.0`, `n_repeats = 8`, `camera_read_noise_e = 2.0`.
- CI `direct_model_comparison.csv` retains the Round-4 faithfulness floor: `axial_fwhm_relative_error = 7.232718201331539e-05`, `nrmse_3d = 3.851839665003354e-06`, `full_spectral_rci_interpolated = False`.

### Open items carried forward
- None for T5.4.

## T5.5 — Document + collapse common-path pupil identity (N4) (2026-05-09, SHA dc2dc0f)

### What changed
- `cop_oct_sim/oct_forward.py` adds `COMMON_PATH_PUPIL_IDENTITY_CONTRACT = "scalar_low_na_common_path_incident_detection_pupils_identical"`.
- `cop_oct_sim/oct_forward.py` adds `_common_path_rci_from_shared_pupil(...)`, collapsing the scalar common-path `h_RCI = u_i * conj(u_d)` form to `u * conj(u)` for one shared OCT pupil.
- `_through_focus_rci_stack(...)`, `_full_spectral_rci_direct_psf(...)`, `simulate_oct_raw_direct(...)`, and `simulate_oct_psf_direct(...)` now use the shared helper rather than rebuilding identical incident/detection pupils.
- Direct OCT outputs now report `common_path_pupil_identity` and `common_path_pupil_identity_contract`; full-spectral RCI metadata also reports the same contract.
- `tests/test_sanity.py` adds identity-collapse and metadata tests.
- `PHYSICS_CONTRACT.md` documents the N4 common-path pupil identity and states that future vector Debye or separated-pupil work must introduce a new explicit contract.

### Verification
- Local AST parse passed for `cop_oct_sim/oct_forward.py` and `tests/test_sanity.py`.
- Local targeted pytest: `python -m pytest tests/test_sanity.py -q -k "common_path_rci or full_spectral_rci_reports_common_path"` -> 2 passed.
- Local pre/post numerical reference check on a `N=12`, `k_samples=24`, full-spectral RCI run:
  - `psf_max_abs_diff = 0.0`
  - `h_max_abs_diff = 0.0`
- Local full pytest: `python -m pytest tests/test_sanity.py -q` -> 51 passed.
- CI run: https://github.com/dongfanghong656/MIESIM/actions/runs/25604462482 — success.
- CI artifact `validation-outputs-py3.11` includes two `validation_summary.json` files with `verdict.pilot_pass = true` and the existing V-gates still passing.
- CI `direct_model_comparison.csv` retains the Round-4 faithfulness floor exactly at the tracked values: `axial_fwhm_relative_error = 7.232718201331539e-05`, `nrmse_3d = 3.851839665003354e-06`, `full_spectral_rci_interpolated = False`.

### Open items carried forward
- None for T5.5. The scalar common-path identity is now explicit; vector Debye or separated illumination/detection pupils remain Round 6 scope.

## T5.6 — Round-4 response package (2026-05-10, SHA e138504 + c1f96d0)

### What changed
- `scripts/plot_round4_axial_lines.py` generates four on-axis axial magnitude PNGs at `scatterer_z_um = 0, 5, 10, 20`.
- `scripts/build_review_package.py --round 4` builds `review_packages/round4_response_<timestamp>/response.docx` plus lightweight CSV/PNG evidence and a `package_manifest.json`.
- The response DOCX covers the Round-4 P0-by-P0 response, faithfulness number, T1.3 finding, and Round-5 V-gate evidence.
- `.github/workflows/ci.yml` now runs `python scripts/build_review_package.py --round 4` after the validation suites and uploads `round4-response-package-py3.11`.
- `tests/test_sanity.py` adds a DOCX structure test for the Round-4 response writer.
- Follow-up fix `c1f96d0` keeps the response package lightweight by copying only metrics and plots; full H5 outputs remain in the validation artifact.

### Verification
- Local AST parse passed for `scripts/build_review_package.py`, `scripts/plot_round4_axial_lines.py`, and `tests/test_sanity.py`.
- Local plot smoke: `python scripts/plot_round4_axial_lines.py --config configs/config_minimal.yaml --output-root .codex-tmp/round4_axial_lines --N 12 --k-samples 24 --pad-factor 1` generated four PNGs.
- Local package smoke with a fake validation directory produced `response.docx`, metrics evidence, four axial-line PNGs, `SHA256SUMS.txt`, and a zip package.
- Local full pytest: `python -m pytest tests/test_sanity.py -q` -> 52 passed.
- First CI run for `e138504` succeeded, but revealed an oversized `round4-response-package-py3.11` artifact because full H5 validation outputs were copied into the response package.
- Fix CI run: https://github.com/dongfanghong656/MIESIM/actions/runs/25606208085 — success.
- Corrected CI artifacts:
  - `round4-response-package-py3.11`: `168847` bytes, containing one `response.docx`, four PNGs, and `package_manifest.json`.
  - `validation-outputs-py3.11`: includes two validation summaries with `verdict.pilot_pass = true`.
- CI `direct_model_comparison.csv` retains the Round-4 faithfulness floor: `axial_fwhm_relative_error = 7.232718201331539e-05`, `nrmse_3d = 3.851839665003354e-06`, `full_spectral_rci_interpolated = False`.

### Open items carried forward
- The DOCX is generated as a compact standards-compliant OOXML file without embedded PNGs; the PNGs are delivered as sibling evidence files in the CI artifact.
