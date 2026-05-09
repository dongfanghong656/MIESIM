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
