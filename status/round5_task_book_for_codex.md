# Round 5 Task Book — for Codex

Author: Claude (planner)
Date: 2026-05-08
Target repo: `dongfanghong656/MIESIM` branch `main`
Prereq read: [status/round4_verification_log.md](round4_verification_log.md), [PHYSICS_CONTRACT.md](../PHYSICS_CONTRACT.md), PRO Round 3 DOCX reports 01–08.

## 0 Starting state — what Round 4 delivered

All five PRO P0s structurally resolved and a reproducible CI baseline exists.

| PRO Round 3 P0 | Round 4 resolution | Evidence |
|---|---|---|
| 1. Defocus unit conflation | Physical-z and Zernike-OPD paths fully separated | [`cop_oct_sim/oct_forward.py:28-48`](../cop_oct_sim/oct_forward.py), [`propagation.py:30-38`](../cop_oct_sim/propagation.py), [`pupil.py:9-13`](../cop_oct_sim/pupil.py) |
| 2. Hybrid as truth | `OCTConfig.direct_psf_model` default flipped to `"full_spectral_rci"`; hybrid is opt-in | [`cop_oct_sim/config_schema.py:24`](../cop_oct_sim/config_schema.py) |
| 3. Pass-gate semantics | N1: `documented_limitations` + `unmet_requirements`; N2: gates on on-axis FWHM | [`validation.py:59-61,390-392,836-857`](../cop_oct_sim/validation.py) |
| 4. Config relative paths | `_resolve_optional_path` + `resolve_config_paths` | [`config_schema.py:119-158`](../cop_oct_sim/config_schema.py) |
| 5. Central-lobe / on-axis FWHM | Functions exist + wired into gates (N2) | [`metrics.py:34-88`](../cop_oct_sim/metrics.py), [`validation.py:59-61,390-392`](../cop_oct_sim/validation.py) |
| Full-spectral low-res caveat | T1.1 real fix: `max(...,3) → max(...,1)` in `build_direct_model_*_rows` | [`validation.py:358,451`](../cop_oct_sim/validation.py) |

Known-good CI baseline (reference for regression detection):

- Latest green run: [`25555010839`](https://github.com/dongfanghong656/MIESIM/actions/runs/25555010839) @ `c1a4d94`.
- Production smoke at `decimation=1`, default `direct_psf_model=full_spectral_rci`:
  - `verdict.pilot_pass=True`, `blocker_count=0`, `failing_check_count=0`.
  - `direct_model_comparison.csv`: on-axis `axial_fwhm_relative_error = 7.2e-5`; `nrmse_3d = 3.85e-6`; `full_spectral_rci_interpolated = False`; `computed_depth_count == depth_count == 32`.
  - These numbers are the **faithfulness floor**. Any Round 5 change that pushes them above ~1% without physics justification is a regression.

Environment contract (confirmed with user twice, see memory [feedback_codex_ci_pytest.md](../../.claude/projects/c--Users-1-OneDrive---fzu-edu-cn--1--Attachments/memory/feedback_codex_ci_pytest.md)):

- Local (Windows, system Py 3.13) has no `h5py`, no `pytest`; corp PyPI proxy blocks `pip install`. Do **not** work around this — it is by design.
- Local sanity is `python -m py_compile` and `python -c "import ast; ast.parse(...)"` only.
- Validation, pytest, strict smoke all run in GitHub Actions ([`.github/workflows/ci.yml`](../.github/workflows/ci.yml)).
- Push → wait for CI → read artifacts is the standard loop.

## 1 Round 5 goals

1. Close remaining V-gate items PRO flagged in reports 02–08 (Strehl, differential dispersion, sensitivity roll-off, phase-stability realism).
2. Document or refactor the N4 pupil-identity finding (`pupil_i ≡ pupil_d` in `_full_spectral_rci_direct_psf`).
3. Produce the Round-4 response DOCX PRO expects for P0 sign-off.
4. Land everything through the existing `status/roundN_*.md` + `tests/test_sanity.py` + CI-green discipline — no new process.

Anything beyond these four items (vector Debye, Mie internal solver, Zemax POP, FDTD, Köhler separation, pupil-angle-resolved dichroic) is **Round 6**, not Round 5.

## 2 Task breakdown

Each task below is self-contained: scope, files/lines, acceptance criteria. Execute in order (T5.1 → T5.6) to minimize merge friction; each task ends with a commit and a CI run to lock in the baseline before the next task touches overlapping files.

### T5.1 Strehl-ratio V-gate

**Why:** PRO 04 flagged Strehl missing. We already compute aberrated and ideal PSFs; Strehl is a one-liner metric that gates focus quality independently of FWHM.

**Scope:**

- Add `strehl_ratio(psf_actual, psf_ideal) -> float` in [`cop_oct_sim/metrics.py`](../cop_oct_sim/metrics.py). Definition: `|psf_actual(0,0,0)|² / |psf_ideal(0,0,0)|²`, both peak-normalized by `max |.|²` before extraction. Guard against division by zero.
- Add a validation check in [`validation.py`](../cop_oct_sim/validation.py) named `strehl_v_gate` that:
  - Builds an ideal reference PSF by calling `simulate_oct_psf_direct` with `config.objective.defocus_um=0.0` and `config.errors.*` zeroed (use `dataclasses.replace` pattern, same as `build_direct_model_comparison_rows` does for hybrid/full).
  - Computes `strehl = strehl_ratio(actual, ideal)`.
  - Pass if `strehl >= 0.80` for the unaberrated minimal config; keep threshold field-configurable via a new `ValidationConfig` entry `strehl_threshold: float = 0.80`.
  - Record both the raw Strehl value and the pass/fail in `checks["strehl_v_gate"]` so it shows up in `validation_summary.json` and feeds `blocker_count` on failure.
- Register the gate in the `checks` dict aggregation that feeds `all_pass`.

**Tests (append to `tests/test_sanity.py`):**

- `test_strehl_ratio_is_one_for_identical_psfs`: identical inputs → `strehl_ratio == 1.0` within 1e-12.
- `test_strehl_v_gate_passes_for_unaberrated_minimal_config`: use `small_config()`; assert check present and `pass=True`.
- `test_strehl_v_gate_fails_when_defocus_is_heavy`: inject `objective=replace(base.objective, defocus_um=2.0)`; assert check `pass=False` (large defocus collapses the on-axis amplitude).

**Acceptance:**

- Local: `python -m py_compile cop_oct_sim/metrics.py cop_oct_sim/validation.py tests/test_sanity.py` clean.
- CI: green. `validation_summary.json.checks.strehl_v_gate.strehl_ratio` present in minimal + production-smoke artifacts; minimal strehl ≥ 0.80.
- `verdict.pilot_pass` still `True` on both configs (otherwise either the threshold is wrong or the ideal-reference construction is wrong — investigate before lowering threshold).

**Commit message prefix:** `round5 T5.1: strehl-ratio V-gate`.

### T5.2 Differential-dispersion modeling in common-path geometry

**Why:** PRO 02 + 04 noted the current dispersion term is absolute, but in a common-path interferometer only the **mismatch** between reference and sample arm should appear in reconstruction phase. Absolute dispersion is cancelled by definition. The current `dispersion_phase(config, k)` applied identically to both arms would cancel in any faithful implementation — leaving it in place as an absolute spectral phase in `_full_spectral_rci_direct_psf` is a modeling inconsistency, not just a naming issue.

**Scope:**

- Read [`cop_oct_sim/spectrometer.py`](../cop_oct_sim/spectrometer.py) `dispersion_phase`. Rename the existing function to `absolute_dispersion_phase` (keep for diagnostic) and introduce `differential_dispersion_phase(config, k_array)` that returns **only the mismatch** phase.
- Add to [`config_schema.py`](../cop_oct_sim/config_schema.py) `ErrorConfig`:
  - `dispersion_quadratic_rad` (existing) — document it as the differential term.
  - New field `reference_dispersion_quadratic_rad: float = 0.0`.
  - `differential_dispersion_phase` returns the quadratic phase for `(dispersion_quadratic_rad - reference_dispersion_quadratic_rad)`.
- Update both `simulate_oct_psf_direct` (line ~381) and `_full_spectral_rci_direct_psf` (line ~193) to consume `differential_dispersion_phase` instead of `dispersion_phase`.
- Update [`PHYSICS_CONTRACT.md`](../PHYSICS_CONTRACT.md) with the common-path dispersion identity and the differential definition.

**Tests:**

- `test_dispersion_phase_cancels_when_reference_matches_sample`: set both quadratic terms equal; assert `max|differential_dispersion_phase|` < 1e-12 over the k-grid.
- `test_dispersion_mismatch_broadens_axial_fwhm`: set `dispersion_quadratic_rad=5.0`, `reference_dispersion_quadratic_rad=0.0`; compare full_spectral axial FWHM against baseline (both zero); assert broadening > 10 %.

**Acceptance:**

- CI green. On the existing minimal + production-smoke configs (both have `dispersion_quadratic_rad=0`), `axial_fwhm_relative_error` remains at the numerical floor (~1e-4), confirming no regression from the rename/refactor.
- `differential_dispersion_phase` referenced in `_full_spectral_rci_direct_psf` and `simulate_oct_psf_direct` only; `absolute_dispersion_phase` used only in diagnostic reporting.

**Commit message prefix:** `round5 T5.2: differential dispersion in common-path`.

### T5.3 Sensitivity roll-off vs depth-of-focus V-gate

**Why:** PRO 04 listed sensitivity roll-off as a required V-gate. Current `ErrorConfig.rolloff_per_um` exists but is neither tied to a spec nor gated.

**Scope:**

- Add to [`cop_oct_sim/spectrometer.py`](../cop_oct_sim/spectrometer.py) a closed-form `theoretical_sensitivity_rolloff_db_per_mm(config)` derived from spectrometer pixel count and source bandwidth (use the standard formula: roll-off ≈ `-6.7 dB` at `z_max/2` scaled by spectral resolution; derive the concrete form from PRO 04 or cite standard Leitgeb/Yun formula in `PHYSICS_CONTRACT.md`).
- Add a validation check `sensitivity_rolloff_v_gate` in [`validation.py`](../cop_oct_sim/validation.py):
  - Measures effective roll-off from `oct_raw` or `oct_complex` by fitting amplitude decay vs z.
  - Computes `relative_error = |measured - theoretical| / theoretical`.
  - Pass if `relative_error <= 0.15`.
- Gate wired into `checks["sensitivity_rolloff_v_gate"]` following the same pattern as `path_e_measured_predictive`.

**Tests:**

- `test_theoretical_sensitivity_rolloff_matches_leitgeb_form`: unit test on the closed-form, 3 spec points (center, z_max/2, z_max).
- `test_sensitivity_rolloff_v_gate_passes_for_minimal`: `small_config()` passes.
- `test_sensitivity_rolloff_v_gate_fails_when_rolloff_grossly_wrong`: inject `errors.rolloff_per_um = 10 * theoretical_per_um`; assert `pass=False`.

**Acceptance:**

- CI green. Both minimal and production smoke surface a finite, non-NaN `sensitivity_rolloff_db_per_mm` in the summary JSON and pass the gate.

**Commit message prefix:** `round5 T5.3: sensitivity-rolloff V-gate`.

### T5.4 Replace synthetic-noise phase-stability test with propagated reflector signal

**Why:** PRO noted `_phase_stability_noise_gate` injects noise directly into the phase vector rather than propagating a reflector signal through the forward model + reconstruction. The current test is therefore a measurement of `numpy.random` behavior, not of the pipeline.

**Scope:**

- Locate the phase-stability gate in [`validation.py`](../cop_oct_sim/validation.py) (grep `phase_stability`). Replace the synthetic-noise path with:
  1. Simulate a single point reflector via `simulate_oct_raw_direct`.
  2. Add `camera_read_noise_e` via the existing `ErrorConfig` knob.
  3. Reconstruct with `reconstruct_sd_oct`.
  4. Compute phase std over M repeats (`M=8` in CI, parameterize for local quicker runs).
- Keep the acceptance criterion (phase std ≤ spec) but drive the input through the actual forward model.

**Tests:**

- `test_phase_stability_uses_propagated_reflector`: assert the test helper function now calls `simulate_oct_raw_direct` (check via a patch-based spy, or simply by structural assertion that the result dict contains a key only the reflector path produces).
- `test_phase_stability_passes_on_minimal`: end-to-end through the new path.

**Acceptance:**

- CI green. Phase-stability summary JSON entry now records `reflector_z_um`, `n_repeats`, `camera_read_noise_e` alongside `phase_std_rad`.

**Commit message prefix:** `round5 T5.4: phase-stability via propagated reflector`.

### T5.5 N4 — incident/detection pupil identity in `_full_spectral_rci_direct_psf`

**Why:** [`oct_forward.py:224-225`](../cop_oct_sim/oct_forward.py) constructs `pupil_i` and `pupil_d` with identical parameters, so `h_rci = ui * conj(ud)` reduces to `|u|²`. For the scalar low-NA common-path geometry this is physically correct (both arms share the same objective pupil in the same field), but it is currently undocumented. PRO did not flag it; surfacing it prevents a future silent divergence.

**Preferred resolution:** document, do not refactor.

**Scope:**

- Collapse the two pupil calls to one: `pupil = ...; ui = fraunhofer_psf_from_pupil(pupil, ...); ud = ui` (or equivalently `h_rci = np.abs(ui)**2`).
- Extend the same simplification to `_through_focus_rci_stack` [`oct_forward.py:146-147`](../cop_oct_sim/oct_forward.py) for consistency.
- Add a short section to [`PHYSICS_CONTRACT.md`](../PHYSICS_CONTRACT.md) titled "Common-path pupil identity" stating: in this scalar low-NA common-path model the incident and detection optical paths share one physical objective pupil; any future bi-directional geometry (folded vs unfolded, polarization diversity, beam-splitter asymmetry) must split them.
- Do **not** introduce a config flag to re-enable two distinct pupils — YAGNI.

**Tests:**

- `test_h_rci_equals_mod_squared_of_u_in_common_path`: sample one `(z, k)` in `_full_spectral_rci_direct_psf` and assert `np.allclose(h_rci, np.abs(ui)**2)` by computing both and diffing, to pin the assumption.
- No numerical regression should appear in any existing test — this is algebraically identical.

**Acceptance:**

- CI green. No change to any `validation_summary.json` numeric field at the floating-point level (≤ 1e-12 delta).

**Commit message prefix:** `round5 T5.5: document + collapse common-path pupil identity (N4)`.

### T5.6 Round-4 response DOCX

**Why:** PRO expects a ≤6-page DOCX response summarizing how each of their Round 3 P0s was closed, for sign-off.

**Scope (follow Round 3 response template if one exists under `review_packages/`):**

- Cover page: title, date, version, scope (scalar low-NA common-path OCT validation harness).
- Section 1: Changelog table (reuse the table in §0 of this task book).
- Section 2: P0-by-P0 response — one paragraph per P0, evidence as file:line references + CI run URL.
- Section 3: Axial-line plots. Generate via a new script [`scripts/plot_round4_axial_lines.py`](../scripts/plot_round4_axial_lines.py):
  - Takes a `--config` path, runs `simulate_oct_psf_direct` at `scatterer_z_um ∈ {0, 5, 10, 20}` µm, saves 4 PNGs of on-axis axial magnitude profile vs z.
  - Writes to `outputs/round4_axial_lines/`.
  - CI job uploads these as an artifact.
- Section 4: Faithfulness number (on-axis axial `7.2e-5` relative error at decimation=1, with reference to artifact path).
- Section 5: T1.3 finding — scalar low-NA hybrid is faithful; negative test reframed as positive pin.
- Section 6: Known limitations (the 3 entries now in `documented_limitations`).
- Section 7: Newly identified issues (N4, handed off to Round 6).

**Acceptance:**

- DOCX ≤ 6 pages, stored as `review_packages/round4_response_<timestamp>/response.docx` plus the supporting CSV/PNG evidence in a sibling folder.
- Confirm it builds cleanly on CI (new workflow step `Build round4 response package` invoking a new script `scripts/build_review_package.py --round 4`, which may already exist — extend, don't duplicate).

**Commit message prefix:** `round5 T5.6: round-4 response package`.

## 3 Environment + conventions

### Commit / PR discipline

- One commit per T5.X. No bundling across tasks.
- Commit message: subject line `round5 T5.X: <short>`, body explains **why** + **how tested** + **CI run URL** once green.
- No `--no-verify`, no `--force`. If a push is rejected, diagnose (likely: CI failed prior run depends on a file changed here).

### Local sanity gates (what Codex can run locally)

```bash
# Always before commit:
python -m py_compile $(git diff --cached --name-only | grep '\.py$')
python -c "import ast; [ast.parse(open(f).read()) for f in '$(git diff --cached --name-only | grep .py)'.split()]"
python scripts/run_validation_suite.py --config configs/config_minimal.yaml  # only works if h5py present — will fail locally, that is OK
```

### CI gates (authoritative)

- `.github/workflows/ci.yml` steps (in order): checkout, setup-python 3.11, install deps, install editable, compileall, pytest `-q --maxfail=5`, minimal config validation, production-smoke `--strict-pass`, upload artifacts.
- **Treat a single local `py_compile` pass as necessary but not sufficient.** The only green light is a full CI run.

### Polling CI from Codex's environment

- Remote: `https://api.github.com/repos/dongfanghong656/MIESIM/actions/runs?per_page=2`.
- Proxy required for the Windows host: `http://127.0.0.1:7897` (Clash Verge). Without it, GitHub API hangs.
- Auth: GitHub PAT; **never embed the PAT in `git remote` — it will leak in logs.** Use `-H "Authorization: Bearer $GH_PAT"` in curl only.

### Work-log protocol

- After each Tx.y commit + CI green: append a dated block to [`status/round5_work_log.md`](round5_work_log.md) (create on first task).
- Block template:

  ```markdown
  ## T5.X — <title> (YYYY-MM-DD, SHA <7char>)

  ### What changed
  - file:line bullets

  ### Verification
  - CI run URL + step outcome
  - Numeric before/after where relevant

  ### Open items carried forward
  - (none | short bullet)
  ```

- This log is the input to the eventual Round 5 review package. Keep it factual; save opinions for the DOCX.

## 4 Deferred to Round 6 (explicit out-of-scope for Round 5)

- Vector Debye / full polarization Jones through the optical train.
- Internal Mie (current path uses calibrated lookup; keep).
- Zemax POP or FDTD spot checks.
- Closed-loop measurement: finite-size scattering ball with dispersion + spherical aberration; round-trip a measured PSF.
- Köhler illumination: aperture stop vs field stop as independent parameters.
- Pupil-angle-resolved dichroic with AOI table.
- Any change to the `hybrid_rci` code path itself. Hybrid is frozen as a legacy/diagnostic model; only full-spectral evolves.

If Codex identifies a Round 6 candidate while doing T5.X work, capture it as a bullet in the T5.X work-log block under "Open items carried forward"; do not expand scope in-flight.

## 5 Handoff checklist for PRO

Round 5 is ready for PRO review when **all** of:

- [ ] Every T5.X has a green CI run referenced in `status/round5_work_log.md`.
- [ ] `validation_summary.json.verdict.pilot_pass == True` for both minimal and production-smoke configs.
- [ ] `direct_model_comparison.csv` on-axis axial relative error still ≤ 1e-3 at decimation=1.
- [ ] `checks` now includes `strehl_v_gate`, `sensitivity_rolloff_v_gate`, and the propagated-reflector phase-stability check.
- [ ] `PHYSICS_CONTRACT.md` has new sections for differential dispersion and common-path pupil identity.
- [ ] `review_packages/round4_response_<ts>/` exists with DOCX + evidence artifacts.
- [ ] No new `unmet_requirements` entry introduced.

## 6 Non-negotiables

- Do not install packages on the local Windows host; CI is the sole test environment.
- Do not change the `hybrid_rci` code path except as mandated by T5.5 (pupil identity collapse, algebraically identical).
- Do not touch `outputs/`, `review_packages/`, `tmp/`, `logs/` — they are `.gitignore`d for a reason; committing outputs is a hard rule violation.
- Do not reopen any Round 4 decision (default `direct_psf_model`, decimation floor, N1 split, on-axis gate wiring). If a Round 5 task seems to need that, stop and escalate before coding.

— End of Round 5 task book —
