# Round 4 Phase 0 Verification Log

Date: 2026-04-30
Verifier: Claude (independent of Codex)
Target package: `review_packages/common_path_oct_sim_round3_response_20260428-222556`
Reference audit: PRO Round 3 Response audit on `095635` (PRO has NOT seen 222556).

## V0.1 — Defocus split (physical z vs Zernike OPD)

**Status: PASS**

Evidence:
- `cop_oct_sim/propagation.py:30-38` defines `physical_defocus_phase(grid, z_um, wavelength_um, NA)` with formula `phase = π · z_um · NA² · ρ² / wavelength_um`, exactly matching `PHYSICS_CONTRACT.md` line 18-20.
- `cop_oct_sim/pupil.py:9-13` `zernike_defocus(grid, coeff_um, wavelength_um)` is the independent Zernike OPD coefficient path.
- `cop_oct_sim/oct_forward.py:32-47` `_make_oct_pupil` correctly takes `physical_defocus_um` (routed through `physical_defocus_phase`) and `zernike_defocus_opd_um` (routed through `build_shared_pupil(defocus_um=...)`) as two separate parameters.
- All through-focus call sites in `oct_forward.py` (lines 146-147, 224-225, 298-299) pass physical depth as `physical_defocus_um`, not as a Zernike coefficient.

PRO's P0 #1 (defocus unit conflation) is genuinely fixed.

## V0.2 — Gate semantics (pilot_pass / review_required / all_pass coherence)

**Status: PASS, but overcorrected.**

Evidence: `cop_oct_sim/validation.py:842-847`
```python
verdict = {
    "pilot_pass": bool(checks["all_pass"] and blocker_count == 0),
    "review_required": bool(blocker_count > 0),
    "blocker_count": blocker_count,
    "unsupported_claims": unsupported_claims,
}
```
The contradictory `(all_pass=True, review_required=True, pilot_pass=True)` triple is now logically impossible.

**However**, `blocker_count` is computed at line 838-841 as `failing_checks + len(unsupported_claims)`, where `unsupported_claims` is a hardcoded 4-line list at `validation.py:832-836`. Therefore `blocker_count >= 4` always, which forces `pilot_pass = False` always. This is honest but engineering-impractical: no canonical run can ever be "pilot pass" until vector Debye / internal Mie / production-resolution full spectral RCI all ship.

**Recommendation (round-4 small fix):** Split `unsupported_claims` into `documented_limitations` (informational; not blocker) and `unmet_requirements` (counted as blocker). Documented scalar low-NA boundary conditions belong in the former.

## V0.3 — direct_model_comparison gates on FWHM

**Status: PASS.**

Evidence: `cop_oct_sim/validation.py:388-396`
```python
pass_gate = (
    finite
    and computed_count > 0
    and depth_count > 0
    and axial_rel_error <= axial_rel_threshold     # 0.25
    and lateral_rel_error <= lateral_rel_threshold # 0.25
)
```
NRMSE not even in the gate; axial-FWHM relative error is the primary pass criterion. Matches PRO 06 requirement "axial FWHM 相差 >25% 必须 fail".

**Note:** 0.25 is PRO's *boundary*, not target. PRO 01 implies the goal is to make hybrid and full-spectral converge much closer. After T1.1 (production-resolution full spectral truth), tighten to 0.10.

## V0.4 — Config relative-path resolution

**Status: PASS.**

Evidence: `cop_oct_sim/config_schema.py:119-134`
```python
def _resolve_optional_path(value, config_path):
    ...
    candidates = [
        config_dir / path,
        config_dir.parent / path,
        Path.cwd().resolve() / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())
    return str((config_dir / path).resolve())
```
Falls back to config_dir-relative even if no candidate exists. Plus `resolve_config_paths` (line 136-158) eagerly resolves all known optional file fields to absolute paths at load time.

PRO's P0 #4 (config_visible_660 FileNotFoundError when run from project subdir) is fixed.

## V0.5 — Central-lobe / on-axis FWHM metrics

**Status: HALF-FIXED.**

Evidence:
- `cop_oct_sim/metrics.py:34-58` defines `central_lobe_fwhm_1d` (peak-anchored, walks down only while monotonically decreasing).
- `cop_oct_sim/metrics.py:72-80` defines `on_axis_axial_fwhm` (extracts the geometric center line, then applies `central_lobe_fwhm_1d`).
- `cop_oct_sim/metrics.py:82-88` defines `integrated_axial_fwhm`.
- `compare_psf_volumes` (line 116-130) reports all three flavours.

**Gap:** the path-E pass gate at `validation.py:_path_e_gate` (line 52-79) uses `fwhm_ax_direct/converted` (max-over-XY), not `fwhm_ax_on_axis_*`. Same for the direct_model_comparison gate (line 386-396). PRO 06 explicitly requires the axial gate to be on-axis.

**Recommendation:** in `_path_e_gate` and `build_direct_model_comparison_rows`, replace `cmp["fwhm_ax_direct_um"]/cmp["fwhm_ax_converted_um"]` with `cmp["fwhm_ax_on_axis_direct_um"]/cmp["fwhm_ax_on_axis_converted_um"]`. One-liner each.

## V0.6 — pytest run

**Status: ENV ISSUE.**

System python (`C:\Users\1\AppData\Local\Programs\Python\Python313`) has no pytest module. Need either `pip install pytest` in this interpreter or activate the project's venv. Defer to next phase; not a code defect.

## Summary

| PRO P0 | Status |
|---|---|
| 1. Defocus unit | ✅ Fixed (V0.1) |
| 2. Hybrid as truth | Pending T1.1 (full_spectral_rci promotion still flagged "review artifact") |
| 3. Pass gate semantics | ✅ Fixed but overcorrected (V0.2) — see recommendation |
| 4. Config paths | ✅ Fixed (V0.4) |
| 5. Central-lobe / on-axis FWHM | ⚠️ Functions exist, gates not yet using them (V0.5) |

## Two newly-found issues to fold into round 4 work

- **N1**: `unsupported_claims` hardcoded as blockers makes `pilot_pass` impossible by construction. Refactor `unsupported_claims` into `documented_limitations` + `unmet_requirements` (only the latter blocks). `validation.py:832-841`.
- **N2**: Path-E and direct-model gates still use max-over-XY axial FWHM. Switch to `fwhm_ax_on_axis_*`. `validation.py:_path_e_gate` and `build_direct_model_comparison_rows`.

## Next actions

1. Install pytest in project venv; rerun full pytest + `scripts/run_validation_suite.py --config configs/config_minimal.yaml --strict-pass` and `--config configs/config_full_spectral_rci_production_smoke.yaml --strict-pass`.
2. Land N1 + N2 patches before starting Phase 1 (full_spectral_rci promotion).
3. After patches, re-run validation suite; if minimal verdict says `review_required=True` only because of `documented_limitations`, that is acceptable progress.

## Phase 0 closeout — N1 + N2 landed (2026-04-30)

Verifier: Claude. Local env: system Python 3.13 (no pytest / no h5py — corp PyPI proxy ConnectionReset 10054). Per user policy, pytest + full_spectral smoke run by GitHub CI, not locally.

### N2 — gates switched to on-axis axial FWHM

- `cop_oct_sim/validation.py:59-61` — `_path_e_gate` `axial_error` fed from `metrics["fwhm_ax_on_axis_direct_um"]` / `fwhm_ax_on_axis_converted_um` (was max-over-XY pair).
- `cop_oct_sim/validation.py:390-392` — `build_direct_model_comparison_rows` `axial_rel_error` likewise switched. CSV still publishes both max-over-XY and on-axis fields (lines 410-413) for traceability.

### N1 — unsupported_claims split

- `cop_oct_sim/validation.py:836-856` — `unsupported_claims` (4 hardcoded entries, all blockers) replaced by:
  - `documented_limitations` (3: scalar low-NA, Mie lookup, Path E acceptance rule) — informational only.
  - `unmet_requirements` (1: `lowres_full_spectral_rci_direct` is a review artifact) — counted as blocker.
  - `blocker_count = failing_checks + len(unmet_requirements)`.
  - `verdict` exposes new fields plus legacy `unsupported_claims = documented + unmet` for downstream readers.

### Runtime check — minimal config

- `python scripts/run_validation_suite.py --config configs/config_minimal.yaml` → `outputs/validation_minimal/20260430-152213/`.
- `verdict.blocker_count = 4 = failing_checks(3) + unmet_requirements(1)` ✓ N1 split verified at runtime.
- `documented_limitations` len 3, `unmet_requirements` len 1, `pilot_pass=False`, `review_required=True`.
- `direct_model_comparison.csv`:
  - max-over-XY: full_spectral 7.372 µm, hybrid 2.541 µm.
  - on-axis: full_spectral 6.928 µm, hybrid 2.502 µm.
  - `axial_fwhm_relative_error = 0.6388 = |6.928 − 2.502| / 6.928` — confirms gate now consumes on-axis (max-over-XY would yield 0.6555).
  - Threshold 0.25 → still **fail**, consistent with PRO 01: hybrid_rci is not a faithful surrogate for full_spectral_rci.

### Real failing checks now visible

Previously masked by 4 hardcoded "unsupported claims" inflating `blocker_count` to 8. After N1 split, 3 actual failing leaves surface, to address in Phase 1+2:

1. `path_e_measured_predictive` — converter vs direct truth FWHM mismatch.
2. `direct_model_comparison` — hybrid vs full_spectral on-axis axial FWHM err 63.88% (gate 25%).
3. Third leaf — see `validation_summary.json` `checks` for exact name.

Maps directly to PRO Round 3 reports 01 (physics) and 06 (next package criteria). No new physics regression from N1+N2.

### Environment notes

- AST + py_compile of `validation.py` post-edit: PASS.
- pytest + `config_full_spectral_rci_production_smoke.yaml --strict-pass` deferred to GitHub CI (Codex's standard flow — local PyPI proxy is hostile).

## Phase 0 status

| PRO P0 | Status |
|---|---|
| 1. Defocus unit | ✅ Fixed (V0.1) |
| 2. Hybrid as truth | Pending Phase 1 T1.1 |
| 3. Pass gate semantics | ✅ Fixed + N1 corrected the overcorrection |
| 4. Config paths | ✅ Fixed (V0.4) |
| 5. Central-lobe / on-axis FWHM | ✅ Functions exist (V0.5) + N2 wired into gates |

Phase 0 closed. Ready for Phase 1 T1.1 (promote `full_spectral_rci` decimation=1 to production direct truth, demote `hybrid_rci` to diagnostic).

## Phase 1 T1.1 — full_spectral_rci decimation=1 unlocked (2026-04-30)

### Changes (all in `cop_oct_sim/validation.py`)

- **L356** `build_direct_model_comparison_rows`: `depth_decimation = max(int(config.oct.full_spectral_rci_depth_decimation), 3)` → `max(..., 1)`. Removed the hardcoded floor that was overriding any config below 3. **Root cause of the "lowres review artifact" P0 finding:** even when `configs/config_full_spectral_rci_production_smoke.yaml` set `full_spectral_rci_depth_decimation: 1`, this floor silently coerced it back to 3, then `_interpolate_complex_stack` ran phase-mixing linear interp on real/imag separately. PRO 07 risk #1 ("wrong direct truth validates theory by coincidence") was directly enabled by this single line.
- **L841-845** `unmet_requirements` now conditional: only added when `int(config.oct.full_spectral_rci_depth_decimation) > 1`. Wording updated to point at the exact mechanism (non-phase-aware complex linear interpolation across z) and the fix (set decimation=1).

### Local verification

- `configs/config_minimal.yaml` (decimation=3, unchanged): `verdict.blocker_count=4`, `failing=3`, `unmet=1` — unchanged outcome, only the message text changed.
- `configs/config_full_spectral_rci_production_smoke.yaml` (decimation=1): cannot run locally (`io_hdf5` requires `h5py`, not installed; PyPI proxy blocked). **Deferred to GitHub CI.**
  - Expected behavior under T1.1 unlock: `verdict.unmet_requirements` becomes empty; `direct_model_comparison.full_spectral_rci_interpolated == False`; `axial_fwhm_relative_error` recomputed against production-resolution full_spectral truth.

### Newly-found N4 (NOT in T1.1 scope; flag for PRO + future round)

- `cop_oct_sim/oct_forward.py:224-225` inside `_full_spectral_rci_direct_psf` per-k loop:
  ```python
  pupil_i, _ = _make_oct_pupil(config, grid, wl_nm, physical_defocus_um=defocus_um)
  pupil_d, _ = _make_oct_pupil(config, grid, wl_nm, physical_defocus_um=defocus_um)
  ```
  Both incident and detection pupils are constructed with **identical** parameters. Subsequent `h_rci = ui * conj(ud)` therefore reduces to `|u|^2`, not a genuine incident×detection RCI kernel. For the scalar low-NA common-path geometry where both arms truly share the same pupil this may be intentional, but it should be documented in `PHYSICS_CONTRACT.md` as the explicit modeling assumption rather than left implicit. PRO Round 3 audit did not flag this; suggest folding into Round 5 questions.

### Phase 1 status

| Task | Status |
|---|---|
| T1.1 Remove low-resolution caveats | ✅ |
| T1.2 One direct truth API (hybrid as diagnostic) | Partial — hybrid is still in `direct_model_comparison` as the "converted" leg, but with decimation=1 the gate now compares against true full-spectral. Explicit demotion of `hybrid_rci` from `simulate_oct_psf_direct` API surface is **not** done in T1.1; defer to T1.2 proper. |
| T1.3 Negative test (low-NA aberration → hybrid vs full disagree >25%) | Pending GitHub CI smoke + dedicated test file. |

### Next actions

1. Push to GitHub; let CI run pytest + `config_full_spectral_rci_production_smoke.yaml --strict-pass`. Capture verdict JSON for `unmet_requirements` empty + `axial_fwhm_relative_error` value at decimation=1.
2. After CI smoke result is in: decide whether `axial_fwhm_relative_error` at decimation=1 is small enough (\u003c0.25) to flip `direct_model_comparison.pass=True`; if not, that becomes the next physics task.
3. Add T1.2 proper: rename/move hybrid out of the truth path, keep it as a diagnostic-only column in `direct_model_comparison.csv`.
4. Open Round 5 ticket for N4 (incident/detection pupil identity) — requires PRO consultation.
