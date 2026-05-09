# Round 5 Codex Alignment Review

Date: 2026-05-10
Repo: `C:\codex-data\MIESIM`
Task source: `status/round5_task_book_for_codex.md`

## 1. Canonical project

The active Round 5 project is `C:\codex-data\MIESIM`, not the older
`C:\codex-data\OCT_Research_System\common-path-oct-sim` mirror.

Evidence:

- `C:\codex-data\MIESIM` has the live `.git` directory.
- `status/round5_task_book_for_codex.md` names target repo
  `dongfanghong656/MIESIM` branch `main`.
- The same repo contains the current `PHYSICS_CONTRACT.md`,
  `cop_oct_sim/`, `tests/test_sanity.py`, and `status/round4_verification_log.md`.

## 2. Starting state accepted from Claude

Claude's task book is accepted as the Round 5 execution contract. I will not
reopen Round 4 decisions unless a task explicitly requires escalation.

Round 4 hard evidence to preserve:

- Known-good CI run: `25555010839` at commit `c1a4d94`.
- Faithfulness floor:
  - on-axis `axial_fwhm_relative_error = 7.2e-5`
  - `nrmse_3d = 3.85e-6`
- Regression rule: any Round 5 change that pushes those numbers above roughly
  1 percent without physics justification is a regression.

Existing uncommitted local context:

- `status/round4_verification_log.md` already contains Claude's Phase 2 closeout
  block. Treat it as baseline context and do not overwrite it.

## 3. Scope alignment

Round 5 is limited to these four themes:

1. Add remaining V-gates:
   - Strehl ratio
   - differential dispersion
   - sensitivity roll-off
   - propagated-reflector phase stability
2. Document N4: common-path incident/detection pupil identity.
3. Produce the Round 4 response DOCX package.
4. Keep the existing status/test/CI discipline.

Explicitly out of scope for Round 5:

- vector Debye
- Mie kernel
- Zemax
- FDTD
- Kohler illumination
- pupil-angle-resolved dichroic
- broad changes to `hybrid_rci`

## 4. Current implementation status

Local inspection on 2026-05-10:

| Task | Current status | Evidence |
|---|---|---|
| T5.1 Strehl V-gate | Implemented | `strehl_ratio`, `checks["strehl_v_gate"]`, tests, CI run `25597551164` |
| T5.2 differential dispersion | Implemented | `reference_dispersion_quadratic_rad`, `differential_dispersion_phase`, tests, CI run `25603263330` |
| T5.3 sensitivity roll-off V-gate | Implemented | theoretical/measured roll-off gate, tests, CI run `25603613595` |
| T5.4 propagated-reflector phase stability | Implemented | `simulate_oct_raw_direct` -> noise -> `reconstruct_sd_oct`, CI run `25604099700` |
| T5.5 N4 pupil identity | Implemented | common-path pupil identity helper and metadata, CI run `25604462482` |
| T5.6 Round 4 response DOCX | Implemented | `build_review_package.py --round 4`, axial-line PNGs, CI run `25606208085` |

The final work-log commit `451b5d7` also has a green CI run `25606395608`.

## 5. Alignment review findings

### Finding A: Strehl normalization is an intentional contract deviation

Claude's task book specified peak-normalized center-voxel Strehl. The current
direct PSF volumes are already globally peak-normalized, so that literal formula
would make heavy defocus insensitive. The implementation uses integrated
intensity normalization and records `strehl_normalization =
"integrated_intensity"`.

Resolution: accepted as an explicit simulation-harness contract and documented
in `PHYSICS_CONTRACT.md`.

### Finding B: Differential dispersion scope is broader than the task-book shorthand

The task book named `_full_spectral_rci_direct_psf` and
`simulate_oct_psf_direct`. Current code also shares the same differential helper
with `simulate_oct_raw_direct` and `theory_conversion.py`, which keeps raw,
direct, and theory paths consistent.

Resolution: accepted as physically consistent and documented in
`PHYSICS_CONTRACT.md`.

### Finding C: Round-4 response package is lightweight by design

The first T5.6 CI artifact copied full H5 validation outputs and became
oversized. The corrected package copies only lightweight metrics and axial-line
PNGs; full H5 outputs remain in the separate validation artifact.

Resolution: accepted. CI run `25606208085` produced a `168847` byte response
artifact with one DOCX, four PNGs, and manifest/hash files.

## 6. Fresh local review evidence

- `python -c "import ast; ..."` over the touched source/test scripts: passed.
- `python -m pytest tests/test_sanity.py -q`: 52 passed.
- `python scripts/build_review_package.py --round 4` with a fake validation
  directory: produced a valid response package, one DOCX, four PNGs, and a zip
  of approximately 146 KB.

## 7. Remaining scope boundaries

- Do not reopen Round 4 defaults, decimation floor, N1 split, or on-axis gate.
- Do not expand into vector Debye, Mie kernel, Zemax, FDTD, Kohler
  illumination, or pupil-angle-resolved dichroic in Round 5.
- The Round-4 DOCX currently references sibling PNG evidence rather than
  embedding plots directly. This keeps the package compact and avoids adding a
  DOCX rendering dependency to CI.
