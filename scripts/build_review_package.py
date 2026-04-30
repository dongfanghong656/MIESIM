from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

EXCLUDE_DIRS = {
    "__pycache__",
    ".pytest_cache",
    ".test-output",
    "outputs",
    "review_packages",
}
EXCLUDE_SUFFIXES = {".pyc", ".pyo"}

def _copytree(src: Path, dst: Path) -> None:
    def ignore(directory: str, names: list[str]) -> set[str]:
        ignored: set[str] = set()
        for name in names:
            path = Path(directory) / name
            if name in EXCLUDE_DIRS:
                ignored.add(name)
            elif path.suffix in EXCLUDE_SUFFIXES:
                ignored.add(name)
        return ignored

    shutil.copytree(src, dst, ignore=ignore)

def _copy_if_exists(src: Path, dst: Path) -> None:
    if not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

def _latest_validation_dir(outputs_root: Path) -> Path:
    candidates = [p for p in outputs_root.glob("validation_*") if p.is_dir()]
    if not candidates:
        raise FileNotFoundError(f"No validation_* directory under {outputs_root}")
    return max(candidates, key=lambda p: p.stat().st_mtime)

def _git_value(args: list[str]) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=PROJECT_ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unavailable"

def _project_git_status_short() -> str:
    try:
        status = subprocess.check_output(
            ["git", "status", "--short", "--", "."],
            cwd=PROJECT_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unavailable"
    return status or "clean"

def _read_summary(validation_dir: Path) -> dict:
    path = validation_dir / "metrics" / "validation_summary.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))

def _write_readme(path: Path, package_name: str, validation_dir: Path, summary: dict, stamp: str) -> None:
    checks = summary.get("checks", {})
    verdict = summary.get("verdict", {})
    lines = [
        f"# {package_name}",
        "",
        "Round3 response review package for the common-path microscope/OCT scalar low-NA simulation baseline.",
        "",
        "## Environment",
        "",
        f"- Python used by package builder: `{sys.version.split()[0]}`",
        "- Install: `python -m pip install -e .[test]`",
        "- Test command: `python -m pytest`",
        "- Validation command: `python scripts/run_validation_suite.py --config configs/config_minimal.yaml --output-root outputs`",
        "- Strict gate command: `python scripts/run_validation_suite.py --config configs/config_minimal.yaml --output-root outputs --strict-pass`",
        "",
        "## Fresh Validation Snapshot",
        "",
        f"- package timestamp: `{stamp}`",
        f"- validation source: `{validation_dir.name}`",
        f"- all_pass: `{checks.get('all_pass', 'unknown')}`",
        f"- pilot_pass: `{verdict.get('pilot_pass', 'unknown')}`",
        f"- review_required: `{verdict.get('review_required', 'unknown')}`",
        f"- blocker_count: `{verdict.get('blocker_count', 'unknown')}`",
        f"- Path E measured predictive pass: `{checks.get('path_e_measured_predictive', {}).get('pass', 'unknown')}`",
        f"- Path E measured predictive NRMSE: `{checks.get('path_e_measured_predictive', {}).get('nrmse_3d', 'unknown')}`",
        f"- Path E measured predictive axial FWHM relative error: `{checks.get('path_e_measured_predictive', {}).get('axial_fwhm_relative_error', 'unknown')}`",
        f"- Path E ideal upper-bound pass: `{checks.get('path_e_ideal_upper_bound', {}).get('pass', 'unknown')}`",
        f"- direct PSF model: `{summary.get('runtime_parameters', {}).get('direct_psf_model', 'unknown')}`",
        f"- full spectral RCI depth decimation: `{summary.get('runtime_parameters', {}).get('full_spectral_rci_depth_decimation', 'unknown')}`",
        f"- direct PSF model comparison pass: `{checks.get('direct_model_comparison', {}).get('pass', 'unknown')}`",
        f"- direct PSF model comparison max NRMSE: `{checks.get('direct_model_comparison', {}).get('max_nrmse_3d', 'unknown')}`",
        f"- direct model axial profile max on-axis residual: `{checks.get('direct_model_axial_profiles', {}).get('max_on_axis_abs_residual', 'unknown')}`",
        f"- negative controls pass: `{checks.get('negative_controls', {}).get('pass', 'unknown')}`",
        f"- convergence sweeps pass: `{checks.get('convergence_sweeps', {}).get('pass', 'unknown')}`",
        "",
        "## Contents",
        "",
        "- `project/`: source code, tests, scripts, docs, pyproject, requirements",
        "- `configs/`: review configs requested for minimal, visible 660/730 nm, OCT 850 nm, and full-spectral-RCI smoke / production-smoke cases",
        "- `calibration_data/`: placeholder CSV/YAML calibration inputs for source spectrum, pixel-to-lambda, dichroic, camera QE, objective notes",
        "- `validation/`: fresh validation data directory copied from the local run, including `oct_direct/oct_psf_direct_lowres_full_spectral_rci.h5`, `metrics/direct_model_comparison.csv`, and `metrics/direct_model_axial_profiles.csv`",
        "- `logs/`: fresh pytest and validation command logs",
        "- `SHA256SUMS.txt`: package file hashes",
        "",
        "## Claim Boundary",
        "",
        "- Supported: scalar low-NA engineering baseline, no per-plane through-focus normalization, spatial RCI sampling in raw direct, finite sphere projection/volume kernels, path-separated dichroic, measured predictive / ideal upper-bound Path E split, axial/lateral FWHM pass gates, negative controls, convergence sweeps.",
        "- Newly exposed direct option: `oct.direct_psf_model=full_spectral_rci` uses low-resolution h_RCI(x,y,z;k) before k-space FFT reconstruction as the main direct PSF path.",
        "- Scaling control: `oct.full_spectral_rci_depth_decimation` records which depth planes were computed directly and whether interpolation was used.",
        "- Diagnostic evidence: `metrics/direct_model_comparison.csv` compares hybrid RCI and full-spectral RCI at low resolution; it is a model-difference diagnostic, not a proof of production truth.",
        "- Direct-truth evidence: `oct_direct/oct_psf_direct_lowres_full_spectral_rci.h5` now includes sampled complex `spectral_rci_sampled_h_rci[depth,k,y,x]`, k grids, source/window/dispersion vectors, fixed scatterer-depth phase origin, and axial profiles; `metrics/direct_model_axial_profiles.csv` exports the profile comparison in table form.",
        "- Gate semantics: `pilot_pass=True` requires `all_pass=True` and `blocker_count=0`; `review_required=True` cannot coexist with a pilot-pass claim.",
        "- Defocus contract: physical z displacement is propagated with low-NA physical defocus phase; Zernike defocus OPD is reserved for system wavefront aberration.",
        "- Not supported as final truth: full vector Debye, Zemax POP, internal Mie solver, production-resolution full spectral truth, FDTD/COMSOL scatterer truth.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def _write_changelog(path: Path) -> None:
    lines = [
        "# CHANGELOG_ROUND3_RESPONSE",
        "",
        "## Fixed Earlier P0 Findings",
        "",
        "- Through-focus propagation no longer peak-normalizes each z plane by default.",
        "- OCT raw direct samples spatial RCI response at scatterer position instead of using a scalar mean peak.",
        "- Finite bead/sphere handling now includes uniform sphere projection/volume kernels and explicit scatterer config.",
        "",
        "## Round3 Response Changes",
        "",
        "- Path E split into `path_e_measured_predictive` and `path_e_ideal_upper_bound` outputs.",
        "- Pass gate now requires axial and lateral FWHM agreement; low 3D NRMSE alone is not enough.",
        "- Validation preserves visible microscope config values instead of silently forcing the microscope wavelength to OCT center wavelength.",
        "- Added low-resolution full spectral RCI direct PSF artifact before k-space FFT reconstruction.",
        "- Added `oct.direct_psf_model=full_spectral_rci` so the full-spectral-RCI path can be promoted to the main direct truth in smoke runs.",
        "- Added `oct.full_spectral_rci_depth_decimation` plus computed-depth metadata for scalable full-spectral-RCI validation.",
        "- Added `metrics/direct_model_comparison.csv` to report low-resolution hybrid RCI versus full-spectral RCI model differences.",
        "- Added sampled full-spectral `h_RCI[depth,k,y,x]`, k/source/window/dispersion/fixed-depth-phase metadata, and `metrics/direct_model_axial_profiles.csv` for direct-truth review.",
        "- Corrected full-spectral RCI depth phase to use fixed scatterer-depth phase; output depth now comes from FFT sampling plus `h_RCI(z;k)` through-focus variation instead of re-centering the axial gate at every depth plane.",
        "- Split physical defocus phase from Zernike OPD defocus and added regression tests for small physical z behavior.",
        "- Tightened gate semantics: direct model mismatch, convergence review rows, and unsupported blockers prevent pilot pass.",
        "- Resolved config-relative calibration paths to absolute paths at load time.",
        "- Added central-lobe and on-axis axial FWHM metrics.",
        "- Illumination NA / aperture / pupil-fill changes now alter epi scatterer response shape instead of only scaling amplitude.",
        "- Path-separated dichroic transmission/reflection coefficients and CSV table loader.",
        "- Source spectrum and pixel-to-lambda CSV loaders.",
        "- Microscope illumination pupil/PSF/envelope output.",
        "- Camera QE, flat field, dark frame, saturation, shot/read noise chain.",
        "- Negative controls for chromatic mismatch, wrong dichroic path, and bad k calibration.",
        "- Convergence sweeps for N, pad factor, k samples, z step, window, bead size, and camera pixel.",
        "- Automated figures including XY/XZ/YZ, A-scan, residual, and sweep trend.",
        "",
        "## Remaining Boundaries",
        "",
        "- Scalar low-NA only.",
        "- Low-resolution full spectral RCI is a review artifact and is not yet the production direct-truth path.",
        "- `rayleigh_mie_lookup` uses supplied lookup tables; it is not an internal Mie solver.",
        "- Vector Debye, Zemax POP, FDTD/COMSOL spot checks remain future work.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def _write_sha256sums(root: Path) -> None:
    entries = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p.name != "SHA256SUMS.txt"):
        rel = path.relative_to(root).as_posix()
        entries.append(f"{_sha256_file(path)}  {rel}")
    (root / "SHA256SUMS.txt").write_text("\n".join(entries) + "\n", encoding="utf-8")

def build_package(validation_dir: Path | None, logs_dir: Path | None, output_root: Path) -> tuple[Path, Path]:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    package_name = f"common_path_oct_sim_round3_response_{stamp}"
    package_dir = output_root / package_name
    if package_dir.exists():
        shutil.rmtree(package_dir)
    package_dir.mkdir(parents=True)

    validation_dir = validation_dir or _latest_validation_dir(PROJECT_ROOT / "outputs")
    validation_dir = validation_dir.resolve()
    logs_dir = logs_dir.resolve() if logs_dir else PROJECT_ROOT / "logs"

    project_dir = package_dir / "project"
    project_dir.mkdir()
    for name in [
        "pyproject.toml",
        "requirements.txt",
        "README.md",
        "AGENTS.md",
        "PHYSICS_CONTRACT.md",
    ]:
        _copy_if_exists(PROJECT_ROOT / name, project_dir / name)
    for name in ["cop_oct_sim", "tests", "scripts", "docs"]:
        src = PROJECT_ROOT / name
        if src.exists():
            _copytree(src, project_dir / name)

    for name in ["configs", "calibration_data"]:
        src = PROJECT_ROOT / name
        if src.exists():
            _copytree(src, package_dir / name)

    validation_target = package_dir / "validation" / validation_dir.name
    _copytree(validation_dir, validation_target)

    logs_target = package_dir / "logs"
    logs_target.mkdir()
    for name in ["pytest.log", "run_validation.log"]:
        _copy_if_exists(logs_dir / name, logs_target / name)

    summary = _read_summary(validation_dir)
    meta = {
        "package_name": package_name,
        "created_at": stamp,
        "source_project": str(PROJECT_ROOT),
        "validation_dir": str(validation_dir),
        "git_commit": _git_value(["rev-parse", "HEAD"]),
        "git_status_short": _project_git_status_short(),
    }
    (package_dir / "package_manifest.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    _write_readme(package_dir / "README_REVIEW.md", package_name, validation_dir, summary, stamp)
    _write_changelog(package_dir / "CHANGELOG_ROUND3_RESPONSE.md")
    _write_sha256sums(package_dir)

    zip_base = output_root / package_name
    zip_path = Path(shutil.make_archive(str(zip_base), "zip", root_dir=output_root, base_dir=package_name))
    return package_dir, zip_path

def main() -> int:
    parser = argparse.ArgumentParser(description="Build Round3 response review package.")
    parser.add_argument("--validation-dir", type=Path, default=None)
    parser.add_argument("--logs-dir", type=Path, default=PROJECT_ROOT / "logs")
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT / "review_packages")
    args = parser.parse_args()

    package_dir, zip_path = build_package(args.validation_dir, args.logs_dir, args.output_root)
    print(json.dumps({"package_dir": str(package_dir), "zip_path": str(zip_path)}, ensure_ascii=False))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
