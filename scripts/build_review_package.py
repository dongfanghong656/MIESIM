from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from xml.sax.saxutils import escape

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.plot_round4_axial_lines import generate_round4_axial_line_plots

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

def _docx_paragraph(text: str, style: str | None = None) -> str:
    style_xml = f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>' if style else ""
    preserve = ' xml:space="preserve"' if text.startswith(" ") or text.endswith(" ") else ""
    return f"<w:p>{style_xml}<w:r><w:t{preserve}>{escape(text)}</w:t></w:r></w:p>"

def _docx_table(rows: list[list[str]]) -> str:
    row_xml = []
    for row in rows:
        cells = "".join(
            "<w:tc><w:tcPr><w:tcW w:w=\"3000\" w:type=\"dxa\"/></w:tcPr>"
            f"{_docx_paragraph(cell)}</w:tc>"
            for cell in row
        )
        row_xml.append(f"<w:tr>{cells}</w:tr>")
    return (
        "<w:tbl><w:tblPr><w:tblW w:w=\"0\" w:type=\"auto\"/>"
        "<w:tblBorders>"
        "<w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>"
        "<w:left w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>"
        "<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>"
        "<w:right w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>"
        "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>"
        "<w:insideV w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"auto\"/>"
        "</w:tblBorders></w:tblPr>"
        + "".join(row_xml)
        + "</w:tbl>"
    )

def _write_minimal_docx(path: Path, body_xml: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    document_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        f"<w:body>{body_xml}<w:sectPr><w:pgSz w:w=\"12240\" w:h=\"15840\"/>"
        "<w:pgMar w:top=\"720\" w:right=\"720\" w:bottom=\"720\" w:left=\"720\" "
        "w:header=\"360\" w:footer=\"360\" w:gutter=\"0\"/></w:sectPr></w:body></w:document>"
    )
    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        '<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/>'
        '<w:rPr><w:sz w:val="21"/></w:rPr></w:style>'
        '<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/>'
        '<w:rPr><w:b/><w:sz w:val="32"/></w:rPr></w:style>'
        '<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/>'
        '<w:rPr><w:b/><w:sz w:val="26"/></w:rPr></w:style>'
        '<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/>'
        '<w:rPr><w:b/><w:sz w:val="23"/></w:rPr></w:style>'
        "</w:styles>"
    )
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
        '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>'
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
        '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
        "</Types>"
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
        "</Relationships>"
    )
    document_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"/>'
    )
    core = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        "<dc:title>Round-4 Response Package</dc:title>"
        "<dc:creator>Codex</dc:creator>"
        "</cp:coreProperties>"
    )
    app = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        "<Application>Codex OOXML builder</Application></Properties>"
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("word/_rels/document.xml.rels", document_rels)
        zf.writestr("word/document.xml", document_xml)
        zf.writestr("word/styles.xml", styles_xml)
        zf.writestr("docProps/core.xml", core)
        zf.writestr("docProps/app.xml", app)

def _read_first_direct_model_row(validation_dir: Path) -> dict[str, str]:
    path = validation_dir / "metrics" / "direct_model_comparison.csv"
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return rows[0] if rows else {}

def _write_round4_response_docx(path: Path, *, validation_dir: Path, summary: dict, plot_manifest: dict) -> None:
    row = _read_first_direct_model_row(validation_dir)
    checks = summary.get("checks", {})
    verdict = summary.get("verdict", {})
    run_url = "https://github.com/dongfanghong656/MIESIM/actions"
    changelog_rows = [
        ["PRO Round 3 P0", "Round 4/Round 5 response", "Evidence"],
        ["Defocus unit conflation", "Physical-z and Zernike-OPD paths separated", "PHYSICS_CONTRACT.md; oct_forward.py; propagation.py"],
        ["Hybrid as truth", "full_spectral_rci remains default; hybrid is opt-in", "config_schema.py; validation_summary.json"],
        ["Pass-gate semantics", "all_pass/pilot_pass/review_required gates wired to blockers", "validation.py; validation_summary.json"],
        ["Config relative paths", "config paths resolve from config directory and project root", "config_schema.py"],
        ["Central-lobe/on-axis FWHM", "metrics and validation gates report central/on-axis behavior", "metrics.py; validation.py"],
        ["Round 5 V-gates", "Strehl, differential dispersion, rolloff, and propagated-reflector phase stability added", "validation_summary.json checks"],
    ]
    plot_names = [
        Path(item["path"]).name
        for item in plot_manifest.get("plots", [])
        if isinstance(item, dict) and "path" in item
    ]
    paragraphs = [
        _docx_paragraph("Round-4 Response for PRO Sign-off", "Title"),
        _docx_paragraph(f"Date: {datetime.now().strftime('%Y-%m-%d')}"),
        _docx_paragraph(f"Version: { _git_value(['rev-parse', '--short', 'HEAD']) }"),
        _docx_paragraph("Scope: scalar low-NA common-path microscope/OCT validation harness. This is not a final vector-Debye, Mie, Zemax, or FDTD truth simulator."),
        _docx_paragraph("1. Changelog", "Heading1"),
        _docx_table(changelog_rows),
        _docx_paragraph("2. P0-by-P0 response", "Heading1"),
        _docx_paragraph("P0-1 defocus: physical sample displacement now propagates through the low-NA physical defocus phase, while objective.defocus_um remains a Zernike OPD coefficient."),
        _docx_paragraph("P0-2 hybrid truth: the default direct model is full_spectral_rci; hybrid_rci is retained as an explicit diagnostic surrogate and is guarded by faithfulness tests."),
        _docx_paragraph("P0-3 pass gates: pilot_pass is now tied to all_pass and blocker_count, with unmet_requirements reported instead of hidden."),
        _docx_paragraph("P0-4 config paths: calibration and config-relative paths are resolved deterministically for clean checkout reproduction."),
        _docx_paragraph("P0-5 FWHM metrics: central-lobe, on-axis, and integrated axial metrics are reported to avoid mirror-peak or off-axis masking."),
        _docx_paragraph("Round 5 closeout: Strehl, differential dispersion, sensitivity rolloff, propagated-reflector phase stability, and N4 common-path pupil identity are now explicit contracts."),
        _docx_paragraph("3. Axial-line plots", "Heading1"),
        _docx_paragraph("Generated on-axis axial magnitude plots are included as sibling PNG evidence in evidence/round4_axial_lines/."),
        *[_docx_paragraph(f"- {name}") for name in plot_names],
        _docx_paragraph("4. Faithfulness number", "Heading1"),
        _docx_paragraph(
            "CI artifact direct_model_comparison.csv reports "
            f"axial_fwhm_relative_error = {row.get('axial_fwhm_relative_error', 'unknown')}, "
            f"nrmse_3d = {row.get('nrmse_3d', 'unknown')}, "
            f"full_spectral_rci_interpolated = {row.get('full_spectral_rci_interpolated', 'unknown')}."
        ),
        _docx_paragraph(f"Validation verdict: pilot_pass = {verdict.get('pilot_pass', 'unknown')}; blocker_count = {verdict.get('blocker_count', 'unknown')}."),
        _docx_paragraph("5. T1.3 finding", "Heading1"),
        _docx_paragraph("The scalar low-NA hybrid_rci path is treated as a faithful diagnostic surrogate only within the pinned low-NA regime. It is not the default truth model and not a claim about high-NA vector or particle-scattering truth."),
        _docx_paragraph("6. V-gate evidence", "Heading1"),
        _docx_paragraph(f"strehl_v_gate pass = {checks.get('strehl_v_gate', {}).get('pass', 'unknown')}; sensitivity_rolloff_v_gate pass = {checks.get('sensitivity_rolloff_v_gate', {}).get('pass', 'unknown')}; phase_stability_v_gate pass = {checks.get('phase_stability_v_gate', {}).get('pass', 'unknown')}."),
        _docx_paragraph(f"CI evidence root: {run_url}. The response package also includes copied validation_summary.json, direct_model_comparison.csv, and plot manifest files."),
    ]
    _write_minimal_docx(path, "".join(paragraphs))

def build_round4_response_package(
    *,
    validation_dir: Path | None,
    logs_dir: Path | None,
    output_root: Path,
    config_path: Path,
) -> tuple[Path, Path]:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    package_dir = output_root / f"round4_response_{stamp}"
    if package_dir.exists():
        shutil.rmtree(package_dir)
    package_dir.mkdir(parents=True)

    validation_dir = (validation_dir or _latest_validation_dir(PROJECT_ROOT / "outputs")).resolve()
    logs_dir = logs_dir.resolve() if logs_dir else PROJECT_ROOT / "logs"
    evidence_dir = package_dir / "evidence"
    evidence_dir.mkdir()

    validation_target = evidence_dir / "validation" / validation_dir.name
    _copytree(validation_dir, validation_target)
    for relative in [
        Path("metrics") / "validation_summary.json",
        Path("metrics") / "direct_model_comparison.csv",
        Path("metrics") / "direct_model_axial_profiles.csv",
    ]:
        _copy_if_exists(validation_dir / relative, evidence_dir / relative)

    plots_dir = evidence_dir / "round4_axial_lines"
    plot_manifest = generate_round4_axial_line_plots(
        config_path=config_path,
        output_root=plots_dir,
        N=24,
        k_samples=64,
        pad_factor=1,
    )

    logs_target = evidence_dir / "logs"
    logs_target.mkdir()
    for name in ["pytest.log", "run_validation.log"]:
        _copy_if_exists(logs_dir / name, logs_target / name)

    summary = _read_summary(validation_dir)
    _write_round4_response_docx(package_dir / "response.docx", validation_dir=validation_dir, summary=summary, plot_manifest=plot_manifest)
    meta = {
        "package_name": package_dir.name,
        "created_at": stamp,
        "round": 4,
        "git_commit": _git_value(["rev-parse", "HEAD"]),
        "git_status_short": _project_git_status_short(),
        "validation_dir": str(validation_dir),
        "response_docx": str(package_dir / "response.docx"),
        "plot_manifest": str(plots_dir / "round4_axial_lines_manifest.json"),
    }
    (package_dir / "package_manifest.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    _write_sha256sums(package_dir)
    zip_path = Path(shutil.make_archive(str(package_dir), "zip", root_dir=output_root, base_dir=package_dir.name))
    return package_dir, zip_path

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
    parser = argparse.ArgumentParser(description="Build response review package.")
    parser.add_argument("--round", type=int, choices=[3, 4], default=3)
    parser.add_argument("--validation-dir", type=Path, default=None)
    parser.add_argument("--logs-dir", type=Path, default=PROJECT_ROOT / "logs")
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT / "review_packages")
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "configs" / "config_minimal.yaml")
    args = parser.parse_args()

    if args.round == 4:
        package_dir, zip_path = build_round4_response_package(
            validation_dir=args.validation_dir,
            logs_dir=args.logs_dir,
            output_root=args.output_root,
            config_path=args.config,
        )
    else:
        package_dir, zip_path = build_package(args.validation_dir, args.logs_dir, args.output_root)
    print(json.dumps({"package_dir": str(package_dir), "zip_path": str(zip_path)}, ensure_ascii=False))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
