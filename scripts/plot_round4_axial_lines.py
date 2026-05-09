from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cop_oct_sim.config_schema import load_config
from cop_oct_sim.oct_forward import simulate_oct_psf_direct

ROUND4_SCATTERER_Z_UM = (0.0, 5.0, 10.0, 20.0)


def _safe_depth_label(z_um: float) -> str:
    if float(z_um).is_integer():
        return f"{int(z_um)}um"
    return f"{z_um:g}um".replace(".", "p")


def generate_round4_axial_line_plots(
    *,
    config_path: Path,
    output_root: Path,
    N: int = 24,
    k_samples: int = 64,
    pad_factor: int = 1,
) -> dict:
    config = load_config(config_path)
    output_root.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {
        "config_path": str(config_path),
        "output_root": str(output_root),
        "scatterer_z_um": list(ROUND4_SCATTERER_Z_UM),
        "N": int(N),
        "k_samples": int(k_samples),
        "pad_factor": int(pad_factor),
        "plots": [],
    }

    full_spectral = config.oct.direct_psf_model == "full_spectral_rci"
    for z_um in ROUND4_SCATTERER_Z_UM:
        direct = simulate_oct_psf_direct(
            config,
            N=N,
            pad_factor=pad_factor,
            k_samples=k_samples,
            scatterer_z_um=float(z_um),
            full_spectral_rci=full_spectral,
        )
        psf = np.asarray(direct["psf"])
        depth_um = np.asarray(direct["depth_um"], dtype=float)
        cy = psf.shape[1] // 2
        cx = psf.shape[2] // 2
        axial = np.abs(psf[:, cy, cx])
        peak = float(np.max(axial))
        if peak > 0:
            axial = axial / peak

        fig, ax = plt.subplots(figsize=(5.0, 3.2), dpi=150)
        ax.plot(depth_um, axial, color="#1f4e79", linewidth=1.8)
        ax.axvline(float(z_um), color="#b03a2e", linestyle="--", linewidth=1.0, label="scatterer z")
        ax.set_title(f"Round-4 on-axis axial line, scatterer z = {z_um:g} um")
        ax.set_xlabel("Depth (um)")
        ax.set_ylabel("Normalized on-axis magnitude")
        ax.set_ylim(bottom=0.0, top=max(1.05, float(np.max(axial)) * 1.05))
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()

        out_path = output_root / f"round4_axial_line_z_{_safe_depth_label(float(z_um))}.png"
        fig.savefig(out_path)
        plt.close(fig)
        manifest["plots"].append(
            {
                "scatterer_z_um": float(z_um),
                "path": str(out_path),
                "depth_count": int(depth_um.size),
                "peak": peak,
            }
        )

    manifest_path = output_root / "round4_axial_lines_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Round-4 axial-line response plots.")
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "configs" / "config_minimal.yaml")
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT / "outputs" / "round4_axial_lines")
    parser.add_argument("--N", type=int, default=24)
    parser.add_argument("--k-samples", type=int, default=64)
    parser.add_argument("--pad-factor", type=int, default=1)
    args = parser.parse_args()

    manifest = generate_round4_axial_line_plots(
        config_path=args.config,
        output_root=args.output_root,
        N=args.N,
        k_samples=args.k_samples,
        pad_factor=args.pad_factor,
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
