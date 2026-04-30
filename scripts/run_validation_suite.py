from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cop_oct_sim.config_schema import load_config
from cop_oct_sim.validation import run_validation_suite

def main() -> int:
    parser = argparse.ArgumentParser(description="Run common-path OCT validation suite.")
    parser.add_argument("--config", type=Path, default=None, help="Optional YAML config file.")
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT / "outputs")
    parser.add_argument("--N", type=int, default=48, help="Pupil grid size for the main validation run.")
    parser.add_argument("--k-samples", type=int, default=160, help="Spectral samples for the main validation run.")
    parser.add_argument("--pad-factor", type=int, default=2, help="Fourier-plane zero-padding factor.")
    parser.add_argument(
        "--strict-pass",
        action="store_true",
        help="Exit non-zero when any physics/contract gate is marked review.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config) if args.config else None
    summary = run_validation_suite(args.output_root, cfg, N=args.N, k_samples=args.k_samples, pad_factor=args.pad_factor)
    print(
        json.dumps(
            {
                "output_dir": summary["output_dir"],
                "all_pass": summary["checks"]["all_pass"],
                "pilot_pass": summary["verdict"]["pilot_pass"],
                "review_required": summary["verdict"]["review_required"],
                "blocker_count": summary["verdict"]["blocker_count"],
            },
            ensure_ascii=False,
        )
    )
    return 1 if args.strict_pass and not summary["verdict"]["pilot_pass"] else 0

if __name__ == "__main__":
    raise SystemExit(main())
