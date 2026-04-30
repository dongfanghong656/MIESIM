from __future__ import annotations
from pathlib import Path

import numpy as np

def _normalize_image(arr: np.ndarray) -> np.ndarray:
    arr = np.abs(np.asarray(arr)).astype(np.float64)
    peak = float(np.max(arr)) if arr.size else 0.0
    return arr / peak if peak > 0 else arr

def _save_image(path: Path, image: np.ndarray, title: str, xlabel: str = "x (um)", ylabel: str = "y (um)") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4.2, 3.6), constrained_layout=True)
    im = ax.imshow(_normalize_image(image), cmap="magma", origin="lower", aspect="auto")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(path, dpi=160)
    plt.close(fig)

def _save_line(path: Path, x: np.ndarray, y: np.ndarray, title: str, xlabel: str, ylabel: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4.8, 3.2), constrained_layout=True)
    ax.plot(x, _normalize_image(y), linewidth=1.5)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=160)
    plt.close(fig)

def generate_validation_figures(
    out_dir: str | Path,
    mic: dict,
    reconstruction: dict,
    direct_psf: np.ndarray,
    converted_psf: np.ndarray,
) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    mic_stack = np.asarray(mic.get("corrected", mic["measured"]))
    z_index = int(np.argmax(np.max(np.abs(mic_stack), axis=(1, 2))))
    y_index = mic_stack.shape[1] // 2
    x_index = mic_stack.shape[2] // 2
    _save_image(out / "mic_central_xy.png", mic_stack[z_index], "Microscope central XY")
    _save_image(out / "mic_xz.png", mic_stack[:, y_index, :], "Microscope XZ", xlabel="x index", ylabel="z index")
    _save_image(out / "mic_yz.png", mic_stack[:, :, x_index], "Microscope YZ", xlabel="y index", ylabel="z index")

    _save_line(
        out / "oct_a_scan.png",
        np.asarray(reconstruction["depth_um"]),
        np.asarray(reconstruction["magnitude"]),
        "OCT reconstructed A-scan",
        "depth (um)",
        "normalized magnitude",
    )

    direct = _normalize_image(direct_psf)
    converted = _normalize_image(converted_psf)
    z = int(np.argmax(np.max(direct, axis=(1, 2))))
    residual = np.abs(direct[z] - converted[z])
    _save_image(out / "path_e_residual_xy.png", residual, "Path E residual XY")

def generate_sweep_trend_figure(out_dir: str | Path, rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    labels = [str(row["run_id"]) for row in rows]
    y = np.array([float(row["nrmse_3d"]) for row in rows], dtype=np.float64)
    x = np.arange(len(rows), dtype=np.float64)

    fig, ax = plt.subplots(figsize=(8.5, 3.6), constrained_layout=True)
    ax.plot(x, y, marker="o", linewidth=1.2, markersize=3.2)
    ax.set_title("Sweep trend")
    ax.set_ylabel("Path E NRMSE")
    ax.set_xlabel("sweep case")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=65, ha="right", fontsize=6)
    ax.grid(True, alpha=0.25)
    fig.savefig(out / "sweep_trend.png", dpi=160)
    plt.close(fig)
