"""Render the report's display equations as high-resolution transparent PNGs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ros_spin_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EQUATIONS = {
    "projectors": r"$P_D=\frac{0.5I-\mathbf{s}\cdot\mathbf{S}}{1.5},\qquad P_Q=\frac{\mathbf{s}\cdot\mathbf{S}+I}{1.5}$",
    "hamiltonian": r"$\frac{H}{\hbar}=\beta_e\mathbf{B}\cdot(g_{\mathrm{SQ}}\mathbf{s}+g_{\mathrm{O_2}}\mathbf{S})+J\mathbf{s}\cdot\mathbf{S}+\mathbf{s}\cdot\mathbf{D}_{\mathrm{dip}}\cdot\mathbf{S}+\mathbf{S}\cdot\mathbf{D}_{\mathrm{O_2}}\cdot\mathbf{S}+\boldsymbol{\Omega}_{\mathrm{local}}\cdot\mathbf{s}$",
    "master_equation": r"$\frac{d\rho}{dt}=-i[H/\hbar,\rho]-\frac{1}{2}\left\{k_DP_D+k_QP_Q+k_{\mathrm{escape}}I,\rho\right\}+\mathcal{L}_{\mathrm{relax}}(\rho)$",
    "reaction_yields": r"$Y_D(t)=\int_0^t k_D\operatorname{Tr}[P_D\rho(\tau)]\,d\tau,\quad Y_Q(t)=\int_0^t k_Q\operatorname{Tr}[P_Q\rho(\tau)]\,d\tau,\quad Y_{\mathrm{primary}}=Y_D+Y_Q$",
    "probability_balance": r"$Y_D+Y_Q+Y_{\mathrm{escape}}+\operatorname{Tr}(\rho)=1$",
    "h2o2_stoichiometry": r"$2\ \mathrm{radical\ equivalents}\longrightarrow 1\ \mathrm{H_2O_2}$",
}


def render(output_dir: Path, dpi: int = 300) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Serif",
        "mathtext.fontset": "dejavuserif",
        "figure.facecolor": "none",
        "savefig.facecolor": "none",
    })
    written = []
    for name, equation in EQUATIONS.items():
        figure = plt.figure(figsize=(12, 0.85))
        figure.text(0.5, 0.5, equation, ha="center", va="center", fontsize=15)
        output = output_dir / f"{name}.png"
        figure.savefig(
            output, dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0.06
        )
        plt.close(figure)
        written.append(output)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", nargs="?", type=Path, default=Path("assets/equations"))
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if not 150 <= args.dpi <= 1200:
        raise SystemExit("--dpi must be between 150 and 1200")
    files = render(args.output_dir.resolve(), args.dpi)
    print(json.dumps({"status": "complete", "files": [str(path) for path in files]}))


if __name__ == "__main__":
    main()
