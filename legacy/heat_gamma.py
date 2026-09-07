# LEGACY SEMICONDUCTOR SUPPORT CODE; NOT ANTHRACYCLINE CHEMISTRY.
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import ROS_Util as rc


def _parse_values(arg: str) -> list[float]:
    """Return a list of floats from ``arg`` which may be ``a,b,c`` or ``start:stop:step``."""
    if ":" in arg:
        parts = arg.split(":")
        if len(parts) != 3:
            raise argparse.ArgumentTypeError("Range must be start:stop:step")
        start, stop, step = map(float, parts)
        n = int(round((stop - start) / step)) + 1
        return [start + i * step for i in range(n)]
    return [float(x) for x in arg.split(",")]


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Plot effective gamma over B_rms and tau")
    p.add_argument("--B_rms", required=True,
                   help="Range or comma list of magnetic field RMS values")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--tau",
                   help="Range or comma list of dephasing times")
    group.add_argument("--beta",
                       help="Range or comma list of dimensionless beta = gamma_base * tau")
    p.add_argument("--weights", type=Path, default=Path("dataset/fig06"),
                   help="Directory containing weighting CSV files")
    p.add_argument("--protocols", default="f_ow,f_fb",
                   help="Comma-separated protocol column names to use from weighting data")
    p.add_argument("--dt", type=float,
                   help="Time step; defaults to value from dataset/fig03")
    p.add_argument("--gamma_k", type=float, default=1.0e4,
                   help=("Scaling coefficient k for rc.gamma_base;"
                        " gamma = min(0.25, k * B_rms * dt)"))
    p.add_argument("--csv_out", type=Path,
                   help="Optional path to save the gamma table as CSV")
    p.add_argument("--figure_out", type=Path,
                   help="Path to save the generated figure")
    p.add_argument("--show", action="store_true",
                   help="Display the generated figure")
    args = p.parse_args(argv)

    if not args.weights.exists():
        sys.exit(f"Weighting data directory not found: {args.weights}")

    B_vals = _parse_values(args.B_rms)
    if args.tau:
        y_vals = _parse_values(args.tau)
        y_label = "tau"
    else:
        y_vals = _parse_values(args.beta)
        y_label = "beta"

    if args.dt is None:
        axis = Path("dataset/fig03/Sim07_20240405_stochastic_field_axis_3.json")
        vals = Path("dataset/fig03/Sim07_20240405_stochastic_field_vals_3.csv")
        if not axis.exists():
            sys.exit(f"Field axis file not found: {axis}")
        if not vals.exists():
            sys.exit(f"Field values file not found: {vals}")
        dt, _ = rc.load_field(axis, vals)
    else:
        dt = args.dt

    grids: dict[str, np.ndarray] = {}
    protocols = [p.strip() for p in args.protocols.split(",") if p.strip()]
    for protocol in protocols:
        data = np.zeros((len(y_vals), len(B_vals)), dtype=float)
        for i, y in enumerate(y_vals):
            for j, B in enumerate(B_vals):
                g_base = rc.gamma_base(B, dt, args.gamma_k)
                if args.tau:
                    tau = y
                    beta = g_base * tau
                else:
                    beta = y
                    tau = beta / g_base
                weight = rc.weight_factor(beta, protocol, args.weights)
                data[i, j] = g_base * weight
        grids[protocol] = data
    
    fig, axes = plt.subplots(1, len(protocols), squeeze=False, figsize=(6 * len(protocols), 4))
    for ax, (protocol, data) in zip(axes.flat, grids.items()):
        mesh = ax.pcolormesh(B_vals, y_vals, data, shading="auto", cmap="viridis")
        ax.set_xlabel("B_rms")
        ax.set_ylabel(y_label)
        ax.set_title(protocol)
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label("gamma_eff")
    plt.tight_layout()
    if args.figure_out:
        plt.savefig(args.figure_out)
    if args.show:
        plt.show()

    if args.csv_out:
        for protocol, data in grids.items():
            df = pd.DataFrame({
                "B_rms": np.tile(B_vals, len(y_vals)),
                y_label: np.repeat(y_vals, len(B_vals)),
                "gamma_eff": data.ravel(),
            })
            out = args.csv_out.with_name(f"{args.csv_out.stem}_{protocol}{args.csv_out.suffix}")
            df.to_csv(out, index=False)


if __name__ == "__main__":
    main()

