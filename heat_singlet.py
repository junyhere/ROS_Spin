from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ROS import build_rp_circuit, simulate, counts_to_ros
from noise import noise_mod
import ROS_Util as rc
import error


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
    p = argparse.ArgumentParser(description="Simulate singlet ratio over B_rms and tau/beta grid")
    p.add_argument("--B_rms", required=True, help="Range or comma list of B_rms values")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--tau", help="Range or comma list of dephasing times")
    group.add_argument("--beta",help="Range or comma list of dimensionless beta = gamma_base * tau")
    p.add_argument("--weights", type=Path, default=Path("dataset/fig06"),
                   help="Directory containing weighting CSV files")
    p.add_argument("--protocols", default="f_ow,f_fb",
                   help="Comma-separated protocol column names to use from weighting data")
    p.add_argument("--delay", type=int, default=4, help="Number of idle gates")
    p.add_argument("--trotter", type=int, help="Optional number of CZ steps")
    p.add_argument("--shots", type=int, default=10000)
    p.add_argument("--seed", type=int,
                   help="Random seed for deterministic runs")
    p.add_argument("--phi_frac", type=float, default=0.0)
    p.add_argument("--gamma_k", type=float, default=1.0e4,
                   help="Scaling coefficient k for rc.gamma_base;"
                        " gamma = min(0.25, k * B_rms * dt)")
    p.add_argument("--error_prefix")
    p.add_argument("--error_method", choices=["MC", "NI"])
    p.add_argument("--target_error", type=float,
                   help="Desired accuracy for parse_fig10.recommend_N")
    p.add_argument("--csv_out", type=Path, help="Optional path to save the table as CSV")
    p.add_argument("--figure_out", type=Path,
                   help="Path to save the generated figure")
    p.add_argument("--show", action="store_true",
                   help="Display the generated figure")
    args = p.parse_args(argv)

    if not args.weights.exists():
        sys.exit(f"Weighting data directory not found: {args.weights}")

    B_vals = np.array(_parse_values(args.B_rms))
    if args.tau:
        y_vals = np.array(_parse_values(args.tau))
        y_label = "tau"
    else:
        y_vals = np.array(_parse_values(args.beta))
        y_label = "beta"
    protocols = [p.strip() for p in args.protocols.split(",") if p.strip()]

    axis = Path("dataset/fig03/Sim07_20240405_stochastic_field_axis_3.json")
    vals = Path("dataset/fig03/Sim07_20240405_stochastic_field_vals_3.csv")
    if not axis.exists():
        sys.exit(f"Field axis file not found: {axis}")
    if not vals.exists():
        sys.exit(f"Field values file not found: {vals}")
    dt, _ = rc.load_field(axis, vals)
    
    # Precompute gamma_base for each magnetic field value once
    gamma_base_vals = np.array([rc.gamma_base(b, dt, args.gamma_k) for b in B_vals])

    grids: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    if args.error_prefix and args.error_method and args.target_error is not None and args.trotter is None:
        try:
            df_err = error.load_error(args.error_prefix, args.error_method)
        except FileNotFoundError as e:
            sys.exit(f"Error data file not found: {e}")
        except ValueError as e:
            sys.exit(str(e))
        rec_n = error.recommend_N(df_err, args.target_error)
        args.trotter = rec_n
        print(f"Recommended N from {args.error_prefix}_{args.error_method}: {rec_n}")

    qc = build_rp_circuit(delay_ids=args.delay, trotter=args.trotter)

    for protocol in protocols:
        g_eff_grid = np.zeros((len(y_vals), len(B_vals)), dtype=float)
        for i, y in enumerate(y_vals):
            for j, g_base in enumerate(gamma_base_vals):
                beta = g_base * y if args.tau else y
                weight = rc.weight_factor(beta, protocol, args.weights)
                g_eff_grid[i, j] = g_base * weight

        g_eff_flat = g_eff_grid.ravel()
        data_flat = np.empty_like(g_eff_flat)
        for idx, g_eff in enumerate(g_eff_flat):
            noise = noise_mod(g_eff, args.phi_frac)
            s, _ = counts_to_ros(simulate(qc, noise, args.shots, seed=args.seed))
            data_flat[idx] = s
        data = data_flat.reshape(g_eff_grid.shape)
        grids[protocol] = (g_eff_grid, data)

    fig, axes = plt.subplots(1, len(protocols), squeeze=False, figsize=(6 * len(protocols), 4))
    for ax, protocol in zip(axes.flat, protocols):
        g_eff_grid, data = grids[protocol]
        mesh = ax.pcolormesh(B_vals, y_vals, data, shading="auto", cmap="viridis")
        ax.set_xlabel("B_rms")
        ax.set_ylabel(y_label)
        ax.set_title(protocol)
        fig.colorbar(mesh, ax=ax, label="singlet ratio")
    plt.tight_layout()
    if args.figure_out:
        plt.savefig(args.figure_out)
    if args.show:
        plt.show()

    if args.csv_out:
        for protocol, (g_eff_grid, data) in grids.items():
            df = pd.DataFrame({
                "B_rms": np.tile(B_vals, len(y_vals)),
                y_label: np.repeat(y_vals, len(B_vals)),
                "gamma_eff": g_eff_grid.ravel(),
                "singlet_ratio": data.ravel(),
            })
            out = args.csv_out.with_name(f"{args.csv_out.stem}_{protocol}{args.csv_out.suffix}")
            df.to_csv(out, index=False)
            
if __name__ == "__main__":
    main()





