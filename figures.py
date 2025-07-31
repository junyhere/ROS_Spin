from __future__ import annotations
import argparse
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
    p.add_argument("--tau", required=True,
                   help="Range or comma list of dephasing times")
    p.add_argument("--weights", type=Path, default=Path("dataset/fig06"),
                   help="Directory containing weighting CSV files")
    p.add_argument("--protocol", default="f_ow",
                   help="Protocol column name to use from weighting data")
    p.add_argument("--dt", type=float,
                   help="Time step; defaults to value from dataset/fig03")
    p.add_argument("--csv_out", type=Path,
                   help="Optional path to save the gamma table as CSV")
    args = p.parse_args(argv)

    B_vals = _parse_values(args.B_rms)
    tau_vals = _parse_values(args.tau)

    if args.dt is None:
        axis = Path("dataset/fig03/Sim07_20240405_stochastic_field_axis_3.json")
        vals = Path("dataset/fig03/Sim07_20240405_stochastic_field_vals_3.csv")
        dt, _ = rc.load_field(axis, vals)
    else:
        dt = args.dt

    data = np.zeros((len(tau_vals), len(B_vals)), dtype=float)
    for i, tau in enumerate(tau_vals):
        for j, B in enumerate(B_vals):
            g_base = rc.gamma_base(B, dt)
            beta = g_base * tau
            weight = rc.weight_factor(beta, args.protocol, args.weights)
            data[i, j] = g_base * weight

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(B_vals, tau_vals, data, shading="auto", cmap="viridis")
    ax.set_xlabel("B_rms")
    ax.set_ylabel("tau")
    fig.colorbar(mesh, ax=ax, label="gamma")
    plt.show()

    if args.csv_out:
        df = pd.DataFrame(data, index=tau_vals, columns=B_vals)
        df.index.name = "tau"
        df.to_csv(args.csv_out)


if __name__ == "__main__":
    main()
