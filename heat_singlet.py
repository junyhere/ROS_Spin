from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ROS import build_rp_circuit, simulate, counts_to_ros
from noise import noise_mod


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
    p = argparse.ArgumentParser(description="Simulate singlet ratio over B_rms and gamma grid")
    p.add_argument("--B_rms", required=True, help="Range or comma list of B_rms values")
    p.add_argument("--gamma", required=True, help="Range or comma list of gamma values")
    p.add_argument("--delay", type=int, default=4, help="Number of idle gates")
    p.add_argument("--trotter", type=int, help="Optional number of CZ steps")
    p.add_argument("--shots", type=int, default=5000)
    p.add_argument("--phi_frac", type=float, default=0.0)
    p.add_argument("--csv_out", type=Path, help="Optional path to save the table as CSV")
    args = p.parse_args(argv)

    B_vals = _parse_values(args.B_rms)
    g_vals = _parse_values(args.gamma)

    data = np.zeros((len(g_vals), len(B_vals)), dtype=float)

    qc = build_rp_circuit(delay_ids=args.delay, trotter=args.trotter)

    for i, g in enumerate(g_vals):
        for j, _B in enumerate(B_vals):
            noise = noise_mod(g, args.phi_frac)
            s, _ = counts_to_ros(simulate(qc, noise, args.shots))
            data[i, j] = s

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(B_vals, g_vals, data, shading="auto", cmap="viridis")
    ax.set_xlabel("B_rms")
    ax.set_ylabel("gamma")
    fig.colorbar(mesh, ax=ax, label="singlet ratio")
    plt.show()

    if args.csv_out:
        df = pd.DataFrame(data, index=g_vals, columns=B_vals)
        df.index.name = "gamma"
        df.to_csv(args.csv_out)


if __name__ == "__main__":
    main()
