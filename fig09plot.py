from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def _parse_values(arg: str | None) -> list[float]:
    """Parse comma lists or start:stop:step ranges."""
    if arg is None:
        return []
    if ":" in arg:
        parts = arg.split(":")
        if len(parts) != 3:
            raise ValueError(
                f"Expected 'start:stop:step' with three values, got {arg!r}"
            )
        try:
            start, stop, step = map(float, parts)
        except ValueError as exc:
            raise ValueError(
                f"Could not parse start, stop, and step from {arg!r}"
            ) from exc
        n = int(round((stop - start) / step)) + 1
        return [start + i * step for i in range(n)]
    try:
        return [float(x) for x in arg.split(",")]
    except ValueError as exc:
        raise ValueError(
            f"Could not parse comma-separated values from {arg!r}"
        ) from exc


def load_fig09(directory: Path, tau_filter: list[float] | None):
    for file in sorted(directory.glob("realistic_tau=*.csv")):
        tau = float(file.stem.split("=")[1])
        if tau_filter and tau not in tau_filter:
            continue
        df = pd.read_csv(file)
        yield tau, df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description="Plot velocity vs numerical infidelity for fig09 datasets"
    )
    p.add_argument("--directory", type=Path, default=Path("dataset/fig09"),
                   help="Directory containing fig09 CSV files")
    p.add_argument("--tau", help="Comma list or start:stop:step of tau values to include")
    p.add_argument("--x-label", default="v", help="Label for velocity axis")
    p.add_argument("--y-label", default=r"$\chi_{\mathrm{NI}}$", help="Label for infidelity axis")
    p.add_argument("--logx", action="store_true", help="Log-scale the x-axis")
    p.add_argument("--logy", action="store_true", help="Log-scale the y-axis")
    p.add_argument("--legend-loc", default="best", help="Legend location")
    p.add_argument("--title", help="Figure title")
    p.add_argument("--save", type=Path, help="Path to save the image")
    p.add_argument("--dpi", type=int, default=300, help="Resolution when saving")
    p.add_argument("--show", action="store_true", help="Display the figure")
    args = p.parse_args(argv)

    tau_filter = _parse_values(args.tau)
    fig, ax = plt.subplots()

    for tau, df in load_fig09(args.directory, tau_filter or None):
        ax.plot(df["v"], df["chi_ni"], marker="o", label=f"tau={tau}")

    ax.set_xlabel(args.x_label)
    ax.set_ylabel(args.y_label)
    if args.title:
        ax.set_title(args.title)
    if args.logx:
        ax.set_xscale("log")
    if args.logy:
        ax.set_yscale("log")
    ax.grid(True, which="both", ls=":")
    ax.legend(loc=args.legend_loc)
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=args.dpi)
    if args.show or not args.save:
        plt.show()


if __name__ == "__main__":
    main()
