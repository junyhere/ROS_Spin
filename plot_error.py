from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from error import load_error


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description="Plot numerical error scaling from fig10-12 tables"
    )
    p.add_argument("prefix", nargs="+",
                   help="Dataset prefix without _MC.csv or _NI.csv")
    p.add_argument("--method", choices=["MC", "NI"], default="MC",
                   help="Dataset type to load")
    p.add_argument("--directory", type=Path, default=Path("dataset/fig10-12"),
                   help="Directory containing tables")
    p.add_argument("--x", choices=["N", "cpu_time", "ram_bytes"], default="N",
                   help="Column for x-axis")
    p.add_argument("--y", choices=["epsilon", "delta", "ratio"], default="ratio",
                   help="Quantity for y-axis; 'ratio' plots epsilon/delta")
    p.add_argument("--logx", action="store_true", help="Log-scale the x-axis")
    p.add_argument("--logy", action="store_true", help="Log-scale the y-axis")
    p.add_argument("--markers", action="store_true", help="Draw markers on lines")
    p.add_argument("--legend-loc", default="best", help="Legend location")
    p.add_argument("--title", help="Figure title")
    p.add_argument("--save", type=Path, help="Path to save the image")
    p.add_argument("--dpi", type=int, default=300, help="Resolution when saving")
    p.add_argument("--show", action="store_true", help="Display the figure")
    args = p.parse_args(argv)

    fig, ax = plt.subplots()
    for pref in args.prefix:
        df = load_error(pref, args.method, args.directory)
        x = df[args.x]
        if args.y == "ratio":
            if "delta" not in df.columns:
                raise ValueError(f"{pref} lacks 'delta' column for epsilon/delta")
            y = df["epsilon"] / df["delta"]
            y_label = r"$\\epsilon/\\delta$"
        else:
            y = df[args.y]
            y_label = args.y
        ax.plot(x, y, marker="o" if args.markers else None, label=pref)

    ax.set_xlabel(args.x)
    ax.set_ylabel(y_label)
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
