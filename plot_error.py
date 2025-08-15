"""Plot N versus epsilon/delta from fig10-12 tables."""

from __future__ import annotations

from pathlib import Path
import argparse

import matplotlib.pyplot as plt

from error import load_error


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot N vs epsilon/delta from numerical error tables")
    parser.add_argument(
        "prefix",
        nargs="+",
        help="Dataset prefix without _MC.csv or _NI.csv",
    )
    parser.add_argument(
        "-m",
        "--method",
        choices=["MC", "NI"],
        default="MC",
        help="Dataset type to load",
    )
    parser.add_argument(
        "-d",
        "--directory",
        type=Path,
        default=Path("dataset/fig10-12"),
        help="Directory containing tables",
    )
    parser.add_argument(
        "-s", "--save", type=Path, help="Path to save image instead of showing"
    )
    args = parser.parse_args()

    fig, ax = plt.subplots()
    for pref in args.prefix:
        df = load_error(pref, args.method, args.directory)
        if "delta" not in df.columns:
            raise ValueError(
                f"{pref} does not contain 'delta' column required for epsilon/delta"
            )
        ratio = df["epsilon"] / df["delta"]
        ax.loglog(df["N"], ratio, marker="o", label=pref)
    ax.set_xlabel("N")
    ax.set_ylabel(r"$\epsilon/\delta$")
    ax.legend()
    ax.grid(True, which="both", ls=":")
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save)
    else:
        plt.show()


if __name__ == "__main__":
    main()