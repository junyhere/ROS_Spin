from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd
import matplotlib.pyplot as plt


def load_fig09(directory: Path) -> list[tuple[float, pd.DataFrame]]:
    """Return list of (tau, dataframe) loaded from fig09 CSV files."""
    files = sorted(directory.glob("realistic_tau=*.csv"))
    data: list[tuple[float, pd.DataFrame]] = []
    for file in files:
        try:
            tau = float(file.stem.split("=")[1])
        except (IndexError, ValueError):
            continue
        df = pd.read_csv(file)
        if "v" not in df.columns or "chi_ni" not in df.columns:
            raise ValueError(f"{file} missing required columns 'v' and 'chi_ni'")
        data.append((tau, df[["v", "chi_ni"]]))
    data.sort(key=lambda x: x[0])
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot v vs chi_ni for fig09 datasets")
    parser.add_argument(
        "-d",
        "--directory",
        type=Path,
        default=Path("dataset/fig09"),
        help="Directory containing fig09 CSV files",
    )
    parser.add_argument(
        "-s", "--save", type=Path, help="Path to save image instead of showing"
    )
    args = parser.parse_args()

    fig, ax = plt.subplots()
    for tau, df in load_fig09(args.directory):
        ax.plot(df["v"], df["chi_ni"], label=f"tau={tau}")
    ax.set_xlabel("v")
    ax.set_ylabel(r"$\chi_{NI}$")
    ax.legend()
    ax.grid(True, which="both")
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save)
    else:
        plt.show()


if __name__ == "__main__":
    main()