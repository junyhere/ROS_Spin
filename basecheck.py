from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd


def load_pair(tau: str, directory: Path = Path("dataset/fig09")) -> tuple[pd.DataFrame, dict]:
    """Return the CSV/JSON data for the given ``tau`` value."""
    base = directory / f"realistic_tau={tau}"
    df = pd.read_csv(base.with_suffix(".csv"))
    with base.with_suffix(".json").open("r", encoding="utf-8") as f:
        meta = json.load(f)
    return df, meta


def baseline_gamma(meta: dict, k: float = 1.0e4) -> float:
    """Estimate ``gamma_base`` using the noise ``sigma`` and first time step."""
    sigma = meta.get("B", {}).get("sigma")
    if sigma is None:
        sigma = meta.get("B", {}).get("σ")
    if sigma is None:
        raise KeyError("sigma not found in metadata")
    n0 = float(meta.get("N", [1])[0])
    T0 = float(meta.get("T0", 1.0))
    dt = T0 / n0
    return min(0.25, k * float(sigma) * dt)


def _cli() -> None:
    p = argparse.ArgumentParser(description="Show how fig09 data can inform the ROS model")
    p.add_argument("tau", help="Value used in the dataset file names, e.g. 0.1")
    args = p.parse_args()

    df, meta = load_pair(args.tau)
    g = baseline_gamma(meta)
    print(f"Loaded {len(df)} rows for tau={args.tau}")
    print(f"sigma = {meta['B'].get('sigma', meta['B'].get('σ'))}")
    print(f"gamma_base = {g:.6f}")
    print(df.head().to_string(index=False))


if __name__ == "__main__":
    _cli()
