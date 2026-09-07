# LEGACY SEMICONDUCTOR SUPPORT CODE; NOT ANTHRACYCLINE CHEMISTRY.
from __future__ import annotations

from pathlib import Path
import json
import pandas as pd


def load_error(prefix: str, method: str, directory: Path = Path("dataset/fig10-12")) -> pd.DataFrame:
    """Load numerical error CSV and attach metadata from the matching JSON.

    Parameters
    ----------
    prefix:
        File prefix shared by the CSV and JSON files.
    method:
        Either ``"MC"`` or ``"NI"`` selecting the Monte Carlo or numerical
        integration table.
    directory:
        Folder containing the dataset files.
    """
    method = method.upper()
    if method not in {"MC", "NI"}:
        raise ValueError("method must be 'MC' or 'NI'")

    csv_file = directory / f"{prefix}_{method}.csv"
    if not csv_file.exists():
        raise FileNotFoundError(csv_file)
    df = pd.read_csv(csv_file)

    json_file = directory / f"{prefix}.json"
    if json_file.exists():
        with json_file.open("r", encoding="utf-8") as f:
            meta = json.load(f)
        for k, v in meta.items():
            df[k] = v
    return df


def recommend_N(df: pd.DataFrame, epsilon: float) -> int:
    """Return the smallest ``N`` for which the error is below ``epsilon``."""
    metric = "delta" if "delta" in df.columns else "epsilon"
    df_sorted = df.sort_values("N")
    mask = df_sorted[metric] <= epsilon
    if mask.any():
        return int(df_sorted.loc[mask, "N"].iloc[0])
    return int(df_sorted["N"].iloc[-1])
