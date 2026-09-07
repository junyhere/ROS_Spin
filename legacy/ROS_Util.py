# LEGACY SEMICONDUCTOR SUPPORT CODE; NOT ANTHRACYCLINE CHEMISTRY.
from pathlib import Path
import json
import warnings
import numpy as np, pandas as pd
from scipy.interpolate import interp1d

"""#Absolute paths used, reproducing could dictate altering"""
AXIS_JSON = Path("/dataset/fig03/Sim07_20240405_stochastic_field_axis_3.json")
VALS_CSV  = Path("/dataset/fig03/Sim07_20240405_stochastic_field_vals_3.csv")

def load_field( axis_path: Path = AXIS_JSON, vals_path: Path = VALS_CSV,) -> tuple[float, float]:
    """Return the mean time step and RMS field magnitude from dataset files."""
    with axis_path.open("r", encoding="utf-8") as f:
        t = np.asarray(json.load(f)["t"], dtype=float)
    if t.size < 2:
        raise ValueError("Need at least two time points to compute change in time")
    dt = float(np.diff(t).mean())
    
    field = pd.read_csv(vals_path, header=None).to_numpy(dtype=float)
    B_rms = float(np.sqrt(np.mean(field**2)))

    return dt, B_rms

def gamma_base(B_rms: float, dt: float, k: float = 1.0e4) -> float:
    """Compute the baseline decoherence rate for a given field RMS and step."""
    return min(0.25, k * B_rms * dt)


_weight_cache: dict[Path, pd.DataFrame] = {}

def weight_factor(beta: float, protocol: str, directory: Path, *, out_of_range: str = "error",) -> float:
    """Interpolate protocol-specific weighting factors from fig06 CSV files."""
    if directory not in _weight_cache:
        files = sorted(directory.glob("Sim08*_OU_*.csv"))
        if not files:
            raise FileNotFoundError(f"No weighting CSVs found in {directory}")
        df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
        df = df.drop_duplicates(subset="beta").sort_values("beta")
        _weight_cache[directory] = df
    tbl = _weight_cache[directory]
    if protocol not in tbl.columns:
        raise KeyError(f"{protocol} not found in weighting data")
    betas = tbl["beta"]
    min_beta = float(betas.iloc[0])
    max_beta = float(betas.iloc[-1])
    if beta < min_beta or beta > max_beta:
        if out_of_range == "error":
            raise ValueError(
                f"beta={beta} outside weighting data range [{min_beta}, {max_beta}]"
            )
        elif out_of_range == "clamp":
            clamped = min(max(beta, min_beta), max_beta)
            warnings.warn(
                f"beta={beta} outside weighting data range [{min_beta}, {max_beta}], clamping to {clamped}",
                RuntimeWarning,
            )
            beta = clamped
        else:
            raise ValueError("out_of_range must be 'error' or 'clamp'")
    interp = interp1d(betas, tbl[protocol])
    return float(interp(beta))
