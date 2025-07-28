from pathlib import Path
import json
import numpy as np, pandas as pd
from scipy.interpolate import interp1d

"""#Absolute paths used, reproducing could dictate altering"""
AXIS_JSON = Path("/dataset/Fig03/Sim07_20240405_stochastic_field_axis_3.json")
VALS_CSV  = Path("/dataset/Fig03/Sim07_20240405_stochastic_field_vals_3.csv")

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

def weight_factor(beta: float, protocol: str, directory: Path) -> float:
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
    interp = interp1d(tbl["beta"], tbl[protocol], fill_value="extrapolate")
    return float(interp(beta))
