from pathlib import Path
import json, numpy as np, pandas as pd
from scipy.interpolate import interp1d

#Absolute paths used, reproducing could dictate altering
AXIS_JSON = Path("/DATASET/Fig03/Sim07_20240405_stochastic_field_axis_3.json")
VALS_CSV  = Path("/DATASET/Fig03/Sim07_20240405_stochastic_field_vals_3.csv")

def load_field( axis_path: Path = AXIS_JSON, vals_path: Path = VALS_CSV,) -> tuple[float, float]:
    #Finding average time step
    with axis_path.open("r", encoding="utf-8") as f:
        t = np.asarray(json.load(f)["t"], dtype=float)
    if t.size < 2:
        raise ValueError("Need at least two time points to compute change in time")
    dt = float(np.diff(t).mean())
    
    #Field Root Mean Square calculation 
    field = pd.read_csv(vals_path, header=None).to_numpy(dtype=float)
    B_rms = float(np.sqrt(np.mean(field**2)))

    return dt, B_rms

#Base decoherence
def gamma_base(B_rms: float, dt: float, k: float = 1.0e4) -> float:
    return min(0.25, k * B_rms * dt)

