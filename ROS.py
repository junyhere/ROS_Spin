from pathlib import Path
import argparse, textwrap, sys
import pandas as pd
import matplotlib.pyplot as plt

import error

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from noise import noise_mod
import ROS_Util as rc

def build_rp_circuit(delay_ids: int = 4, trotter: int | None = None) -> QuantumCircuit:
    """Create the two-qubit circuit used in the simulation."""
    qc = QuantumCircuit(2, 2)
    qc.x(1); qc.h(0); qc.cx(0,1); qc.z(0); qc.x([0,1])
    if trotter:
        for _ in range(trotter): qc.cz(0,1)
    else:
        for _ in range(delay_ids): qc.id([0,1])
    qc.measure([0,1],[0,1])
    return qc

def simulate(qc: QuantumCircuit, noise, shots: int):
    """Run the AerSimulator with ``noise`` and return raw counts."""
    backend = AerSimulator(noise_model=noise)
    tcirc   = transpile(qc, backend, optimization_level=0)
    job     = backend.run(tcirc, shots=shots)
    return job.result().get_counts()

def counts_to_ros(c: dict[str, int]):
    """Convert singlet/triplet counts into ROS fractions."""
    tot = sum(c.values())
    trip = c.get("00", 0) + c.get("11", 0)
    return (tot - trip) / tot, trip / tot

p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent("""
      Example (OW+FB heat‑map across delays):
        for d in 1 2 3 4 5; do
          python3 ROS.py --json axis.json --csv vals.csv --tau 0.1 \\
                         --protocols f_ow,f_fb --delay $d --heatmap done"""))
p.add_argument("--json"); p.add_argument("--csv")            #Sim07 file calls
p.add_argument("--surface_vals"); p.add_argument("--surface_axis")  #Sim06 file calls
p.add_argument("--tau_index", type=int); p.add_argument("--T0_index", type=int)
p.add_argument("--tau", type=float)
p.add_argument("--protocols", default="f_ow")
p.add_argument("--weights", type=Path, default=Path("dataset/fig06"))
p.add_argument("--phi_frac", type=float, default=0.0)
p.add_argument("--delay", type=int, default=4); p.add_argument("--trotter", type=int)
p.add_argument("--shots", type=int, default=10000)
p.add_argument("--error_prefix")
p.add_argument("--error_method", choices=["MC", "NI"])
p.add_argument("--target_error", type=float,
               help="Desired accuracy for parse_fig10.recommend_N")
p.add_argument("--plot", action="store_true")
p.add_argument("--heatmap", action="store_true",
               help="Accumulate multi‑protocol, multi‑delay results into a heat‑map")
p.add_argument("--csv_out")
a = p.parse_args()

if a.error_prefix and a.error_method and a.target_error is not None:
    df_err = error.load_error(a.error_prefix, a.error_method)
    rec_n = error.recommend_N(df_err, a.target_error)
    if a.trotter is None:
        a.trotter = rec_n
    print(f"Recommended N from {a.error_prefix}_{a.error_method}: {rec_n}")

"""Effective gamma for decoherence"""
records=[]
tau_val = None
T0_val = None
if a.surface_vals:
    if not a.surface_axis:
        sys.exit("--surface_axis required when using --surface_vals")
    vals = pd.read_csv(a.surface_vals, header=None).values
    axis = pd.read_csv(a.surface_axis)
    tau_axis = axis.iloc[:,0].to_numpy()
    T0_axis = axis.iloc[:,1].to_numpy()
    try:
        tau_val = float(tau_axis[a.tau_index])
        T0_val = float(T0_axis[a.T0_index])
    except IndexError as e:
        sys.exit(f"Index out of range: {e}")
    records.append(("surface", float(vals[a.tau_index, a.T0_index])))
else:
    if not (a.json and a.csv):
        sys.exit("Need --json/--csv or surface files")
    dt, B = rc.load_field(Path(a.json), Path(a.csv))
    g_base = rc.gamma_base(B, dt)
    for proto in a.protocols.split(","):
        beta = g_base * a.tau if a.tau else 0
        f = rc.weight_factor(beta, proto, a.weights) if a.tau else 1.0
        records.append((proto, g_base * f))

"""Simulation using obtained effective gamma"""
rows=[]
for proto, g_eff in records:
    qc = build_rp_circuit(delay_ids=a.delay, trotter=a.trotter)
    noise = noise_mod(g_eff, a.phi_frac)
    s,t = counts_to_ros(simulate(qc, noise, a.shots))
    row = {"delay":a.delay,"protocol":proto,
           "gamma":g_eff,"singlet":s,"triplet":t}
    if tau_val is not None:
        row["tau"] = tau_val
    if T0_val is not None:
        row["T0"] = T0_val
    rows.append(row)

df = pd.DataFrame(rows)

if a.heatmap:
    out = Path(a.csv_out or "heatmap.csv")
    df.to_csv(out, mode="a", index=False, header=not out.exists())
    if a.plot:
        data = pd.read_csv(out)
        pivot = data.pivot(index="delay", columns="protocol", values="singlet")
        plt.imshow(pivot, origin="lower", aspect="auto", cmap="viridis")
        plt.xlabel("protocol")
        plt.ylabel("delay")
        plt.colorbar(label="singlet")
        plt.show()
else:
    print(df.to_string(index=False, float_format="%.4f"))
    if a.csv_out:
        df.to_csv(a.csv_out, index=False)
