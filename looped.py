from pathlib import Path
import argparse, subprocess, sys


def main() -> None:
    """Batch script to run ROS.py displaying results from different delays and protocols."""
    p = argparse.ArgumentParser()
    p.add_argument("--delays", default="1,2,3,4,5")
    p.add_argument("--json", default="dataset/fig03/Sim07_20240405_stochastic_field_axis_3.json")
    p.add_argument("--csv", default="dataset/fig03/Sim07_20240405_stochastic_field_vals_3.csv")
    p.add_argument("--tau", default="0.1")
    p.add_argument("--protocols", default="f_ow,f_fb")
    p.add_argument("--out", default="heatmap.csv")
    args = p.parse_args()

    for d in args.delays.split(","):
        cmd = [
            sys.executable,
            "ROS.py",
            "--json", args.json,
            "--csv", args.csv,
            "--tau", args.tau,
            "--protocols", args.protocols,
            "--delay", d,
            "--heatmap",
            "--csv_out", args.out,
        ]
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
