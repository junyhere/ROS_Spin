import argparse
import itertools
import subprocess
import sys


def _parse_values(text: str) -> list[str]:
    """Return a list parsed from ``text`` as CSV or ``start:stop:step``."""
    if ":" in text:
        start, stop, step = (float(x) for x in text.split(":"))
        vals = []
        v = start
        # ensure inclusive range accounting for floating point error
        while v <= stop + 1e-12:
            vals.append(str(int(v)) if v.is_integer() else str(v))
            v += step
        return vals
    return [t for t in text.split(",") if t]

def main() -> None:
    """Batch utility to run ROS.py displaying results from different delays and protocols."""
    p = argparse.ArgumentParser()
    p.add_argument("--delays", default="1,2,3,4,5")
    p.add_argument("--tau", default="0.1")
    p.add_argument("--protocols", default="f_ow,f_fb")
    p.add_argument("--phi_frac", default="0.0")
    p.add_argument("--trotter", default="0")
    p.add_argument("--out", default="table.csv")
    p.add_argument("--trotter", nargs="?",
                   help="Optional CZ step count; empty to auto‑recommend")
    args, ros_args = p.parse_known_args()

    delays = _parse_values(args.delays)
    taus = _parse_values(args.tau)
    protos = _parse_values(args.protocols)
    phis = _parse_values(args.phi_frac)
    trotters = _parse_values(args.trotter)

    
    for d, tau, proto, phi, N in itertools.product(
        delays, taus, protos, phis, trotters
    ):
        cmd = [
            sys.executable,
            "ROS.py",
            *ros_args,
            "--delay",
            d,
            "--tau",
            tau,
            "--protocols",
            proto,
            "--phi_frac",
            phi,
            "--trotter",
            N,
            "--table",
            "--csv_out",
            args.out,
        ]
        if args.trotter:
            cmd.extend(["--trotter", args.trotter])
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
