"""Runtime boundary for retained semiconductor/two-qubit entry points."""
from __future__ import annotations

from pathlib import Path


MESSAGE = (
    "LEGACY SEMICONDUCTOR/TWO-QUBIT MATERIAL: this entry point does not produce "
    "results for the doxorubicin semiquinone/O2 model. The historical source is "
    "retained under legacy/. Use `python3 paper_analysis.py --quick --output-dir "
    "results/paper` for the active chemical sensitivity study."
)


def retired_entrypoint(script_name: str) -> None:
    historical = Path(__file__).resolve().parent / "legacy" / script_name
    raise SystemExit(f"{MESSAGE}\nHistorical source: {historical}")
