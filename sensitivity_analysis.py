"""Small, dimensionless, non-predictive encounter sensitivity sweep."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from spin_chemistry import (
    EncounterParameters,
    propagate_encounter_reference,
    validate_authority_bundle,
)


ROOT = Path(__file__).resolve().parent
LIMITATIONS = (
    "Illustrative one-factor-at-a-time grid; bounds are not measured ranges or "
    "probability priors. Results are populations and per-encounter yields, not "
    "concentrations or biological fluxes. No encounter yield is converted to "
    "molarity, and no continuous biological ROS source is modeled."
)

INITIAL_STATES = (
    ("unpolarized", None, "unpolarized"),
    ("doublet", None, "doublet-manifold mixture"),
    ("quartet", None, "quartet-manifold mixture"),
    ("mixture", 0.25, "pD=0.25 manifold mixture"),
    ("mixture", 0.75, "pD=0.75 manifold mixture"),
)


def normalized_axis(size: int, low: float = 1e-2, high: float = 1e2) -> np.ndarray:
    """Return zero plus a bounded logarithmic axis including exactly one.

    The values are illustrative computational coordinates, not measured ranges,
    confidence intervals, or priors.
    """
    if isinstance(size, bool) or not isinstance(size, int) or not 4 <= size <= 41:
        raise ValueError("grid size must be an integer from 4 through 41")
    if not (np.isfinite(low) and np.isfinite(high) and 0 < low < 1 < high):
        raise ValueError("normalized nonzero bounds must satisfy 0 < low < 1 < high")
    nonzero = np.geomspace(low, high, size - 1)
    nonzero[np.argmin(np.abs(np.log(nonzero)))] = 1.0
    return np.concatenate(([0.0], nonzero))


def selectivity_axis(size: int) -> np.ndarray:
    """Return a kQ/kD axis containing the zero and unity controls exactly."""
    return normalized_axis(size, 1e-2, 1e1)


def encounter_result_row(
    *,
    scenario_id: str,
    reference_rate_s: float,
    mixing_over_reference: float,
    radical_relaxation_over_reference: float,
    oxygen_relaxation_over_reference: float,
    escape_over_reference: float,
    kq_over_kd: float,
    kd_over_reference: float = 1.0,
    initial_state: str = "unpolarized",
    p_doublet: float | None = None,
    duration_over_reference: float = 8.0,
    samples: int = 2,
) -> dict:
    """Resolve and execute one fully recorded non-predictive encounter case."""
    if not np.isfinite(reference_rate_s) or reference_rate_s <= 0:
        raise ValueError("reference_rate_s must be finite and positive")
    coordinates = {
        "mixing_over_reference": mixing_over_reference,
        "radical_relaxation_over_reference": radical_relaxation_over_reference,
        "oxygen_relaxation_over_reference": oxygen_relaxation_over_reference,
        "escape_over_reference": escape_over_reference,
        "kq_over_kd": kq_over_kd,
        "kd_over_reference": kd_over_reference,
        "duration_over_reference": duration_over_reference,
    }
    for name, value in coordinates.items():
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    kd = kd_over_reference * reference_rate_s
    params = EncounterParameters(
        local_field_proxy_rad_s=(mixing_over_reference * reference_rate_s, 0.0, 0.0),
        radical_relaxation_s=radical_relaxation_over_reference * reference_rate_s,
        oxygen_relaxation_s=oxygen_relaxation_over_reference * reference_rate_s,
        k_doublet_s=kd,
        k_quartet_s=kq_over_kd * kd,
        k_escape_s=escape_over_reference * reference_rate_s,
    )
    duration_s = duration_over_reference / reference_rate_s
    result = propagate_encounter_reference(
        params, duration_s, initial_state, samples, p_doublet
    )
    initial_definitions = {
        "unpolarized": "unpolarized I6/6",
        "doublet": "doublet-manifold mixture PD/2",
        "quartet": "quartet-manifold mixture PQ/4",
        "mixture": f"manifold mixture pD={p_doublet}",
    }
    return {
        "scenario_id": scenario_id,
        "initial_state": initial_state,
        "initial_state_definition": initial_definitions[initial_state],
        "p_doublet_requested": "" if p_doublet is None else p_doublet,
        "p_doublet_initial": result["p_doublet_initial"],
        "reference_rate_s^-1": reference_rate_s,
        "duration_over_reference": duration_over_reference,
        "duration_s": duration_s,
        "mixing_over_reference": mixing_over_reference,
        "local_mixing_proxy_rad_s": mixing_over_reference * reference_rate_s,
        "radical_relaxation_over_reference": radical_relaxation_over_reference,
        "radical_relaxation_s^-1": radical_relaxation_over_reference * reference_rate_s,
        "oxygen_relaxation_over_reference": oxygen_relaxation_over_reference,
        "oxygen_relaxation_s^-1": oxygen_relaxation_over_reference * reference_rate_s,
        "escape_over_reference": escape_over_reference,
        "k_escape_s^-1": escape_over_reference * reference_rate_s,
        "kd_over_reference": kd_over_reference,
        "k_doublet_s^-1": kd,
        "kq_over_kd": kq_over_kd,
        "k_quartet_s^-1": kq_over_kd * kd,
        "final_doublet_population": float(result["p_doublet"][-1]),
        "final_quartet_population": float(result["p_quartet"][-1]),
        "final_total_survival": float(result["survival"][-1]),
        "doublet_reaction_yield_per_encounter": result["doublet_reaction_yield"],
        "quartet_reaction_yield_per_encounter": result["quartet_reaction_yield"],
        "primary_superoxide_yield_per_encounter": result["primary_superoxide_yield"],
        "escape_yield_per_encounter": result["escape_yield"],
        "unresolved_probability": result["unresolved_probability"],
        "probability_balance": result["probability_balance"],
        "probability_balance_error": result["probability_balance_error"],
        "solver": result["numerical_method"],
        "solver_samples": samples,
        "output_class": "dimensionless populations/per-encounter yields",
        "units": "dimensionless per encounter unless field name gives s^-1 or s",
        "non_predictive": True,
        "limitations": LIMITATIONS,
    }


def run_mixing_escape_sweep(
    reference_rate_s: float, grid_size: int, kq_over_kd_values=(0.0, 0.1, 1.0, 10.0)
) -> list[dict]:
    """Bounded two-dimensional mixing/escape sweep for an unpolarized input."""
    axis = normalized_axis(grid_size)
    rows = []
    for ratio in kq_over_kd_values:
        for mixing in axis:
            for escape in axis:
                rows.append(encounter_result_row(
                    scenario_id=f"mix_escape_kqkd_{ratio:g}",
                    reference_rate_s=reference_rate_s,
                    mixing_over_reference=float(mixing),
                    radical_relaxation_over_reference=0.1,
                    oxygen_relaxation_over_reference=0.1,
                    escape_over_reference=float(escape),
                    kq_over_kd=float(ratio),
                ))
    return rows


def run_mixing_relaxation_sweep(
    reference_rate_s: float,
    grid_size: int,
    kq_over_kd_values=(0.0, 0.1, 1.0, 10.0),
) -> list[dict]:
    """Bounded mixing/relaxation sweep across the selectivity controls.

    The encounter starts from the unpolarized benchmark, with kD/kref=1 and
    kescape/kref=1 fixed.  The four kQ/kD panels are sensitivity coordinates;
    equality is the spin-independent null condition.
    """
    axis = normalized_axis(grid_size)
    rows = []
    for ratio in kq_over_kd_values:
        for mixing in axis:
            for relaxation in axis:
                rows.append(encounter_result_row(
                    scenario_id=f"mix_relax_kqkd_{ratio:g}",
                    reference_rate_s=reference_rate_s,
                    mixing_over_reference=float(mixing),
                    radical_relaxation_over_reference=float(relaxation),
                    oxygen_relaxation_over_reference=float(relaxation),
                    escape_over_reference=1.0,
                    kq_over_kd=float(ratio),
                    initial_state="unpolarized",
                ))
    return rows


def run_selectivity_sweep(reference_rate_s: float, grid_size: int) -> list[dict]:
    """Sweep kQ/kD for initial states, mixing strengths, and escape ratios."""
    ratios = selectivity_axis(grid_size)
    rows = []
    for initial_state, p_doublet, label in INITIAL_STATES:
        for mixing in (0.0, 1.0, 10.0):
            for escape in (0.1, 1.0, 10.0):
                for ratio in ratios:
                    rows.append(encounter_result_row(
                        scenario_id=f"selectivity_{label.replace(' ', '_')}",
                        reference_rate_s=reference_rate_s,
                        mixing_over_reference=mixing,
                        radical_relaxation_over_reference=0.1,
                        oxygen_relaxation_over_reference=0.1,
                        escape_over_reference=escape,
                        kq_over_kd=float(ratio),
                        initial_state=initial_state,
                        p_doublet=p_doublet,
                    ))
    return rows


def _metadata(config: Path, provenance: Path) -> dict:
    digest = hashlib.sha256()
    for path in (ROOT / "spin_chemistry.py", Path(__file__), config, provenance):
        try:
            identity = str(path.resolve().relative_to(ROOT))
        except ValueError:
            identity = str(path.resolve())
        digest.update(identity.encode())
        digest.update(path.read_bytes())
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        text=True, capture_output=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"], cwd=ROOT, check=True,
        text=True, capture_output=True,
    ).stdout.splitlines()
    return {
        "commit": commit,
        "dirty_tree": bool(dirty),
        "dirty_entry_count": len(dirty),
        "source_fingerprint_sha256": digest.hexdigest(),
        "config_identity": (
            str(config.resolve().relative_to(ROOT))
            if config.resolve().is_relative_to(ROOT)
            else str(config.resolve())
        ),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "provenance_identity": (
            str(provenance.resolve().relative_to(ROOT))
            if provenance.resolve().is_relative_to(ROOT)
            else str(provenance.resolve())
        ),
        "provenance_sha256": hashlib.sha256(provenance.read_bytes()).hexdigest(),
    }


def _scenarios() -> list[dict[str, float | str]]:
    baseline = {
        "mixing_over_reference": 1.0,
        "relaxation_over_reference": 0.1,
        "escape_over_reference": 1.0,
        "kq_over_kd": 0.1,
    }
    definitions: list[dict[str, float | str]] = [
        {"scenario": "baseline", **baseline}
    ]
    grids = {
        "mixing_over_reference": [0.0, 0.1, 10.0],
        "relaxation_over_reference": [0.0, 1.0, 100.0],
        "escape_over_reference": [0.1, 10.0],
        "kq_over_kd": [0.0, 1.0, 10.0],
    }
    for coordinate, values in grids.items():
        for value in values:
            item = dict(baseline)
            item[coordinate] = value
            item["scenario"] = f"vary_{coordinate}:{value:g}"
            definitions.append(item)
    return definitions


def run_sweep(reference_rate_s: float, samples: int) -> list[dict]:
    if not np.isfinite(reference_rate_s) or reference_rate_s <= 0:
        raise ValueError("reference_rate_s must be finite and positive")
    if isinstance(samples, bool) or not isinstance(samples, int) or samples < 2:
        raise ValueError("samples must be an integer >= 2")
    states = [
        ("unpolarized", None), ("doublet", None), ("quartet", None),
        ("mixture", 0.25), ("mixture", 0.75),
    ]
    rows = []
    for scenario in _scenarios():
        for state, p_doublet in states:
            row = encounter_result_row(
                scenario_id=str(scenario["scenario"]),
                reference_rate_s=reference_rate_s,
                mixing_over_reference=float(scenario["mixing_over_reference"]),
                radical_relaxation_over_reference=float(scenario["relaxation_over_reference"]),
                oxygen_relaxation_over_reference=float(scenario["relaxation_over_reference"]),
                escape_over_reference=float(scenario["escape_over_reference"]),
                kq_over_kd=float(scenario["kq_over_kd"]),
                initial_state=state,
                p_doublet=p_doublet,
                samples=samples,
            )
            row.update({
                "scenario": scenario["scenario"],
                "relaxation_over_reference": scenario["relaxation_over_reference"],
                "samples": samples,
                "survival_probability": row["unresolved_probability"],
            })
            rows.append(row)
    return rows


def _write_svg(rows: list[dict], path: Path) -> None:
    """Write a dependency-free summary figure for the four coordinates."""
    width, height = 960, 720
    panels = [
        "mixing_over_reference", "relaxation_over_reference",
        "escape_over_reference", "kq_over_kd",
    ]
    selected = [row for row in rows if row["initial_state"] == "unpolarized"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<title>Dimensionless encounter sensitivity summary</title>',
        '<desc>Four one-factor plots of unpolarized primary superoxide yield per encounter.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="30" y="30" font-family="sans-serif" font-size="18">'
        'Non-predictive one-factor sensitivity: unpolarized per-encounter yield</text>',
    ]
    for panel_index, coordinate in enumerate(panels):
        x0 = 60 + (panel_index % 2) * 470
        y0 = 70 + (panel_index // 2) * 320
        plot_w, plot_h = 390, 235
        varied = [
            row for row in selected
            if row["scenario"] == "baseline"
            or str(row["scenario"]).startswith(f"vary_{coordinate}:")
        ]
        unique = {}
        for row in varied:
            unique[float(row[coordinate])] = float(
                row["primary_superoxide_yield_per_encounter"]
            )
        points = sorted(unique.items())
        transformed = [np.log10(max(x, 1e-3)) for x, _ in points]
        xmin, xmax = min(transformed), max(transformed)
        if xmax == xmin:
            xmax = xmin + 1
        parts.extend([
            f'<line x1="{x0}" y1="{y0 + plot_h}" x2="{x0 + plot_w}" '
            f'y2="{y0 + plot_h}" stroke="black"/>',
            f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0 + plot_h}" stroke="black"/>',
            f'<text x="{x0}" y="{y0 - 12}" font-family="sans-serif" font-size="14">'
            f'{coordinate} (log axis; zero shown at 10^-3)</text>',
        ])
        coordinates = []
        for transformed_x, (raw_x, y) in zip(transformed, points):
            px = x0 + (transformed_x - xmin) / (xmax - xmin) * plot_w
            py = y0 + (1 - y) * plot_h
            coordinates.append(f"{px:.2f},{py:.2f}")
            parts.append(
                f'<circle cx="{px:.2f}" cy="{py:.2f}" r="4" fill="#2356a8">'
                f'<title>{raw_x:g}: {y:.6g}</title></circle>'
            )
            parts.append(
                f'<text x="{px:.2f}" y="{y0 + plot_h + 18}" text-anchor="middle" '
                f'font-family="sans-serif" font-size="11">{raw_x:g}</text>'
            )
        parts.append(
            f'<polyline points="{" ".join(coordinates)}" fill="none" '
            'stroke="#2356a8" stroke-width="2"/>'
        )
        for tick in (0.0, 0.5, 1.0):
            py = y0 + (1 - tick) * plot_h
            parts.append(
                f'<text x="{x0 - 42}" y="{py + 4}" font-family="sans-serif" '
                f'font-size="11">{tick:.1f}</text>'
            )
        parts.append(
            f'<text x="{x0 - 48}" y="{y0 + plot_h / 2}" text-anchor="middle" '
            f'transform="rotate(-90 {x0 - 48} {y0 + plot_h / 2})" '
            'font-family="sans-serif" font-size="11">Yield / encounter</text>'
        )
    parts.extend([
        '<text x="30" y="700" font-family="sans-serif" font-size="12">'
        'Illustrative bounds, not measured ranges or priors; yields are not concentrations or fluxes.</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-rate-s", type=float, required=True)
    parser.add_argument("--samples", type=int, default=81)
    parser.add_argument(
        "--config", type=Path,
        default=ROOT / "configs" / "doxorubicin_parameters.json",
    )
    parser.add_argument("--provenance", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = args.config.resolve()
    provenance = (
        args.provenance.resolve()
        if args.provenance
        else config.with_name("parameter_provenance.csv")
    )
    validate_authority_bundle(config, provenance)
    metadata = _metadata(config, provenance)
    rows = run_sweep(args.reference_rate_s, args.samples)
    for row in rows:
        row.update({
            "commit": metadata["commit"],
            "dirty_tree": metadata["dirty_tree"],
            "source_fingerprint_sha256": metadata["source_fingerprint_sha256"],
            "config_identity": metadata["config_identity"],
            "config_sha256": metadata["config_sha256"],
            "provenance_identity": metadata["provenance_identity"],
            "provenance_sha256": metadata["provenance_sha256"],
            "limitations": LIMITATIONS,
        })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "dimensionless_sensitivity.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    svg_path = args.output_dir / "dimensionless_sensitivity.svg"
    _write_svg(rows, svg_path)
    summary = {
        "output_label": "NON-PREDICTIVE DIMENSIONLESS SENSITIVITY OUTPUT",
        "scenario": "one-factor-at-a-time dimensionless encounter sweep",
        "normalization": {
            "reference_rate_s^-1": args.reference_rate_s,
            "convention": (
                "kD=k_ref; mixing angular frequency, both local relaxation rates, "
                "and escape are divided by k_ref; kQ/kD is dimensionless; "
                "duration=8/k_ref"
            ),
        },
        "resolved_inputs": {
            "coordinates": _scenarios(),
            "initial_states": ["unpolarized", "doublet", "quartet", "mixture(0.25)", "mixture(0.75)"],
            "zero_mixing_control": True,
            "fast_relaxation_control": "relaxation/reference=100",
            "spin_independent_null": "kQ/kD=1",
        },
        "units": {
            "rates_and_frequencies": "normalized by positive k_ref",
            "duration": "s",
            "outputs": "dimensionless populations and per-encounter yields",
        },
        "numerical_settings": {
            "solver": "six-state constant-generator matrix exponential Pade(13)",
            "samples": args.samples,
            "row_count": len(rows),
        },
        "limitations": LIMITATIONS,
        "metadata": metadata,
        "files": [csv_path.name, svg_path.name],
    }
    summary_path = args.output_dir / "dimensionless_sensitivity_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"csv": str(csv_path), "figure": str(svg_path), "summary": str(summary_path), "rows": len(rows)}))


if __name__ == "__main__":
    main()
