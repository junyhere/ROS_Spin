"""Shared, deterministic plotting for the active ROS_Spin paper pipeline."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ros_spin_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np


COLORS = {
    "doublet": "#0072B2",
    "quartet": "#E69F00",
    "survival": "#222222",
    "reaction": "#009E73",
    "superoxide": "#CC79A7",
    "escape": "#D55E00",
    "unresolved": "#999999",
    "ho2": "#D55E00",
    "o2minus": "#0072B2",
    "h2o2": "#009E73",
    "loss": "#CC79A7",
}

STATE_COLORS = {
    "unpolarized I6/6": "#222222",
    "doublet-manifold mixture PD/2": "#0072B2",
    "quartet-manifold mixture PQ/4": "#E69F00",
    "manifold mixture pD=0.25": "#009E73",
    "manifold mixture pD=0.75": "#CC79A7",
}

SENSITIVITY_WARNING = "NON-PREDICTIVE DIMENSIONLESS SENSITIVITY ANALYSIS"


def _add_sensitivity_warning(fig, *, y: float = 0.006) -> None:
    """Add the required visible warning to a sensitivity-result figure."""
    fig.text(
        0.5,
        y,
        SENSITIVITY_WARNING,
        ha="center",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        color="#8B1A1A",
    )


def configure_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "svg.hashsalt": "ROS_Spin",
        "pdf.compression": 9,
    })


def save_figure(fig, stem: Path, formats: Iterable[str], dpi: int) -> list[Path]:
    """Save one figure in requested formats and close it unconditionally."""
    written = []
    try:
        for extension in formats:
            path = stem.with_suffix(f".{extension}")
            fig.savefig(path, dpi=dpi if extension == "png" else None)
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"figure was not written correctly: {path}")
            written.append(path)
    finally:
        plt.close(fig)
    return written


def plot_model_overview(rows: list[dict]):
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.set_xlim(0, len(rows) + 0.2)
    ax.set_ylim(-0.6, 2.1)
    ax.axis("off")
    boundary_colors = {"upstream": "#E1F5FE", "encounter": "#FFF3E0", "downstream": "#E8F5E9"}
    for index, row in enumerate(rows):
        x = index + 0.15
        box = FancyBboxPatch(
            (x, 0.35), 0.83, 0.75,
            boxstyle="round,pad=0.04,rounding_size=0.05",
            edgecolor="#444444", facecolor=boundary_colors[row["stage"]], linewidth=1.1,
        )
        ax.add_patch(box)
        ax.text(x + 0.415, 0.73, row["label"], ha="center", va="center", fontsize=8.2, wrap=True)
        ax.text(x + 0.415, 0.1, row["stage"], ha="center", va="center", fontsize=7.5, color="#555555")
        if index < len(rows) - 1:
            ax.add_patch(FancyArrowPatch((x + 0.85, 0.73), (x + 1.12, 0.73), arrowstyle="->", mutation_scale=12, color="#555555"))
    ax.text(
        0.2, 1.62,
        "Preparation is an upstream condition; coherence and D/Q preparation are not experimentally established",
        ha="left", va="center", fontsize=9.5, fontweight="bold",
    )
    ax.text(
        0.2, -0.35,
        "Conceptual workflow created by the authors; not a simulation output.",
        ha="left", va="center", fontsize=8.5, color="#555555",
    )
    return fig


def plot_benchmark_trajectories(rows: list[dict]):
    labels = []
    for row in rows:
        if row["initial_state_definition"] not in labels:
            labels.append(row["initial_state_definition"])
    fig, axes = plt.subplots(3, 2, figsize=(11, 11), sharex=True, sharey=True)
    axes = axes.flat
    for ax, label in zip(axes, labels):
        subset = [r for r in rows if r["initial_state_definition"] == label]
        x = np.array([float(r["normalized_time"]) for r in subset])
        series = (
            ("surviving_doublet_population", "Surviving D", COLORS["doublet"], "-"),
            ("surviving_quartet_population", "Surviving Q", COLORS["quartet"], "-"),
            ("total_survival", "Total survival", COLORS["survival"], "-"),
            ("cumulative_doublet_reaction_yield", "D reaction", COLORS["doublet"], "--"),
            ("cumulative_quartet_reaction_yield", "Q reaction", COLORS["quartet"], "--"),
            ("cumulative_escape_yield", "Escape", COLORS["escape"], ":"),
            ("cumulative_primary_superoxide_yield", "Primary O$_2^{\u2022-}$", COLORS["superoxide"], "-."),
        )
        for key, name, color, linestyle in series:
            ax.plot(x, [float(r[key]) for r in subset], label=name, color=color, linestyle=linestyle, linewidth=1.5)
        ax.set_title(label)
        ax.set_ylim(-0.025, 1.025)
        ax.grid(alpha=0.2)
    axes[-1].axis("off")
    handles, legend_labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(
        handles, legend_labels, loc="upper center", frameon=False, ncol=2,
        bbox_to_anchor=(0.5, 0.88),
    )
    axes[-1].text(
        0.5, 0.25,
        r"$Y_D + Y_Q + Y_{escape} + P_{survival} = 1$",
        transform=axes[-1].transAxes, ha="center", va="center",
        fontsize=10, fontweight="bold",
    )
    for ax in axes[:4]:
        ax.set_ylabel("Probability or cumulative yield")
    for ax in axes[4:]:
        ax.set_xlabel(r"Normalized time, $t k_{ref}$")
    fig.suptitle("Illustrative benchmark encounter trajectories", y=0.995, fontsize=12)
    _add_sensitivity_warning(fig)
    fig.tight_layout(rect=(0, 0.025, 1, 0.98))
    return fig


def plot_initial_state_comparison(rows: list[dict]):
    labels = [r["initial_state_definition"] for r in rows]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    fields = [
        ("doublet_reaction_yield_per_encounter", "D reaction", COLORS["doublet"]),
        ("quartet_reaction_yield_per_encounter", "Q reaction", COLORS["quartet"]),
        ("escape_yield_per_encounter", "Escape", COLORS["escape"]),
        ("unresolved_probability", "Unresolved", COLORS["unresolved"]),
    ]
    bottom = np.zeros(len(rows))
    for field, label, color in fields:
        values = np.array([float(r[field]) for r in rows])
        ax.bar(x, values, bottom=bottom, label=label, color=color, width=0.7)
        bottom += values
    primary = [float(r["primary_superoxide_yield_per_encounter"]) for r in rows]
    ax.scatter(x, primary, color=COLORS["superoxide"], marker="D", s=35, label="Total primary O$_2^{\u2022-}$")
    ax.set_xticks(x, labels, rotation=18, ha="right")
    ax.set_ylabel("Final probability or yield per encounter")
    ax.set_ylim(0, 1.08)
    ax.legend(ncol=5, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    ax.grid(axis="y", alpha=0.2)
    _add_sensitivity_warning(fig)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    return fig


def _heatmap_grid(rows: list[dict], x_field: str, y_field: str):
    xs = np.array(sorted({float(r[x_field]) for r in rows}))
    ys = np.array(sorted({float(r[y_field]) for r in rows}))
    lookup = {(float(r[x_field]), float(r[y_field])): float(r["primary_superoxide_yield_per_encounter"]) for r in rows}
    z = np.array([[lookup[(x, y)] for x in xs] for y in ys])
    return xs, ys, z


def _set_normalized_symlog(ax, x_label: str, y_label: str) -> None:
    ax.set_xscale("symlog", linthresh=1e-2, linscale=0.5)
    ax.set_yscale("symlog", linthresh=1e-2, linscale=0.5)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def plot_mixing_escape_heatmaps(rows: list[dict]):
    ratios = sorted({float(r["kq_over_kd"]) for r in rows})
    vmin = min(float(r["primary_superoxide_yield_per_encounter"]) for r in rows)
    vmax = max(float(r["primary_superoxide_yield_per_encounter"]) for r in rows)
    fig, axes = plt.subplots(
        2, 2, figsize=(10.5, 8.5), sharex=True, sharey=True,
        layout="constrained",
    )
    mesh = None
    for ax, ratio in zip(axes.flat, ratios):
        subset = [r for r in rows if float(r["kq_over_kd"]) == ratio]
        xs, ys, z = _heatmap_grid(subset, "mixing_over_reference", "escape_over_reference")
        mesh = ax.pcolormesh(xs, ys, z, shading="nearest", cmap="viridis", vmin=vmin, vmax=vmax, rasterized=True)
        _set_normalized_symlog(ax, r"$\omega_{local}/k_{ref}$", r"$k_{escape}/k_{ref}$")
        ax.set_title(rf"$k_Q/k_D={ratio:g}$" + (" (null)" if ratio == 1 else ""))
    fig.colorbar(
        mesh, ax=axes.ravel().tolist(),
        label="Primary-superoxide yield per encounter", shrink=0.82, pad=0.03,
    )
    _add_sensitivity_warning(fig, y=-0.022)
    return fig


def plot_mixing_relaxation_heatmaps(rows: list[dict]):
    ratios = sorted({float(r["kq_over_kd"]) for r in rows})
    vmin = min(float(r["primary_superoxide_yield_per_encounter"]) for r in rows)
    vmax = max(float(r["primary_superoxide_yield_per_encounter"]) for r in rows)
    fig, axes = plt.subplots(
        2, 2, figsize=(10.5, 8.5), sharex=True, sharey=True,
        layout="constrained",
    )
    mesh = None
    for ax, ratio in zip(axes.flat, ratios):
        subset = [r for r in rows if float(r["kq_over_kd"]) == ratio]
        xs, ys, z = _heatmap_grid(subset, "mixing_over_reference", "radical_relaxation_over_reference")
        mesh = ax.pcolormesh(xs, ys, z, shading="nearest", cmap="cividis", vmin=vmin, vmax=vmax, rasterized=True)
        _set_normalized_symlog(ax, r"$\omega_{local}/k_{ref}$", r"$\gamma_R/k_{ref}=\gamma_O/k_{ref}$")
        ax.set_title(rf"$k_Q/k_D={ratio:g}$" + (" (null)" if ratio == 1 else ""))
    fig.colorbar(
        mesh, ax=axes.ravel().tolist(),
        label="Primary-superoxide yield per encounter", shrink=0.85, pad=0.03,
    )
    _add_sensitivity_warning(fig, y=-0.022)
    return fig


def plot_reaction_selectivity(rows: list[dict]):
    mixings = sorted({float(r["mixing_over_reference"]) for r in rows})
    escapes = sorted({float(r["escape_over_reference"]) for r in rows})
    states = []
    for row in rows:
        if row["initial_state_definition"] not in states:
            states.append(row["initial_state_definition"])
    linestyles = {0.1: "-", 1.0: "--", 10.0: ":"}
    fig, axes = plt.subplots(len(mixings), len(escapes), figsize=(12, 9), sharex=True, sharey=True)
    for i, mixing in enumerate(mixings):
        for j, escape in enumerate(escapes):
            ax = axes[i, j]
            subset = [r for r in rows if float(r["mixing_over_reference"]) == mixing and float(r["escape_over_reference"]) == escape]
            for state in states:
                series = [r for r in subset if r["initial_state_definition"] == state]
                series.sort(key=lambda r: float(r["kq_over_kd"]))
                x = np.array([float(r["kq_over_kd"]) for r in series])
                y = np.array([float(r["primary_superoxide_yield_per_encounter"]) for r in series])
                ax.plot(x, y, color=STATE_COLORS[state], linewidth=1.4, label=state)
            ax.axvline(1, color="#777777", linestyle="--", linewidth=1)
            ax.set_xscale("symlog", linthresh=1e-2, linscale=0.5)
            ax.grid(alpha=0.18)
            ax.set_title(rf"$\omega/k_{{ref}}={mixing:g}$; $k_{{esc}}/k_{{ref}}={escape:g}$", fontsize=8.5)
            if i == len(mixings) - 1:
                ax.set_xlabel(r"$k_Q/k_D$ (1 = null)")
            if j == 0:
                ax.set_ylabel("Primary O$_2^{\u2022-}$ yield / encounter")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.035))
    _add_sensitivity_warning(fig, y=-0.006)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    return fig


def plot_controls(rows: list[dict]):
    labels = [r["control_label"] for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.4), gridspec_kw={"height_ratios": [2.2, 1]})
    bottom = np.zeros(len(rows))
    for field, label, color in (
        ("primary_superoxide_yield_per_encounter", "Primary O$_2^{\u2022-}$", COLORS["reaction"]),
        ("escape_yield_per_encounter", "Escape", COLORS["escape"]),
        ("unresolved_probability", "Unresolved", COLORS["unresolved"]),
    ):
        values = np.array([float(r[field]) for r in rows])
        axes[0].bar(x, values, bottom=bottom, color=color, label=label)
        bottom += values
    axes[0].set_ylim(0, 1.05)
    axes[0].set_xticks([])
    axes[0].set_ylabel("Final probability")
    axes[0].legend(
        ncol=3,
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        borderaxespad=0.0,
    )
    axes[0].grid(axis="y", alpha=0.2)
    commutators = np.array([float(r["dq_commutator_frobenius"]) for r in rows])
    axes[1].bar(x, commutators, color=[COLORS["doublet"] if v > 1e-12 else "#BBBBBB" for v in commutators])
    axes[1].axhline(1e-12, color="#777777", linestyle="--", linewidth=1)
    axes[1].set_yscale("symlog", linthresh=1e-14)
    axes[1].set_ylabel(r"$\|[H,P_D]\|_F/k_{ref}$")
    axes[1].set_xticks(x, labels, rotation=24, ha="right")
    axes[1].grid(axis="y", alpha=0.2)
    _add_sensitivity_warning(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    return fig


def plot_solver_validation(rows: list[dict], tolerances: dict):
    x = np.array([float(r["normalized_time"]) for r in rows])
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), sharex=True)
    floor = 1e-17
    axes[0, 0].plot(x, np.maximum([float(r["max_density_matrix_abs_error"]) for r in rows], floor), color="#222222")
    axes[0, 0].axhline(tolerances["density"], color="#D55E00", linestyle="--", label="declared tolerance")
    axes[0, 0].set_ylabel("Max |density-matrix difference|")
    axes[0, 0].legend(frameon=False)
    for key, label, color in (
        ("doublet_population_abs_difference", "D population", COLORS["doublet"]),
        ("quartet_population_abs_difference", "Q population", COLORS["quartet"]),
        ("survival_abs_difference", "Survival", COLORS["survival"]),
    ):
        axes[0, 1].plot(x, np.maximum([float(r[key]) for r in rows], floor), label=label, color=color)
    axes[0, 1].axhline(tolerances["observable"], color="#D55E00", linestyle="--")
    axes[0, 1].set_ylabel("Population absolute difference")
    axes[0, 1].legend(
        frameon=True, facecolor="white", framealpha=1.0, edgecolor="none"
    )
    for key, label, color in (
        ("doublet_reaction_yield_abs_difference", "D reaction", COLORS["doublet"]),
        ("quartet_reaction_yield_abs_difference", "Q reaction", COLORS["quartet"]),
        ("primary_superoxide_yield_abs_difference", "Total primary reaction", COLORS["superoxide"]),
        ("escape_yield_abs_difference", "Escape", COLORS["escape"]),
    ):
        axes[1, 0].plot(x, np.maximum([float(r[key]) for r in rows], floor), label=label, color=color)
    axes[1, 0].axhline(tolerances["observable"], color="#D55E00", linestyle="--")
    axes[1, 0].set_ylabel("Cumulative-yield absolute difference")
    axes[1, 0].legend(frameon=False)
    axes[1, 1].plot(x, np.maximum([float(r["reference_probability_accounting_error"]) for r in rows], floor), label="Matrix exponential", color="#0072B2")
    axes[1, 1].plot(x, np.maximum([float(r["independent_probability_accounting_error"]) for r in rows], floor), label="Dormand–Prince", color="#E69F00")
    axes[1, 1].axhline(tolerances["balance"], color="#D55E00", linestyle="--")
    axes[1, 1].set_ylabel("Probability-accounting error")
    axes[1, 1].legend(frameon=False)
    for ax in axes.flat:
        ax.set_yscale("log")
        ax.grid(alpha=0.2)
    for ax in axes[1]:
        ax.set_xlabel(r"Normalized time, $t k_{ref}$")
    _add_sensitivity_warning(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    return fig


def plot_downstream(rows: list[dict]):
    scenarios = []
    for row in rows:
        if row["scenario_id"] not in scenarios:
            scenarios.append(row["scenario_id"])
    fig, axes = plt.subplots(1, len(scenarios), figsize=(13, 4.3), sharex=True, sharey=True)
    for ax, scenario in zip(axes, scenarios):
        subset = [r for r in rows if r["scenario_id"] == scenario]
        x = np.array([float(r["time_s"]) * 1e3 for r in subset])
        for key, label, color, style in (
            ("radical_pool_m", "Total radical pool", COLORS["survival"], "-"),
            ("ho2_m", "HO$_2^{\u2022}$", COLORS["ho2"], "--"),
            ("o2minus_m", "O$_2^{\u2022-}$", COLORS["o2minus"], "--"),
            ("hydrogen_peroxide_m", "H$_2$O$_2$", COLORS["h2o2"], "-"),
            ("accumulated_hydrogen_peroxide_loss_m", "Accumulated H$_2$O$_2$ loss", COLORS["loss"], ":"),
        ):
            ax.plot(x, np.array([float(r[key]) for r in subset]) * 1e6, label=label, color=color, linestyle=style)
        ax.set_title(subset[0]["scenario_label"])
        ax.set_xlabel("Time (ms)")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Concentration (µM; illustrative pulse)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.04))
    _add_sensitivity_warning(fig, y=-0.006)
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    return fig


def plot_bulk_rates(rows: list[dict]):
    import textwrap

    fig, (ax, detail) = plt.subplots(
        1, 2, figsize=(14, 5.4), gridspec_kw={"width_ratios": [1.0, 1.65]},
        layout="constrained",
    )
    y = np.arange(len(rows))
    values = np.array([float(r["value_m^-1_s^-1"]) for r in rows])
    errors = np.array([float(r["uncertainty_m^-1_s^-1"]) for r in rows])
    labels = [
        f"{'1981' if r['pH'] == 'unknown' else 'pH ' + str(r['pH'])}\n"
        f"{r['temperature']}"
        for r in rows
    ]
    complete = np.array([
        r["pH"] != "unknown" and r["temperature"] != "unknown"
        and str(r.get("condition_matched_allowed", "False")) == "True"
        for r in rows
    ])
    for mask, label, color, marker, fill in (
        (complete, "Complete matched-condition record", "#0072B2", "o", "#0072B2"),
        (~complete, "Incomplete or validation-only conditions", "#D55E00", "^", "white"),
    ):
        if np.any(mask):
            ax.errorbar(
                values[mask], y[mask], xerr=errors[mask], fmt=marker,
                color=color, markerfacecolor=fill, markeredgecolor=color,
                ecolor="#777777", capsize=4, label=label,
            )
    ax.set_xscale("log")
    ax.set_yticks(y, labels)
    ax.set_xlabel(r"Documented bulk total rate constant (M$^{-1}$ s$^{-1}$)")
    ax.set_title("Condition-specific literature records (not a universal fitted rate)")
    ax.grid(axis="x", which="both", alpha=0.2)
    ax.legend(frameon=False, fontsize=7.2, loc="best")
    detail.axis("off")
    table_rows = []
    for row in rows:
        condition = f"pH {row['pH']}; {row['temperature']}; {row['environment']}"
        chemistry = f"{row['species']}; protonation: {row['protonation']}"
        source_use = f"{row['source']}; allowed: {row['allowed_use']}"
        table_rows.append([
            row["parameter_id"].replace("doxorubicin_semiquinone_plus_oxygen_", ""),
            textwrap.fill(condition, 34),
            textwrap.fill(chemistry, 38),
            textwrap.fill(source_use, 38),
        ])
    table = detail.table(
        cellText=table_rows,
        colLabels=["Record", "Conditions/environment", "Species/protonation", "Source/allowed use"],
        cellLoc="left", colLoc="left", loc="center",
        colWidths=[0.12, 0.27, 0.30, 0.31],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.2)
    # The longest source entry spans seven lines.  Give every row enough
    # vertical room so text cannot cross the horizontal cell rules in any
    # output backend.
    table.scale(1, 3.5)
    detail.set_title("Records remain conditionally distinct", fontsize=10)
    return fig


def plot_circuit_validation(rows: list[dict]):
    fig, ax = plt.subplots(figsize=(9, 4.8))
    numeric = [r for r in rows if r.get("value") not in (None, "")]
    labels = [r["metric"].replace("_", "\n") for r in numeric]
    values = np.array([float(r["value"]) for r in numeric])
    tolerances = np.array([float(r["tolerance"]) for r in numeric])
    x = np.arange(len(numeric))
    display = np.maximum(values, 1e-18)
    ax.bar(x, display, color="#0072B2", label="Observed error (zeros shown at 10$^{-18}$)")
    ax.scatter(x, tolerances, color="#D55E00", marker="_", s=500, linewidths=2, label="Tolerance")
    ax.set_yscale("log")
    ax.set_xticks(x, labels, rotation=22, ha="right")
    ax.set_ylabel("Absolute error or probability")
    ax.set_title(f"Three-qubit dense-unitary consistency: {rows[0]['status']}")
    ax.legend(frameon=False)
    ax.grid(axis="y", which="both", alpha=0.2)
    fig.tight_layout()
    return fig


configure_style()
