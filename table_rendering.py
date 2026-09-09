"""Render publication tables from their authoritative CSV files.

The renderer intentionally uses a curated presentation schema while preserving
every source row and every source column in the CSV.  Each panel repeats all
rows, long panels are paginated with repeated headings, and the returned
linkage metadata records both displayed and source-only columns.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import textwrap
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.transforms import Bbox


PAGE_WIDTH_IN = 11.0
PAGE_HEIGHT_IN = 8.5
PREVIEW_DPI_MINIMUM = 300
MISSING_DISPLAY = "Not available"


def read_source_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one nonempty CSV and retain its exact column order."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not fields:
        raise ValueError(f"table CSV has no columns: {path}")
    if not rows:
        raise ValueError(f"table CSV has no rows: {path}")
    return fields, rows


def _ascii_hyphens(value: str) -> str:
    return (
        value.replace("\u2010", "-")
        .replace("\u2011", "-")
        .replace("\u2012", "-")
        .replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u2212", "-")
    )


def _format_number(field: str, number: float) -> str:
    if number == 0:
        return "0"
    absolute = abs(number)
    if (
        any(token in field for token in ("probability", "population"))
        or field.startswith("p_doublet")
    ) and absolute < 1e-15:
        return "0 (roundoff)"
    if any(token in field for token in ("error", "tolerance")):
        return f"{number:.3e}"
    if field.endswith("_bytes") or field in {"repetitions", "qubits", "seed"}:
        return str(int(round(number)))
    if "rss_mib" in field or "time_ms" in field:
        return f"{number:.4f}"
    if any(token in field for token in ("yield", "probability", "population")):
        return f"{number:.3e}" if absolute < 1e-5 else f"{number:.6f}"
    if absolute >= 1e4 or absolute < 1e-3:
        return f"{number:.3e}"
    return f"{number:.5g}"


def format_presentation_cell(field: str, value: Any) -> str:
    """Format a displayed cell without altering the source CSV precision."""
    if value is None:
        return MISSING_DISPLAY
    text = _ascii_hyphens(str(value).replace("\r", " ").replace("\n", " ").strip())
    if not text:
        return MISSING_DISPLAY
    if text == "True":
        return "Yes"
    if text == "False":
        return "No"
    if field in {
        "reported_result", "reported_values_and_uncertainties", "does_not_validate",
    } and text[:1] in {"{", "["}:
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            pass
        else:
            if isinstance(parsed, dict):
                text = "; ".join(
                    f"{key.replace('_', ' ')}: "
                    + (", ".join(str(item) for item in item_value)
                       if isinstance(item_value, list) else str(item_value))
                    for key, item_value in parsed.items()
                )
            elif isinstance(parsed, list):
                text = "; ".join(str(item) for item in parsed)
    if field in {
        "status", "evidence_status", "source_category", "mode",
        "problem_group", "phase", "implementation", "disposition",
        "permitted_use", "summary_type",
    }:
        text = text.replace("_", " ")
    try:
        number = float(text)
    except ValueError:
        return text
    return _format_number(field, number)


def _column_weights(
    rows: list[dict[str, str]], fields: list[str], labels: dict[str, str]
) -> list[float]:
    narrative_tokens = (
        "meaning", "definition", "implementation", "consequence", "limitation",
        "conditions", "source", "scope", "evidence", "criticism", "attribution",
        "treatment", "prohibited", "method", "reason", "does_not_validate",
    )
    compact_tokens = (
        "status", "passes", "unit", "symbol", "number", "id", "pH",
        "repetitions", "qubits", "seed", "value", "error", "tolerance",
        "probability", "yield", "time_ms", "rss_mib", "bytes",
    )
    weights: list[float] = []
    for field in fields:
        values = [format_presentation_cell(field, row.get(field)) for row in rows]
        typical = sum(min(len(value), 180) for value in values) / max(len(values), 1)
        heading = labels.get(field, field.replace("_", " "))
        if any(token in field for token in narrative_tokens):
            weight = 2.0 if typical < 70 else 2.7
        elif any(token in field for token in compact_tokens):
            weight = 0.85 if typical < 25 else 1.2
        elif typical > 65:
            weight = 2.2
        elif typical > 30:
            weight = 1.55
        else:
            weight = 1.05
        weights.append(max(weight, min(1.7, len(heading) / 24)))
    return weights


def _wrap(value: str, width_inches: float, *, heading: bool = False) -> str:
    characters_per_inch = 12 if heading else 14
    width = max(8, int(width_inches * characters_per_inch))
    parts = textwrap.wrap(
        value,
        width=width,
        break_long_words=True,
        break_on_hyphens=True,
        replace_whitespace=True,
        drop_whitespace=True,
    )
    return "\n".join(parts) if parts else MISSING_DISPLAY


def build_table_pages(
    rows: list[dict[str, str]],
    spec: dict[str, Any],
    column_labels: dict[str, str],
) -> list[dict[str, Any]]:
    """Create a deterministic pagination plan covering every row in every panel."""
    pages: list[dict[str, Any]] = []
    panel_titles = spec.get("panel_titles", [])
    available_width = 10.15
    available_height = 5.55
    for panel_index, fields in enumerate(spec["panels"]):
        missing = [field for field in fields if field not in rows[0]]
        if missing:
            raise ValueError(
                f"presentation fields missing from table {spec['number']}: {missing}"
            )
        weights = _column_weights(rows, fields, column_labels)
        widths = [available_width * weight / sum(weights) for weight in weights]
        headers = [
            _wrap(column_labels.get(field, field.replace("_", " ").capitalize()), width, heading=True)
            for field, width in zip(fields, widths)
        ]
        formatted_rows = []
        row_heights = []
        for source_index, row in enumerate(rows):
            values = [
                _wrap(format_presentation_cell(field, row.get(field)), width)
                for field, width in zip(fields, widths)
            ]
            line_count = max(value.count("\n") + 1 for value in values)
            height = max(0.31, 0.145 * line_count + 0.10)
            formatted_rows.append((source_index, values))
            row_heights.append(height)

        start = 0
        while start < len(rows):
            used = 0.44
            stop = start
            while stop < len(rows):
                proposed = row_heights[stop]
                if stop > start and used + proposed > available_height:
                    break
                used += proposed
                stop += 1
            page_rows = formatted_rows[start:stop]
            pages.append({
                "panel_index": panel_index,
                "panel_title": (
                    panel_titles[panel_index]
                    if panel_index < len(panel_titles)
                    else spec["title"]
                ),
                "fields": list(fields),
                "headers": headers,
                "widths_in": widths,
                "source_row_indices": [item[0] for item in page_rows],
                "cell_rows": [item[1] for item in page_rows],
                "row_heights_in": [0.44, *row_heights[start:stop]],
            })
            start = stop
    return pages


def _draw_page(
    page: dict[str, Any],
    spec: dict[str, Any],
    source_csv: Path,
    page_number: int,
    total_pages: int,
) -> plt.Figure:
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, PAGE_HEIGHT_IN), facecolor="white")
    axis = fig.add_axes([0, 0, 1, 1])
    axis.axis("off")
    panel_letter = chr(65 + page["panel_index"])
    continuation = " (continued)" if page["source_row_indices"][0] else ""
    fig.text(
        0.035,
        0.952,
        f"Table {spec['number']}. {spec['title']}",
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
        color="#111111",
    )
    fig.text(
        0.035,
        0.910,
        f"Panel {panel_letter}. {page['panel_title']}{continuation}",
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
        color="#222222",
    )

    total_height = sum(page["row_heights_in"])
    height = min(0.68, total_height / PAGE_HEIGHT_IN)
    top = 0.865
    bbox = [0.035, top - height, 0.93, height]
    table = axis.table(
        cellText=page["cell_rows"],
        colLabels=page["headers"],
        colWidths=[width / sum(page["widths_in"]) for width in page["widths_in"]],
        cellLoc="left",
        loc="upper left",
        bbox=bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.0)
    height_scale = bbox[3] / sum(page["row_heights_in"])
    narrative_columns = {
        index for index, field in enumerate(page["fields"])
        if any(token in field for token in (
            "meaning", "definition", "implementation", "consequence", "limitation",
            "conditions", "source", "scope", "evidence", "criticism", "attribution",
            "treatment", "prohibited", "method", "reason", "does_not_validate",
        ))
    }
    for (row_index, column_index), cell in table.get_celld().items():
        cell.set_edgecolor("#C9CED4")
        cell.set_linewidth(0.45)
        cell.PAD = 0.035
        cell.set_height(page["row_heights_in"][row_index] * height_scale)
        text = cell.get_text()
        text.set_fontfamily("DejaVu Sans")
        text.set_verticalalignment("center")
        if row_index == 0:
            cell.set_facecolor("#244A6B")
            text.set_color("white")
            text.set_fontweight("bold")
            text.set_ha("center")
            text.set_fontsize(7.8)
        else:
            cell.set_facecolor("#F2F6F9" if row_index % 2 == 0 else "white")
            text.set_color("#111111")
            text.set_ha("left" if column_index in narrative_columns else "center")
            text.set_fontsize(8.0)

    source_only = [
        field for field in spec.get("source_columns", [])
        if field not in {item for panel in spec["panels"] for item in panel}
    ]
    note = _ascii_hyphens(spec["note"])
    completeness = (
        f"Presentation values are rounded. All {spec['source_row_count']} source rows and "
        f"all {spec['source_column_count']} source columns remain in {source_csv.name}; "
        f"{len(source_only)} provenance/support columns are source-only. "
        f"'{MISSING_DISPLAY}' is distinct from numeric zero."
    )
    footer_lines = textwrap.wrap(note + " " + completeness, width=190)
    fig.text(
        0.035,
        0.095,
        "\n".join(footer_lines),
        ha="left",
        va="bottom",
        fontsize=7.1,
        color="#333333",
        linespacing=1.25,
    )
    fig.text(
        0.965,
        0.035,
        f"Page {page_number} of {total_pages}",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#444444",
    )
    return fig


def render_table_from_csv(
    source_csv: Path,
    output_pdf: Path,
    spec: dict[str, Any],
    column_labels: dict[str, str],
    *,
    requested_dpi: int,
) -> dict[str, Any]:
    """Render one CSV into one multipage PDF and numbered PNG previews."""
    source_fields, rows = read_source_csv(source_csv)
    displayed_fields = [field for panel in spec["panels"] for field in panel]
    spec_for_render = dict(spec)
    spec_for_render.update({
        "source_columns": source_fields,
        "source_row_count": len(rows),
        "source_column_count": len(source_fields),
    })
    pages = build_table_pages(rows, spec_for_render, column_labels)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    preview_dpi = max(int(requested_dpi), PREVIEW_DPI_MINIMUM)
    preview_paths: list[Path] = []
    page_bbox = Bbox.from_bounds(0, 0, PAGE_WIDTH_IN, PAGE_HEIGHT_IN)
    with PdfPages(output_pdf) as pdf:
        for page_number, page in enumerate(pages, start=1):
            figure = _draw_page(
                page, spec_for_render, source_csv, page_number, len(pages)
            )
            pdf.savefig(
                figure, facecolor="white", bbox_inches=page_bbox, pad_inches=0
            )
            preview = output_pdf.with_name(
                f"{output_pdf.stem}_page_{page_number:02d}.png"
            )
            figure.savefig(
                preview,
                dpi=preview_dpi,
                facecolor="white",
                bbox_inches=page_bbox,
                pad_inches=0,
            )
            preview_paths.append(preview)
            plt.close(figure)
    return {
        "source_csv": str(source_csv.name),
        "source_row_count": len(rows),
        "source_column_count": len(source_fields),
        "presentation_columns": list(dict.fromkeys(displayed_fields)),
        "source_only_columns": [
            field for field in source_fields if field not in displayed_fields
        ],
        "presented_rows_per_panel": {
            chr(65 + index): len(rows) for index in range(len(spec["panels"]))
        },
        "pdf": output_pdf.name,
        "png_previews": [path.name for path in preview_paths],
        "preview_dpi": preview_dpi,
        "rendered_page_count": len(pages),
    }
