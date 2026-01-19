"""Output writers for analysis results."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from epitype.core.types import InterfaceMetrics


def write_json(
    metrics: InterfaceMetrics,
    filepath: str | Path,
    structure_name: str | None = None,
    indent: int = 2,
) -> None:
    """
    Write metrics to JSON file.

    Args:
        metrics: Analysis results
        filepath: Output file path
        structure_name: Optional structure identifier
        indent: JSON indentation
    """
    filepath = Path(filepath)

    data: dict[str, Any] = {
        "structure": structure_name,
        "metrics": metrics.to_dict(),
    }

    with open(filepath, "w") as f:
        json.dump(data, f, indent=indent)


def write_csv(
    metrics: InterfaceMetrics,
    filepath: str | Path,
    structure_name: str | None = None,
    append: bool = False,
) -> None:
    """
    Write metrics to CSV file.

    Args:
        metrics: Analysis results
        filepath: Output file path
        structure_name: Optional structure identifier
        append: If True, append to existing file
    """
    filepath = Path(filepath)

    fieldnames = [
        "structure",
        "dG_separated",
        "dSASA_int",
        "dG_per_dSASA",
        "sc_value",
        "packstat",
        "delta_unsat_hbonds",
        "hbonds_int",
        "sasa_complex",
        "sasa_separated",
    ]

    row = {
        "structure": structure_name or "",
        **{k: v for k, v in metrics.to_dict().items() if k in fieldnames},
    }

    mode = "a" if append and filepath.exists() else "w"
    write_header = not (append and filepath.exists())

    with open(filepath, mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def format_metrics_table(
    metrics: InterfaceMetrics,
    structure_name: str | None = None,
) -> str:
    """
    Format metrics as human-readable table.

    Args:
        metrics: Analysis results
        structure_name: Optional structure identifier

    Returns:
        Formatted string
    """
    lines = []

    if structure_name:
        lines.append(f"Structure: {structure_name}")
        lines.append("=" * 50)

    lines.append(f"{'Metric':<25} {'Value':>15} {'Quality':>10}")
    lines.append("-" * 50)

    # dG_separated
    quality = "good" if metrics.dG_separated < 0 else "poor"
    lines.append(f"{'dG_separated (REU)':<25} {metrics.dG_separated:>15.2f} {quality:>10}")

    # dSASA
    lines.append(f"{'dSASA_int (A^2)':<25} {metrics.dSASA_int:>15.1f}")

    # dG/dSASA
    if metrics.dG_per_dSASA < -1.5:
        quality = "good"
    elif metrics.dG_per_dSASA < 0:
        quality = "moderate"
    else:
        quality = "poor"
    lines.append(f"{'dG/dSASA (REU/A^2)':<25} {metrics.dG_per_dSASA:>15.3f} {quality:>10}")

    # Shape complementarity
    if metrics.sc_value > 0.65:
        quality = "good"
    elif metrics.sc_value > 0.5:
        quality = "moderate"
    else:
        quality = "poor"
    lines.append(f"{'Shape complementarity':<25} {metrics.sc_value:>15.3f} {quality:>10}")

    # Packstat
    if metrics.packstat > 0.65:
        quality = "good"
    elif metrics.packstat > 0.5:
        quality = "moderate"
    else:
        quality = "poor"
    lines.append(f"{'Packstat':<25} {metrics.packstat:>15.3f} {quality:>10}")

    # H-bonds
    lines.append(f"{'Unsat. H-bonds':<25} {metrics.delta_unsat_hbonds:>15d}")
    lines.append(f"{'Interface H-bonds':<25} {metrics.hbonds_int:>15d}")

    return "\n".join(lines)
