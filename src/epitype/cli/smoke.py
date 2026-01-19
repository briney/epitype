"""Smoke test command for epitype."""
from __future__ import annotations

from importlib.resources import as_file, files
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from epitype.analysis.pipeline import AnalysisConfig, analyze_interface
from epitype.io.writers import format_metrics_table, write_json

console = Console()

# Smoke test configurations
SMOKE_TESTS = [
    {"pdb": "1YY9.pdb", "chains": "CD_A", "name": "1YY9"},
    {"pdb": "4FQI.pdb", "chains": "HL_A", "name": "4FQI"},
]

# Expected metric ranges for validation
EXPECTED_RANGES: dict[str, dict[str, tuple[float, float]]] = {
    "1YY9": {
        "dSASA_int": (2000.0, 4000.0),
        "hbonds_int": (0.0, 50.0),
        "packstat": (0.0, 1.0),
    },
    "4FQI": {
        "dSASA_int": (3000.0, 8000.0),
        "hbonds_int": (0.0, 50.0),
        "packstat": (0.0, 1.0),
    },
}


def _validate_metrics(
    metrics_dict: dict,
    expected: dict[str, tuple[float, float]],
) -> list[tuple[str, float, tuple[float, float], bool]]:
    """
    Validate metrics against expected ranges.

    Returns list of (metric_name, value, expected_range, passed).
    """
    results = []
    for metric_name, (low, high) in expected.items():
        value = metrics_dict.get(metric_name)
        if value is not None:
            passed = low <= value <= high
            results.append((metric_name, value, (low, high), passed))
    return results


def _print_validation(
    validations: list[tuple[str, float, tuple[float, float], bool]],
) -> int:
    """Print validation results and return number of failures."""
    failures = 0
    console.print("\n[bold]Validation:[/bold]")
    for metric_name, value, (low, high), passed in validations:
        if passed:
            console.print(
                f"  [green]\u2713[/green] {metric_name}: {value:.2f} "
                f"(expected: {low:.0f}-{high:.0f})"
            )
        else:
            console.print(
                f"  [red]\u2717[/red] {metric_name}: [red]{value:.2f}[/red] "
                f"(expected: {low:.0f}-{high:.0f}) [red]\u2190 FAIL[/red]"
            )
            failures += 1
    return failures


def smoke_command(
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Directory to write JSON results (optional)",
    ),
    no_energy: bool = typer.Option(
        False,
        "--no-energy",
        help="Skip energy calculation (faster)",
    ),
    no_shape: bool = typer.Option(
        False,
        "--no-shape",
        help="Skip shape complementarity (faster, no NanoShaper required)",
    ),
    no_packstat: bool = typer.Option(
        False,
        "--no-packstat",
        help="Skip PackStat calculation (faster)",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output",
    ),
):
    """
    Run smoke tests on bundled test PDB files.

    Analyzes the test structures (1YY9 and 4FQI) and validates that
    computed metrics fall within expected ranges. Useful for verifying
    that epitype is installed and working correctly.

    Example:
        epitype smoke
        epitype smoke --no-energy --no-shape
        epitype smoke --output-dir ./results
    """
    # Get path to bundled data files
    data_files = files("epitype.data")

    # Track overall results
    results_summary: list[tuple[str, int, int]] = []  # (name, passed, total)

    for test in SMOKE_TESTS:
        pdb_name = test["pdb"]
        chains = test["chains"]
        name = test["name"]

        console.print()
        console.print(f"[bold cyan]\u2550\u2550\u2550 {name} ({chains}) \u2550\u2550\u2550[/bold cyan]")

        # Access the PDB file from package data
        with as_file(data_files.joinpath(pdb_name)) as pdb_path:
            config = AnalysisConfig(
                chain_spec=chains,
                interface_cutoff=8.0,
                compute_energy=not no_energy,
                compute_shape=not no_shape,
                compute_packstat=not no_packstat,
            )

            try:
                if not quiet:
                    with Progress(
                        SpinnerColumn(),
                        TextColumn("[progress.description]{task.description}"),
                        console=console,
                    ) as progress:
                        task = progress.add_task(f"Analyzing {name}...", total=None)

                        def progress_callback(step: str, frac: float):
                            progress.update(task, description=f"{step}...")

                        metrics = analyze_interface(pdb_path, config, progress_callback)
                else:
                    metrics = analyze_interface(pdb_path, config)

            except (ValueError, FileNotFoundError) as e:
                console.print(f"[red]Error analyzing {name}: {e}[/red]")
                results_summary.append((name, 0, 0))
                continue

        # Print metrics table
        console.print()
        console.print(format_metrics_table(metrics, structure_name=name))

        # Validate metrics
        metrics_dict = metrics.to_dict()
        expected = EXPECTED_RANGES.get(name, {})
        validations = _validate_metrics(metrics_dict, expected)
        failures = _print_validation(validations)

        passed = len(validations) - failures
        results_summary.append((name, passed, len(validations)))

        # Write output if requested
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = output_dir / f"{name}_metrics.json"
            write_json(metrics, output_file, structure_name=name)
            console.print(f"\n[dim]Results written to {output_file}[/dim]")

    # Print summary
    console.print()
    console.print("[bold cyan]\u2550\u2550\u2550 Summary \u2550\u2550\u2550[/bold cyan]")

    all_passed = True
    for name, passed, total in results_summary:
        if total == 0:
            console.print(f"  {name}: [red]ERROR[/red]")
            all_passed = False
        elif passed == total:
            console.print(f"  {name}: [green]PASS[/green] ({passed}/{total} metrics)")
        else:
            console.print(f"  {name}: [red]FAIL[/red] ({passed}/{total} metrics)")
            all_passed = False

    console.print()
    if all_passed:
        console.print("[green bold]All smoke tests passed![/green bold]")
    else:
        console.print("[red bold]Some smoke tests failed.[/red bold]")
        raise typer.Exit(1)
