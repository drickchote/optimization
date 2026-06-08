#!/usr/bin/env python3
"""Compute hypervolume indicators for algorithm runs in src/runs/."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable

# Table 2 from Anagnostopoulos & Mamanis (2011), MVCCPO paper.
# Objective 1: variance (risk), minimize.
# Objective 2: expected return, maximize.
OBJECTIVE_BOUNDS = {
    "port1": (0.000578, 0.005253, 0.00234, 0.011950),
    "port2": (0.000130, 0.003120, 0.00140, 0.010800),
    "port3": (0.000185, 0.001668, 0.00211, 0.009030),
    "port4": (0.000120, 0.003233, 0.00156, 0.010000),
    "port5": (0.000270, 0.001800, -0.00034, 0.004370),
}

REFERENCE_POINT = (1.0, 0.0)

ALGORITHM_LABELS = {
    "nsgaii": "nsgaii",
    "bb": "branch_and_bound",
    "branch_and_bound": "branch_and_bound",
}

RUN_FILENAME_PATTERN = re.compile(r"^(port\d+)_k(\d+)\.csv$", re.IGNORECASE)
HEADER_PATTERN = re.compile(
    r"^Port\s+(\d+)\s+-\s+K\s*=\s*(\d+)\s+-\s+(.+)$",
    re.IGNORECASE,
)


def script_runs_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "runs"


def normalize_objective(value: float, minimum: float, maximum: float) -> float:
    return (value - minimum) / (maximum - minimum)


def normalize_point(
    risk: float,
    expected_return: float,
    bounds: tuple[float, float, float, float],
) -> tuple[float, float]:
    f_min_1, f_max_1, f_min_2, f_max_2 = bounds
    risk_norm = normalize_objective(risk, f_min_1, f_max_1)
    return_norm = normalize_objective(expected_return, f_min_2, f_max_2)
    return risk_norm, return_norm


def to_minimization(point: tuple[float, float]) -> tuple[float, float]:
    risk_norm, return_norm = point
    return risk_norm, -return_norm


def nondominated_min(points: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    unique_points = list(dict.fromkeys(points))
    nondominated: list[tuple[float, float]] = []

    for x1, y1 in unique_points:
        dominated = False
        for x2, y2 in unique_points:
            if x2 <= x1 and y2 <= y1 and (x2 < x1 or y2 < y1):
                dominated = True
                break
        if not dominated:
            nondominated.append((x1, y1))

    return nondominated


def hypervolume_2d(
    points: Iterable[tuple[float, float]],
    reference: tuple[float, float] = REFERENCE_POINT,
) -> float:
    nondominated = nondominated_min(points)
    if not nondominated:
        return 0.0

    nondominated.sort(key=lambda point: (point[0], point[1]))

    filtered: list[tuple[float, float]] = []
    best_y = float("inf")
    for x, y in nondominated:
        if y < best_y:
            filtered.append((x, y))
            best_y = y

    hypervolume = 0.0
    for index, (x, y) in enumerate(filtered):
        x_next = filtered[index + 1][0] if index + 1 < len(filtered) else reference[0]
        if x_next <= x or y > reference[1]:
            continue
        hypervolume += (x_next - x) * (reference[1] - y)

    return hypervolume


def portfolio_key_from_filename(filename: str) -> str | None:
    match = RUN_FILENAME_PATTERN.match(filename)
    if not match:
        return None
    return match.group(1).lower()


def run_key_from_filename(filename: str) -> str | None:
    match = RUN_FILENAME_PATTERN.match(filename)
    if not match:
        return None
    return f"{match.group(1).lower()}_k{match.group(2)}"


def algorithm_label(directory_name: str) -> str:
    return ALGORITHM_LABELS.get(directory_name.lower(), directory_name.lower())


def parse_run_file(path: Path) -> tuple[list[tuple[float, float]], str | None]:
    points: list[tuple[float, float]] = []
    portfolio_key: str | None = portfolio_key_from_filename(path.name)

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle):
            line = raw_line.strip()
            if not line:
                continue

            if line_number == 0:
                header_match = HEADER_PATTERN.match(line)
                if header_match:
                    portfolio_key = f"port{header_match.group(1)}"
                continue

            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue

            risk = float(parts[0])
            expected_return = float(parts[1])
            points.append((risk, expected_return))

    return points, portfolio_key


def compute_hypervolume_for_run(path: Path) -> float:
    points, portfolio_key = parse_run_file(path)
    if not points:
        return 0.0

    if portfolio_key is None:
        raise ValueError(f"Could not infer portfolio for run file: {path}")

    bounds = OBJECTIVE_BOUNDS.get(portfolio_key)
    if bounds is None:
        raise ValueError(
            f"No objective bounds configured for {portfolio_key}. "
            "Add Table 2 values from the MVCCPO paper to OBJECTIVE_BOUNDS."
        )

    normalized_min_points = [
        to_minimization(normalize_point(risk, expected_return, bounds))
        for risk, expected_return in points
    ]

    return hypervolume_2d(normalized_min_points, REFERENCE_POINT)


def collect_run_files(runs_dir: Path) -> dict[str, dict[str, Path]]:
    runs: dict[str, dict[str, Path]] = {}

    if not runs_dir.exists():
        return runs

    for algorithm_dir in sorted(runs_dir.iterdir()):
        if not algorithm_dir.is_dir():
            continue

        algorithm = algorithm_label(algorithm_dir.name)
        runs.setdefault(algorithm, {})

        for csv_path in sorted(algorithm_dir.glob("*.csv")):
            run_key = run_key_from_filename(csv_path.name)
            if run_key is None:
                continue
            runs[algorithm][run_key] = csv_path

    return runs


def generate_hypervolume_csv(
    runs_dir: Path | None = None,
    output_path: Path | None = None,
) -> Path:
    runs_dir = runs_dir or script_runs_dir()
    output_path = output_path or (runs_dir / "hypervolume.csv")

    collected_runs = collect_run_files(runs_dir)
    run_keys = sorted(
        {run_key for algorithm_runs in collected_runs.values() for run_key in algorithm_runs}
    )
    algorithms = sorted(collected_runs)

    rows: list[dict[str, str | float]] = []
    for algorithm in algorithms:
        row: dict[str, str | float] = {"algorithm": algorithm}
        for run_key in run_keys:
            run_path = collected_runs.get(algorithm, {}).get(run_key)
            row[run_key] = (
                compute_hypervolume_for_run(run_path)
                if run_path is not None
                else ""
            )
        rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["algorithm", *run_keys])
        writer.writeheader()
        writer.writerows(rows)

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a CSV with hypervolume values for all runs."
    )
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=script_runs_dir(),
        help="Directory containing algorithm run folders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path (default: <runs-dir>/hypervolume.csv).",
    )
    args = parser.parse_args()

    output_path = generate_hypervolume_csv(args.runs_dir, args.output)
    print(f"Wrote hypervolume summary to {output_path}")


if __name__ == "__main__":
    main()
