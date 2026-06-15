#!/usr/bin/env python3
"""Run portfolio/K experiment sweeps defined in experiments.json."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent
RUNS_DIR = SRC_DIR / "runs"
LOGS_DIR = RUNS_DIR / "logs"

ALGORITHM_DIRS = {
    "nsgaii": "nsgaii",
    "branch_and_bound": "bb",
}

sys.path.insert(0, str(SCRIPT_DIR))
from indicators import generate_hypervolume_csv  # noqa: E402


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def portfolio_path(config: dict, portfolio_name: str) -> Path:
    inputs_dir = Path(config["inputs_dir"])
    if not inputs_dir.is_absolute():
        inputs_dir = (SRC_DIR / inputs_dir).resolve()
    return inputs_dir / f"{portfolio_name}.txt"


def expected_output_path(algorithm: str, portfolio_name: str, k: int) -> Path:
    algorithm_dir = ALGORITHM_DIRS[algorithm]
    return RUNS_DIR / algorithm_dir / f"{portfolio_name}_k{k}.csv"


def build_command(config: dict, algorithm: str) -> list[str]:
    build = config["build"]
    source_dir = Path(build.get("source_dir", "."))
    if not source_dir.is_absolute():
        source_dir = (SRC_DIR / source_dir).resolve()

    if algorithm == "nsgaii":
        sources = build["common_sources"] + build["nsgaii_sources"]
        binary = build["nsgaii_binary"]
        flags = build["nsgaii_flags"]
    elif algorithm == "branch_and_bound":
        sources = build["common_sources"] + build["branch_and_bound_sources"]
        binary = build["branch_and_bound_binary"]
        flags = build["branch_and_bound_flags"]
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")

    output_path = source_dir / binary
    command = ["g++", *flags, *[str(source_dir / source) for source in sources], "-o", str(output_path)]
    return command


def run_command(command: list[str], log_path: Path, cwd: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"command: {' '.join(command)}\n")
        log_handle.write(f"cwd: {cwd}\n")
        log_handle.write(f"started: {datetime.now().isoformat()}\n\n")
        log_handle.flush()

        process = subprocess.run(
            command,
            cwd=cwd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )

        log_handle.write(f"\nfinished: {datetime.now().isoformat()}\n")
        log_handle.write(f"exit_code: {process.returncode}\n")

    return process.returncode


def algorithm_command(config: dict, algorithm: str, portfolio_file: Path, k: int) -> list[str]:
    build = config["build"]
    source_dir = Path(build.get("source_dir", "."))
    if not source_dir.is_absolute():
        source_dir = (SRC_DIR / source_dir).resolve()

    portfolio_arg = str(portfolio_file)

    if algorithm == "nsgaii":
        binary = source_dir / build["nsgaii_binary"]
        return [str(binary), "--portfolio", portfolio_arg, "--k", str(k)]

    if algorithm == "branch_and_bound":
        binary = source_dir / build["branch_and_bound_binary"]
        return [str(binary), "--portfolio", portfolio_arg, "--k", str(k), "-", "-"]

    raise ValueError(f"Unknown algorithm: {algorithm}")


def build_binaries(config: dict, algorithms: list[str]) -> None:
    for algorithm in sorted(set(algorithms)):
        print(f"Building {algorithm}...")
        command = build_command(config, algorithm)
        log_path = LOGS_DIR / f"build_{algorithm}.log"
        exit_code = run_command(command, log_path, SRC_DIR)
        if exit_code != 0:
            raise RuntimeError(
                f"Build failed for {algorithm} with exit code {exit_code}. See {log_path}"
            )


def run_experiments(
    config_path: Path | None = None,
    force: bool = False,
    skip_build: bool = False,
    algorithms: list[str] | None = None,
) -> None:
    config_path = config_path or (SRC_DIR / "experiments.json")
    config = load_config(config_path)

    selected_algorithms = algorithms or config["algorithms"]
    portfolios = config["portfolios"]
    k_values = config["k_values"]

    if not skip_build:
        build_binaries(config, selected_algorithms)

    failures: list[str] = []

    for algorithm in selected_algorithms:
        for portfolio_name in portfolios:
            for k in k_values:
                output_path = expected_output_path(algorithm, portfolio_name, k)
                run_label = f"{algorithm}/{portfolio_name}_k{k}"

                if output_path.exists() and not force:
                    print(f"Skipping existing run: {run_label}")
                    continue

                portfolio_file = portfolio_path(config, portfolio_name)
                if not portfolio_file.exists():
                    failures.append(f"{run_label}: missing portfolio file {portfolio_file}")
                    continue

                command = algorithm_command(config, algorithm, portfolio_file, k)
                log_path = LOGS_DIR / f"{algorithm}_{portfolio_name}_k{k}.log"

                print(f"Running {run_label}...")
                exit_code = run_command(command, log_path, SRC_DIR)
                if exit_code != 0:
                    failures.append(f"{run_label}: exit code {exit_code} (see {log_path})")

    hypervolume_path = generate_hypervolume_csv(RUNS_DIR)
    print(f"Wrote hypervolume summary to {hypervolume_path}")

    if failures:
        print("\nFailures:")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run portfolio/K experiment sweeps.")
    parser.add_argument(
        "--config",
        type=Path,
        default=SRC_DIR / "experiments.json",
        help="Path to experiments.json",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run experiments even if output CSV already exists",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip compiling binaries before running",
    )
    parser.add_argument(
        "--algorithms",
        nargs="+",
        choices=["nsgaii", "branch_and_bound"],
        help="Run only the selected algorithms",
    )
    args = parser.parse_args()

    run_experiments(
        config_path=args.config,
        force=args.force,
        skip_build=args.skip_build,
        algorithms=args.algorithms,
    )


if __name__ == "__main__":
    main()
