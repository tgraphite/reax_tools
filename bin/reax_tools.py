#!/usr/bin/env python3
"""Unified ReaxTools command router."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = ROOT / "src" / "python"
CORE = ROOT / "bin" / "reax_tools_core"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/reax_tools_matplotlib")

# If someone runs reax_tools.py directly with a newer system Python and the
# installer created a project-local venv, switch to the venv automatically.
VENV_PYTHON = ROOT / "venv" / "bin" / "python"
if VENV_PYTHON.exists() and Path(sys.executable).resolve() != VENV_PYTHON.resolve():
    os.execv(str(VENV_PYTHON), [str(VENV_PYTHON), str(Path(__file__).resolve()), *sys.argv[1:]])

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))


def run_core(args: list[str]) -> int:
    if not CORE.exists():
        print(f"error: C++ core binary not found: {CORE}", file=sys.stderr)
        print("hint: run `bash install_reax_tools.sh` first", file=sys.stderr)
        return 127
    return subprocess.call([str(CORE), *args])


def run_python(args: list[str]) -> int:
    from reax_tools_viz.cli import main

    old_argv = sys.argv
    try:
        sys.argv = ["reax_tools", *args]
        return main()
    finally:
        sys.argv = old_argv


def print_help() -> None:
    print(
        "ReaxTools\n\n"
        "Version: 2.3\n\n"
        "Analyze mode (C++ core):\n"
        "  reax_tools analyze -f <trajectory> [analysis options]\n"
        "  reax_tools -f <trajectory> [analysis options]\n"
        "\n"
        "If the first argument is an option such as -f, ReaxTools routes to analyze mode.\n"
        "\n"
        "Common analyze options:\n"
        "  -f, --traj <file>          input trajectory (.xyz/.lammpstrj)\n"
        "  -o, --output <dir>         output directory\n"
        "  -t, --types C,H,O,N        element types for LAMMPS dump files\n"
        "  -r, --rescale-vdw <value>  vdw radius scale factor\n"
        "  -tr, --type-radius N:1.5   override element radius\n"
        "  -tv, --type-valence N:4    override element valence\n"
        "  --no-rings                 disable ring detection\n"
        "  --no-reactions             disable transfer-flow analysis\n"
        "  --no-track-reactions       disable reaction-event tracking\n"
        "  --stable-time <frames>     molecule stability threshold\n"
        "\n"
        "Python post-processing commands:\n"
        "  reax_tools plot -f <output_dir>                 default plots + key molecules\n"
        "  reax_tools counts -f <count.csv|output_dir>     count line plots\n"
        "  reax_tools network -f <flow.csv|output_dir>     transfer network graph\n"
        "  reax_tools flow -f <flow.csv|output_dir>        experimental Sankey flow\n"
        "  reax_tools focus -f <flow.csv|output_dir>       centered Sankey views\n"
        "  reax_tools events -f <events.csv|output_dir>    reaction-event plot\n"
        "  reax_tools molecules -f <csv|output_dir>        SMILES molecule drawings\n"
        "  reax_tools summary -f <output_dir>              output validation summary\n"
        "\n"
        "Use `reax_tools analyze --help` for the full C++ option table.\n"
        "Use `reax_tools <command> --help` for a Python command.\n",
        file=sys.stdout,
    )


def main() -> int:
    args = sys.argv[1:]
    if not args:
        print_help()
        return 0

    command = args[0]
    if command == "analyze":
        return run_core(args[1:])

    if command in {"plot", "network", "flow", "focus", "events", "counts", "summary", "molecules"}:
        return run_python(args)

    if command in {"-h", "--help", "help"}:
        print_help()
        return 0

    # Backward-compatible shortcut: `reax_tools -f traj.xyz ...`.
    if command.startswith("-"):
        return run_core(args)

    print(f"error: unknown command: {command}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
