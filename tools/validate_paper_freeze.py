#!/usr/bin/env python3
"""Validate the metadata and integrity of doc/paper-freeze.

This checker never runs MD and never treats pending scientific gates as
completed. It is intentionally small so it can be used before manuscript
writing or after moving the archive to another checkout.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "doc" / "paper-freeze"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_yaml(path: Path):
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - environment diagnostic
        raise RuntimeError("PyYAML is required to validate MANIFEST.yaml") from exc
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def check_csv(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.reader(handle) if row]
    if not rows or not rows[0]:
        raise ValueError(f"empty CSV: {path}")
    width = len(rows[0])
    if any(len(row) != width for row in rows[1:]):
        raise ValueError(f"inconsistent CSV columns: {path}")


def main() -> int:
    errors: list[str] = []
    required = [
        "README.md",
        "MANIFEST.yaml",
        "build-environment.txt",
        "commands.txt",
        "DECISIONS.md",
        "LICENSE_DATA.md",
        "checksums.sha256",
        "verification-report.md",
        "output-index.csv",
        "validation/README.md",
        "validation/event-labels.template.csv",
        "validation/synthetic-expected.json",
        "benchmarks/README.md",
        "benchmarks/benchmark-matrix.csv",
        "figures/README.md",
        "figures/source-map.csv",
        "literature/references.bib",
    ]
    for relative in required:
        if not (FREEZE / relative).is_file():
            errors.append(f"missing freeze artifact: {relative}")

    manifest = {}
    if not errors:
        try:
            manifest = load_yaml(FREEZE / "MANIFEST.yaml")
        except Exception as exc:  # pragma: no cover - diagnostic branch
            errors.append(f"invalid MANIFEST.yaml: {exc}")

    if manifest:
        if manifest.get("manuscript_written") is not False:
            errors.append("MANIFEST.yaml manuscript_written must remain false")
        if manifest.get("nhn_conversion_equation_required") is not False:
            errors.append("NHN conversion equation must not be a freeze gate")
        artifacts = manifest.get("artifacts", {})
        flattened: list[str] = []
        for key, value in artifacts.items():
            if isinstance(value, str):
                flattened.append(value)
            elif key == "scripts" and isinstance(value, list):
                flattened.extend(value)
        for relative in flattened:
            if not (ROOT / relative).is_file() and not (FREEZE / relative).is_file():
                errors.append(f"manifest artifact missing: {relative}")

    # Verify the static checksum records. External paths are accepted when the
    # checkout is moved away from the author's local data volumes.
    checksum_file = FREEZE / "checksums.sha256"
    checked = 0
    external_missing = 0
    if checksum_file.is_file():
        for line_number, line in enumerate(checksum_file.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            match = re.match(r"^([0-9a-f]{64})  (.+)$", line)
            if not match:
                errors.append(f"malformed checksum line {line_number}")
                continue
            expected, raw_path = match.groups()
            path = Path(raw_path)
            if not path.is_absolute():
                path = ROOT / path
            if not path.is_file():
                if raw_path.startswith(("/mnt/f/", "/mnt/d/CurrentProjects/")):
                    external_missing += 1
                    continue
                errors.append(f"checksum target missing: {raw_path}")
                continue
            checked += 1
            actual = sha256(path)
            if actual != expected:
                errors.append(f"checksum mismatch: {raw_path}")

    for relative in [
        "validation/synthetic-expected.json",
    ]:
        try:
            json.loads((FREEZE / relative).read_text(encoding="utf-8"))
        except Exception as exc:
            errors.append(f"invalid JSON {relative}: {exc}")

    for relative in [
        "output-index.csv",
        "validation/event-labels.template.csv",
        "benchmarks/benchmark-matrix.csv",
        "figures/source-map.csv",
    ]:
        try:
            check_csv(FREEZE / relative)
        except Exception as exc:
            errors.append(str(exc))

    if errors:
        print("paper-freeze validation FAILED")
        for error in errors:
            print(f"- {error}")
        return 1

    pending = sum(bool(value) for value in manifest.get("pending_gates", {}).values())
    print("paper-freeze validation PASSED")
    print(f"  checksum records checked: {checked}")
    print(f"  external records unavailable: {external_missing}")
    print(f"  pending scientific/publication gates: {pending}")
    print("  manuscript_written: false")
    print("  nhn_conversion_equation_required: false")
    return 0


if __name__ == "__main__":
    sys.exit(main())
