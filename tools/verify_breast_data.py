#!/usr/bin/env python3
"""Verify the SHA-256 hashes of the raw breast-cancer inputs used by the paper."""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "cancer_data" / "BREAST_SHA256SUMS.txt"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def read_manifest() -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    for line_number, line in enumerate(MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            raise ValueError(f"invalid manifest line {line_number}")
        expected, relative = parts
        path = PurePosixPath(relative)
        if not SHA256_RE.fullmatch(expected):
            raise ValueError(f"invalid SHA-256 on line {line_number}")
        if path.is_absolute() or ".." in path.parts or "\\" in relative:
            raise ValueError(f"unsafe relative path on line {line_number}")
        entries.append((expected, relative))
    if not entries:
        raise ValueError("manifest contains no input files")
    return entries


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default=os.environ.get("LROT_DATA_DIR"),
        help="folder containing the two H5 files and their section subfolders "
             "(defaults to LROT_DATA_DIR, then repository/cancer_data)",
    )
    args = parser.parse_args()
    data_dir = Path(args.data_dir).expanduser() if args.data_dir else ROOT / "cancer_data"

    try:
        entries = read_manifest()
    except (OSError, ValueError) as error:
        print(f"Breast input check FAILED: {error}", file=sys.stderr)
        return 1

    errors = []
    print(f"Checking breast raw inputs under: {data_dir}")
    for expected, relative in entries:
        path = data_dir.joinpath(*PurePosixPath(relative).parts)
        if not path.is_file():
            errors.append(f"missing: {relative}")
            print(f"MISSING  {relative}")
            continue
        actual = sha256_file(path)
        if actual != expected:
            errors.append(f"SHA-256 mismatch: {relative}")
            print(f"MISMATCH {relative}\n  expected {expected}\n  actual   {actual}")
        else:
            print(f"OK       {relative}")

    if errors:
        print("Breast input check FAILED:\n- " + "\n- ".join(errors), file=sys.stderr)
        return 1
    print(f"Breast input check passed: all {len(entries)} files match their recorded SHA-256.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
