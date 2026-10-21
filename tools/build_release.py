#!/usr/bin/env python3
"""Build an LROT release directory from the reviewed RELEASE_FILES.txt allowlist.

The destination must not exist. This tool never deletes or overwrites a release.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path, PurePosixPath
import shutil

EXCLUDED = {"MANIFEST.sha256", "MANIFEST_size.txt"}


def release_paths(root: Path) -> list[str]:
    path = root / "RELEASE_FILES.txt"
    entries = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
               if line.strip() and not line.startswith("#")]
    if entries != sorted(set(entries)) or "RELEASE_FILES.txt" not in entries:
        raise ValueError("RELEASE_FILES.txt must contain sorted, unique paths including itself")
    for name in entries:
        pure = PurePosixPath(name)
        if pure.is_absolute() or ".." in pure.parts or "\\" in name or name in EXCLUDED:
            raise ValueError(f"invalid release path: {name}")
        if not (root / name).is_file():
            raise FileNotFoundError(root / name)
    return entries


def write_manifests(root: Path, entries: list[str]) -> None:
    hashes, sizes, total = [], [], 0
    for name in entries:
        data = (root / name).read_bytes()
        total += len(data)
        hashes.append(f"{hashlib.sha256(data).hexdigest()}  {name}")
        sizes.append(f"{len(data) / 1024:10.1f} KB  {name}")
    (root / "MANIFEST.sha256").write_bytes(("\n".join(hashes) + "\n").encode("utf-8"))
    heading = f"LROT 代码包文件清单（{len(entries)} 个文件，{total / 1024**2:.1f} MB）"
    (root / "MANIFEST_size.txt").write_bytes(
        ("\n".join([heading, "", *sizes]) + "\n").encode("utf-8"))
    print(f"Built {len(entries)} listed files, {total / 1024**2:.1f} MiB")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path,
                        help="new destination directory; it must not exist")
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]
    destination = args.output.resolve()
    if destination == source or source in destination.parents:
        raise ValueError("destination must be outside the source package")
    entries = release_paths(source)
    destination.mkdir(parents=True, exist_ok=False)
    for name in entries:
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / name, target)
    write_manifests(destination, entries)
    print(f"Release candidate: {destination}")


if __name__ == "__main__":
    main()
