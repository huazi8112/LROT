#!/usr/bin/env python3
"""Verify release-file hashes and displayed sizes against the reviewed allowlist.

``--write`` rebuilds inventories from ``RELEASE_FILES.txt`` without requiring Git.
``--strict`` also rejects files that are not in the release list (except manifests).
Use strict mode on an extracted release package; a development checkout may keep
additional historical files outside the release list.
"""
from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
SHA_FILE = ROOT / "MANIFEST.sha256"
SIZE_FILE = ROOT / "MANIFEST_size.txt"
LIST_FILE = ROOT / "RELEASE_FILES.txt"
MANIFESTS = {"MANIFEST.sha256", "MANIFEST_size.txt"}
SIZE_PATTERN = re.compile(r"^\s*(\d+(?:\.\d+)?) KB\s{2}(.+)$")


def release_names() -> list[str]:
    names = [line.strip() for line in LIST_FILE.read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.startswith("#")]
    if names != sorted(set(names)) or "RELEASE_FILES.txt" not in names:
        raise ValueError("RELEASE_FILES.txt must contain sorted unique paths including itself")
    for name in names:
        pure = PurePosixPath(name)
        if pure.is_absolute() or ".." in pure.parts or "\\" in name or name in MANIFESTS:
            raise ValueError(f"invalid release path: {name}")
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError(f"release path escapes package root: {name}")
    return names


def read_manifests() -> tuple[dict[str, str], dict[str, float], int | None]:
    hashes: dict[str, str] = {}
    for line in SHA_FILE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, rel = line.split("  ", 1)
        if rel in hashes:
            raise ValueError(f"duplicate SHA entry: {rel}")
        hashes[rel] = digest
    size_lines = SIZE_FILE.read_text(encoding="utf-8").splitlines()
    header = re.search(r"(\d+) 个文件", size_lines[0]) if size_lines else None
    header_count = int(header.group(1)) if header else None
    sizes: dict[str, float] = {}
    for line in size_lines[2:]:
        if not line.strip():
            continue
        match = SIZE_PATTERN.match(line)
        if not match:
            raise ValueError(f"invalid size entry: {line}")
        size_kb, rel = match.groups()
        if rel in sizes:
            raise ValueError(f"duplicate size entry: {rel}")
        sizes[rel] = float(size_kb)
    return hashes, sizes, header_count


def write_manifests(names: list[str]) -> None:
    hashes, sizes, total = [], [], 0
    for name in names:
        data = (ROOT / name).read_bytes()
        total += len(data)
        hashes.append(f"{hashlib.sha256(data).hexdigest()}  {name}")
        sizes.append(f"{len(data) / 1024:10.1f} KB  {name}")
    SHA_FILE.write_bytes(("\n".join(hashes) + "\n").encode("utf-8"))
    heading = f"LROT 代码包文件清单（{len(names)} 个文件，{total / 1024**2:.1f} MB）"
    SIZE_FILE.write_bytes(("\n".join([heading, "", *sizes]) + "\n").encode("utf-8"))
    print(f"Rebuilt SHA-256 and size inventories: {len(names)} files, {total / 1024**2:.1f} MB")


def verify(strict: bool) -> int:
    try:
        names = release_names()
        hashes, sizes, header_count = read_manifests()
    except (OSError, ValueError) as error:
        print(f"Manifest verification FAILED: {error}")
        return 1
    errors = []
    if set(names) != set(hashes) or set(names) != set(sizes):
        errors.append("release list, SHA and size entries differ")
    if header_count is not None and header_count != len(names):
        errors.append(f"size header lists {header_count} files, expected {len(names)}")
    for name in names:
        path = ROOT / name
        if not path.is_file():
            errors.append(f"missing file: {name}")
            continue
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != hashes.get(name):
            errors.append(f"SHA-256 mismatch: {name}")
        if round(len(data) / 1024, 1) != sizes.get(name):
            errors.append(f"size mismatch: {name}")
    if strict:
        actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
                  if p.is_file() and ".git" not in p.relative_to(ROOT).parts}
        extras = sorted(actual - set(names) - MANIFESTS)
        if extras:
            errors.append(f"{len(extras)} unlisted files: {', '.join(extras[:10])}")
    if errors:
        print("Manifest verification FAILED:\n" + "\n".join(f"- {e}" for e in errors))
        return 1
    print(f"Manifest verification passed: {len(names)} files; SHA-256 and sizes match.")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="rebuild inventories from RELEASE_FILES.txt")
    parser.add_argument("--strict", action="store_true", help="reject unlisted files")
    args = parser.parse_args()
    if args.write:
        write_manifests(release_names())
    raise SystemExit(verify(args.strict))
