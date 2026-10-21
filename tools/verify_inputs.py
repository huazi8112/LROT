"""Check separately supplied input files without downloading or modifying data.

Checksums identify bytes only. See data/input_files.json for the distinction
between previously bound inputs and the current packaging inventory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ("breast", "mouse", "dlpfc", "go", "cellchat")


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(rows: list[dict], datasets: list[str], data_root: Path | None,
           resource_dir: Path) -> list[dict]:
    results = []
    for row in rows:
        if row["dataset"] not in datasets:
            continue
        relative = PurePosixPath(row["path"])
        if relative.is_absolute() or ".." in relative.parts or "\\" in row["path"]:
            raise ValueError("Unsafe input-manifest path: " + row["path"])
        if row["dataset"] == "cellchat":
            if len(relative.parts) != 2 or relative.parts[0] != "data":
                raise ValueError("Invalid resource location: " + row["path"])
            path = resource_dir / relative.name
        else:
            if data_root is None:
                raise ValueError("--data-root is required for raw datasets")
            path = data_root.joinpath(*relative.parts)
        result = {"dataset": row["dataset"], "path": row["path"]}
        try:
            if not path.is_file():
                result["status"] = "missing"
            elif path.stat().st_size != row["size_bytes"]:
                result["status"] = "size_mismatch"
            elif checksum(path) != row["sha256"]:
                result["status"] = "sha256_mismatch"
            else:
                result["status"] = "matched"
        except OSError as error:
            result.update(status="unreadable", reason=type(error).__name__)
        results.append(result)
    for dataset in datasets:
        if not any(row["dataset"] == dataset for row in results):
            raise ValueError("No manifest entries for " + dataset)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list required files; no input reads")
    parser.add_argument("--dataset", action="append", choices=DATASETS)
    parser.add_argument("--data-root", type=Path,
                        help="external parent of cancer_data/, dlpfc_data/ and data/go/")
    parser.add_argument("--resource-dir", type=Path, default=ROOT / "data",
                        help="directory containing separately supplied CellChatDB files")
    parser.add_argument("--json", action="store_true", help="print machine-readable results")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "data/input_files.json").read_text(encoding="utf-8"))
    datasets = list(dict.fromkeys(args.dataset or DATASETS))
    rows = manifest["files"]
    if args.list:
        for row in rows:
            if row["dataset"] in datasets:
                print(row["dataset"], row["path"], row["sha256"], row["identity_basis"], sep="\t")
        return 0
    if not args.dataset:
        parser.error("Select one or more --dataset values, or use --list.")
    if any(dataset != "cellchat" for dataset in datasets) and args.data_root is None:
        parser.error("--data-root is required for raw datasets")
    try:
        results = verify(rows, datasets, args.data_root, args.resource_dir)
    except (ValueError, KeyError) as error:
        parser.error(str(error))
    passed = all(row["status"] == "matched" for row in results)
    if args.json:
        print(json.dumps({"passed": passed, "files": results,
                          "scope": "Byte identity only; no download, input edits or experiment execution."}, indent=2))
    else:
        for row in results:
            print(row["status"], row["path"], sep="\t")
        print("Input check passed." if passed else "Input check failed; do not substitute different inputs.")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
