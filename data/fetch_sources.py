#!/usr/bin/env python3
"""Verify bundled public GEO sources; restore missing files by their recorded URL."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


def verify_sources(data_dir: Path, offline: bool = True) -> list[dict]:
    """Return provenance after checking every source byte against sources.json."""
    data_dir = Path(data_dir).resolve()
    manifest = json.loads((data_dir / "sources.json").read_text())
    verified = []
    for record in manifest["files"]:
        relative = Path(record["relative_path"])
        path = (data_dir / relative).resolve()
        if not path.is_relative_to(data_dir):
            raise ValueError(f"Source path escapes data directory: {relative}")
        if not path.exists():
            if offline:
                raise FileNotFoundError(f"Missing {path}; rerun fetch_sources.py without --offline")
            request = urllib.request.Request(
                record["requested_url"], headers={"User-Agent": "ExtremaRank reproducible source fetch"}
            )
            with urllib.request.urlopen(request, timeout=90) as response:
                payload = response.read()
            if len(payload) != record["bytes"] or hashlib.sha256(payload).hexdigest() != record["sha256"]:
                raise ValueError(f"Remote source differs from pinned snapshot: {relative}")
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix + ".partial")
            temporary.write_bytes(payload)
            temporary.replace(path)
        payload = path.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if len(payload) != record["bytes"] or digest != record["sha256"]:
            raise ValueError(f"Source checksum mismatch: {relative}")
        verified.append({"relative_path": str(relative), "bytes": len(payload), "sha256": digest})
    return verified


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--offline", action="store_true", help="Verify only; make no network requests")
    args = parser.parse_args()
    records = verify_sources(args.data_dir, offline=args.offline)
    print(json.dumps({"verified_source_files": len(records), "files": records}, indent=2))


if __name__ == "__main__":
    main()
