"""Download the pinned source objects used by prepare_hypertension.R.

Run from the main project folder:
    python scripts/python/fetch_hypertension_sources.py
Then run the base-R preparation helper. Ready-to-use TSVs are already supplied.
Existing source files must match the expected SHA-256; they are never overwritten.
"""
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/provenance/hypertension"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    manifest = json.loads((SOURCE / "download-manifest.json").read_text())
    for record in manifest:
        target = SOURCE / record["filename"]
        if target.exists():
            if digest(target) != record["sha256"]:
                raise SystemExit(f"Source checksum mismatch: {target.name}; preserve and investigate it.")
            print(f"Verified existing source: {target.name}")
            continue
        partial = target.with_suffix(target.suffix + ".partial")
        try:
            with urllib.request.urlopen(record["url"], timeout=120) as response, partial.open("wb") as handle:
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    handle.write(chunk)
            if digest(partial) != record["sha256"]:
                raise ValueError(f"Downloaded source checksum mismatch: {target.name}")
            partial.replace(target)
        finally:
            partial.unlink(missing_ok=True)
        print(f"Downloaded and verified: {target.name}")


if __name__ == "__main__":
    main()
