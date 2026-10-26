#!/usr/bin/env python3
"""Check delivered file bytes against manifest.json; this is not a proof check."""
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import sys

ROOT = Path(__file__).resolve().parent


def main():
    try:
        manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
        if manifest.get("schema_version") != 1 or manifest.get("algorithm") != "sha256":
            raise ValueError("Unsupported manifest format")
        files = manifest["files"]
        if not isinstance(files, dict) or not files:
            raise ValueError("The manifest must contain a nonempty file map")
        failures = []
        for name, expected in sorted(files.items()):
            relative = PurePosixPath(name)
            if (relative.is_absolute() or ".." in relative.parts
                    or chr(92) in name or ":" in name or str(relative) != name):
                raise ValueError("Invalid manifest path: " + repr(name))
            path = ROOT.joinpath(*relative.parts)
            if path.is_symlink() or not path.is_file():
                failures.append(name + ": missing, not a regular file, or a symbolic link")
                continue
            if not path.resolve().is_relative_to(ROOT):
                raise ValueError("Manifest path leaves the repository: " + name)
            digest = sha256()
            size = 0
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    size += len(chunk)
                    digest.update(chunk)
            if size != expected["bytes"] or digest.hexdigest() != expected["sha256"]:
                failures.append(name + ": size or SHA-256 mismatch")
        if failures:
            for failure in failures:
                print("FAIL:", failure, file=sys.stderr)
            print("A verification rerun may legitimately replace recorded logs and summaries.",
                  file=sys.stderr)
            return 1
        print(f"PASS: all {len(files)} manifest-listed files match their supplied SHA-256 hashes.")
        print("Integrity only: this does not check the mathematics or rerun interval calculations.")
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        print("FAIL:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
