#!/usr/bin/env python3
"""Rebuild every active certificate and then audit the certified range."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent


def main() -> None:
    subprocess.run(
        [sys.executable, str(HERE / "verify_certified_results.py")],
        cwd=HERE.parent,
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(HERE / "n9_to_14" / "run_verification.py"),
            "--jobs",
            "4",
            "--seconds",
            "240",
        ],
        cwd=HERE / "n9_to_14",
        check=True,
    )
    subprocess.run(
        [sys.executable, str(HERE / "verify_coverage.py")],
        cwd=HERE.parent,
        check=True,
    )


if __name__ == "__main__":
    main()
