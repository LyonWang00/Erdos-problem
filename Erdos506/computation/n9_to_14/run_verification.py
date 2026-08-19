#!/usr/bin/env python3
"""Reproduce the final n=9,...,14 verification chain.

All generated certificates are written below ``certificates``.  The runner
accepts a stage only when the producing program terminates successfully; its
final audit additionally checks every terminal status, the unique n=10
residual signature, and SHA-256 digests of all source and certificate files.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
CERTIFICATES = HERE / "certificates"
RANGES = {
    9: (17, 24),
    10: (25, 32),
    11: (33, 40),
    12: (41, 50),
}
EXPECTED_N10 = [10, 0, 0, 0, 10, 20, 2, 0, 15]


def run(script: str, *arguments: object) -> None:
    command = [sys.executable, str(HERE / script)]
    command.extend(map(str, arguments))
    subprocess.run(command, cwd=HERE, check=True)


def read(name: str) -> dict[str, object]:
    return json.loads((CERTIFICATES / name).read_text(encoding="utf-8"))


def retained_top_signatures(document: dict[str, object]) -> list[list[int]]:
    return [
        record["signature"]
        for record in document["levels"][str(document["n"])]
        if not str(record["status"]).startswith("INFEASIBLE")
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--seconds", type=float, default=240.0)
    args = parser.parse_args()
    CERTIFICATES.mkdir(parents=True, exist_ok=True)

    run("verify_largest_block_reduction.py", "--quiet")

    for n, (first, last) in RANGES.items():
        prefix = f"n{n}"
        local = CERTIFICATES / f"{prefix}_01_signature_filter.json"
        summed = CERTIFICATES / f"{prefix}_02_augmented_summed.json"
        rows = CERTIFICATES / f"{prefix}_03_augmented_point_rows.json"
        run(
            "verify_local_type_filter.py",
            "--n", n,
            "--first-circle-count", first,
            "--last-circle-count", last,
            "--output", local,
            "--quiet",
        )
        run(
            "verify_augmented_line_melchior.py",
            "--n", n, "--input-json", local,
            "--summed-only", "--jobs", args.jobs,
            "--seconds", args.seconds, "--output", summed,
        )
        run(
            "verify_augmented_line_melchior.py",
            "--n", n, "--input-json", summed,
            "--jobs", args.jobs, "--seconds", args.seconds,
            "--output", rows,
        )

        if n == 9:
            run(
                "verify_recursive_inversion_rows.py",
                "--n", 9, "--input-json", rows,
                "--jobs", args.jobs, "--seconds", args.seconds,
                "--output", CERTIFICATES / "n9_04_exact_support.json",
            )
            continue

        endpoint = CERTIFICATES / f"{prefix}_04_endpoint_balance.json"
        moments = CERTIFICATES / f"{prefix}_05_pair_moments.json"
        run(
            "filter_by_classified_split_endpoints.py",
            "--n", n, "--input-json", rows,
            "--no-zero-catalogue", "--no-combined-catalogue",
            "--augmented-arrangement-types",
            "--jobs", args.jobs, "--seconds", args.seconds,
            "--output", endpoint,
        )
        run(
            "verify_pair_moment_filter.py",
            "--n", n, "--input", endpoint,
            "--jobs", args.jobs, "--seconds", args.seconds,
            "--output", moments, "--quiet",
        )

        recursive_input = moments
        recursive_number = 6
        if n == 10:
            blocks = CERTIFICATES / "n10_06_block_rows.json"
            run(
                "block_row_existence_filter.py",
                "--n", 10, "--input", moments,
                "--jobs", args.jobs,
                "--point-row-capacities", "--point-row-pair-sums",
                "--deletion-bounds", "--subset-inheritance-bounds",
                "--conditioned-line-bounds", "--output", blocks,
            )
            recursive_input = blocks
            recursive_number = 7
        run(
            "verify_recursive_inversion_rows.py",
            "--n", n, "--input-json", recursive_input,
            "--jobs", args.jobs, "--seconds", args.seconds,
            "--output",
            CERTIFICATES / f"{prefix}_{recursive_number:02d}_recursive.json",
        )

    run("verify_n10_final_radical_axis.py")
    run(
        "verify_n13_n14_global_inequalities.py",
        "--output", CERTIFICATES / "n13_n14_global_inequalities.json",
    )

    largest_block = read("largest_block_reduction.json")
    n9 = read("n9_04_exact_support.json")
    n10 = read("n10_07_recursive.json")
    n10_geometry = read("n10_08_final_geometry.json")
    n11 = read("n11_06_recursive.json")
    n12 = read("n12_06_recursive.json")
    n13_n14 = read("n13_n14_global_inequalities.json")
    checks = {
        "largest_block_reduction": largest_block["status"] == "PASS"
            and [record["n"] for record in largest_block["records"]]
            == list(range(7, 16)),
        "n9_closed": n9["status"] == "PASS"
            and n9["retained_signature_count"] == 0,
        "n10_unique_residual": retained_top_signatures(n10)
            == [EXPECTED_N10],
        "n10_geometry_closed": n10_geometry["status"] == "PASS",
        "n11_closed": n11["status"] == "PASS"
            and n11["retained_signature_count"] == 0,
        "n12_closed": n12["status"] == "PASS"
            and n12["retained_signature_count"] == 0,
        "n13_n14_closed": n13_n14["status"] == "PASS"
            and all(record["retained_signature_count"] == 0
                    for record in n13_n14["records"]),
    }
    files = sorted(
        path for path in HERE.glob("*.py")
    ) + sorted(path for path in CERTIFICATES.glob("*.json")
               if path.name != "manifest.json")
    manifest = {
        "schema_version": 1,
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "files": {
            str(path.relative_to(HERE)).replace("\\", "/"): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path.read_bytes()).hexdigest(),
            }
            for path in files
        },
    }
    (CERTIFICATES / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": manifest["status"],
        "checks": checks,
        "manifest": str(CERTIFICATES / "manifest.json"),
    }, indent=2))
    if manifest["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
