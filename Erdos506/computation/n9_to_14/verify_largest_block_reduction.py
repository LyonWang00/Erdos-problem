#!/usr/bin/env python3
"""Verify the closed-form largest-block reduction used for 7 <= n <= 15."""

from __future__ import annotations

import argparse
from math import comb
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "certificates" / "largest_block_reduction.json"

# If the largest maximal line or circle has at least the indicated size, the
# two elementary counting bounds below already reach F(n).  The remaining
# sizes are handled by the arguments specific to the corresponding order.
FIRST_DIRECT_SIZE = {
    7: 5,
    8: 5,
    9: 6,
    10: 7,
    11: 7,
    12: 7,
    13: 7,
    14: 8,
    15: 8,
}


def extremal_value(n: int) -> int:
    return 1 + comb(n - 1, 2) - (n - 1) // 2


def line_bound(n: int, block_size: int) -> int:
    outside = n - block_size
    return outside * comb(block_size, 2) - comb(outside, 2) * (block_size // 2)


def circle_bound(n: int, block_size: int) -> int:
    outside = n - block_size
    return (
        1
        + outside * (comb(block_size, 2) - block_size // 2)
        - comb(outside, 2) * (block_size // 2)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, object]] = []
    checked_branches = 0
    for n, first_size in FIRST_DIRECT_SIZE.items():
        target = extremal_value(n)
        values = []
        for block_size in range(first_size, n):
            line = line_bound(n, block_size)
            circle = circle_bound(n, block_size)
            outside = n - block_size
            assert line - circle == outside * (block_size // 2) - 1
            assert circle >= target
            values.append(
                {
                    "block_size": block_size,
                    "line_bound": line,
                    "circle_bound": circle,
                    "margin": circle - target,
                }
            )
            checked_branches += 1
        records.append(
            {
                "n": n,
                "target": target,
                "all_block_sizes_at_least": first_size,
                "remaining_block_sizes": [3, first_size - 1],
                "bounds": values,
            }
        )

    payload = {
        "schema_version": 2,
        "status": "PASS",
        "proof_type": "closed-form largest-block counting",
        "coordinate_or_orbit_search": False,
        "checked_branches": checked_branches,
        "records": records,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if not args.quiet:
        print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
