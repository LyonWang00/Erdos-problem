#!/usr/bin/env python3
"""Fast logical audit of the current certificate map and final certificates."""

from __future__ import annotations

from math import comb
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(relative: str) -> dict:
    return json.loads((HERE / relative).read_text(encoding="utf-8"))


def f(n: int) -> int:
    return 1 + comb(n - 1, 2) - (n - 1) // 2


def main() -> None:
    index = read("certificate_map.json")
    assert index["schema_version"] == 6
    assert index["certified_statement"]["large_range"] == (
        "c(n)=F(n) for n>=9"
    )
    assert index["certified_statement"]["small_values"] == {
        "4": 3, "5": 5, "6": 8, "7": 11, "8": 17,
    }
    assert index["certified_statement"]["no_three_collinear"] == (
        "c_nc(n)=1+C(n-1,2) for n>=4, except c_nc(8)=20"
    )

    unified = read("n9_to_14/certificates/manifest.json")
    assert unified["status"] == "PASS"
    assert all(unified["checks"].values())

    largest = read(
        "n9_to_14/certificates/largest_block_reduction.json"
    )
    assert largest["status"] == "PASS"
    assert [row["n"] for row in largest["records"]] == list(range(7, 16))

    n8 = read("n4_to_8/n8/boundary_layer_classification.json")
    assert n8["status"] == "PASS"
    assert [row["circle_count"] for row in n8["layers"]] == list(range(12, 17))
    assert [row["class_count"] for row in n8["layers"]] == [1, 1, 1, 3, 3]

    assert [f(n) for n in range(4, 9)] == [3, 5, 9, 13, 19]
    assert index["certified_statement"]["small_values"]["6"] == f(6) - 1
    assert index["certified_statement"]["small_values"]["7"] == f(7) - 2
    assert index["certified_statement"]["small_values"]["8"] == f(8) - 2

    print(json.dumps({
        "status": "PASS",
        "certified_range": "all n>=4",
        "formula_exceptions": {"6": 8, "7": 11, "8": 17},
        "no_three_collinear_exception": {"8": 20},
        "n9_n14_unknown": 0,
        "n8_boundary_class_counts": [1, 1, 1, 3, 3],
    }, indent=2))


if __name__ == "__main__":
    main()
