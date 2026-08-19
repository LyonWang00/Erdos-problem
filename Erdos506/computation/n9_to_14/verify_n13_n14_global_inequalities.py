#!/usr/bin/env python3
"""Exact global-signature verification for n=13 and n=14.

No point-support or orbit search is used.  Candidate signatures are generated
from the pair/triple identities and the documented classical line-arrangement,
deletion, and subset inequalities in ``SignatureFilter``.  They are then
excluded by the augmented-line Melchior inequality, the elementary union
capacity of selected rich blocks, and (only for n=13) two local inequalities
for arrangements of twelve projective lines.

For a twelve-line local arrangement write delta for its Melchior defect and
t4,t5 for the numbers of vertices of multiplicity four and five.  The two
pointwise inequalities are

    -2*delta + t4 + 3*t5 <= 6,
    -delta + t5 <= 1.

The first follows from Melchior and the Bojanowski inequality, integrality,
and the classified zero-defect equality case.  For the second, four
five-fold vertices violate the union bound; the cases delta=0 and
(delta,t5)=(1,3) are removed respectively by the zero-defect classification
and by the pair identity plus the same union bound.  The paper gives the
full derivation.  This program checks their summed applications exactly.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
COMPUTATIONS = HERE.parent
AUTHORITATIVE = HERE
sys.path.insert(0, str(AUTHORITATIVE))

from signature_filters import SignatureFilter, block_family_convexity_witness
from projected_arrangement_model import melchior_defect, projected_local_types
from verify_augmented_line_melchior import summed_gap, total_defect


CIRCLE_RANGES = {13: range(51, 61), 14: range(61, 73)}
EXPECTED = {
    13: {
        "raw": 3962,
        "augmented": 137,
        "block_capacity": 73,
        "local_first": 71,
        "local_second": 2,
    },
    14: {
        "raw": 8981,
        "augmented": 1,
        "block_capacity": 0,
        "local_first": 0,
        "local_second": 0,
    },
}


def all_signatures(n: int) -> tuple[tuple[int, ...], ...]:
    generator = SignatureFilter(n)
    return tuple(
        signature
        for circle_count in CIRCLE_RANGES[n]
        for signature in generator.signatures(circle_count)
    )


def local_inequality_values(
    n: int, signature: tuple[int, ...],
) -> tuple[int, int]:
    dimension = len(signature) // 2
    blocks = tuple(
        signature[index] + signature[dimension + index]
        for index in range(dimension)
    )
    defect = total_defect(n, signature)
    # Sum_x t4(x)=5 B5 and sum_x t5(x)=6 B6.
    first = -2 * defect + 5 * blocks[2] + 18 * blocks[3]
    second = -defect + 6 * blocks[3]
    return first, second


def verify(n: int) -> dict[str, object]:
    generator = SignatureFilter(n)
    raw = all_signatures(n)
    after_augmented = tuple(
        signature for signature in raw if summed_gap(n, signature) <= 0
    )
    block_records = []
    after_block_capacity = []
    for signature in after_augmented:
        witness = block_family_convexity_witness(
            n, generator.sizes, signature
        )
        block_records.append({
            "signature": list(signature),
            "augmented_summed_gap": summed_gap(n, signature),
            "block_capacity_witness": witness,
        })
        if witness is None:
            after_block_capacity.append(signature)

    first_count = second_count = 0
    local_records = []
    for signature in after_block_capacity:
        first, second = local_inequality_values(n, signature)
        if n == 13 and first > 6 * n:
            reason = "minus_2_delta_plus_t4_plus_3_t5"
            first_count += 1
        elif n == 13 and second > n:
            reason = "minus_delta_plus_t5"
            second_count += 1
        else:
            reason = "NOT_EXCLUDED"
        local_records.append({
            "signature": list(signature),
            "summed_first_left": first,
            "summed_first_right": 6 * n,
            "summed_second_left": second,
            "summed_second_right": n,
            "reason": reason,
        })

    counts = {
        "raw": len(raw),
        "augmented": len(after_augmented),
        "block_capacity": len(after_block_capacity),
        "local_first": first_count,
        "local_second": second_count,
    }
    retained = sum(
        record["reason"] == "NOT_EXCLUDED" for record in local_records
    )
    passed = counts == EXPECTED[n] and retained == 0
    return {
        "n": n,
        "status": "PASS" if passed else "FAIL",
        "claim": (
            "Complete exact arithmetic reduction at the global-signature "
            "level; no labelled supports or geometric orbits."
        ),
        "counts": counts,
        "retained_signature_count": retained,
        "augmented_survivors": block_records,
        "post_capacity_records": local_records,
    }


def audit_twelve_line_local_inequalities() -> dict[str, object]:
    local_types = projected_local_types(13, 6)
    values = []
    for local_type in local_types:
        defect = melchior_defect(local_type)
        t4, t5 = local_type[2], local_type[3]
        values.append((-2 * defect + t4 + 3 * t5, -defect + t5))
    first_maximum = max(first for first, _second in values)
    second_maximum = max(second for _first, second in values)
    return {
        "number_of_retained_numerical_types": len(local_types),
        "first_left_maximum": first_maximum,
        "first_right": 6,
        "second_left_maximum": second_maximum,
        "second_right": 1,
        "status": (
            "PASS" if first_maximum <= 6 and second_maximum <= 1
            else "FAIL"
        ),
    }
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = [verify(13), verify(14)]
    local_audit = audit_twelve_line_local_inequalities()
    document = {
        "schema_version": 1,
        "status": (
            "PASS" if (
                all(record["status"] == "PASS" for record in records)
                and local_audit["status"] == "PASS"
            )
            else "FAIL"
        ),
        "claim": (
            "The n=13 and n=14 counterexample ranges are excluded by exact "
            "global-signature inequalities."
        ),
        "twelve_line_local_inequality_audit": local_audit,
        "records": records,
    }
    canonical = json.dumps(document, sort_keys=True).encode()
    document["payload_sha256"] = sha256(canonical).hexdigest()
    args.output.write_text(json.dumps(document, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "status": document["status"],
        "counts": {str(record["n"]): record["counts"] for record in records},
        "output": str(args.output),
    }, indent=2))
    if document["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
