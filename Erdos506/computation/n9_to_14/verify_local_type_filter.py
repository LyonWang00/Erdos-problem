#!/usr/bin/env python3
"""Uniform aggregate local-arrangement filter for 9 <= n <= 14.

For every global signature, sum the Melchior defect and the incidences of
local vertices of multiplicity at least four over all choices of the marked
point.  A small integer multiplicity model tests each distinct required
total.  No point labels, block supports, or block orbits are used.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

COMPUTATIONS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(COMPUTATIONS))
from runtime_dependencies import activate_ortools  # noqa: E402

activate_ortools()
import ortools  # noqa: E402

from projected_arrangement_model import melchior_defect, projected_local_types
from signature_filters import (
    SignatureFilter,
    block_family_convexity_witness,
    summed_local_line_witness,
)

ORTOOLS_VERSION = ortools.__version__
if ORTOOLS_VERSION != "9.15.6755":
    raise RuntimeError(
        f"expected OR-Tools 9.15.6755, found {ORTOOLS_VERSION}"
    )
from ortools.sat.python import cp_model  # noqa: E402


def target(
    n: int,
    sizes: tuple[int, ...],
    signature: tuple[int, ...],
) -> tuple[int, ...]:
    dimension = len(sizes)
    blocks = tuple(
        signature[index] + signature[dimension + index]
        for index in range(dimension)
    )
    defect = (
        3 * blocks[0] - 3 * n
        - sum(size * (size - 4) * number
              for size, number in zip(sizes, blocks) if size >= 5)
    )
    high_incidences = tuple(
        size * blocks[size - 3] for size in sizes if size >= 5
    )
    return (defect, *high_incidences)


def row_feature(row: tuple[int, ...]) -> tuple[int, ...]:
    return (melchior_defect(row), *row[2:])


def reachable_targets(
    n: int,
    rows: tuple[tuple[int, ...], ...],
    targets: tuple[tuple[int, ...], ...],
) -> tuple[set[tuple[int, ...]], list[dict[str, object]]]:
    """Test each distinct global target by a tiny multiplicity model."""
    features = tuple(row_feature(row) for row in rows)
    reachable = set()
    audits = []
    for target_value in targets:
        model = cp_model.CpModel()
        multiplicities = tuple(
            model.NewIntVar(0, n, f"row_{index}")
            for index in range(len(features))
        )
        model.Add(sum(multiplicities) == n)
        for coordinate, required in enumerate(target_value):
            model.Add(sum(
                multiplicities[index] * feature[coordinate]
                for index, feature in enumerate(features)
            ) == required)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 5.0
        solver.parameters.num_search_workers = 1
        status = solver.Solve(model)
        status_name = solver.StatusName(status)
        if status_name not in ("INFEASIBLE", "FEASIBLE", "OPTIMAL"):
            raise RuntimeError(
                f"unresolved local-target solver status {status_name}: "
                f"{target_value}"
            )
        if status_name in ("FEASIBLE", "OPTIMAL"):
            reachable.add(target_value)
        audits.append({
            "target": list(target_value),
            "status": status_name,
            "wall_time_seconds": solver.WallTime(),
        })
    return reachable, audits


def verify(n: int, first: int, last: int) -> dict[str, object]:
    filt = SignatureFilter(n)
    sizes = filt.sizes
    prepared = []
    all_targets = []
    for circle_count in range(first, last + 1):
        raw = filt.signatures(circle_count)
        convexity = [
            signature for signature in raw
            if block_family_convexity_witness(n, sizes, signature) is None
        ]
        summed = (
            [signature for signature in convexity
             if summed_local_line_witness(n, sizes, signature) is None]
            if n >= 11 else convexity
        )
        prepared.append((circle_count, raw, convexity, summed))
        all_targets.extend(target(n, sizes, signature) for signature in summed)
    feature_dimension = len(sizes) - 1
    bounds = tuple(
        max(values) if values else 0
        for values in zip(*all_targets)
    ) if all_targets else (0,) * feature_dimension
    rows = tuple(
        row for row in projected_local_types(n, sizes[-1])
        if all(value <= bound
               for value, bound in zip(row_feature(row), bounds))
    )
    distinct_targets = tuple(sorted(set(all_targets)))
    states, target_audits = reachable_targets(n, rows, distinct_targets)
    layers = []
    for circle_count, raw, convexity, summed in prepared:
        survivors = []
        for signature in summed:
            totals = target(n, sizes, signature)
            if totals in states:
                survivors.append({
                    "signature": list(signature),
                    "target": list(totals),
                    "local_type_sum_reachable": True,
                })
        layers.append({
            "circle_count": circle_count,
            "raw_signature_count": len(raw),
            "after_block_family_convexity": len(convexity),
            "after_summed_local_line_inequalities": len(summed),
            "after_local_type_sum": len(survivors),
            "survivors": survivors,
        })
    return {
        "schema_version": 1,
        "status": "PASS",
        "claim": "Necessary aggregate local-arrangement filter.",
        "n": n,
        "sizes": list(sizes),
        "proof_scope": (
            "global signature inequalities, block-family convexity, summed "
            "local line inequalities, and target-wise integer balance of "
            "local types; no labelled supports"
        ),
        "local_type_count": len(rows),
        "reachable_target_count": len(states),
        "distinct_target_count": len(distinct_targets),
        "target_audits": target_audits,
        "ortools_version": ORTOOLS_VERSION,
        "coordinate_bounds": list(bounds),
        "layers": layers,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n", type=int, choices=(9, 10, 11, 12, 13, 14), required=True
    )
    parser.add_argument("--first-circle-count", type=int, required=True)
    parser.add_argument("--last-circle-count", type=int, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    result = verify(args.n, args.first_circle_count, args.last_circle_count)
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    if not args.quiet:
        print(rendered, end="")


if __name__ == "__main__":
    main()
