#!/usr/bin/env python3
"""Filter signatures by existence of every required block row."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

from block_intersection_rows import intersection_rows

HERE = Path(__file__).resolve().parent
AUTHORITATIVE = HERE
sys.path.insert(0, str(AUTHORITATIVE))
from verify_pair_moment_filter import conditioned_rows  # noqa: E402


# A preceding necessary-condition model may return UNKNOWN after its time
# limit.  Such a record has not been excluded and must therefore be carried
# into every later filter.  SURVIVES is used by exact, solver-free filters.
RETAINED_SOURCE_STATUSES = frozenset({
    "OPTIMAL", "FEASIBLE", "UNKNOWN", "SURVIVES",
})


def point_incidence_inequalities(
    n: int,
    signature: tuple[int, ...],
    gamma: int,
    rows: tuple[tuple[int, ...], ...],
) -> tuple[tuple[tuple[int, ...], int, int], ...]:
    """Return one- and two-coordinate bounds for points on one block."""
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    block_size = sizes[gamma % dimension]
    reduced = tuple(
        tuple(row[alpha] - int(alpha == gamma)
              for alpha in range(2 * dimension))
        for row in rows if row[gamma] > 0
    )
    if not reduced:
        return ()
    weights = []
    for alpha in range(2 * dimension):
        weight = [0] * (2 * dimension)
        weight[alpha] = 1
        weights.append(tuple(weight))
    for alpha in range(2 * dimension):
        for beta in range(alpha + 1, 2 * dimension):
            total = [0] * (2 * dimension)
            total[alpha] = total[beta] = 1
            weights.append(tuple(total))
            difference = [0] * (2 * dimension)
            difference[alpha] = 1
            difference[beta] = -1
            weights.append(tuple(difference))
    answer = []
    for weight in weights:
        values = tuple(
            sum(coefficient * value
                for coefficient, value in zip(weight, row))
            for row in reduced
        )
        answer.append((
            weight,
            block_size * min(values),
            block_size * max(values),
        ))
    return tuple(answer)


def point_incidence_pair_sets(
    signature: tuple[int, ...],
    gamma: int,
    rows: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, int, frozenset[tuple[int, int]]], ...]:
    """Exact two-coordinate sums of the local rows on one block."""
    dimension = len(signature) // 2
    block_size = gamma % dimension + 3
    reduced = tuple(
        tuple(row[alpha] - int(alpha == gamma)
              for alpha in range(2 * dimension))
        for row in rows if row[gamma] > 0
    )
    answer = []
    for alpha in range(2 * dimension):
        for beta in range(alpha + 1, 2 * dimension):
            options = {
                (row[alpha], row[beta]) for row in reduced
            }
            attainable = {(0, 0)}
            for _ in range(block_size):
                attainable = {
                    (left + add_left, right + add_right)
                    for left, right in attainable
                    for add_left, add_right in options
                }
            answer.append((alpha, beta, frozenset(attainable)))
    return tuple(answer)


def test(task: tuple[
    int, tuple[int, ...], bool, bool, bool, bool, bool,
]) -> dict[str, object]:
    (n, signature, use_point_capacities, use_point_pair_sums,
     use_deletion_bounds, use_subset_inheritance_bounds,
     use_conditioned_line_bounds) = task
    started = perf_counter()
    witnesses = []
    row_counts = {}
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    local_rows = (
        conditioned_rows(n, sizes, signature)
        if (use_point_capacities or use_point_pair_sums) else ()
    )
    for gamma, number in enumerate(signature):
        if number == 0:
            continue
        inequalities = (
            point_incidence_inequalities(
                n, signature, gamma, local_rows
            ) if use_point_capacities else ()
        )
        pair_sets = (
            point_incidence_pair_sets(signature, gamma, local_rows)
            if use_point_pair_sums else ()
        )
        if (use_point_capacities or use_point_pair_sums) and not any(
                row[gamma] > 0 for row in local_rows):
            count = 0
        else:
            count = len(intersection_rows(
                n, signature, gamma, 1, inequalities, pair_sets,
                use_deletion_bounds,
                use_subset_inheritance_bounds,
                use_conditioned_line_bounds,
            ))
        row_counts[str(gamma)] = count
        if count == 0:
            witnesses.append(gamma)
            break
    return {
        "signature": list(signature),
        "status": "INFEASIBLE" if witnesses else "SURVIVES",
        "empty_block_categories": witnesses,
        "tested_row_counts": row_counts,
        "wall_time_seconds": perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=range(9, 15), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--point-row-capacities", action="store_true")
    parser.add_argument("--point-row-pair-sums", action="store_true")
    parser.add_argument("--deletion-bounds", action="store_true")
    parser.add_argument("--subset-inheritance-bounds", action="store_true")
    parser.add_argument("--conditioned-line-bounds", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    source_records = (
        source["records"]
        if "records" in source
        else [
            record
            for layer in source.get("layers", [])
            for record in layer["records"]
        ]
    )
    source_status_counts = Counter(
        record["status"] for record in source_records
    )
    unexpected_statuses = set(source_status_counts) - (
        RETAINED_SOURCE_STATUSES
        | {"INFEASIBLE", "INFEASIBLE_NO_POINT_ROW"}
    )
    if unexpected_statuses:
        raise ValueError(
            "unrecognised source statuses: "
            + ", ".join(sorted(unexpected_statuses))
        )
    signatures = tuple(
        tuple(record["signature"])
        for record in source_records
        if record["status"] in RETAINED_SOURCE_STATUSES
    )
    tasks = tuple(
        (args.n, signature, args.point_row_capacities,
         args.point_row_pair_sums,
         args.deletion_bounds, args.subset_inheritance_bounds,
         args.conditioned_line_bounds)
        for signature in signatures
    )
    if args.jobs == 1:
        records = list(map(test, tasks))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            records = list(executor.map(test, tasks, chunksize=1))
    result = {
        "schema_version": 1,
        "status": "DIAGNOSTIC_ONLY",
        "n": args.n,
        "source_file": args.input.name,
        "source_sha256": sha256(args.input.read_bytes()).hexdigest(),
        "point_row_capacities": args.point_row_capacities,
        "point_row_pair_sums": args.point_row_pair_sums,
        "deletion_bounds": args.deletion_bounds,
        "subset_inheritance_bounds": args.subset_inheritance_bounds,
        "conditioned_line_bounds": args.conditioned_line_bounds,
        "retained_source_statuses": sorted(RETAINED_SOURCE_STATUSES),
        "source_status_counts": dict(source_status_counts),
        "claim": (
            "Every represented rich block must have a nonnegative integral "
            "intersection row satisfying the three exact triple strata and "
            "the cross-block matching inequality."
        ),
        "input_count": len(records),
        "status_counts": dict(Counter(r["status"] for r in records)),
        "wall_time_seconds_sum": sum(
            float(r["wall_time_seconds"]) for r in records
        ),
        "records": records,
    }
    args.output.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        key: result[key] for key in (
            "n", "input_count", "status_counts", "wall_time_seconds_sum"
        )
    }, indent=2))


if __name__ == "__main__":
    main()
