#!/usr/bin/env python3
"""Recursive inversion-row verifier with an exact nine-point base.

The verifier never enumerates global block orbits at orders ten and eleven.
For a point x, inversion gives a complete signature on one fewer point:
blocks through x become lines and blocks avoiding x become circles.  A point
row is retained only if that child signature survives the next level.  At
order nine the necessary abstract block system is checked exactly on nine
labelled points.  UNKNOWN outcomes are retained at every level.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import time

from low_defect_point_support_filter import (
    augmented_line_arrangement_row_ok,
    conditioned_rows,
    retained_records,
)
from ortools.sat.python import cp_model


SIZES = (3, 4, 5, 6)
KNOWN = {4: 3, 5: 5, 6: 8, 7: 11, 8: 17}
INFEASIBLE = frozenset({"INFEASIBLE", "INFEASIBLE_NO_POINT_ROW"})


def ordinary_line_count(n: int, rich_signature: tuple[int, ...]) -> int:
    line = rich_signature[:4]
    return comb(n, 2) - sum(
        comb(size, 2) * number for size, number in zip(SIZES, line)
    )


def padded_rich_signature(signature: tuple[int, ...]) -> tuple[int, ...]:
    """Pad line and circle coordinates to the common sizes 3,4,5,6."""
    dimension = len(signature) // 2
    if len(signature) != 2 * dimension or dimension > len(SIZES):
        raise ValueError("invalid rich-block signature")
    padding = (0,) * (len(SIZES) - dimension)
    return signature[:dimension] + padding + signature[dimension:] + padding


def child_signature(
    n: int, rich_signature: tuple[int, ...], row: tuple[int, ...],
) -> tuple[int, ...]:
    """Signature of the inverted (n-1)-point configuration."""
    totals = tuple(
        rich_signature[index] + rich_signature[4 + index]
        for index in range(4)
    )
    through = tuple(row[index] + row[4 + index] for index in range(4))
    child_lines = through[1:] + (0,)
    child_circles = tuple(
        totals[index] - through[index] for index in range(4)
    )
    return child_lines + child_circles + (through[0],)


def candidate_rows(n: int, full_signature: tuple[int, ...]):
    rich_signature = full_signature[:8]
    rows = []
    children = []
    for row in conditioned_rows(n, SIZES, rich_signature):
        if not augmented_line_arrangement_row_ok(n, row):
            continue
        rows.append(row)
        children.append(child_signature(n, rich_signature, row))
    return tuple(rows), tuple(children)


def exact_nine_point_signature(full_signature: tuple[int, ...], seconds: float):
    start = time.perf_counter()
    rich_signature = full_signature[:8]
    recorded_l2 = full_signature[8]
    if ordinary_line_count(9, rich_signature) != recorded_l2:
        return {
            "signature": list(full_signature), "status": "INFEASIBLE",
            "reason": "ordinary_line_count", "wall_time_seconds": 0.0,
            "variables": 0, "constraints": 0,
        }
    line_counts = rich_signature[:4]
    circle_counts = rich_signature[4:]
    points = tuple(range(9))
    by_size = {size: tuple(combinations(points, size)) for size in SIZES}
    blocks = tuple(block for size in SIZES for block in by_size[size])
    block_sets = {block: frozenset(block) for block in blocks}
    model = cp_model.CpModel()
    line = {block: model.NewBoolVar("L_" + "_".join(map(str, block)))
            for block in blocks}
    circle = {block: model.NewBoolVar("C_" + "_".join(map(str, block)))
              for block in blocks}
    for triple in combinations(points, 3):
        triple_set = frozenset(triple)
        model.Add(sum(
            line[block] + circle[block] for block in blocks
            if triple_set.issubset(block_sets[block])
        ) == 1)
    for size, number in zip(SIZES, line_counts):
        model.Add(sum(line[block] for block in by_size[size]) == number)
    for size, number in zip(SIZES, circle_counts):
        model.Add(sum(circle[block] for block in by_size[size]) == number)
    for pair in combinations(points, 2):
        pair_set = frozenset(pair)
        model.Add(sum(
            line[block] for block in blocks
            if pair_set.issubset(block_sets[block])
        ) <= 1)
        model.Add(sum(
            (len(block) - 2) * (line[block] + circle[block])
            for block in blocks if pair_set.issubset(block_sets[block])
        ) == 7)

    # Every subset inherits the known lower bound unless it is itself
    # collinear or concyclic.  For a concyclic exceptional subset its defining
    # circle remains, which gives the coefficient KNOWN[m]-1.
    for subset_size in range(4, 9):
        for subset in combinations(points, subset_size):
            subset_set = frozenset(subset)
            inherited = sum(
                circle[block] for block in blocks
                if len(block_sets[block] & subset_set) >= 3
            )
            if subset_size >= 7:
                model.Add(inherited >= KNOWN[subset_size])
            else:
                supersets = [
                    block for block in blocks
                    if subset_set.issubset(block_sets[block])
                ]
                model.Add(
                    inherited >= KNOWN[subset_size]
                    - KNOWN[subset_size]
                    * sum(line[block] for block in supersets)
                    - (KNOWN[subset_size] - 1)
                    * sum(circle[block] for block in supersets)
                )

    anchor_kind = None
    anchor_size = None
    for size in reversed(SIZES):
        if line_counts[size - 3]:
            anchor_kind, anchor_size = "line", size
            break
        if circle_counts[size - 3]:
            anchor_kind, anchor_size = "circle", size
            break
    if anchor_size is None:
        return {
            "signature": list(full_signature), "status": "INFEASIBLE",
            "reason": "no_rich_block", "wall_time_seconds": 0.0,
            "variables": len(model.Proto().variables),
            "constraints": len(model.Proto().constraints),
        }
    anchor = tuple(range(anchor_size))
    model.Add((line if anchor_kind == "line" else circle)[anchor] == 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.symmetry_level = 2
    status = solver.Solve(model)
    return {
        "signature": list(full_signature),
        "status": solver.StatusName(status),
        "reason": "exact_abstract_support",
        "wall_time_seconds": time.perf_counter() - start,
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "anchor": {"kind": anchor_kind, "size": anchor_size},
    }


def solve_base_task(task):
    return exact_nine_point_signature(*task)


def balance_signature(
    n: int,
    full_signature: tuple[int, ...],
    child_status: dict[tuple[int, ...], str],
    seconds: float,
):
    rows, children = candidate_rows(n, full_signature)
    retained = [
        (row, child) for row, child in zip(rows, children)
        if child_status.get(child, "UNKNOWN") not in INFEASIBLE
    ]
    if not retained:
        return {
            "signature": list(full_signature),
            "status": "INFEASIBLE_NO_POINT_ROW",
            "candidate_row_count": len(rows), "retained_row_count": 0,
            "wall_time_seconds": 0.0,
        }
    model = cp_model.CpModel()
    variables = tuple(
        model.NewIntVar(0, n, f"row_{index}")
        for index in range(len(retained))
    )
    model.Add(sum(variables) == n)
    rich_signature = full_signature[:8]
    for coordinate, (size, number) in enumerate(
            zip(SIZES + SIZES, rich_signature)):
        model.Add(sum(
            variable * row[coordinate]
            for variable, (row, _child) in zip(variables, retained)
        ) == size * number)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    status = solver.Solve(model)
    return {
        "signature": list(full_signature),
        "status": solver.StatusName(status),
        "candidate_row_count": len(rows),
        "retained_row_count": len(retained),
        "distinct_child_count": len(set(child for _row, child in retained)),
        "wall_time_seconds": solver.WallTime(),
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(9, 10, 11, 12), required=True)
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = retained_records(args.n, args.input_json)
    current = set()
    for record in source:
        rich_signature = padded_rich_signature(tuple(record["signature"]))
        current.add(
            rich_signature
            + (ordinary_line_count(args.n, rich_signature),)
        )
    levels = {args.n: current}
    row_maps = {}
    prebalance_counts = {}
    for order in range(args.n, 9, -1):
        children = set()
        counts = Counter()
        for signature in levels[order]:
            preliminary = balance_signature(
                order, signature, {}, args.seconds
            )
            counts[preliminary["status"]] += 1
            if preliminary["status"] in INFEASIBLE:
                row_maps[order, signature] = ((), ())
                continue
            rows, child_signatures = candidate_rows(order, signature)
            row_maps[order, signature] = (rows, child_signatures)
            children.update(child_signatures)
        levels[order - 1] = children
        prebalance_counts[order] = dict(counts)

    base_tasks = [(signature, args.seconds) for signature in sorted(levels[9])]
    if args.jobs == 1:
        base_records = list(map(solve_base_task, base_tasks))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            base_records = list(executor.map(
                solve_base_task, base_tasks, chunksize=1
            ))
    statuses = {
        tuple(record["signature"]): record["status"]
        for record in base_records
    }
    level_records = {9: base_records}
    for order in range(10, args.n + 1):
        records = [
            balance_signature(order, signature, statuses, args.seconds)
            for signature in sorted(levels[order])
        ]
        statuses = {
            tuple(record["signature"]): record["status"]
            for record in records
        }
        level_records[order] = records

    source_records = level_records[args.n]
    retained_final = [record for record in source_records
                      if record["status"] not in INFEASIBLE]
    payload = {
        "schema_version": 1,
        "status": "PASS" if not retained_final else "INCOMPLETE",
        "claim": (
            "Recursive inversion-row balance down to an exact nine-point "
            "abstract-support base; only INFEASIBLE statuses exclude and "
            "UNKNOWN is retained."
        ),
        "n": args.n,
        "source_file": args.input_json.name,
        "source_sha256": sha256(args.input_json.read_bytes()).hexdigest(),
        "input_signature_count": len(source),
        "seconds_per_base_signature": args.seconds,
        "level_signature_counts": {
            str(order): len(signatures) for order, signatures in levels.items()
        },
        "level_status_counts": {
            str(order): dict(Counter(record["status"] for record in records))
            for order, records in level_records.items()
        },
        "pre_recursion_balance_status_counts": {
            str(order): counts for order, counts in prebalance_counts.items()
        },
        "retained_signature_count": len(retained_final),
        "levels": {
            str(order): records for order, records in level_records.items()
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "input_signature_count": len(source),
        "level_signature_counts": payload["level_signature_counts"],
        "level_status_counts": payload["level_status_counts"],
        "retained_signature_count": len(retained_final),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
