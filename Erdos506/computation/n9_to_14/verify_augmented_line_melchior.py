#!/usr/bin/env python3
"""Verify the augmented-line Melchior inequality and point-row balance.

For every point p, inversion followed by projective duality gives an
essential arrangement of n-1 lines.  Adjoining the line dual to the inversion
centre turns every original rich line through p into an old vertex on the
new line.  Melchior's inequality for the enlarged arrangement gives

    sum_{k>=3} k d^L_k(p) <= n-1 + delta(p).

The first stage below checks the sum of these inequalities directly from a
global signature.  The optional second stage asks only whether n admissible
point rows can realize the global incidence margins.  It never constructs
labelled blocks or geometric orbits.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path

from low_defect_point_support_filter import (
    augmented_line_melchior_row_ok,
    conditioned_rows,
    melchior_defect,
    retained_records,
)
from ortools.sat.python import cp_model


def total_defect(n: int, signature: tuple[int, ...]) -> int:
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    line = signature[:dimension]
    circle = signature[dimension:]
    return (
        3 * (line[0] + circle[0]) - 3 * n
        - sum(
            (size - 4) * size * (line[index] + circle[index])
            for index, size in enumerate(sizes) if size >= 5
        )
    )


def summed_gap(n: int, signature: tuple[int, ...]) -> int:
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    line = signature[:dimension]
    return (
        sum(size * size * number for size, number in zip(sizes, line))
        - n * (n - 1) - total_defect(n, signature)
    )


def solve(task):
    n, signature, seconds, point_rows = task
    gap = summed_gap(n, signature)
    if gap > 0:
        return {
            "signature": list(signature),
            "status": "INFEASIBLE_SUMMED_INEQUALITY",
            "summed_gap": gap,
            "row_count": 0,
            "wall_time_seconds": 0.0,
        }
    if not point_rows:
        return {
            "signature": list(signature),
            "status": "SURVIVES_SUMMED_INEQUALITY",
            "summed_gap": gap,
            "row_count": None,
            "wall_time_seconds": 0.0,
        }

    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    rows = tuple(
        row for row in conditioned_rows(n, sizes, signature)
        if augmented_line_melchior_row_ok(n, row)
    )
    if not rows:
        return {
            "signature": list(signature),
            "status": "INFEASIBLE_NO_POINT_ROW",
            "summed_gap": gap,
            "row_count": 0,
            "wall_time_seconds": 0.0,
        }
    model = cp_model.CpModel()
    multiplicities = tuple(
        model.NewIntVar(0, n, f"row_{index}")
        for index in range(len(rows))
    )
    model.Add(sum(multiplicities) == n)
    for coordinate, (size, number) in enumerate(zip(sizes + sizes, signature)):
        model.Add(sum(
            multiplicity * row[coordinate]
            for multiplicity, row in zip(multiplicities, rows)
        ) == size * number)
    model.Add(sum(
        multiplicity * melchior_defect(tuple(
            row[index] + row[dimension + index]
            for index in range(dimension)
        ))
        for multiplicity, row in zip(multiplicities, rows)
    ) == total_defect(n, signature))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    status = solver.Solve(model)
    result = {
        "signature": list(signature),
        "status": solver.StatusName(status),
        "summed_gap": gap,
        "row_count": len(rows),
        "wall_time_seconds": solver.WallTime(),
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        result["point_row_multiplicities"] = [
            {"row": list(row), "multiplicity": solver.Value(variable)}
            for row, variable in zip(rows, multiplicities)
            if solver.Value(variable)
        ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=5.0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--summed-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source_records = retained_records(args.n, args.input_json)
    tasks = [
        (args.n, tuple(record["signature"]), args.seconds,
         not args.summed_only)
        for record in source_records
    ]
    if args.jobs == 1:
        outcomes = list(map(solve, tasks))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            outcomes = list(executor.map(solve, tasks, chunksize=1))
    excluded = {
        "INFEASIBLE", "INFEASIBLE_NO_POINT_ROW",
        "INFEASIBLE_SUMMED_INEQUALITY",
    }
    retained = [record for record in outcomes
                if record["status"] not in excluded]
    payload = {
        "schema_version": 1,
        "status": "PASS",
        "claim": (
            "The summed augmented-line Melchior inequality is checked "
            "exactly; when requested, the remaining signatures are checked "
            "by integer point-row balance only.  UNKNOWN is retained."
        ),
        "n": args.n,
        "source_file": args.input_json.name,
        "source_sha256": sha256(args.input_json.read_bytes()).hexdigest(),
        "input_signature_count": len(outcomes),
        "summed_only": args.summed_only,
        "status_counts": dict(Counter(record["status"] for record in outcomes)),
        "retained_signature_count": len(retained),
        "records": outcomes,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "input_signature_count": len(outcomes),
        "status_counts": payload["status_counts"],
        "retained_signature_count": len(retained),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
