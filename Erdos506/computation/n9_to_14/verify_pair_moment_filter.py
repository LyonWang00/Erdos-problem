#!/usr/bin/env python3
"""Uniform support-free point/pair moment filter for small n.

For every global line/circle signature, the verifier balances multiplicities
of conditioned point rows and exact point-pair profiles.  Their first two
intersection moments are coupled.  For every distinguished block family it
then counts, in two ways, triples disjoint from a distinguished block.  The
model has no labelled point or block support variables.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from hashlib import sha256
from math import comb
import json
from pathlib import Path
import sys

import numpy as np

COMPUTATIONS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(COMPUTATIONS))
from runtime_dependencies import activate_ortools  # noqa: E402

activate_ortools()
import ortools  # noqa: E402

from pair_profile_generator import profiles as pair_profiles
from projected_arrangement_model import category_point_rows
from signature_filters import (
    KNOWN,
    conditioned_line_subset_ok,
    conditioned_point_row_ok,
    conditioned_subset_coefficients,
)

ORTOOLS_VERSION = ortools.__version__
if ORTOOLS_VERSION != "9.15.6755":
    raise RuntimeError(
        f"expected OR-Tools 9.15.6755, found {ORTOOLS_VERSION}"
    )
from ortools.sat.python import cp_model  # noqa: E402


def choose(n: int, r: int) -> int:
    return comb(n, r) if 0 <= r <= n else 0


def moment_column(
    rows: tuple[tuple[int, ...], ...],
    alpha: int,
    beta: int,
) -> tuple[int, ...]:
    """Return the coefficient column of an unordered incidence moment."""
    if alpha == beta:
        return tuple(choose(row[alpha], 2) for row in rows)
    return tuple(row[alpha] * row[beta] for row in rows)


def outside_triple_coefficients(
    signature: tuple[int, ...],
    category_sizes: tuple[int, ...],
    point_rows: tuple[tuple[int, ...], ...],
    pair_rows: tuple[tuple[int, ...], ...],
    gamma: int,
) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    """Exact coefficient vector for triples disjoint from a gamma block.

    For a category ``alpha`` of block size ``b``, its contribution is

        C(b,3) Q - C(b-1,2) M_1 + (b-2) M_2.

    When ``alpha == gamma`` the moment functions count unordered pairs of
    distinct blocks, whereas the displayed identity is for ordered pairs;
    hence both moment columns have factor two.  Returning plain integer
    vectors avoids aliasing mutable OR-Tools ``LinearExpr`` objects.
    """
    point_coefficients = [0] * len(point_rows)
    pair_coefficients = [0] * len(pair_rows)
    constant = 0
    for alpha, size in enumerate(category_sizes):
        if alpha == gamma:
            number = signature[gamma] * (signature[gamma] - 1)
            factor = 2
        else:
            number = signature[gamma] * signature[alpha]
            factor = 1
        constant += choose(size, 3) * number
        point_column = moment_column(point_rows, gamma, alpha)
        pair_column = moment_column(pair_rows, gamma, alpha)
        point_weight = -choose(size - 1, 2) * factor
        pair_weight = (size - 2) * factor
        for index, coefficient in enumerate(point_column):
            point_coefficients[index] += point_weight * coefficient
        for index, coefficient in enumerate(pair_column):
            pair_coefficients[index] += pair_weight * coefficient
    return constant, tuple(point_coefficients), tuple(pair_coefficients)


@lru_cache(maxsize=None)
def point_row_array(
    n: int, maximum: int,
) -> tuple[tuple[tuple[int, ...], ...], np.ndarray]:
    rows = category_point_rows(n, maximum)
    matrix = np.asarray(rows, dtype=np.int64)
    matrix.setflags(write=False)
    return rows, matrix


def conditioned_rows(
    n: int,
    sizes: tuple[int, ...],
    signature: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    """Vectorized form of the exact scalar point-row filter.

    The two child-signature membership tests specific to n=12,13 are still
    checked by the scalar reference function on the already reduced rows.
    """
    dimension = len(sizes)
    rows, matrix = point_row_array(n, sizes[-1])
    signature_array = np.asarray(signature, dtype=np.int64)
    lines = signature_array[:dimension]
    circles = signature_array[dimension:]
    local_lines = matrix[:, :dimension]
    local_circles = matrix[:, dimension:]
    mask = np.all(matrix <= signature_array, axis=1)

    ordinary_lines = choose(n, 2) - sum(
        choose(size, 2) * int(number)
        for size, number in zip(sizes, lines)
    )
    ordinary_through = n - 1 - local_lines @ np.asarray(
        tuple(size - 1 for size in sizes), dtype=np.int64
    )
    mask &= (ordinary_through >= 0) & (ordinary_through <= ordinary_lines)

    previous = KNOWN[n - 1]
    inverse_delete = KNOWN[n - 2]
    circle_count = int(circles.sum())
    block_count = int(signature_array.sum())
    local_block_count = matrix.sum(axis=1)
    mask &= local_circles[:, 0] <= circle_count - previous
    mask &= block_count - local_block_count >= previous
    global_three = int(lines[0] + circles[0])
    local_three = local_lines[:, 0] + local_circles[:, 0]
    mask &= (
        (n - 1) * (block_count - local_block_count)
        - 3 * (global_three - local_three)
        >= (n - 1) * inverse_delete
    )

    indices = np.flatnonzero(mask)
    if len(indices):
        coefficients = conditioned_subset_coefficients(
            n, sizes, tuple(sorted(KNOWN.items()))
        )
        local_coefficients = np.asarray(
            [line + circle for _constant, _global_line, _global_circle,
             line, circle in coefficients],
            dtype=np.int64,
        )
        thresholds = np.asarray([
            constant
            - sum(int(value) * coefficient
                  for value, coefficient in zip(lines, global_line))
            - sum(int(value) * coefficient
                  for value, coefficient in zip(circles, global_circle))
            for constant, global_line, global_circle,
            _local_line, _local_circle in coefficients
        ], dtype=np.int64)
        subset_left = matrix[indices] @ local_coefficients.T
        indices = indices[np.all(subset_left >= thresholds, axis=1)]

    selected = tuple(rows[int(index)] for index in indices)
    selected = tuple(
        row for row in selected
        if conditioned_line_subset_ok(n, sizes, signature, row)
    )
    if n in (12, 13):
        selected = tuple(
            row for row in selected
            if conditioned_point_row_ok(n, sizes, signature, row, KNOWN)
        )
    return selected


def audit_signature(
    n: int,
    signature: tuple[int, ...],
    seconds: float,
) -> dict[str, object]:
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    maximum = sizes[-1]
    rows = conditioned_rows(n, sizes, signature)
    if not rows:
        return {
            "signature": list(signature),
            "point_row_count": 0,
            "status": "INFEASIBLE_NO_POINT_ROW",
            "wall_time_seconds": 0.0,
        }
    pair_types = pair_profiles(n, maximum)
    model = cp_model.CpModel()
    z = tuple(model.NewIntVar(0, n, f"point_{index}")
              for index in range(len(rows)))
    model.Add(sum(z) == n)
    category_sizes = sizes + sizes
    for coordinate, (size, number) in enumerate(
            zip(category_sizes, signature)):
        model.Add(sum(z[index] * row[coordinate]
                      for index, row in enumerate(rows)) == size * number)

    y = tuple(model.NewIntVar(0, choose(n, 2), f"pair_{index}")
              for index in range(len(pair_types)))
    model.Add(sum(y) == choose(n, 2))
    for coordinate, (size, number) in enumerate(
            zip(category_sizes, signature)):
        model.Add(sum(y[index] * row[coordinate]
                      for index, row in enumerate(pair_types))
                  == choose(size, 2) * number)

    def point_moment(alpha: int, beta: int):
        if alpha == beta:
            return sum(z[index] * choose(row[alpha], 2)
                       for index, row in enumerate(rows))
        return sum(z[index] * row[alpha] * row[beta]
                   for index, row in enumerate(rows))

    def pair_moment(alpha: int, beta: int):
        if alpha == beta:
            return sum(y[index] * choose(row[alpha], 2)
                       for index, row in enumerate(pair_types))
        return sum(y[index] * row[alpha] * row[beta]
                   for index, row in enumerate(pair_types))

    coupling_constraint_count = 0
    for alpha in range(2 * dimension):
        for beta in range(alpha, 2 * dimension):
            first = point_moment(alpha, beta)
            second = pair_moment(alpha, beta)
            block_pairs = (
                choose(signature[alpha], 2) if alpha == beta
                else signature[alpha] * signature[beta]
            )
            model.Add(2 * second <= first)
            model.Add(first <= second + block_pairs)
            coupling_constraint_count += 2

    outside_triple_constraint_count = 0
    for gamma, distinguished_size in enumerate(category_sizes):
        if signature[gamma] == 0:
            continue
        constant, point_coefficients, pair_coefficients = (
            outside_triple_coefficients(
                signature, category_sizes, rows, pair_types, gamma
            )
        )
        outside_triples = (
            constant
            + cp_model.LinearExpr.weighted_sum(z, point_coefficients)
            + cp_model.LinearExpr.weighted_sum(y, pair_coefficients)
        )
        model.Add(
            outside_triples
            == signature[gamma] * choose(n - distinguished_size, 3)
        )
        outside_triple_constraint_count += 1

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    status = solver.Solve(model)
    return {
        "signature": list(signature),
        "point_row_count": len(rows),
        "pair_profile_count": len(pair_types),
        "status": solver.StatusName(status),
        "coupling_constraint_count": coupling_constraint_count,
        "outside_triple_constraint_count": outside_triple_constraint_count,
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
        "wall_time_seconds": solver.WallTime(),
    }


def audit_signature_task(
    task: tuple[int, tuple[int, ...], float],
) -> dict[str, object]:
    return audit_signature(*task)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--first-circle-count", type=int)
    parser.add_argument("--last-circle-count", type=int)
    parser.add_argument("--seconds", type=float, default=5.0)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    if source.get("status") != "PASS":
        raise ValueError("the input signature filter did not pass its audit")
    first = args.first_circle_count
    last = args.last_circle_count
    source_layers = source.get("layers")
    if source_layers is None and isinstance(source.get("records"), list):
        dimension = len(source["records"][0]["signature"]) // 2 if source["records"] else 0
        grouped = {}
        for record in source["records"]:
            if str(record.get("status", "SURVIVES")).startswith("INFEASIBLE"):
                continue
            circle_count = sum(record["signature"][dimension:])
            grouped.setdefault(circle_count, []).append(record)
        source_layers = [
            {"circle_count": circle_count, "survivors": records}
            for circle_count, records in sorted(grouped.items())
        ]
    if not isinstance(source_layers, list):
        raise ValueError("the input has neither layers nor records")
    selected_layers = []
    flat_tasks = []
    for layer in source_layers:
        circle_count = layer["circle_count"]
        if first is not None and circle_count < first:
            continue
        if last is not None and circle_count > last:
            continue
        tasks = [
            (args.n, tuple(entry["signature"]), args.seconds)
            for entry in layer["survivors"]
        ]
        selected_layers.append((circle_count, len(tasks)))
        flat_tasks.extend(tasks)
    if args.jobs == 1:
        flat_records = [audit_signature_task(task) for task in flat_tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            flat_records = list(executor.map(
                audit_signature_task, flat_tasks, chunksize=1
            ))

    layers = []
    offset = 0
    for circle_count, task_count in selected_layers:
        records = flat_records[offset:offset + task_count]
        offset += task_count
        status_counts = Counter(record["status"] for record in records)
        layers.append({
            "circle_count": circle_count,
            "input_signature_count": len(records),
            "status_counts": dict(status_counts),
            "survivor_count": sum(
                status in ("FEASIBLE", "OPTIMAL", "UNKNOWN")
                for status in (record["status"] for record in records)
            ),
            "survivors": [
                record for record in records
                if record["status"] in ("FEASIBLE", "OPTIMAL", "UNKNOWN")
            ],
            "records": records,
        })
    if offset != len(flat_records):
        raise AssertionError((offset, len(flat_records)))
    unknown = sum(
        layer["status_counts"].get("UNKNOWN", 0) for layer in layers
    )
    result = {
        "schema_version": 1,
        "status": "PASS" if unknown == 0 else "INCOMPLETE",
        "claim": (
            "Necessary support-free point/pair moment and outside-triple "
            "filter."
        ),
        "n": args.n,
        "source_file": args.input.name,
        "source_sha256": sha256(args.input.read_bytes()).hexdigest(),
        "input_signature_count": len(flat_records),
        "proof_scope": (
            "conditioned point rows, exact pair profiles, coupled first two "
            "intersection moments, and the exact count of triples disjoint "
            "from a distinguished block; no labelled supports"
        ),
        "ortools_version": ORTOOLS_VERSION,
        "layers": layers,
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    if not args.quiet:
        print(rendered, end="")


if __name__ == "__main__":
    main()
