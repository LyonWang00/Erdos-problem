#!/usr/bin/env python3
"""Unified orbit classification for the eight-point boundary layers 12--16.

The largest-block estimate treats every block of size at least five before
this script is invoked.  Hence the labelled model needs only three- and
four-point lines and four-point circles.  Each circle count is split into the
two exhaustive branches: a four-point line exists (and is relabelled 0123),
or no four-point line exists.
"""

from __future__ import annotations

import argparse
from itertools import combinations, permutations
import json
from pathlib import Path
import sys
import time


HERE = Path(__file__).resolve().parent
COMPUTATIONS = HERE.parents[1]
sys.path.insert(0, str(COMPUTATIONS))
from runtime_dependencies import activate_ortools  # noqa: E402

activate_ortools()
from ortools.sat.python import cp_model  # noqa: E402


N = 8
POINTS = tuple(range(N))
TRIPLES = tuple(combinations(POINTS, 3))
FOURS = tuple(combinations(POINTS, 4))
PERMUTATIONS = tuple(permutations(POINTS))
FIXED_LINE = (0, 1, 2, 3)
EXPECTED_CLASS_COUNTS = {
    12: 1,
    13: 1,
    14: 1,
    15: 3,
    16: 3,
}


def block_string(block) -> str:
    return "".join(map(str, block))


def parse_blocks(strings):
    return tuple(tuple(int(character) for character in text)
                 for text in strings)


def transform(block, permutation):
    return tuple(sorted(permutation[index] for index in block))


def canonical(lines, circles):
    answer = None
    for permutation in PERMUTATIONS:
        row = (
            tuple(sorted("L" + block_string(transform(block, permutation))
                         for block in lines))
            + ("|",)
            + tuple(sorted("C" + block_string(transform(block, permutation))
                           for block in circles))
        )
        if answer is None or row < answer:
            answer = row
    return answer


def build_model(circle_count: int, branch: str):
    model = cp_model.CpModel()
    line3 = {
        block: model.NewBoolVar("L3_" + block_string(block))
        for block in TRIPLES
    }
    line4 = {
        block: model.NewBoolVar("L4_" + block_string(block))
        for block in FOURS
    }
    circle4 = {
        block: model.NewBoolVar("C4_" + block_string(block))
        for block in FOURS
    }
    selected = (
        [("L3", block, line3[block], 1) for block in TRIPLES]
        + [("L4", block, line4[block], 4) for block in FOURS]
        + [("C4", block, circle4[block], 4) for block in FOURS]
    )

    for triple in TRIPLES:
        model.Add(sum(
            variable for _kind, block, variable, _weight in selected
            if set(triple).issubset(block)
        ) <= 1)

    line_blocks = (
        [(block, line3[block]) for block in TRIPLES]
        + [(block, line4[block]) for block in FOURS]
    )
    for index, (first, first_variable) in enumerate(line_blocks):
        first_set = set(first)
        for second, second_variable in line_blocks[index + 1:]:
            if len(first_set & set(second)) >= 2:
                model.Add(first_variable + second_variable <= 1)

    l3 = sum(line3.values())
    l4 = sum(line4.values())
    q4 = sum(circle4.values())
    l2 = 28 - 3 * l3 - 6 * l4
    model.Add(l3 <= 7)              # eight-point orchard bound
    model.Add(l2 >= 3 + l4)         # Melchior
    model.Add(56 - l3 - 4 * l4 - 3 * q4 == circle_count)

    if branch == "fixed_four_line":
        model.Add(line4[FIXED_LINE] == 1)
    elif branch == "no_four_line":
        for variable in line4.values():
            model.Add(variable == 0)
    else:
        raise ValueError(branch)
    return model, {"L3": line3, "L4": line4, "C4": circle4}


def extract_record(solver, variables):
    lines = tuple(sorted(
        [block for block, variable in variables["L3"].items()
         if solver.Value(variable)]
        + [block for block, variable in variables["L4"].items()
           if solver.Value(variable)],
        key=lambda block: (len(block), block),
    ))
    circles = tuple(sorted(
        block for block, variable in variables["C4"].items()
        if solver.Value(variable)
    ))
    owned = {
        triple
        for block in lines + circles
        for triple in combinations(block, 3)
    }
    l3 = sum(len(block) == 3 for block in lines)
    l4 = sum(len(block) == 4 for block in lines)
    q4 = len(circles)
    q3 = len(TRIPLES) - len(owned)
    return {
        "stats": {
            "l3": l3,
            "l4": l4,
            "q3": q3,
            "q4": q4,
            "c": q3 + q4,
        },
        "lines": [block_string(block) for block in lines],
        "circle4": [block_string(block) for block in circles],
    }


def verify_record(record, expected_circle_count: int) -> None:
    lines = parse_blocks(record["lines"])
    circles = parse_blocks(record["circle4"])
    owners = {}
    for kind, blocks in (("line", lines), ("circle", circles)):
        for block in blocks:
            for triple in combinations(block, 3):
                if triple in owners:
                    raise AssertionError((triple, owners[triple], kind, block))
                owners[triple] = (kind, block)
    l3 = sum(len(block) == 3 for block in lines)
    l4 = sum(len(block) == 4 for block in lines)
    q4 = len(circles)
    q3 = len(TRIPLES) - len(owners)
    expected = {"l3": l3, "l4": l4, "q3": q3,
                "q4": q4, "c": q3 + q4}
    if expected != record["stats"] or q3 + q4 != expected_circle_count:
        raise AssertionError(record)


def add_orbit_block(model, variables, record, branch: str):
    true_base = {
        "L3": {block for block in parse_blocks(record["lines"])
               if len(block) == 3},
        "L4": {block for block in parse_blocks(record["lines"])
               if len(block) == 4},
        "C4": set(parse_blocks(record["circle4"])),
    }
    universes = {"L3": set(TRIPLES), "L4": set(FOURS), "C4": set(FOURS)}
    seen = set()
    clauses = 0
    literals = 0
    for permutation in PERMUTATIONS:
        true_sets = {
            kind: {transform(block, permutation) for block in blocks}
            for kind, blocks in true_base.items()
        }
        if (branch == "fixed_four_line"
                and FIXED_LINE not in true_sets["L4"]):
            continue
        key = tuple(
            tuple(sorted(block_string(block) for block in true_sets[kind]))
            for kind in ("L3", "L4", "C4")
        )
        if key in seen:
            continue
        seen.add(key)
        clause = []
        for kind in ("L3", "L4", "C4"):
            for block in universes[kind]:
                variable = variables[kind][block]
                clause.append(
                    variable.Not() if block in true_sets[kind] else variable
                )
        model.AddBoolOr(clause)
        clauses += 1
        literals += len(clause)
    return clauses, literals


def classify_branch(circle_count: int, branch: str, seconds: float):
    started = time.time()
    model, variables = build_model(circle_count, branch)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.cp_model_presolve = True
    records = {}
    iterations = []
    while True:
        status = solver.Solve(model)
        status_name = solver.StatusName(status)
        if status_name == "INFEASIBLE":
            break
        if status_name not in ("OPTIMAL", "FEASIBLE"):
            raise RuntimeError((circle_count, branch, status_name))
        record = extract_record(solver, variables)
        verify_record(record, circle_count)
        key = canonical(
            parse_blocks(record["lines"]),
            parse_blocks(record["circle4"]),
        )
        records[key] = record
        clauses, literals = add_orbit_block(
            model, variables, record, branch
        )
        iterations.append({
            "status": status_name,
            "record": record,
            "orbit_blocking_clauses": clauses,
            "orbit_blocking_literals": literals,
        })
    return {
        "branch": branch,
        "final_status": "INFEASIBLE",
        "elapsed_seconds": time.time() - started,
        "iterations": iterations,
        "unique_records": sorted(
            records.values(),
            key=lambda record: (
                record["stats"]["l4"], record["stats"]["l3"],
                record["lines"], record["circle4"],
            ),
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int, default=12)
    parser.add_argument("--last", type=int, default=16)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.first < 12 or args.last > 16 or args.first > args.last:
        raise ValueError("the verified boundary range is 12 <= c <= 16")
    layers = []
    for circle_count in range(args.first, args.last + 1):
        branches = [
            classify_branch(circle_count, branch, args.seconds)
            for branch in ("fixed_four_line", "no_four_line")
        ]
        records = {
            canonical(parse_blocks(record["lines"]),
                      parse_blocks(record["circle4"])): record
            for branch in branches for record in branch["unique_records"]
        }
        layers.append({
            "circle_count": circle_count,
            "branches": branches,
            "class_count": len(records),
            "classes": sorted(
                records.values(),
                key=lambda record: (
                    record["stats"]["l4"], record["stats"]["l3"],
                    record["lines"], record["circle4"],
                ),
            ),
        })
        expected = EXPECTED_CLASS_COUNTS[circle_count]
        if len(records) != expected:
            raise AssertionError(
                (circle_count, len(records), expected)
            )
        print(circle_count, len(records), flush=True)
    result = {
        "schema_version": 1,
        "status": "PASS",
        "scope": "all abstract eight-point boundary layers 12<=c<=16 after the mathematical c>=12 and size-five largest-block reductions",
        "uses_five_point_block_variables": False,
        "expected_class_counts": EXPECTED_CLASS_COUNTS,
        "layers": layers,
    }
    output = args.output or HERE / "boundary_layer_classification.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("wrote", output)


if __name__ == "__main__":
    main()
