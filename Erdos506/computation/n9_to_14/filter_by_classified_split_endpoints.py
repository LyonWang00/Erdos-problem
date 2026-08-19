#!/usr/bin/env python3
"""Balance exact split pair profiles forced at zero-defect points."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from math import comb
import json
from pathlib import Path

from low_defect_point_support_filter import (
    augmented_line_arrangement_row_ok,
    augmented_line_melchior_row_ok,
    classified_zero_split_profile_histograms,
    conditioned_rows,
    melchior_defect,
    retained_records,
)
from ortools.sat.python import cp_model
from pair_profile_generator import profiles as pair_profiles
from projected_arrangement_model import projected_local_types


HERE = Path(__file__).resolve().parent


def low_defect_combined_histograms(n: int):
    if n not in (12, 13):
        return {}
    path = HERE / (
        "n12_low_defect_line_histograms_wiring_retry300.json"
        if n == 12 else "n13_low_defect_line_histograms_pbd.json"
    )
    source = json.loads(path.read_text(encoding="utf-8"))
    if (source.get("status") != "PASS"
            or source.get("materialised_parent_count")
            != source.get("parent_count")):
        raise ValueError(f"incomplete low-defect catalogue: {path}")
    answer = {}
    for record in source["records"]:
        profiles = tuple(tuple(profile) for profile in record["profiles"])
        answer[tuple(record["parent_type"])] = tuple(
            Counter({profiles[index]: multiplicity
                     for index, multiplicity in enumerate(histogram)
                     if multiplicity})
            for histogram in record["histograms"]
        )
    return answer


def solve(task):
    (n, signature, seconds, zero_branch_pairing,
     use_zero_catalogue, use_combined_catalogue,
     use_augmented_arrangement_types) = task
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    zero_catalogue = (
        classified_zero_split_profile_histograms(n - 1, dimension)
        if use_zero_catalogue else {}
    )
    combined_catalogue = (
        low_defect_combined_histograms(n) if use_combined_catalogue else {}
    )
    rows = []
    for row in conditioned_rows(n, sizes, signature):
        if not augmented_line_melchior_row_ok(n, row):
            continue
        if (use_augmented_arrangement_types
                and not augmented_line_arrangement_row_ok(n, row)):
            continue
        parent = tuple(row[index] + row[dimension + index]
                       for index in range(dimension))
        if (use_zero_catalogue and melchior_defect(parent) == 0
                and row not in zero_catalogue):
            continue
        rows.append(row)
    if not rows:
        return {"signature": list(signature), "status": "INFEASIBLE",
                "wall_time_seconds": 0.0, "variables": 0,
                "constraints": 0}

    profiles = pair_profiles(n, sizes[-1])
    profile_index_of = {
        tuple(profile): index for index, profile in enumerate(profiles)
    }
    children = frozenset(projected_local_types(n - 1, sizes[-1]))
    model = cp_model.CpModel()
    row_variables = tuple(
        model.NewIntVar(0, n, f"row_{index}")
        for index in range(len(rows))
    )
    model.Add(sum(row_variables) == n)
    for coordinate, (size, number) in enumerate(zip(sizes + sizes, signature)):
        model.Add(sum(
            row_variables[index] * row[coordinate]
            for index, row in enumerate(rows)
        ) == size * number)

    total_defect = (
        3 * (signature[0] + signature[dimension]) - 3 * n
        - sum((size - 4) * size
              * (signature[index] + signature[dimension + index])
              for index, size in enumerate(sizes) if size >= 5)
    )
    model.Add(sum(
        row_variables[index] * melchior_defect(tuple(
            row[j] + row[dimension + j] for j in range(dimension)
        )) for index, row in enumerate(rows)
    ) == total_defect)

    profile_variables = tuple(
        model.NewIntVar(0, comb(n, 2), f"pair_{index}")
        for index in range(len(profiles))
    )
    model.Add(sum(profile_variables) == comb(n, 2))
    for coordinate, (size, number) in enumerate(zip(sizes + sizes, signature)):
        model.Add(sum(
            profile_variables[index] * profile[coordinate]
            for index, profile in enumerate(profiles)
        ) == comb(size, 2) * number)

    endpoints = {}
    by_profile = [[] for _ in profiles]
    zero_pair_branches = []
    for row_index, row in enumerate(rows):
        combined_row = tuple(
            row[index] + row[dimension + index]
            for index in range(dimension)
        )
        allowed_indices = []
        for profile_index, profile in enumerate(profiles):
            if any(profile[index] > row[index]
                   for index in range(2 * dimension)):
                continue
            combined_profile = tuple(
                profile[index] + profile[dimension + index]
                for index in range(dimension)
            )
            child = tuple(
                combined_row[index] - combined_profile[index]
                + (combined_profile[index + 1]
                   if index + 1 < dimension else 0)
                for index in range(dimension)
            )
            if child not in children:
                continue
            variable = model.NewIntVar(
                0, n * (n - 1), f"endpoint_{row_index}_{profile_index}"
            )
            endpoints[row_index, profile_index] = variable
            by_profile[profile_index].append(variable)
            allowed_indices.append(profile_index)
        model.Add(sum(endpoints[row_index, index]
                      for index in allowed_indices)
                  == (n - 1) * row_variables[row_index])
        for coordinate, size in enumerate(sizes + sizes):
            model.Add(sum(
                endpoints[row_index, index] * profiles[index][coordinate]
                for index in allowed_indices
            ) == (size - 1) * row[coordinate]
                 * row_variables[row_index])

        if row in zero_catalogue:
            histograms = tuple(Counter(histogram)
                               for histogram in zero_catalogue[row])
            choices = tuple(
                model.NewIntVar(0, n, f"choice_{row_index}_{index}")
                for index in range(len(histograms))
            )
            model.Add(sum(choices) == row_variables[row_index])
            zero_pair_branches.extend(
                (choice, histogram)
                for choice, histogram in zip(choices, histograms)
            )
            for profile_index, profile in enumerate(profiles):
                target = sum(
                    choices[index] * histogram.get(tuple(profile), 0)
                    for index, histogram in enumerate(histograms)
                )
                variable = endpoints.get((row_index, profile_index))
                if variable is None:
                    model.Add(target == 0)
                else:
                    model.Add(variable == target)

        if combined_row in combined_catalogue:
            histograms = combined_catalogue[combined_row]
            choices = tuple(
                model.NewIntVar(0, n, f"combined_choice_{row_index}_{index}")
                for index in range(len(histograms))
            )
            model.Add(sum(choices) == row_variables[row_index])
            combined_profile_types = sorted({
                tuple(profile[index] + profile[dimension + index]
                      for index in range(dimension))
                for profile in profiles
            })
            for combined_profile in combined_profile_types:
                left = sum(
                    endpoints[row_index, profile_index]
                    for profile_index, profile in enumerate(profiles)
                    if (row_index, profile_index) in endpoints
                    and tuple(profile[index] + profile[dimension + index]
                              for index in range(dimension))
                    == combined_profile
                )
                right = sum(
                    choices[index] * histogram.get(combined_profile, 0)
                    for index, histogram in enumerate(histograms)
                )
                model.Add(left == right)

            # At the centre, rich original-line blocks are pairwise disjoint
            # outside the centre.  Thus each endpoint of a combined profile
            # is either on no rich original line through the centre or is
            # marked by one unique block size.  For a size-k line, exactly
            # k-1 endpoints receive that mark.  This converts every retained
            # combined-profile histogram into a necessary split-profile
            # histogram without enumerating labelled supports.
            if row not in zero_catalogue:
                generated = {index: [] for index in range(len(profiles))}
                for branch_index, histogram in enumerate(histograms):
                    by_coordinate = [[] for _ in range(dimension)]
                    for combined_profile, count in histogram.items():
                        options = [-1] + [
                            coordinate for coordinate in range(dimension)
                            if combined_profile[coordinate] > 0
                        ]
                        marks = []
                        for option in options:
                            mark = model.NewIntVar(
                                0, (n - 1) * n,
                                "combined_mark_"
                                f"{row_index}_{branch_index}_"
                                f"{'none' if option < 0 else option}_"
                                + "_".join(map(str, combined_profile)),
                            )
                            marks.append(mark)
                            if option < 0:
                                split_profile = ((0,) * dimension
                                                 + combined_profile)
                            else:
                                line_part = tuple(
                                    int(index == option)
                                    for index in range(dimension)
                                )
                                circle_part = tuple(
                                    combined_profile[index]
                                    - int(index == option)
                                    for index in range(dimension)
                                )
                                split_profile = line_part + circle_part
                                by_coordinate[option].append(mark)
                            split_index = profile_index_of.get(split_profile)
                            if split_index is None:
                                model.Add(mark == 0)
                            else:
                                generated[split_index].append(mark)
                        model.Add(sum(marks) == count * choices[branch_index])
                    for coordinate, size in enumerate(sizes):
                        model.Add(
                            sum(by_coordinate[coordinate])
                            == (size - 1) * row[coordinate]
                            * choices[branch_index]
                        )
                for profile_index in range(len(profiles)):
                    variable = endpoints.get((row_index, profile_index))
                    terms = generated[profile_index]
                    if variable is None:
                        model.Add(sum(terms) == 0)
                    else:
                        model.Add(variable == sum(terms))

    for profile_index, variable in enumerate(profile_variables):
        model.Add(sum(by_profile[profile_index]) == 2 * variable)

    if zero_branch_pairing:
        zero_edges_by_profile = [[] for _ in profiles]
        zero_incident = [
            {profile: [] for profile in histogram}
            for _choice, histogram in zero_pair_branches
        ]
        for first in range(len(zero_pair_branches)):
            first_choice, first_histogram = zero_pair_branches[first]
            for second in range(first, len(zero_pair_branches)):
                second_choice, second_histogram = zero_pair_branches[second]
                common = tuple(sorted(
                    set(first_histogram) & set(second_histogram)
                ))
                edge_variables = []
                for profile in common:
                    variable = model.NewIntVar(
                        0, comb(n, 2),
                        f"zero_edge_{first}_{second}_"
                        + "_".join(map(str, profile)),
                    )
                    edge_variables.append(variable)
                    zero_edges_by_profile[profile_index_of[profile]].append(
                        variable
                    )
                    zero_incident[first][profile].append(variable)
                    zero_incident[second][profile].append(variable)
                total = model.NewIntVar(
                    0, comb(n, 2), f"zero_pair_total_{first}_{second}"
                )
                model.Add(total == sum(edge_variables))
                if first == second:
                    model.AddAllowedAssignments(
                        (first_choice, total),
                        tuple((value, comb(value, 2))
                              for value in range(n + 1)),
                    )
                else:
                    model.AddAllowedAssignments(
                        (first_choice, second_choice, total),
                        tuple((left, right, left * right)
                              for left in range(n + 1)
                              for right in range(n + 1)),
                    )
        for branch_index, (choice, histogram) in enumerate(
                zero_pair_branches):
            for profile, count in histogram.items():
                # A same-branch edge was inserted twice into this incidence
                # list, exactly as required for its two endpoints.
                model.Add(
                    sum(zero_incident[branch_index][profile])
                    <= count * choice
                )
        for profile_index, profile in enumerate(profiles):
            zero_supply = sum(
                histogram.get(tuple(profile), 0) * choice
                for choice, histogram in zero_pair_branches
            )
            zero_to_positive = model.NewIntVar(
                0, n * (n - 1), f"zero_positive_{profile_index}"
            )
            model.Add(
                zero_to_positive
                == zero_supply - 2 * sum(zero_edges_by_profile[profile_index])
            )
            positive_supply = sum(by_profile[profile_index]) - zero_supply
            model.Add(positive_supply >= zero_to_positive)
            positive_positive_endpoints = model.NewIntVar(
                0, n * (n - 1), f"positive_positive_{profile_index}"
            )
            model.Add(
                positive_positive_endpoints
                == positive_supply - zero_to_positive
            )
            parity = model.NewIntVar(0, 1, f"positive_parity_{profile_index}")
            model.AddModuloEquality(parity, positive_positive_endpoints, 2)
            model.Add(parity == 0)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    status = solver.Solve(model)
    return {
        "signature": list(signature),
        "status": solver.StatusName(status),
        "wall_time_seconds": solver.WallTime(),
        "row_count": len(rows),
        "profile_count": len(profiles),
        "zero_row_count": sum(row in zero_catalogue for row in rows),
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
        "zero_branch_pairing": zero_branch_pairing,
        "zero_catalogue": use_zero_catalogue,
        "combined_catalogue": use_combined_catalogue,
        "augmented_arrangement_types": use_augmented_arrangement_types,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--only-source-status", action="append",
        help="retain only source records with one of these status strings",
    )
    parser.add_argument("--zero-branch-pairing", action="store_true")
    parser.add_argument("--no-zero-catalogue", action="store_true")
    parser.add_argument("--no-combined-catalogue", action="store_true")
    parser.add_argument("--augmented-arrangement-types", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = retained_records(args.n, args.input_json)
    if args.only_source_status:
        wanted = set(args.only_source_status)
        records = tuple(
            record for record in records
            if str(record.get("status", "SURVIVES")) in wanted
        )
    if args.limit is not None:
        records = records[:args.limit]
    tasks = [(args.n, tuple(record["signature"]), args.seconds,
              args.zero_branch_pairing, not args.no_zero_catalogue,
              not args.no_combined_catalogue,
              args.augmented_arrangement_types)
             for record in records]
    if args.jobs == 1:
        outcomes = list(map(solve, tasks))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            outcomes = list(executor.map(solve, tasks, chunksize=1))
    retained = [record for record in outcomes
                if record["status"] != "INFEASIBLE"]
    payload = {
        "schema_version": 1,
        "status": "PASS",
        "claim": (
            "Exact line/circle-split endpoint balance with child-arrangement "
            "constraints.  Optional zero-defect and low-defect catalogues are "
            "used only when their corresponding switches are enabled; "
            "disabling a catalogue retains rather than deletes its rows.  "
            "UNKNOWN outcomes are retained."
        ),
        "n": args.n,
        "source_file": args.input_json.name,
        "source_sha256": sha256(args.input_json.read_bytes()).hexdigest(),
        "input_signature_count": len(outcomes),
        "status_counts": dict(Counter(record["status"] for record in outcomes)),
        "retained_signature_count": len(retained),
        "zero_branch_pairing": args.zero_branch_pairing,
        "zero_catalogue": not args.no_zero_catalogue,
        "combined_catalogue": not args.no_combined_catalogue,
        "augmented_arrangement_types": args.augmented_arrangement_types,
        "records": outcomes,
    }
    args.output.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "PASS",
        "input_signature_count": len(outcomes),
        "status_counts": payload["status_counts"],
        "retained_signature_count": len(retained),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
