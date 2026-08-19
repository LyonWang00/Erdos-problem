#!/usr/bin/env python3
"""Exact rich-block support model anchored at an arbitrary low-defect point.

The total Melchior defect guarantees a point whose local arrangement has
small defect.  Relabel it as 0.  For each complete retained line-profile
histogram, relabel the other points so equal profiles occur consecutively.
The model then selects every maximal line and circle support and enforces the
profile of each pair {0,p} exactly.  This covers all incidence realizations of
the histogram; no representative PBD is fixed.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
COMPUTATIONS = HERE.parent
AUTHORITATIVE = HERE
sys.path.insert(0, str(COMPUTATIONS))
from runtime_dependencies import activate_ortools  # noqa: E402

activate_ortools()
sys.path.insert(0, str(AUTHORITATIVE))
from ortools.sat.python import cp_model  # noqa: E402
from signature_filters import KNOWN  # noqa: E402
from verify_pair_moment_filter import conditioned_rows  # noqa: E402
from projected_arrangement_model import projected_local_types  # noqa: E402


INPUTS = {
    12: "n12_sparse_split_transport_from_correct_156.json",
    13: "n13_block_row_strengthened_corrected_complete.json",
    14: "n14_summed_local_inequalities.json",
}
RETAINED = frozenset({"OPTIMAL", "FEASIBLE", "UNKNOWN", "SURVIVES"})
SIMPLICIAL_NORMALS = HERE / "simplicial_arrangement_normals_q11_q13.json"
PSEUDOLINE_TRANSVERSALS = (
    HERE / "classified_pseudoline_transversal_families.json"
)

# Exact generalized-allowable-sequence exclusions for all retained numerical
# types of eleven projective lines.  A normalized copy and independent rerun
# certificate are generated alongside the final verifier.
INFEASIBLE_LOCAL_TYPES_11 = frozenset({
    (6, 11, 1, 1),
    (7, 10, 3, 0), (8, 9, 0, 2), (9, 6, 3, 1), (10, 3, 6, 0),
    (7, 12, 2, 0),
    (7, 14, 1, 0), (11, 6, 1, 2), (12, 3, 4, 1),
})


def melchior_defect(parent):
    return (
        parent[0] - 3
        - sum((r - 3) * parent[r - 2]
              for r in range(4, len(parent) + 2))
    )


def augmented_line_melchior_row_ok(n: int, row: tuple[int, ...]) -> bool:
    """Check Melchior after adjoining the dual of the inversion centre.

    If a rich original line of size k passes through the centre, the added
    dual line passes through an existing vertex of multiplicity k-1.  Adding
    that line changes the local Melchior defect by

        (n-1) - sum(k * d_k^line).

    The enlarged essential real projective arrangement has nonnegative
    defect, which gives the tested inequality.
    """
    dimension = len(row) // 2
    parent = tuple(
        row[index] + row[dimension + index]
        for index in range(dimension)
    )
    return sum(
        size * row[index]
        for index, size in enumerate(range(3, 3 + dimension))
    ) <= n - 1 + melchior_defect(parent)


def augmented_line_arrangement_vector(
    n: int, row: tuple[int, ...],
) -> tuple[int, ...]:
    """Intersection vector after adjoining the dual centre line."""
    dimension = len(row) // 2
    parent = tuple(
        row[index] + row[dimension + index]
        for index in range(dimension)
    )
    line_degrees = row[:dimension]
    sizes = tuple(range(3, 3 + dimension))
    generic_double_points = n - 1 - sum(
        (size - 1) * line_degrees[index]
        for index, size in enumerate(sizes)
    )
    answer = [parent[0] - line_degrees[0] + generic_double_points]
    for index in range(1, dimension + 1):
        answer.append(
            (parent[index] if index < dimension else 0)
            - (line_degrees[index] if index < dimension else 0)
            + line_degrees[index - 1]
        )
    return tuple(answer)


@lru_cache(maxsize=None)
def _augmented_line_arrangement_types(n: int, maximum_global_size: int):
    # projected_local_types(m,K) describes arrangements of m-1 lines with
    # multiplicities through K-1.  The augmented arrangement has n lines and
    # may increase the former maximum multiplicity by one.
    return frozenset(projected_local_types(n + 1, maximum_global_size + 1))


def augmented_line_arrangement_row_ok(n: int, row: tuple[int, ...]) -> bool:
    dimension = len(row) // 2
    maximum_global_size = dimension + 2
    return augmented_line_arrangement_vector(n, row) in (
        _augmented_line_arrangement_types(n, maximum_global_size)
    )


def retained_records(n: int, path: Path):
    source = json.loads(path.read_text(encoding="utf-8"))
    source_records = source.get("records")
    if source_records is None and isinstance(source.get("layers"), list):
        source_records = [
            record for layer in source["layers"]
            for record in layer.get("survivors", [])
        ]
    if not isinstance(source_records, list):
        raise ValueError(f"no record list in {path}")
    answer = []
    for record in source_records:
        if record.get("excluded") is True:
            continue
        status = str(record.get("status", "SURVIVES"))
        if status.startswith("INFEASIBLE"):
            continue
        answer.append(record)
    signatures = [tuple(record["signature"]) for record in answer]
    if len(signatures) != len(set(signatures)):
        raise ValueError("duplicate retained signatures")
    return tuple(answer)


def load_anchor_profiles(n: int, maximum_block_size: int):
    path = HERE / (
        "n12_low_defect_line_histograms_wiring_retry300.json"
        if n == 12 else f"n{n}_low_defect_line_histograms_pbd.json"
    )
    source = json.loads(path.read_text(encoding="utf-8"))
    if (source.get("status") != "PASS"
            or source.get("n") != n
            or source.get("maximum_block_size") != maximum_block_size
            or source.get("materialised_parent_count")
            != source.get("parent_count")):
        raise ValueError(f"incomplete anchor catalogue: {path}")
    branches = []
    for record in source["records"]:
        profiles = tuple(tuple(profile) for profile in record["profiles"])
        for histogram in record["histograms"]:
            expanded = tuple(sorted(
                profile
                for profile, multiplicity in zip(profiles, histogram)
                for _ in range(multiplicity)
            ))
            if len(expanded) != n - 1:
                raise AssertionError((n, len(expanded)))
            branches.append((tuple(record["parent_type"]), expanded))
    zero_profiles = classified_zero_profile_histograms(
        n - 1, maximum_block_size - 2
    )
    branches = [
        branch for branch in branches
        if (melchior_defect(branch[0]) != 0
            or branch[1] in zero_profiles.get(branch[0], ()))
    ]
    if not branches:
        raise ValueError("anchor catalogue has no retained branch")
    return path, tuple(branches)


def load_anchor_parent_types(n: int, maximum_block_size: int):
    """Load only local multiplicity vectors, without profile histograms."""
    candidates = [
        HERE / f"n{n}_low_defect_parent_types.json",
        HERE / f"n{n}_low_defect_line_histograms_corrected.json",
        HERE / f"n{n}_low_defect_line_histograms.json",
        HERE / f"n{n}_low_defect_line_histograms_pbd.json",
    ]
    path = next((candidate for candidate in candidates
                 if candidate.is_file()), None)
    if path is None:
        raise FileNotFoundError(f"no low-defect parent catalogue for n={n}")
    source = json.loads(path.read_text(encoding="utf-8"))
    if (source.get("status") != "PASS"
            or source.get("n") != n
            or source.get("maximum_block_size") != maximum_block_size
            or source.get("materialised_parent_count")
            != source.get("parent_count")):
        raise ValueError(f"incomplete parent catalogue: {path}")
    parents = tuple(sorted({
        tuple(record["parent_type"]) for record in source["records"]
        if n != 12
        or tuple(record["parent_type"]) not in INFEASIBLE_LOCAL_TYPES_11
    }))
    if not parents:
        raise ValueError("parent catalogue is empty")
    return path, tuple((parent, None) for parent in parents)


def quadratic(value):
    """Return a+b*tau as (a,b), where tau^2=tau+1."""
    return (int(value), 0) if isinstance(value, int) else tuple(value)


def qadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def qneg(value):
    return -value[0], -value[1]


def qmul(left, right):
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def determinant_zero(first, second, third):
    u = tuple(quadratic(value) for value in first)
    v = tuple(quadratic(value) for value in second)
    w = tuple(quadratic(value) for value in third)

    def minor(a, b, c, d):
        return qadd(qmul(a, d), qneg(qmul(b, c)))

    determinant = qadd(
        qadd(qmul(u[0], minor(v[1], v[2], w[1], w[2])),
             qneg(qmul(u[1], minor(v[0], v[2], w[0], w[2])))),
        qmul(u[2], minor(v[0], v[1], w[0], w[1])),
    )
    return determinant == (0, 0)


def simplicial_incidence_catalogue(number_of_lines: int, dimension: int):
    """Reconstruct and audit every relevant zero-defect incidence table."""
    source = json.loads(SIMPLICIAL_NORMALS.read_text(encoding="utf-8"))
    if source.get("status") != "PASS":
        raise ValueError("simplicial normal-vector catalogue is not complete")
    answer = []
    for entry in source["arrangements"]:
        normals = entry["normals"]
        if len(normals) != number_of_lines:
            continue
        blocks = set()
        for first, second in combinations(range(number_of_lines), 2):
            block = tuple(
                point for point in range(number_of_lines)
                if determinant_zero(
                    normals[first], normals[second], normals[point]
                )
            )
            if len(block) < 2:
                raise AssertionError("a pair did not determine a vertex")
            blocks.add(block)
        blocks = tuple(sorted(blocks, key=lambda block: (len(block), block)))
        covered_pairs = Counter(
            pair for block in blocks for pair in combinations(block, 2)
        )
        if (len(covered_pairs) != number_of_lines * (number_of_lines - 1) // 2
                or set(covered_pairs.values()) != {1}):
            raise ValueError(f"invalid incidence table for {entry['name']}")
        full_vector = tuple(
            sum(len(block) == multiplicity for block in blocks)
            for multiplicity in range(2, number_of_lines + 1)
        )
        recorded = tuple(entry["t_vector"])
        recorded += (0,) * (len(full_vector) - len(recorded))
        if full_vector != recorded:
            raise ValueError(f"incorrect t-vector for {entry['name']}")
        if full_vector[0] - 3 - sum(
                (multiplicity - 3) * full_vector[multiplicity - 2]
                for multiplicity in range(4, number_of_lines + 1)) != 0:
            raise ValueError(f"nonzero defect for {entry['name']}")
        if any(full_vector[dimension:]):
            continue
        answer.append({
            "name": entry["name"],
            "parent": full_vector[:dimension],
            "blocks": blocks,
            "normals": tuple(
                tuple(quadratic(value) for value in normal)
                for normal in normals
            ),
        })
    return tuple(answer)


def classified_zero_profile_histograms(number_of_lines: int, dimension: int):
    """Return the point-profile histograms from the complete classification."""
    answer = {}
    for entry in simplicial_incidence_catalogue(number_of_lines, dimension):
        histogram = tuple(sorted(
            tuple(
                sum(point in block and len(block) == multiplicity
                    for block in entry["blocks"])
                for multiplicity in range(2, dimension + 2)
            )
            for point in range(number_of_lines)
        ))
        answer.setdefault(entry["parent"], set()).add(histogram)
    return {parent: frozenset(histograms)
            for parent, histograms in answer.items()}


def classified_zero_pair_statistics(number_of_lines: int, dimension: int):
    """Count endpoint-profile pairs by their common vertex multiplicity."""
    answer = {}
    for entry in simplicial_incidence_catalogue(number_of_lines, dimension):
        blocks = entry["blocks"]
        profiles = {
            point: tuple(
                sum(point in block and len(block) == multiplicity
                    for block in blocks)
                for multiplicity in range(2, dimension + 2)
            )
            for point in range(number_of_lines)
        }
        statistics = Counter()
        for first, second in combinations(range(number_of_lines), 2):
            multiplicities = [
                len(block) for block in blocks
                if first in block and second in block
            ]
            if len(multiplicities) != 1:
                raise ValueError("classified PBD does not own each pair once")
            endpoints = tuple(sorted((profiles[first], profiles[second])))
            statistics[endpoints + (multiplicities[0],)] += 1
        answer[entry["parent"]] = statistics
    return answer


def classified_zero_line_families(number_of_lines: int, dimension: int):
    """Return a necessary catalogue for blocks that remain original lines.

    Invert about the chosen point and dualise the resulting point
    configuration.  Blocks that were lines before inversion become vertices
    on one additional projective line.  A lift of that line crosses every old
    arrangement line exactly once.  The companion verifier enumerates every
    such monotone path in the tope graph of each classified arrangement.

    The catalogue deliberately allows pseudoline paths, not only straight
    lines in the displayed coordinate realization.  It is therefore a
    realization-independent necessary condition and does not assume
    projective rigidity of the classified arrangements.
    """
    source = json.loads(PSEUDOLINE_TRANSVERSALS.read_text(encoding="utf-8"))
    if source.get("status") != "PASS":
        raise ValueError("pseudoline-transversal catalogue is not complete")
    by_name = {entry["name"]: entry for entry in source["arrangements"]}
    answer = {}
    for arrangement in simplicial_incidence_catalogue(
            number_of_lines, dimension):
        entry = by_name.get(arrangement["name"])
        if entry is None:
            raise ValueError(
                f"missing transversal catalogue for {arrangement['name']}"
            )
        if (entry["number_of_lines"] != number_of_lines
                or tuple(entry["parent"]) != arrangement["parent"]):
            raise ValueError("transversal catalogue metadata mismatch")
        families = tuple(
            frozenset(int(index) for index in family)
            for family in entry["families"]
        )
        if len(families) != entry["transversal_family_count"]:
            raise ValueError("transversal-family count mismatch")
        disjoint = set(pairwise_disjoint_subfamilies(arrangement["blocks"]))
        if not set(families).issubset(disjoint):
            raise ValueError("a transversal family is not pairwise disjoint")
        answer[arrangement["parent"]] = families
    return answer


def classified_zero_line_degree_vectors(number_of_lines: int, dimension: int):
    """Enumerate size counts of classified pseudoline-transversal families."""
    answer = {}
    for entry in simplicial_incidence_catalogue(number_of_lines, dimension):
        blocks = entry["blocks"]
        counts = {
            tuple(
                sum(index in family and len(block) == multiplicity
                    for index, block in enumerate(blocks))
                for multiplicity in range(2, dimension + 2)
            )
            for family in classified_zero_line_families(
                number_of_lines, dimension
            )[entry["parent"]]
        }
        answer[entry["parent"]] = tuple(sorted(counts))
    return answer


def classified_zero_split_profile_histograms(
        number_of_lines: int, dimension: int):
    """Classify line/circle-split endpoint profiles at a zero-defect point.

    After inversion at the point, vertices assigned to original lines form a
    pairwise-disjoint subfamily of the classified arrangement vertices.  Each
    such colouring fixes both the line/circle degree row at the point and the
    complete multiset of split profiles of its pairs with the other points.
    """
    answer = {}
    for entry in simplicial_incidence_catalogue(number_of_lines, dimension):
        blocks = entry["blocks"]
        parent = entry["parent"]
        combined_profiles = tuple(
            tuple(
                sum(point in block and len(block) == multiplicity
                    for block in blocks)
                for multiplicity in range(2, dimension + 2)
            )
            for point in range(number_of_lines)
        )
        for subfamily in classified_zero_line_families(
                number_of_lines, dimension)[parent]:
            line_degrees = tuple(
                sum(index in subfamily and len(block) == multiplicity
                    for index, block in enumerate(blocks))
                for multiplicity in range(2, dimension + 2)
            )
            row = line_degrees + tuple(
                parent[index] - line_degrees[index]
                for index in range(dimension)
            )
            histogram = tuple(sorted(
                tuple(
                    sum(
                        block_index in subfamily
                        and point in blocks[block_index]
                        and len(blocks[block_index]) == multiplicity
                        for block_index in range(len(blocks))
                    )
                    for multiplicity in range(2, dimension + 2)
                )
                + tuple(
                    combined_profiles[point][index]
                    - sum(
                        block_index in subfamily
                        and point in blocks[block_index]
                        and len(blocks[block_index]) == index + 2
                        for block_index in range(len(blocks))
                    )
                    for index in range(dimension)
                )
                for point in range(number_of_lines)
            ))
            answer.setdefault(row, set()).add(histogram)
    return {row: frozenset(histograms) for row, histograms in answer.items()}


def pairwise_disjoint_subfamilies(blocks):
    """Enumerate all subfamilies whose blocks are pairwise disjoint."""
    answer = []

    def search(index, used_points, selected):
        if index == len(blocks):
            answer.append(frozenset(selected))
            return
        search(index + 1, used_points, selected)
        block_points = set(blocks[index])
        if used_points.isdisjoint(block_points):
            search(index + 1, used_points | block_points,
                   selected + [index])

    search(0, set(), [])
    return tuple(answer)


def solve_signature(
    n: int,
    signature: tuple[int, ...],
    seconds: float,
    anchor_branch_index: int | None = None,
    maximum_deletion_size: int = 2,
    solver_workers: int = 1,
    random_seed: int = 0,
    parent_only: bool = False,
    anchor_row_index: int | None = None,
    zero_pair_statistics_enabled: bool = True,
    anchor_line_family_index: int | None = None,
    zero_split_histogram_enabled: bool = True,
):
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    line_numbers = signature[:dimension]
    circle_numbers = signature[dimension:]
    circle_count = sum(circle_numbers)
    rows = conditioned_rows(n, sizes, signature)
    rows = tuple(
        row for row in rows if augmented_line_melchior_row_ok(n, row)
    )
    if n == 12:
        rows = tuple(
            row for row in rows
            if tuple(
                row[j] + row[dimension + j]
                for j in range(dimension)
            ) not in INFEASIBLE_LOCAL_TYPES_11
        )
    if not rows:
        return {
            "signature": list(signature), "circle_count": circle_count,
            "status": "INFEASIBLE_NO_POINT_ROW", "wall_time_seconds": 0.0,
        }
    loader = load_anchor_parent_types if parent_only else load_anchor_profiles
    _catalogue_path, catalogue_branches = loader(n, sizes[-1])
    anchor_branches = catalogue_branches
    catalogue_source = json.loads(_catalogue_path.read_text(encoding="utf-8"))
    maximum_defect = int(catalogue_source["maximum_defect"])
    total_defect = (
        3 * (line_numbers[0] + circle_numbers[0]) - 3 * n
        - sum(
            (size - 4) * size * (line_number + circle_number)
            for size, line_number, circle_number in zip(
                sizes[2:], line_numbers[2:], circle_numbers[2:]
            )
        )
    )
    if not 0 <= total_defect < n * (maximum_defect + 1):
        raise ValueError(
            "catalogue defect threshold does not guarantee an anchor: "
            f"total={total_defect}, threshold={maximum_defect}"
        )
    anchor_defect_bound = total_defect // n

    if anchor_branch_index is not None:
        if not 0 <= anchor_branch_index < len(anchor_branches):
            raise ValueError("anchor branch index out of range")
        selected_anchor = anchor_branches[anchor_branch_index]
        if melchior_defect(selected_anchor[0]) > anchor_defect_bound:
            return {
                "signature": list(signature),
                "circle_count": circle_count,
                "anchor_branch_count": 1,
                "anchor_branch_index": anchor_branch_index,
                "anchor_defect_bound": anchor_defect_bound,
                "status": "INFEASIBLE_BY_DEFECT_AVERAGE",
                "wall_time_seconds": 0.0,
                "variables": 0,
                "constraints": 0,
            }
        anchor_branches = (selected_anchor,)
    else:
        anchor_branches = tuple(
            branch for branch in anchor_branches
            if melchior_defect(branch[0]) <= anchor_defect_bound
        )
        if not anchor_branches:
            raise AssertionError("defect average left no anchor branch")

    anchor_row_options = None
    if anchor_row_index is not None:
        if anchor_branch_index is None:
            raise ValueError("an anchor row requires a fixed anchor branch")
        parent = anchor_branches[0][0]
        anchor_row_options = tuple(
            row for row in rows
            if tuple(
                row[coordinate] + row[dimension + coordinate]
                for coordinate in range(dimension)
            ) == parent
        )
        if not 0 <= anchor_row_index < len(anchor_row_options):
            raise ValueError("anchor row index out of range")
    local_catalogue_branches = catalogue_branches
    if parent_only and anchor_branch_index is not None:
        anchor_parent = anchor_branches[0][0]
        anchor_defect = melchior_defect(anchor_parent)
        local_catalogue_branches = tuple(
            branch for branch in catalogue_branches
            if not (
                melchior_defect(branch[0]) == anchor_defect
                and branch[0] < anchor_parent
            )
        )
    elif anchor_branch_index is not None:
        anchor_defect = melchior_defect(anchor_branches[0][0])
        local_catalogue_branches = tuple(
            branch for index, branch in enumerate(catalogue_branches)
            if (melchior_defect(branch[0]) > anchor_defect
                or (melchior_defect(branch[0]) == anchor_defect
                    and index >= anchor_branch_index))
        )

    points = tuple(range(n))
    pairs = tuple(combinations(points, 2))
    triples = tuple(combinations(points, 3))
    subsets = {
        size: tuple(combinations(points, size))
        for size, line_number, circle_number in zip(
            sizes, line_numbers, circle_numbers
        ) if line_number or circle_number
    }
    model = cp_model.CpModel()
    line = {}
    circle = {}
    for coordinate, size in enumerate(sizes):
        candidates = subsets.get(size, ())
        if line_numbers[coordinate]:
            for block in candidates:
                line[size, block] = model.NewBoolVar(
                    f"L{size}_" + "_".join(map(str, block))
                )
            model.Add(sum(line[size, block] for block in candidates)
                      == line_numbers[coordinate])
        if circle_numbers[coordinate]:
            for block in candidates:
                circle[size, block] = model.NewBoolVar(
                    f"C{size}_" + "_".join(map(str, block))
                )
            model.Add(sum(circle[size, block] for block in candidates)
                      == circle_numbers[coordinate])

    # Every triple is contained in its unique maximal line or circle.
    for triple in triples:
        triple_set = set(triple)
        owners = []
        for family in (line, circle):
            owners.extend(
                selected for (size, block), selected in family.items()
                if triple_set.issubset(block)
            )
        model.Add(sum(owners) == 1)

    # Two distinct maximal connecting lines cannot contain the same pair.
    for pair in pairs:
        pair_set = set(pair)
        model.Add(sum(
            selected for (_size, block), selected in line.items()
            if pair_set.issubset(block)
        ) <= 1)

    degree = {}
    for point in points:
        point_row = []
        for family_index, family in enumerate((line, circle)):
            for coordinate, size in enumerate(sizes):
                variable = model.NewIntVar(
                    0, n, f"degree_{point}_{family_index}_{size}"
                )
                degree[point, family_index, coordinate] = variable
                model.Add(variable == sum(
                    selected for (block_size, block), selected
                    in family.items()
                    if block_size == size and point in block
                ))
                point_row.append(variable)
        model.AddAllowedAssignments(
            point_row,
            (anchor_row_options[anchor_row_index],)
            if point == 0 and anchor_row_options is not None else rows,
        )
        # Distinct maximal lines through one point are disjoint away from
        # that point.  This follows from the pair constraints below, but the
        # summed form is a much stronger propagation constraint.
        model.Add(sum(
            (size - 1) * degree[point, 0, coordinate]
            for coordinate, size in enumerate(sizes)
        ) <= n - 1)
        model.Add(
            circle_count - degree[point, 1, 0] >= KNOWN[n - 1]
        )

    # Exact two-point deletion inequality.
    c3_index = 0
    c4_index = 1
    loss_bound = circle_count - KNOWN[n - 2]
    for first, second in pairs:
        c3_pair = sum(
            selected for (size, block), selected in circle.items()
            if size == 3 and first in block and second in block
        )
        c4_pair = sum(
            selected for (size, block), selected in circle.items()
            if size == 4 and first in block and second in block
        )
        model.Add(
            degree[first, 1, c3_index]
            + degree[second, 1, c3_index]
            - c3_pair + c4_pair <= loss_bound
        )

    # For every larger deleted set D, all original circles retaining at
    # least three points remain distinct circles of P\D.  Hence the number
    # destroyed is at most c(P)-c(n-|D|).  This is the exact support form of
    # the one- and two-point deletion inequalities above.
    maximum_deletion_size = min(maximum_deletion_size, n - 4)
    for deletion_size in range(3, maximum_deletion_size + 1):
        survivor_bound = KNOWN[n - deletion_size]
        for deleted in combinations(points, deletion_size):
            deleted_set = set(deleted)
            destroyed = sum(
                selected
                for (size, block), selected in circle.items()
                if len(deleted_set.intersection(block)) > size - 3
            )
            model.Add(circle_count - destroyed >= survivor_bound)

    # At least one point has low defect.  Choose it as 0 and relabel the
    # remaining points according to one complete profile histogram.
    branch_variables = []
    for branch_index, (_parent, profiles) in enumerate(anchor_branches):
        selected_branch = (
            None if anchor_branch_index is not None
            else model.NewBoolVar(f"anchor_branch_{branch_index}")
        )
        if selected_branch is not None:
            branch_variables.append(selected_branch)
        for coordinate in range(dimension):
            constraint = model.Add(
                degree[0, 0, coordinate] + degree[0, 1, coordinate]
                == _parent[coordinate]
            )
            if selected_branch is not None:
                constraint.OnlyEnforceIf(selected_branch)
        if profiles is not None:
            for other_point, profile in zip(points[1:], profiles):
                pair_set = {0, other_point}
                for coordinate, size in enumerate(sizes):
                    incidence = sum(
                        selected for family in (line, circle)
                        for (block_size, block), selected in family.items()
                        if block_size == size and pair_set.issubset(block)
                    )
                    constraint = model.Add(incidence == profile[coordinate])
                    if selected_branch is not None:
                        constraint.OnlyEnforceIf(selected_branch)
    if branch_variables:
        model.Add(sum(branch_variables) == 1)

    # If the fixed anchor has zero local Melchior defect, its projected line
    # arrangement is simplicial.  Cuntz's complete classification through
    # thirteen lines, together with exact normal vectors, supplies every
    # possible incidence table.  The multiplicity vector identifies a unique
    # table in the orders used here.  We fix only its supports; every support
    # remains free to be a line or a circle in the original configuration.
    anchor_line_family_count = None
    if (anchor_branch_index is not None
            and melchior_defect(anchor_branches[0][0]) == 0):
        zero_entries = [
            entry for entry in simplicial_incidence_catalogue(
                n - 1, dimension
            )
            if entry["parent"] == anchor_branches[0][0]
        ]
        if len(zero_entries) != 1:
            raise ValueError(
                "the zero-defect parent does not identify one classified "
                "simplicial incidence table"
            )
        local_blocks = zero_entries[0]["blocks"]
        # The anchor constraints above label points 1,...,11 in
        # lexicographic order of their pair profiles.  The published/catalogue
        # labels of A(11,1) are arbitrary, so canonicalise them by the same
        # profiles before fixing its supports.  Ties may be broken arbitrarily:
        # the anchor constraints do not distinguish equal-profile points.
        local_profile = {
            point: tuple(
                sum(point in block and len(block) == multiplicity
                    for block in local_blocks)
                for multiplicity in range(2, len(sizes) + 2)
            )
            for point in range(n - 1)
        }
        old_to_canonical = {
            old: new for new, old in enumerate(sorted(
                range(n - 1), key=lambda point: (local_profile[point], point)
            ))
        }
        if (anchor_branches[0][1] is not None
                and tuple(sorted(local_profile.values()))
                != anchor_branches[0][1]):
            raise ValueError(
                "simplicial incidence table and anchor histogram disagree"
            )
        canonical_supports = tuple(
            tuple(sorted(
                (0,) + tuple(old_to_canonical[point] + 1 for point in block)
            ))
            for block in local_blocks
        )
        expected_supports = set(canonical_supports)
        for size in sizes:
            for block in subsets.get(size, ()):
                if 0 not in block:
                    continue
                selected = []
                if (size, block) in line:
                    selected.append(line[size, block])
                if (size, block) in circle:
                    selected.append(circle[size, block])
                model.Add(sum(selected) == int(block in expected_supports))
        # A support through the anchor is a line exactly when its local
        # vertex block lies on the additional projective line obtained by
        # inversion and duality.  The classified tope-graph catalogue is a
        # stronger necessary condition than pairwise disjointness alone.
        line_positions = [
            index for index, support in enumerate(canonical_supports)
            if (len(support), support) in line
        ]
        if line_positions:
            allowed_line_vectors = sorted({
                tuple(int(index in subfamily) for index in line_positions)
                for subfamily in classified_zero_line_families(
                    n - 1, dimension
                )[anchor_branches[0][0]]
                if all(
                    index not in subfamily
                    for index, support in enumerate(canonical_supports)
                    if (len(support), support) not in line
                )
            })
            anchor_line_family_count = len(allowed_line_vectors)
            if anchor_line_family_index is not None:
                if not 0 <= anchor_line_family_index < len(
                        allowed_line_vectors):
                    raise ValueError("anchor line-family index out of range")
                allowed_line_vectors = (
                    allowed_line_vectors[anchor_line_family_index],
                )
            model.AddAllowedAssignments(
                [line[len(canonical_supports[index]),
                      canonical_supports[index]]
                 for index in line_positions],
                tuple(allowed_line_vectors),
            )
    elif anchor_line_family_index is not None:
        raise ValueError(
            "an anchor line family requires a fixed zero-defect branch"
        )

    # Apply the same complete histogram catalogue to every low-defect point,
    # not only to the symmetry-breaking anchor.  For each ordered pair (x,y)
    # a profile indicator records the numbers of size-k blocks containing the
    # pair.  If x has defect at most the catalogue threshold, these indicators
    # must have one of the retained histogram multisets for x's parent type.
    histogram_branches = []
    all_line_profiles = set()
    if parent_only:
        local_parent_branches = local_catalogue_branches
        zero_histograms = classified_zero_profile_histograms(
            n - 1, dimension
        )
        allowed_local_parents = {
            parent for parent, _unused in local_parent_branches
        }
        zero_histograms = {
            parent: histograms for parent, histograms
            in zero_histograms.items() if parent in allowed_local_parents
        }
        zero_histogram_by_parent = {}
        for parent, histograms in zero_histograms.items():
            if len(histograms) != 1:
                raise ValueError(
                    "a zero-defect parent has multiple profile histograms"
                )
            histogram = Counter(next(iter(histograms)))
            zero_histogram_by_parent[parent] = histogram
            all_line_profiles.update(histogram)
        all_line_profiles = tuple(sorted(all_line_profiles))
        zero_pair_statistics = classified_zero_pair_statistics(
            n - 1, dimension
        )
        zero_pair_statistics = {
            parent: statistics for parent, statistics
            in zero_pair_statistics.items() if parent in allowed_local_parents
        }
        zero_pair_keys = tuple(sorted({
            key for statistics in zero_pair_statistics.values()
            for key in statistics
        }))
        zero_line_degree_vectors = classified_zero_line_degree_vectors(
            n - 1, dimension
        )
        zero_line_degree_vectors = {
            parent: vectors for parent, vectors
            in zero_line_degree_vectors.items()
            if parent in allowed_local_parents
        }
        zero_split_branches = []
        all_split_profiles = set()
        if zero_split_histogram_enabled:
            zero_split_catalogue = classified_zero_split_profile_histograms(
                n - 1, dimension
            )
            for row, histograms in zero_split_catalogue.items():
                parent = tuple(
                    row[index] + row[dimension + index]
                    for index in range(dimension)
                )
                if parent not in allowed_local_parents:
                    continue
                for histogram in histograms:
                    counter = Counter(histogram)
                    zero_split_branches.append((row, counter))
                    all_split_profiles.update(counter)
            all_split_profiles = tuple(sorted(all_split_profiles))
    else:
        for parent, expanded_profiles in local_catalogue_branches:
            histogram = Counter(expanded_profiles)
            histogram_branches.append((parent, histogram))
            all_line_profiles.update(histogram)
        all_line_profiles = tuple(sorted(all_line_profiles))
        local_parent_branches = histogram_branches
    point_defects = []
    fixed_minimum_defect = (
        melchior_defect(anchor_branches[0][0])
        if anchor_branch_index is not None else None
    )
    for centre in points:
        defect_expression = (
            degree[centre, 0, 0] + degree[centre, 1, 0] - 3
            - sum(
                (size - 4) * (
                    degree[centre, 0, coordinate]
                    + degree[centre, 1, coordinate]
                )
                for coordinate, size in enumerate(sizes)
                if size >= 5
            )
        )
        defect = model.NewIntVar(
            -n * n, n * n, f"melchior_defect_{centre}"
        )
        model.Add(defect == defect_expression)
        point_defects.append(defect)
        model.Add(sum(
            size * degree[centre, 0, coordinate]
            for coordinate, size in enumerate(sizes)
        ) <= n - 1 + defect)
        # In a split branch, point 0 is chosen among the points of minimum
        # local Melchior defect.  Relabelling any genuine configuration in
        # this way loses no cases and removes duplicate branches in which a
        # smaller-defect point was left elsewhere.
        if anchor_branch_index is not None:
            model.Add(defect >= fixed_minimum_defect)
            model.Add(
                defect <= total_defect - fixed_minimum_defect * (n - 1)
            )
        is_low = model.NewBoolVar(f"low_defect_{centre}")
        model.Add(defect <= maximum_defect).OnlyEnforceIf(is_low)
        model.Add(defect >= maximum_defect + 1).OnlyEnforceIf(is_low.Not())
        local_choices = []
        for branch_index, (parent, _histogram) in enumerate(
                local_parent_branches):
            choice = model.NewBoolVar(
                f"local_histogram_{centre}_{branch_index}"
            )
            local_choices.append(choice)
            for coordinate in range(dimension):
                model.Add(
                    degree[centre, 0, coordinate]
                    + degree[centre, 1, coordinate]
                    == parent[coordinate]
                ).OnlyEnforceIf(choice)
            if parent_only and melchior_defect(parent) == 0:
                model.AddAllowedAssignments(
                    [degree[centre, 0, coordinate]
                     for coordinate in range(dimension)],
                    zero_line_degree_vectors[parent],
                ).OnlyEnforceIf(choice)
            if (parent_only and anchor_branch_index is not None
                    and melchior_defect(parent) == fixed_minimum_defect
                    and parent < anchor_branches[0][0]):
                # Among all minimum-defect points, choose the anchor with
                # lexicographically least local multiplicity vector.
                model.Add(choice == 0)
        model.Add(sum(local_choices) == is_low)
        if parent_only:
            zero_choice_indices = [
                index for index, (parent, _unused)
                in enumerate(local_parent_branches)
                if melchior_defect(parent) == 0
            ]
            is_zero = sum(
                local_choices[index] for index in zero_choice_indices
            )
            profile_indicators = {
                (other, profile_index): model.NewBoolVar(
                    f"zero_pair_profile_{centre}_{other}_{profile_index}"
                )
                for other in points if other != centre
                for profile_index in range(len(all_line_profiles))
            }
            for other in points:
                if other == centre:
                    continue
                model.Add(sum(
                    profile_indicators[other, profile_index]
                    for profile_index in range(len(all_line_profiles))
                ) == is_zero)
                pair_set = {centre, other}
                for profile_index, profile in enumerate(all_line_profiles):
                    indicator = profile_indicators[other, profile_index]
                    for coordinate, size in enumerate(sizes):
                        incidence = sum(
                            selected for family in (line, circle)
                            for (block_size, block), selected
                            in family.items()
                            if block_size == size
                            and pair_set.issubset(block)
                        )
                        model.Add(
                            incidence == profile[coordinate]
                        ).OnlyEnforceIf(indicator)
            for profile_index, profile in enumerate(all_line_profiles):
                model.Add(sum(
                    profile_indicators[other, profile_index]
                    for other in points if other != centre
                ) == sum(
                    local_choices[index]
                    * zero_histogram_by_parent[parent].get(profile, 0)
                    for index, (parent, _unused)
                    in enumerate(local_parent_branches)
                    if melchior_defect(parent) == 0
                ))
            if zero_split_histogram_enabled:
                split_choices = tuple(
                    model.NewBoolVar(
                        f"zero_split_choice_{centre}_{branch_index}"
                    )
                    for branch_index in range(len(zero_split_branches))
                )
                model.Add(sum(split_choices) == is_zero)
                for branch_index, (row, _histogram) in enumerate(
                        zero_split_branches):
                    for coordinate in range(dimension):
                        model.Add(
                            degree[centre, 0, coordinate]
                            == row[coordinate]
                        ).OnlyEnforceIf(split_choices[branch_index])
                        model.Add(
                            degree[centre, 1, coordinate]
                            == row[dimension + coordinate]
                        ).OnlyEnforceIf(split_choices[branch_index])
                split_indicators = {
                    (other, profile_index): model.NewBoolVar(
                        f"zero_split_profile_{centre}_{other}_{profile_index}"
                    )
                    for other in points if other != centre
                    for profile_index in range(len(all_split_profiles))
                }
                for other in points:
                    if other == centre:
                        continue
                    model.Add(sum(
                        split_indicators[other, profile_index]
                        for profile_index in range(len(all_split_profiles))
                    ) == is_zero)
                    pair_set = {centre, other}
                    for profile_index, profile in enumerate(all_split_profiles):
                        indicator = split_indicators[other, profile_index]
                        for family_index, family in enumerate((line, circle)):
                            for coordinate, size in enumerate(sizes):
                                incidence = sum(
                                    selected for (block_size, block), selected
                                    in family.items()
                                    if block_size == size
                                    and pair_set.issubset(block)
                                )
                                model.Add(
                                    incidence
                                    == profile[
                                        family_index * dimension + coordinate
                                    ]
                                ).OnlyEnforceIf(indicator)
                for profile_index, profile in enumerate(all_split_profiles):
                    model.Add(sum(
                        split_indicators[other, profile_index]
                        for other in points if other != centre
                    ) == sum(
                        split_choices[branch_index]
                        * histogram.get(profile, 0)
                        for branch_index, (_row, histogram)
                        in enumerate(zero_split_branches)
                    ))
            if not zero_pair_statistics_enabled:
                continue

            profile_index_of = {
                profile: index for index, profile
                in enumerate(all_line_profiles)
            }
            pair_key_variables = {
                key: [] for key in zero_pair_keys
            }
            other_points = tuple(
                point for point in points if point != centre
            )
            for first, second in combinations(other_points, 2):
                pair_variables = []
                triple_set = {centre, first, second}
                for low_profile, high_profile, multiplicity in zero_pair_keys:
                    orientations = [(low_profile, high_profile)]
                    if low_profile != high_profile:
                        orientations.append((high_profile, low_profile))
                    relation = sum(
                        selected for family in (line, circle)
                        for (block_size, block), selected in family.items()
                        if block_size == multiplicity + 1
                        and triple_set.issubset(block)
                    )
                    for first_profile, second_profile in orientations:
                        variable = model.NewBoolVar(
                            "zero_pair_stat_"
                            f"{centre}_{first}_{second}_"
                            f"{profile_index_of[first_profile]}_"
                            f"{profile_index_of[second_profile]}_"
                            f"{multiplicity}"
                        )
                        first_indicator = profile_indicators[
                            first, profile_index_of[first_profile]
                        ]
                        second_indicator = profile_indicators[
                            second, profile_index_of[second_profile]
                        ]
                        model.Add(variable <= first_indicator)
                        model.Add(variable <= second_indicator)
                        model.Add(variable <= relation)
                        model.Add(
                            variable >= first_indicator
                            + second_indicator + relation - 2
                        )
                        pair_variables.append(variable)
                        pair_key_variables[
                            low_profile, high_profile, multiplicity
                        ].append(variable)
                model.Add(sum(pair_variables) == is_zero)
            for key, variables in pair_key_variables.items():
                model.Add(sum(variables) == sum(
                    local_choices[index]
                    * zero_pair_statistics[parent].get(key, 0)
                    for index, (parent, _unused)
                    in enumerate(local_parent_branches)
                    if melchior_defect(parent) == 0
                ))
            continue

        profile_indicators = {
            (other, profile_index): model.NewBoolVar(
                f"pair_profile_{centre}_{other}_{profile_index}"
            )
            for other in points if other != centre
            for profile_index in range(len(all_line_profiles))
        }
        for other in points:
            if other == centre:
                continue
            model.Add(sum(
                profile_indicators[other, profile_index]
                for profile_index in range(len(all_line_profiles))
            ) == is_low)
            pair_set = {centre, other}
            for profile_index, profile in enumerate(all_line_profiles):
                indicator = profile_indicators[other, profile_index]
                for coordinate, size in enumerate(sizes):
                    incidence = sum(
                        selected for family in (line, circle)
                        for (block_size, block), selected in family.items()
                        if block_size == size and pair_set.issubset(block)
                    )
                    model.Add(
                        incidence == profile[coordinate]
                    ).OnlyEnforceIf(indicator)
        for profile_index, profile in enumerate(all_line_profiles):
            model.Add(sum(
                profile_indicators[other, profile_index]
                for other in points if other != centre
            ) == sum(
                local_choices[branch_index]
                * histogram.get(profile, 0)
                for branch_index, (_parent, histogram)
                in enumerate(histogram_branches)
            ))
    # This equality follows by summing the local incidence degrees, but
    # recording it explicitly greatly strengthens propagation in the small
    # positive-defect branches.
    model.Add(sum(point_defects) == total_defect)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = solver_workers
    solver.parameters.random_seed = random_seed
    status = solver.Solve(model)
    return {
        "signature": list(signature),
        "circle_count": circle_count,
        "anchor_branch_count": len(anchor_branches),
        "anchor_branch_index": anchor_branch_index,
        "anchor_defect_bound": anchor_defect_bound,
        "maximum_deletion_size": maximum_deletion_size,
        "solver_workers": solver_workers,
        "random_seed": random_seed,
        "parent_only": parent_only,
        "anchor_row_index": anchor_row_index,
        "anchor_line_family_index": anchor_line_family_index,
        "anchor_line_family_count": anchor_line_family_count,
        "zero_pair_statistics_enabled": zero_pair_statistics_enabled,
        "zero_split_histogram_enabled": zero_split_histogram_enabled,
        "status": solver.StatusName(status),
        "wall_time_seconds": solver.WallTime(),
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
    }


def solve_task(task):
    return solve_signature(*task)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=tuple(INPUTS), required=True)
    parser.add_argument("--input-json", type=Path)
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--anchor-branch", type=int)
    parser.add_argument("--anchor-row", type=int)
    parser.add_argument("--anchor-line-family", type=int)
    parser.add_argument("--split-anchor-branches", action="store_true")
    parser.add_argument("--maximum-deletion-size", type=int, default=2)
    parser.add_argument("--solver-workers", type=int, default=1)
    parser.add_argument("--random-seed", type=int, default=0)
    parser.add_argument(
        "--parent-only", action="store_true",
        help="use only local multiplicity vectors, omitting pair-profile histograms",
    )
    parser.add_argument("--no-zero-pair-statistics", action="store_true")
    parser.add_argument("--no-zero-split-histogram", action="store_true")
    parser.add_argument("--maximum-anchor-defect-bound", type=int)
    parser.add_argument("--minimum-anchor-defect-bound", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.input_json or HERE / INPUTS[args.n]
    records = retained_records(args.n, source)
    if args.maximum_anchor_defect_bound is not None:
        filtered = []
        for record in records:
            signature = tuple(record["signature"])
            dimension = len(signature) // 2
            sizes = tuple(range(3, 3 + dimension))
            line_numbers = signature[:dimension]
            circle_numbers = signature[dimension:]
            total_defect = (
                3 * (line_numbers[0] + circle_numbers[0]) - 3 * args.n
                - sum(
                    (size - 4) * size * (line_number + circle_number)
                    for size, line_number, circle_number in zip(
                        sizes[2:], line_numbers[2:], circle_numbers[2:]
                    )
                )
            )
            if (total_defect // args.n
                    <= args.maximum_anchor_defect_bound
                    and (args.minimum_anchor_defect_bound is None
                         or total_defect // args.n
                         >= args.minimum_anchor_defect_bound)):
                filtered.append(record)
        records = tuple(filtered)
    elif args.minimum_anchor_defect_bound is not None:
        raise ValueError(
            "--minimum-anchor-defect-bound requires the maximum bound"
        )
    if args.limit is not None:
        records = records[:args.limit]
    if args.split_anchor_branches and args.anchor_branch is not None:
        raise ValueError("choose only one anchor-branch mode")
    if args.anchor_row is not None and args.anchor_branch is None:
        raise ValueError("--anchor-row requires --anchor-branch")
    if args.anchor_line_family is not None and args.anchor_branch is None:
        raise ValueError(
            "--anchor-line-family requires --anchor-branch"
        )
    if args.split_anchor_branches:
        dimension = len(records[0]["signature"]) // 2 if records else 0
        maximum_block_size = 2 + dimension
        branch_loader = (
            load_anchor_parent_types if args.parent_only
            else load_anchor_profiles
        )
        _path, branches = branch_loader(args.n, maximum_block_size)
        tasks = [
            (args.n, tuple(record["signature"]), args.seconds, branch_index,
             args.maximum_deletion_size, args.solver_workers)
            + (args.random_seed, args.parent_only, None,
               not args.no_zero_pair_statistics, None,
               not args.no_zero_split_histogram)
            for record in records
            for branch_index in range(len(branches))
        ]
    else:
        tasks = [
            (args.n, tuple(record["signature"]), args.seconds,
             args.anchor_branch, args.maximum_deletion_size,
             args.solver_workers, args.random_seed, args.parent_only,
             args.anchor_row, not args.no_zero_pair_statistics,
             args.anchor_line_family, not args.no_zero_split_histogram)
            for record in records
        ]
    if args.jobs == 1:
        results = list(map(solve_task, tasks))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            results = list(executor.map(solve_task, tasks, chunksize=1))
    counts = Counter(record["status"] for record in results)
    if args.split_anchor_branches:
        by_signature = {}
        for result in results:
            by_signature.setdefault(tuple(result["signature"]), []).append(
                result
            )
        signature_records = []
        for signature, branch_results in by_signature.items():
            statuses = Counter(result["status"] for result in branch_results)
            if any(status in ("OPTIMAL", "FEASIBLE") for status in statuses):
                status = "FEASIBLE_BRANCH"
            elif any(not status.startswith("INFEASIBLE")
                     for status in statuses):
                status = "UNKNOWN"
            else:
                status = "INFEASIBLE"
            signature_records.append({
                "signature": list(signature),
                "status": status,
                "branch_status_counts": dict(statuses),
            })
        unclosed = [
            record for record in signature_records
            if record["status"] != "INFEASIBLE"
        ]
    else:
        signature_records = results
        unclosed = [
            record for record in results
            if not record["status"].startswith("INFEASIBLE")
        ]
    if args.parent_only:
        catalogue_path = load_anchor_parent_types(
            args.n, 2 + len(records[0]["signature"]) // 2
            if records else {12: 6, 13: 6, 14: 7}[args.n]
        )[0]
    else:
        catalogue_path = HERE / (
            "n12_low_defect_line_histograms_wiring_retry300.json"
            if args.n == 12
            else f"n{args.n}_low_defect_line_histograms_pbd.json"
        )
    payload = {
        "schema_version": 1,
        "status": "PASS" if not unclosed else "INCOMPLETE",
        "claim": (
            "Exact support filter anchored at a low-defect point; every "
            + ("retained local multiplicity type"
               if args.parent_only else "retained profile histogram")
            + " and every block support is covered."
        ),
        "n": args.n,
        "source_file": source.name,
        "source_sha256": sha256(source.read_bytes()).hexdigest(),
        "catalogue_file": catalogue_path.name,
        "catalogue_sha256": sha256(catalogue_path.read_bytes()).hexdigest(),
        "input_signature_count": len(records),
        "branch_task_count": len(tasks),
        "split_anchor_branches": args.split_anchor_branches,
        "parent_only": args.parent_only,
        "zero_pair_statistics_enabled": (
            not args.no_zero_pair_statistics
        ),
        "zero_split_histogram_enabled": (
            not args.no_zero_split_histogram
        ),
        "seconds_per_signature": args.seconds,
        "status_counts": dict(counts),
        "unclosed_count": len(unclosed),
        "unclosed": unclosed,
        "signature_records": signature_records,
        "records": results,
    }
    args.output.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": payload["status"],
        "input_signature_count": len(records),
        "branch_task_count": len(tasks),
        "status_counts": dict(counts),
        "unclosed_count": len(unclosed),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
