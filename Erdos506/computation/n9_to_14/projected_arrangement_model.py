#!/usr/bin/env python3
"""Shared projected-line-arrangement constraints for 9 <= n <= 15.

For a point x of an admissible configuration P, invert about x and dualize
the remaining n-1 points.  A maximal k-point line or circle through x then
becomes a vertex of multiplicity k-1 in an arrangement of n-1 real
projective lines.  This module generates necessary local multiplicity types
and their splits into line-block and circle-block incidences.

The module is deliberately independent of the historical node tree and of
the corrected support-free signature diagnostics.  It contains no theorem
claim and writes no output file.  Dedicated verifiers decide which uses of
the generated constraints are complete certificates.
"""

from __future__ import annotations

from functools import lru_cache
from math import ceil, comb


MAXIMUM_BLOCK_SIZE = {
    9: 5,
    10: 6,
    11: 6,
    12: 6,
    13: 6,
    14: 7,
    15: 7,
}


# Invariant vectors (t_2,t_3,...) of the stretchable simplicial
# arrangements in Cuntz's complete table, with the near-pencil omitted.
# A near-pencil has a vertex of multiplicity q-1 and is automatically absent
# in every branch covered by MAXIMUM_BLOCK_SIZE.
CUNTZ_SIMPLICIAL_TYPES = {
    8: {(4, 6, 1)},
    9: {(6, 4, 3)},
    10: {(5, 10, 0, 1), (6, 7, 3)},
    11: {(7, 8, 4)},
    12: {(6, 15, 0, 0, 1), (8, 10, 3, 1), (9, 7, 6)},
    13: {(9, 12, 3, 0, 1), (12, 4, 9),
         (10, 10, 3, 2), (6, 18, 3)},
    14: {(7, 21, 0, 0, 0, 1), (11, 12, 4, 2),
         (9, 16, 4, 1), (10, 14, 4, 0, 1)},
}


# Classical orchard upper bounds for the number of lines containing exactly
# three points.  Burr--Gruenbaum--Sloane's definition does not exclude other
# lines containing four or more points.  The entries through twelve points
# are exact; their Table I gives upper bounds 24 and 27 at thirteen and
# fourteen points, respectively.
ORCHARD_UPPER_BOUNDS = {
    4: 1,
    5: 2,
    6: 4,
    7: 6,
    8: 7,
    9: 10,
    10: 12,
    11: 16,
    12: 19,
    13: 24,
    14: 27,
}


def choose(n: int, r: int) -> int:
    return comb(n, r) if 0 <= r <= n else 0


def ordinary_line_lower_bound(number_of_points: int) -> int:
    """The finite bounds used for the dual point arrangements."""
    if number_of_points >= 9:
        return ceil(6 * number_of_points / 13)
    return {8: 4, 7: 3}.get(number_of_points, 2)


def melchior_defect(local_type: tuple[int, ...]) -> int:
    """Return t_2-3-sum_{r>=4}(r-3)t_r."""
    return (local_type[0] - 3
            - sum((r - 3) * local_type[r - 2]
                  for r in range(4, len(local_type) + 2)))


def edge_count(local_type: tuple[int, ...]) -> int:
    """Number of projective edges, sum r t_r."""
    return sum(r * local_type[r - 2]
               for r in range(2, len(local_type) + 2))


def _standard_line_inequalities(
    number_of_lines: int,
    local_type: tuple[int, ...],
    use_elliott_line_count_bound: bool = True,
) -> bool:
    """Classical necessary inequalities for a real line arrangement."""
    multiplicity = {
        r: local_type[r - 2] if r <= len(local_type) + 1 else 0
        for r in range(2, number_of_lines)
    }
    t2 = multiplicity[2]
    t3 = multiplicity.get(3, 0)
    if t2 < ordinary_line_lower_bound(number_of_lines):
        return False
    if 4 * t2 + 3 * t3 < (
            4 * number_of_lines
            + 4 * sum((2 * r - 9) * multiplicity[r]
                      for r in multiplicity if r >= 5)):
        return False
    maximum_multiplicity = max(
        (r for r, number in multiplicity.items() if number), default=0
    )
    if (3 * maximum_multiplicity < 2 * number_of_lines
            and 4 * t2 + 3 * t3 < (
                4 * number_of_lines
                + sum(r * (r - 4) * multiplicity[r]
                      for r in multiplicity if r >= 5))):
        return False
    if 2 * t2 + 3 * t3 < (
            16 + sum((4 * r - 15) * multiplicity[r]
                     for r in multiplicity if r >= 4)):
        return False
    # Elliott, Theorem 3: q points with no q-1 collinear determine at least
    # 2q-4 connecting lines when q>=10.  After point-line duality this is a
    # lower bound for the number of vertices.  The maximum-multiplicity
    # condition below is exactly the missing-hyperplane hypothesis.
    if (use_elliott_line_count_bound
            and number_of_lines >= 10
            and max((r for r, number in multiplicity.items() if number),
                    default=0) <= number_of_lines - 2
            and sum(local_type) < 2 * number_of_lines - 4):
        return False
    return True


def _block_capacity_inequalities(
    number_of_points: int,
    local_type: tuple[int, ...],
) -> bool:
    """Necessary incidence bounds for all maximal connecting lines.

    Regard the dual arrangement as a set of ``number_of_points`` points and
    each multiplicity-r vertex as its maximal r-point connecting line.  The
    union of j distinct lines has at least the sum of their sizes minus
    C(j,2) points.  Moreover, a fixed r-point line needs r(q-r) incidences
    with the other maximal lines, and two fixed lines of sizes a and b need
    at least (a-1)(b-1) further connecting lines.
    """
    sizes = tuple(range(2, len(local_type) + 2))
    line_sizes = sorted(
        (size for size, number in zip(sizes, local_type)
         for _ in range(number)),
        reverse=True,
    )
    partial = 0
    for number, size in enumerate(line_sizes, 1):
        partial += size
        if partial - choose(number, 2) > number_of_points:
            return False

    weighted = sum(
        (size - 1) * number
        for size, number in zip(sizes, local_type)
    )
    line_count = sum(local_type)
    for size, number in zip(sizes, local_type):
        if number and weighted - (size - 1) < size * (
                number_of_points - size):
            return False
    for index, first_size in enumerate(sizes):
        for second_index in range(index, len(sizes)):
            second_size = sizes[second_index]
            if index == second_index:
                present = local_type[index] >= 2
            else:
                present = bool(
                    local_type[index] and local_type[second_index]
                )
            if present and line_count - 2 < (
                    first_size - 1) * (second_size - 1):
                return False
    return True


def _subset_orchard_inequalities(
    number_of_points: int,
    local_type: tuple[int, ...],
) -> bool:
    """Average the orchard bound over every subset of each size."""
    sizes = tuple(range(2, len(local_type) + 2))
    for subset_size in range(4, number_of_points + 1):
        inherited_three_point_lines = sum(
            number * choose(size, 3)
            * choose(number_of_points - size, subset_size - 3)
            for size, number in zip(sizes, local_type)
        )
        if inherited_three_point_lines > (
                choose(number_of_points, subset_size)
                * ORCHARD_UPPER_BOUNDS[subset_size]):
            return False
    return True


def _elliott_avoiding_point_subset_inequalities(
    number_of_image_points: int,
    circle_part: tuple[int, ...],
) -> bool:
    """Sum Elliott's external-point line lemma over image-point subsets.

    A circle of size k through the inversion centre becomes a connecting
    line on r=k-1 image points which avoids that centre.  Its contribution
    to the sum over m-subsets is the number of subsets meeting it at least
    twice.  For m larger than the maximum image-line size every subset is
    noncollinear, so Elliott's Lemma 3 gives m-1 avoiding lines.  The
    strengthening immediately following the lemma gives m unless m-1 of the
    subset are collinear.
    """
    image_line_sizes = tuple(range(2, len(circle_part) + 2))
    maximum_image_line_size = image_line_sizes[-1]
    for subset_size in range(
            maximum_image_line_size + 1, number_of_image_points + 1):
        inherited_avoiding_lines = sum(
            number * (
                choose(number_of_image_points, subset_size)
                - choose(number_of_image_points - size, subset_size)
                - size * choose(
                    number_of_image_points - size, subset_size - 1
                )
            )
            for size, number in zip(image_line_sizes, circle_part)
        )
        per_subset = (
            subset_size
            if maximum_image_line_size <= subset_size - 2
            else subset_size - 1
        )
        if inherited_avoiding_lines < (
                per_subset
                * choose(number_of_image_points, subset_size)):
            return False
    return True


def _summed_subset_line_inequalities(
    number_of_points: int,
    local_type: tuple[int, ...],
    use_elliott_line_count_bound: bool = True,
) -> bool:
    """Sum standard line inequalities over all noncollinear subsets."""
    maximum_line_size = len(local_type) + 1
    sizes = tuple(range(2, maximum_line_size + 1))
    for subset_size in range(maximum_line_size + 1, number_of_points):
        subset_count = choose(number_of_points, subset_size)
        ordinary = defect = bojanowski = shnurnikov = 0
        for size, number in zip(sizes, local_type):
            for intersection in range(2, min(size, subset_size) + 1):
                multiplicity = (
                    number * choose(size, intersection)
                    * choose(
                        number_of_points - size,
                        subset_size - intersection,
                    )
                )
                if intersection == 2:
                    ordinary += multiplicity
                defect += multiplicity * (
                    1 if intersection == 2 else
                    (-(intersection - 3) if intersection >= 4 else 0)
                )
                bojanowski += multiplicity * (
                    4 if intersection == 2 else
                    (3 if intersection == 3 else
                     (4 * intersection - intersection * intersection
                      if intersection >= 5 else 0))
                )
                shnurnikov += multiplicity * (
                    2 if intersection == 2 else
                    (3 if intersection == 3 else
                     (15 - 4 * intersection
                     if intersection >= 4 else 0))
                )
        if ordinary < (
                ordinary_line_lower_bound(subset_size) * subset_count):
            return False
        if defect < 3 * subset_count:
            return False
        if (3 * maximum_line_size < 2 * subset_size
                and bojanowski < 4 * subset_size * subset_count):
            return False
        if (subset_size >= maximum_line_size + 3
                and shnurnikov < 16 * subset_count):
            return False
        # Sum Elliott's Theorem 3 over all subsets for which its two
        # hypotheses are forced: subset_size>=10 and no subset_size-1
        # collinear points.
        inherited_line_count = sum(
            number * (
                subset_count
                - choose(number_of_points - size, subset_size)
                - size * choose(
                    number_of_points - size, subset_size - 1
                )
            )
            for size, number in zip(sizes, local_type)
        )
        if (use_elliott_line_count_bound
                and subset_size >= 10
                and maximum_line_size <= subset_size - 2
                and inherited_line_count
                < (2 * subset_size - 4) * subset_count):
            return False
    return True


def _cuntz_zero_defect_type(
    number_of_lines: int,
    local_type: tuple[int, ...],
) -> bool:
    padded = local_type + (0,) * (number_of_lines - 1 - len(local_type))
    return padded in {
        row + (0,) * (number_of_lines - 1 - len(row))
        for row in CUNTZ_SIMPLICIAL_TYPES[number_of_lines]
    }


def _projective_face_congruence(
    number_of_lines: int,
    local_type: tuple[int, ...],
) -> bool:
    """Euler/face congruence, with the even-line two-colour gap.

    In every projective arrangement delta is congruent to twice the edge
    count modulo three.  If the number of lines is even, the product of the
    homogeneous line equations colours the projective faces.  The two colour
    classes have the same edge residue; consequently a positive defect is at
    least 2 when E=1 mod 3 and at least 4 when E=2 mod 3.
    """
    defect = melchior_defect(local_type)
    residue = edge_count(local_type) % 3
    if defect % 3 != (2 * residue) % 3:
        return False
    if number_of_lines % 2 == 0 and residue and defect < 2 * residue:
        return False
    return True


@lru_cache(maxsize=None)
def projected_local_types(
    n: int,
    maximum_block_size: int | None = None,
    use_elliott_line_count_bound: bool = True,
) -> tuple[tuple[int, ...], ...]:
    """Generate all retained (t_2,...,t_{K-1}) local types."""
    if maximum_block_size is None:
        maximum_block_size = MAXIMUM_BLOCK_SIZE[n]
    number_of_lines = n - 1
    maximum_multiplicity = maximum_block_size - 1
    work = [0] * (maximum_multiplicity - 1)
    answer: list[tuple[int, ...]] = []

    def rec(multiplicity: int, remaining_pairs: int) -> None:
        if multiplicity == 2:
            local_type = (remaining_pairs,) + tuple(
                work[r - 2] for r in range(3, maximum_multiplicity + 1))
            defect = melchior_defect(local_type)
            if defect < 0:
                return
            if not _projective_face_congruence(
                    number_of_lines, local_type):
                return
            if defect == 0 and not _cuntz_zero_defect_type(
                    number_of_lines, local_type):
                return
            # Under point-line duality, t_3 is exactly the number of
            # three-point lines in the corresponding orchard problem.
            t3 = local_type[1] if len(local_type) >= 2 else 0
            if t3 > ORCHARD_UPPER_BOUNDS[number_of_lines]:
                return
            if not _standard_line_inequalities(
                    number_of_lines, local_type,
                    use_elliott_line_count_bound):
                return
            if not _block_capacity_inequalities(
                    number_of_lines, local_type):
                return
            if not _subset_orchard_inequalities(
                    number_of_lines, local_type):
                return
            if not _summed_subset_line_inequalities(
                    number_of_lines, local_type,
                    use_elliott_line_count_bound):
                return
            answer.append(local_type)
            return
        weight = choose(multiplicity, 2)
        for count in range(remaining_pairs // weight + 1):
            work[multiplicity - 2] = count
            rec(multiplicity - 1, remaining_pairs - count * weight)
        work[multiplicity - 2] = 0

    rec(maximum_multiplicity, choose(number_of_lines, 2))
    return tuple(answer)


@lru_cache(maxsize=None)
def category_point_rows(
    n: int,
    maximum_block_size: int | None = None,
) -> tuple[tuple[int, ...], ...]:
    """Split every local type into line-block and circle-block counts."""
    if maximum_block_size is None:
        maximum_block_size = MAXIMUM_BLOCK_SIZE[n]
    local_types = projected_local_types(n, maximum_block_size)
    dimension = maximum_block_size - 2
    line_part = [0] * dimension
    answer: list[tuple[int, ...]] = []

    for local_type in local_types:
        def rec(index: int, used_other_points: int) -> None:
            if index == dimension:
                circle_part = tuple(
                    local_type[j] - line_part[j]
                    for j in range(dimension)
                )
                # Zhang's pointed ordinary-line theorem, applied after
                # inversion about x, gives at least ceil((n-1)/6) ordinary
                # image lines avoiding the inversion centre.  These are
                # exactly the original three-point circles through x.
                if circle_part[0] < ceil((n - 1) / 6):
                    return
                # Elliott, Lemma 3 and the strengthening immediately after
                # it: q>=4 noncollinear points determine at least q lines
                # avoiding a prescribed external point unless q-1 of them
                # are collinear.  Here q=n-1 and a collinear q-1 image set
                # would give an original block of size n-1, excluded by the
                # present maximum-block branch.  These avoiding image lines
                # are exactly the original circles through the centre.
                if sum(circle_part) < n - 1:
                    return
                if not _elliott_avoiding_point_subset_inequalities(
                        n - 1, circle_part):
                    return
                answer.append(
                    tuple(line_part)
                    + circle_part)
                return
            multiplicity = index + 2
            upper = min(
                local_type[index],
                (n - 1 - used_other_points) // multiplicity,
            )
            for count in range(upper + 1):
                line_part[index] = count
                rec(index + 1,
                    used_other_points + multiplicity * count)
            line_part[index] = 0

        rec(0, 0)
    return tuple(answer)


@lru_cache(maxsize=None)
def projected_global_block_totals(
    n: int,
    maximum_block_size: int | None = None,
) -> tuple[tuple[int, ...], ...]:
    """Return block-count totals permitted by all projected point rows.

    If ``B_k`` denotes the total number of selected k-point lines and
    circles, double counting incidences gives

        sum_x t_{k-1}(x) = k B_k.

    The returned rows are the possible vectors ``(B_3,...,B_K)`` obtained by
    summing ``n`` retained local types and applying these divisibilities.  It
    is only an aggregate necessary condition; the support model still imposes
    every point row separately.
    """
    if maximum_block_size is None:
        maximum_block_size = MAXIMUM_BLOCK_SIZE[n]
    local_types = projected_local_types(n, maximum_block_size)
    dimension = maximum_block_size - 2
    sums = {(0,) * dimension}
    for _ in range(n):
        sums = {
            tuple(left + right for left, right in zip(total, local_type))
            for total in sums
            for local_type in local_types
        }
    rows = {
        tuple(total[k - 3] // k
              for k in range(3, maximum_block_size + 1))
        for total in sums
        if all(total[k - 3] % k == 0
               for k in range(3, maximum_block_size + 1))
    }
    return tuple(sorted(rows))


@lru_cache(maxsize=None)
def projected_global_zero_counts(
    n: int,
    maximum_block_size: int | None = None,
) -> dict[tuple[int, ...], tuple[int, ...]]:
    """Possible numbers of zero-defect centres at each block total.

    The key is ``(B_3,...,B_K)`` as in
    :func:`projected_global_block_totals`; the value lists every number of
    local Melchior-defect-zero rows which can occur in an n-row sum.  This is
    a small dynamic programme on local weak vectors, not a support search.
    """
    if maximum_block_size is None:
        maximum_block_size = MAXIMUM_BLOCK_SIZE[n]
    local_types = projected_local_types(n, maximum_block_size)
    dimension = maximum_block_size - 2
    states = {((0,) * dimension, 0)}
    for _ in range(n):
        states = {
            (tuple(left + right
                   for left, right in zip(total, local_type)),
             zero_count + int(melchior_defect(local_type) == 0))
            for total, zero_count in states
            for local_type in local_types
        }
    answer: dict[tuple[int, ...], set[int]] = {}
    for total, zero_count in states:
        if not all(total[k - 3] % k == 0
                   for k in range(3, maximum_block_size + 1)):
            continue
        row = tuple(
            total[k - 3] // k
            for k in range(3, maximum_block_size + 1)
        )
        answer.setdefault(row, set()).add(zero_count)
    return {
        row: tuple(sorted(counts))
        for row, counts in sorted(answer.items())
    }


def n15_local_weight(local_type: tuple[int, ...]) -> int:
    """The local weight 2t_3+5t_4+9t_5+14t_6 for n=15."""
    assert len(local_type) == 5
    return (2 * local_type[1] + 5 * local_type[2]
            + 9 * local_type[3] + 14 * local_type[4])
