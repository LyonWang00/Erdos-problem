#!/usr/bin/env python3
"""Self-contained global signature generator and mathematical prefilters.

Only the general inequalities used by the active n<=15 diagnostics occur in
this module.  It has no solver, finite arrangement catalogue, historical
feature switch, command-line entry point, certificate reader, or output path.
"""

from __future__ import annotations

from functools import lru_cache
from math import ceil, comb


KNOWN = {4: 3, 5: 5, 6: 8, 7: 11, 8: 17,
         9: 25, 10: 33, 11: 41, 12: 51, 13: 61, 14: 73}

# Maximum number of lines containing exactly three points in a real
# m-point set.  The entries through m=12 are exact; the last two are the
# published upper bounds needed here.  Importantly, the orchard problem
# allows other lines to contain four or more points, so these bounds apply
# directly to the maximal-line signature used below.
ORCHARD_UPPER_BOUNDS = {
    4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12,
    11: 16, 12: 19, 13: 24, 14: 27,
}

SUMMED_LOCAL_LINE_INEQUALITIES = {
    11: (((1, 2, 3), 13),),
    12: (((2, 5, 9), 36), ((1, 2, 3), 16), ((1, 3, 6), 21)),
    13: (((1, 3, 5), 26), ((2, 4, 7), 40), ((0, 1, 1), 8)),
    14: (((2, 5, 9, 14), 53), ((8, 23, 42, 77), 242),
         ((1, 3, 5, 10), 31), ((7, 19, 36, 64), 205),
         ((1, 4, 7, 13), 41)),
    15: (((2, 5, 9, 14), 63), ((3, 9, 17, 30), 114),
         ((1, 2, 4, 6), 29), ((1, 3, 5, 9), 37),
         ((2, 6, 11, 19), 75)),
}


def choose(n: int, r: int) -> int:
    return comb(n, r) if 0 <= r <= n else 0


def ordinary_lower(n: int) -> int:
    if n >= 9:
        return ceil(6 * n / 13)
    return {8: 4, 7: 3}.get(n, 2)


def balanced_second_moment(total: int, positions: int) -> int:
    quotient, remainder = divmod(total, positions)
    return positions * choose(quotient, 2) + quotient * remainder


def block_family_convexity_witness(
    n: int,
    sizes: tuple[int, ...],
    signature: tuple[int, ...],
) -> dict[str, int | str] | None:
    """Use the largest numerical subfamilies and block-intersection bounds."""
    d = len(sizes)
    assert len(signature) == 2 * d

    def prefix_sums(counts: tuple[int, ...], pair_level: bool) -> list[int]:
        weights = sorted(
            ((choose(k, 2) if pair_level else k)
             for k, number in zip(sizes, counts) for _ in range(number)),
            reverse=True,
        )
        answer = [0]
        for weight in weights:
            answer.append(answer[-1] + weight)
        return answer

    line_points = prefix_sums(signature[:d], False)
    circle_points = prefix_sums(signature[d:], False)
    line_pairs = prefix_sums(signature[:d], True)
    circle_pairs = prefix_sums(signature[d:], True)
    for lines in range(len(line_points)):
        for circles in range(len(circle_points)):
            if lines + circles < 2:
                continue
            lower = balanced_second_moment(
                line_points[lines] + circle_points[circles], n)
            upper = (choose(lines, 2) + 2 * lines * circles
                     + 2 * choose(circles, 2))
            if lower > upper:
                return {"level": "points", "lines": lines,
                        "circles": circles, "lower": lower, "upper": upper}
            lower = balanced_second_moment(
                line_pairs[lines] + circle_pairs[circles], choose(n, 2))
            upper = lines * circles + choose(circles, 2)
            if lower > upper:
                return {"level": "point_pairs", "lines": lines,
                        "circles": circles, "lower": lower, "upper": upper}
    return None


def summed_local_line_witness(
    n: int,
    sizes: tuple[int, ...],
    signature: tuple[int, ...],
) -> dict[str, object] | None:
    d = len(sizes)
    blocks = tuple(signature[j] + signature[d + j] for j in range(d))
    totals = tuple(k * number for k, number in zip(sizes, blocks))[1:]
    for coefficients, local_bound in SUMMED_LOCAL_LINE_INEQUALITIES[n]:
        value = sum(a * total for a, total in zip(coefficients, totals))
        bound = n * local_bound
        if value > bound:
            return {"coefficients_for_t3_upwards": coefficients,
                    "local_bound": local_bound, "summed_value": value,
                    "summed_bound": bound, "excess": value - bound}
    return None


@lru_cache(maxsize=None)
def candidate_signature_set(
    n: int, circle_count: int,
) -> frozenset[tuple[int, ...]]:
    """All globally permitted signatures at a smaller order.

    The function body is evaluated only after this module has finished
    defining :class:`SignatureFilter`.  It is used for exact one-point
    deletion and inversion transfers, never as a theorem-level assumption.
    """
    return frozenset(SignatureFilter(n).signatures(circle_count))


@lru_cache(maxsize=None)
def conditioned_subset_coefficients(
    n: int,
    sizes: tuple[int, ...],
    known_items: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, tuple[int, ...], tuple[int, ...],
                       tuple[int, ...], tuple[int, ...]], ...]:
    """Expand every retained-subset condition into fixed linear coefficients.

    Each returned row consists of a constant and coefficients for global
    lines, global circles, local lines, and local circles.  Evaluating these
    rows is exactly the two inequalities formerly assembled by nested
    binomial sums for every candidate point row.
    """
    known_bounds = dict(known_items)
    rows = []
    for retained_size in range(4, n):
        if retained_size not in known_bounds:
            continue
        bound = known_bounds[retained_size]
        containing_global_lines = []
        containing_global_circles = []
        containing_local_lines = []
        containing_local_circles = []
        avoiding_global_lines = []
        avoiding_global_circles = []
        avoiding_local_lines = []
        avoiding_local_circles = []
        for block_size in sizes:
            containing_through = sum(
                choose(block_size - 1, j - 1)
                * choose(n - block_size, retained_size - j)
                for j in range(3, min(block_size, retained_size) + 1)
            )
            containing_avoiding = sum(
                choose(block_size, j)
                * choose(n - 1 - block_size, retained_size - 1 - j)
                for j in range(3, min(block_size, retained_size - 1) + 1)
            )
            avoiding_through = sum(
                choose(block_size - 1, j)
                * choose(n - block_size, retained_size - j)
                for j in range(3, min(block_size - 1, retained_size) + 1)
            )
            avoiding_avoiding = sum(
                choose(block_size, j)
                * choose(n - 1 - block_size, retained_size - j)
                for j in range(3, min(block_size, retained_size) + 1)
            )
            containing_bad = choose(block_size - 1, retained_size - 1)
            avoiding_bad_through = choose(block_size - 1, retained_size)
            avoiding_bad_other = choose(block_size, retained_size)

            containing_global_lines.append(0)
            containing_global_circles.append(containing_avoiding)
            containing_local_lines.append(bound * containing_bad)
            containing_local_circles.append(
                containing_through - containing_avoiding
                + (bound - 1) * containing_bad
            )
            avoiding_global_lines.append(bound * avoiding_bad_other)
            avoiding_global_circles.append(
                avoiding_avoiding + (bound - 1) * avoiding_bad_other
            )
            avoiding_local_lines.append(
                bound * (avoiding_bad_through - avoiding_bad_other)
            )
            avoiding_local_circles.append(
                avoiding_through - avoiding_avoiding
                + (bound - 1)
                * (avoiding_bad_through - avoiding_bad_other)
            )

        rows.append((
            choose(n - 1, retained_size - 1) * bound,
            tuple(containing_global_lines),
            tuple(containing_global_circles),
            tuple(containing_local_lines),
            tuple(containing_local_circles),
        ))
        rows.append((
            choose(n - 1, retained_size) * bound,
            tuple(avoiding_global_lines),
            tuple(avoiding_global_circles),
            tuple(avoiding_local_lines),
            tuple(avoiding_local_circles),
        ))
    return tuple(rows)


def conditioned_point_row_ok(
    n: int,
    sizes: tuple[int, ...],
    signature: tuple[int, ...],
    local_row: tuple[int, ...],
    known_bounds: dict[int, int],
) -> bool:
    """Test deletion, inversion, and subset bounds at one marked point.

    ``signature`` and ``local_row`` list the line families first and the
    circle families second.  The subset sums are stratified according as the
    retained subset contains the marked point or avoids it.  Subsets lying
    inside a maximal line or circle are subtracted exactly.
    """
    dimension = len(sizes)
    if len(signature) != 2 * dimension or len(local_row) != 2 * dimension:
        raise ValueError("signature and local row have incompatible lengths")
    lines = signature[:dimension]
    circles = signature[dimension:]
    local_lines = local_row[:dimension]
    local_circles = local_row[dimension:]
    if any(local_lines[index] > lines[index]
           or local_circles[index] > circles[index]
           for index in range(dimension)):
        return False
    ordinary_lines = choose(n, 2) - sum(
        choose(size, 2) * number for size, number in zip(sizes, lines)
    )
    ordinary_through_point = n - 1 - sum(
        (size - 1) * number
        for size, number in zip(sizes, local_lines)
    )
    if not 0 <= ordinary_through_point <= ordinary_lines:
        return False
    previous = known_bounds[n - 1]
    inverse_delete = known_bounds[n - 2]
    circle_count = sum(circles)
    block_count = sum(signature)
    local_block_count = sum(local_row)
    if local_circles[0] > circle_count - previous:
        return False
    if block_count - local_block_count < previous:
        return False
    global_three = lines[0] + circles[0]
    local_three = local_lines[0] + local_circles[0]
    if ((n - 1) * (block_count - local_block_count)
            - 3 * (global_three - local_three)
            < (n - 1) * inverse_delete):
        return False

    if not conditioned_line_subset_ok(n, sizes, signature, local_row):
        return False

    # For n=12,13 the maximum block size does not change after deleting or
    # inverting about the marked point.  Both operations determine the full
    # child signature, which must satisfy the already established global
    # necessary inequalities one order below.
    if n in (12, 13):
        child_n = n - 1
        child_lines = tuple(
            lines[index] - local_lines[index]
            + (local_lines[index + 1]
               if index + 1 < dimension else 0)
            for index in range(dimension)
        )
        child_circles = tuple(
            circles[index] - local_circles[index]
            + (local_circles[index + 1]
               if index + 1 < dimension else 0)
            for index in range(dimension)
        )
        child_circle_count = sum(child_circles)
        if child_lines + child_circles not in candidate_signature_set(
                child_n, child_circle_count):
            return False

        combined = tuple(
            local_lines[index] + local_circles[index]
            for index in range(dimension)
        )
        inverted_lines = tuple(
            combined[index + 1] if index + 1 < dimension else 0
            for index in range(dimension)
        )
        inverted_circles = tuple(
            lines[index] - local_lines[index]
            + circles[index] - local_circles[index]
            for index in range(dimension)
        )
        inverted_circle_count = sum(inverted_circles)
        if inverted_lines + inverted_circles not in candidate_signature_set(
                child_n, inverted_circle_count):
            return False

    coefficients = conditioned_subset_coefficients(
        n, sizes, tuple(sorted(known_bounds.items()))
    )
    for (constant, global_line_coefficients, global_circle_coefficients,
         local_line_coefficients, local_circle_coefficients) in coefficients:
        left = (
            sum(value * coefficient for value, coefficient
                in zip(lines, global_line_coefficients))
            + sum(value * coefficient for value, coefficient
                  in zip(circles, global_circle_coefficients))
            + sum(value * coefficient for value, coefficient
                  in zip(local_lines, local_line_coefficients))
            + sum(value * coefficient for value, coefficient
                  in zip(local_circles, local_circle_coefficients))
        )
        if left < constant:
            return False
    return True


def conditioned_line_subset_ok(
    n: int,
    sizes: tuple[int, ...],
    signature: tuple[int, ...],
    local_row: tuple[int, ...],
) -> bool:
    """Sum line-arrangement inequalities in the two fixed-point strata.

    For every retained subset size, subsets containing the marked point and
    subsets avoiding it are counted separately.  A maximal r-point line
    through the marked point and one avoiding it have different binomial
    multiplicities in the two strata, so this is stronger than the unmarked
    subset average while remaining a closed-form integer calculation.
    """
    dimension = len(sizes)
    lines = signature[:dimension]
    local_lines = local_row[:dimension]
    ordinary_lines = choose(n, 2) - sum(
        choose(size, 2) * number
        for size, number in zip(sizes, lines)
    )
    ordinary_through = n - 1 - sum(
        (size - 1) * number
        for size, number in zip(sizes, local_lines)
    )

    def stratified_count(
        size: int, through: bool, subset_size: int,
        contains_marked: bool, intersection: int,
    ) -> int:
        if through:
            if contains_marked:
                return (choose(size - 1, intersection - 1)
                        * choose(n - size, subset_size - intersection))
            return (choose(size - 1, intersection)
                    * choose(n - size, subset_size - intersection))
        if contains_marked:
            return (choose(size, intersection)
                    * choose(n - 1 - size,
                             subset_size - 1 - intersection))
        return (choose(size, intersection)
                * choose(n - 1 - size,
                         subset_size - intersection))

    maximum = sizes[-1]
    for subset_size in range(maximum + 1, n):
        for contains_marked in (True, False):
            ordinary = defect = bojanowski = shnurnikov = line_number = 0
            data = ((2, ordinary_through,
                     ordinary_lines - ordinary_through),) + tuple(
                (size, local_lines[index],
                 lines[index] - local_lines[index])
                for index, size in enumerate(sizes)
            )
            for size, through_number, avoiding_number in data:
                for through, number in (
                        (True, through_number), (False, avoiding_number)):
                    for intersection in range(
                            2, min(size, subset_size) + 1):
                        multiplicity = number * stratified_count(
                            size, through, subset_size, contains_marked,
                            intersection,
                        )
                        if intersection == 2:
                            ordinary += multiplicity
                        defect += multiplicity * (
                            1 if intersection == 2 else
                            (-(intersection - 3)
                             if intersection >= 4 else 0)
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
                    line_number += number * sum(
                        stratified_count(
                            size, through, subset_size, contains_marked,
                            intersection,
                        )
                        for intersection in range(
                            2, min(size, subset_size) + 1)
                    )
            subset_count = choose(
                n - 1,
                subset_size - 1 if contains_marked else subset_size,
            )
            if ordinary < ordinary_lower(subset_size) * subset_count:
                return False
            if defect < 3 * subset_count:
                return False
            if (3 * maximum < 2 * subset_size
                    and bojanowski < 4 * subset_size * subset_count):
                return False
            if (subset_size >= maximum + 3
                    and shnurnikov < 16 * subset_count):
                return False
            if (subset_size >= 10
                    and maximum <= subset_size - 2
                    and line_number < (2 * subset_size - 4) * subset_count):
                return False
    return True


class SignatureFilter:
    """Generate signatures using only the documented general inequalities."""

    def __init__(self, n: int):
        if n not in (9, 10, 11, 12, 13, 14, 15):
            raise ValueError("implemented only for n=9,...,15")
        self.n = n
        self.kmax = {
            9: 5, 10: 6, 11: 6, 12: 6, 13: 6, 14: 7, 15: 7
        }[n]
        self.sizes = tuple(range(3, self.kmax + 1))

    def line_vector_ok(self, counts: dict[int, int]) -> bool:
        n = self.n
        l2 = counts[2]
        if counts.get(3, 0) > ORCHARD_UPPER_BOUNDS[n]:
            return False
        if l2 < ordinary_lower(n):
            return False
        if l2 < 3 + sum((r - 3) * counts[r]
                        for r in range(4, self.kmax + 1)):
            return False
        if 4 * l2 + 3 * counts.get(3, 0) < (
                4 * n + 4 * sum((2 * r - 9) * counts[r]
                                for r in range(5, self.kmax + 1))):
            return False
        if 4 * l2 + 3 * counts.get(3, 0) < (
                4 * n + sum(r * (r - 4) * counts[r]
                            for r in range(5, self.kmax + 1))):
            return False
        if 2 * l2 + 3 * counts.get(3, 0) < (
                16 + sum((4 * r - 15) * counts[r]
                         for r in range(4, self.kmax + 1))):
            return False
        # Elliott's Theorem 3 applies because n>=10 and the present largest
        # line has size at most kmax<=7<n-1.  It is deliberately not used at
        # n=9, where the theorem's q>=10 hypothesis fails.
        if n >= 10 and sum(counts.values()) < 2 * n - 4:
            return False

        for subset_size in range(self.kmax + 1, n):
            subsets = choose(n, subset_size)
            ordinary = defect = bojanowski = shnurnikov = 0
            orchard_triples = 0
            line_number = 0
            for r, number in counts.items():
                for j in range(2, min(r, subset_size) + 1):
                    multiplicity = (number * choose(r, j)
                                    * choose(n - r, subset_size - j))
                    if j == 2:
                        ordinary += multiplicity
                    defect += multiplicity * (
                        1 if j == 2 else (-(j - 3) if j >= 4 else 0))
                    bojanowski += multiplicity * (
                        4 if j == 2 else
                        (3 if j == 3 else (4 * j - j * j if j >= 5 else 0)))
                    shnurnikov += multiplicity * (
                        2 if j == 2 else
                        (3 if j == 3 else (15 - 4 * j if j >= 4 else 0)))
                if r >= 3:
                    orchard_triples += (
                        number * choose(r, 3)
                        * choose(n - r, subset_size - 3)
                    )
                line_number += number * (
                    subsets - choose(n - r, subset_size)
                    - r * choose(n - r, subset_size - 1)
                )
            if ordinary < ordinary_lower(subset_size) * subsets:
                return False
            if defect < 3 * subsets:
                return False
            if (3 * self.kmax < 2 * subset_size
                    and bojanowski < 4 * subset_size * subsets):
                return False
            if subset_size >= self.kmax + 3 and shnurnikov < 16 * subsets:
                return False
            if orchard_triples > (
                    subsets * ORCHARD_UPPER_BOUNDS[subset_size]):
                return False
            # Sum Elliott's line-number bound over all retained subsets.  Its
            # hypotheses hold only from subset_size=10 onward; the maximum
            # retained collinearity is at most kmax<subset_size-1.
            if (subset_size >= 10
                    and line_number < (2 * subset_size - 4) * subsets):
                return False
        return True

    @lru_cache(maxsize=None)
    def line_vectors(self) -> tuple[tuple[int, ...], ...]:
        answer = []
        work = [0] * (self.kmax + 1)

        def rec(r: int, remaining_pairs: int) -> None:
            if r == 2:
                work[2] = remaining_pairs
                counts = {j: work[j] for j in range(2, self.kmax + 1)}
                if self.line_vector_ok(counts):
                    answer.append(tuple(work[j] for j in self.sizes))
                return
            weight = choose(r, 2)
            for number in range(remaining_pairs // weight + 1):
                work[r] = number
                rec(r - 1, remaining_pairs - weight * number)
            work[r] = 0

        rec(self.kmax, choose(self.n, 2))
        return tuple(answer)

    @lru_cache(maxsize=None)
    def circle_vectors(self, circle_count: int,
                       triple_total: int) -> tuple[tuple[int, ...], ...]:
        answer = []
        work = [0] * (self.kmax + 1)

        def rec(r: int, remaining_count: int, remaining_triples: int) -> None:
            if r == 3:
                if remaining_triples == remaining_count:
                    work[3] = remaining_count
                    answer.append(tuple(work[j] for j in self.sizes))
                return
            weight = choose(r, 3)
            upper = min(remaining_count, remaining_triples // weight)
            for number in range(upper + 1):
                work[r] = number
                rec(r - 1, remaining_count - number,
                    remaining_triples - weight * number)
            work[r] = 0

        rec(self.kmax, circle_count, triple_total)
        return tuple(answer)

    def subset_circle_ok(self, circles: tuple[int, ...],
                         lines: tuple[int, ...]) -> bool:
        for subset_size in range(4, self.n):
            lhs = 0
            for r, number in zip(self.sizes, circles):
                coefficient = sum(
                    choose(r, j) * choose(self.n - r, subset_size - j)
                    for j in range(3, min(r, subset_size) + 1))
                lhs += coefficient * number
            bad_lines = sum(choose(r, subset_size) * number
                            for r, number in zip(self.sizes, lines))
            bad_circles = sum(choose(r, subset_size) * number
                              for r, number in zip(self.sizes, circles))
            rhs = (choose(self.n, subset_size) * KNOWN[subset_size]
                   - KNOWN[subset_size] * bad_lines
                   - (KNOWN[subset_size] - 1) * bad_circles)
            if lhs < rhs:
                return False
        return True

    def signatures(self, circle_count: int) -> list[tuple[int, ...]]:
        answer = []
        for lines in self.line_vectors():
            line_triples = sum(choose(r, 3) * number
                               for r, number in zip(self.sizes, lines))
            remaining = choose(self.n, 3) - line_triples
            if remaining < 0:
                continue
            for circles in self.circle_vectors(circle_count, remaining):
                c3 = circles[0]
                if self.n * circle_count - 3 * c3 < self.n * KNOWN[self.n - 1]:
                    continue
                if 3 * c3 < self.n * ceil((self.n - 1) / 6):
                    continue
                blocks = tuple(lines[i] + circles[i]
                               for i in range(len(self.sizes)))
                if (3 * blocks[0] - 3 * self.n
                        - sum(r * (r - 4) * blocks[r - 3]
                              for r in range(5, self.kmax + 1))) < 0:
                    continue
                if (12 * blocks[0] + 12 * blocks[1]
                        - 4 * self.n * (self.n - 1)
                        - sum(4 * (2 * r - 11) * r * blocks[r - 3]
                              for r in range(6, self.kmax + 1))) < 0:
                    continue
                if (12 * blocks[0] + 12 * blocks[1]
                        - 4 * self.n * (self.n - 1)
                        - sum(r * (r - 1) * (r - 5) * blocks[r - 3]
                              for r in range(6, self.kmax + 1))) < 0:
                    continue
                if (6 * blocks[0] + 12 * blocks[1] - 16 * self.n
                        - sum(r * (4 * r - 19) * blocks[r - 3]
                              for r in range(5, self.kmax + 1))) < 0:
                    continue
                if 3 * blocks[0] < self.n * ordinary_lower(self.n - 1):
                    continue
                if not self.subset_circle_ok(circles, lines):
                    continue
                answer.append(lines + circles)
        return answer
