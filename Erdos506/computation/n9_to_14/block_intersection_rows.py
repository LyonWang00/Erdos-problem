#!/usr/bin/env python3
"""Generate anonymous per-block intersection rows.

This research helper contains no solver and no theorem claim.  For a fixed
global signature and a distinguished block category it enumerates the
possible numbers of other rich line/circle blocks meeting one block in one
or two points, subject to the three exact triple-stratum identities and the
cross-block matching inequality.
"""

from __future__ import annotations

from functools import lru_cache
from math import ceil, comb


KNOWN_CIRCLE_LOWER_BOUNDS = {
    4: 3, 5: 5, 6: 8, 7: 11, 8: 17,
    9: 25, 10: 33, 11: 41, 12: 51, 13: 61, 14: 73,
}
ORCHARD_UPPER_BOUNDS = {
    4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12,
    11: 16, 12: 19, 13: 24, 14: 27,
}


def choose(n: int, r: int) -> int:
    return comb(n, r) if 0 <= r <= n else 0


def inherited_circle_count(
    n: int,
    distinguished_size: int,
    subset_size: int,
    stratum: int,
    circle_size: int,
    intersection_size: int,
) -> int:
    """Count stratum subsets inheriting one fixed maximal circle."""
    total = 0
    for common in range(intersection_size + 1):
        for circle_only in range(circle_size - intersection_size + 1):
            if common + circle_only < 3:
                continue
            total += (
                choose(intersection_size, common)
                * choose(
                    distinguished_size - intersection_size,
                    stratum - common,
                )
                * choose(circle_size - intersection_size, circle_only)
                * choose(
                    n - distinguished_size - circle_size
                    + intersection_size,
                    subset_size - stratum - circle_only,
                )
            )
    return total


def selected_block_count(
    n: int,
    distinguished_size: int,
    subset_size: int,
    stratum: int,
    block_size: int,
    intersection_size: int,
    selected_from_block: int,
) -> int:
    """Count stratum subsets selecting exactly s points of one block."""
    total = 0
    for common in range(intersection_size + 1):
        block_only = selected_from_block - common
        total += (
            choose(intersection_size, common)
            * choose(
                distinguished_size - intersection_size,
                stratum - common,
            )
            * choose(block_size - intersection_size, block_only)
            * choose(
                n - distinguished_size - block_size + intersection_size,
                subset_size - stratum - block_only,
            )
        )
    return total


def ordinary_line_lower_bound(number_of_points: int) -> int:
    if number_of_points >= 9:
        return ceil(6 * number_of_points / 13)
    return {8: 4, 7: 3}.get(number_of_points, 2)


@lru_cache(maxsize=None)
def intersection_rows(
    n: int,
    signature: tuple[int, ...],
    gamma: int,
    limit: int = 0,
    incidence_inequalities: tuple[
        tuple[tuple[int, ...], int, int], ...
    ] = (),
    incidence_pair_sets: tuple[
        tuple[int, int, frozenset[tuple[int, int]]], ...
    ] = (),
    use_deletion_bounds: bool = False,
    use_subset_inheritance_bounds: bool = False,
    use_conditioned_line_bounds: bool = False,
) -> tuple[tuple[int, ...], ...]:
    dimension = len(signature) // 2
    sizes = tuple(range(3, 3 + dimension))
    category_sizes = sizes + sizes
    a = category_sizes[gamma]
    available = tuple(
        number - int(alpha == gamma)
        for alpha, number in enumerate(signature)
    )

    # A 3-circle category can realize all intersection sizes and its three
    # stratum contributions are simply (x_2,x_1,x_0).  Leave one such
    # category to the end, where its values are forced by the residuals.
    pivot = dimension if available[dimension] >= 0 else None
    if pivot is None:
        raise AssertionError("the signature has no 3-circle coordinate")
    order = tuple(
        sorted(
            (alpha for alpha in range(2 * dimension) if alpha != pivot),
            key=lambda alpha: (
                category_sizes[alpha], alpha >= dimension, alpha
            ),
            reverse=True,
        )
    )
    target = (
        choose(a, 2) * (n - a),
        a * choose(n - a, 2),
        choose(n - a, 3),
    )
    matching_bound = (a // 2) * choose(n - a, 2)
    maximum_block_size = max(
        (size for size, number in zip(category_sizes, signature) if number),
        default=2,
    )
    circle_count = sum(signature[dimension:])
    one = [0] * (2 * dimension)
    two = [0] * (2 * dimension)
    rows: list[tuple[int, ...]] = []

    def passes_point_incidence_bounds() -> bool:
        incidences = tuple(
            one[alpha] + 2 * two[alpha]
            for alpha in range(2 * dimension)
        )
        if not all(
            lower <= sum(
                coefficient * value
                for coefficient, value in zip(coefficients, incidences)
            ) <= upper
            for coefficients, lower, upper in incidence_inequalities
        ):
            return False
        return all(
            (incidences[alpha], incidences[beta]) in allowed
            for alpha, beta, allowed in incidence_pair_sets
        )

    def passes_deletion_bounds() -> bool:
        """Check the circle loss after deleting subsets of this one block."""
        for deleted in range(1, a + 1):
            remainder = n - deleted
            if (remainder <= maximum_block_size
                    or remainder not in KNOWN_CIRCLE_LOWER_BOUNDS):
                continue
            loss = (
                choose(a - 1, deleted - 1) * one[dimension]
                + (choose(a, deleted) - choose(a - 2, deleted))
                * two[dimension]
            )
            if dimension >= 2:
                loss += choose(a - 2, deleted - 2) * two[dimension + 1]
            if gamma >= dimension and a - deleted < 3:
                loss += choose(a, deleted)
            if loss > (
                choose(a, deleted)
                * (circle_count - KNOWN_CIRCLE_LOWER_BOUNDS[remainder])
            ):
                return False
        return True

    def passes_subset_inheritance_bounds() -> bool:
        """Sum the known circle lower bound over every block stratum."""
        for subset_size in range(maximum_block_size + 1, n):
            lower_bound = KNOWN_CIRCLE_LOWER_BOUNDS.get(subset_size)
            if lower_bound is None:
                continue
            for stratum in range(a + 1):
                subset_count = (
                    choose(a, stratum)
                    * choose(n - a, subset_size - stratum)
                )
                if subset_count == 0:
                    continue
                inherited = 0
                for alpha in range(dimension, 2 * dimension):
                    b = category_sizes[alpha]
                    zero = available[alpha] - one[alpha] - two[alpha]
                    inherited += (
                        inherited_circle_count(
                            n, a, subset_size, stratum, b, 0
                        ) * zero
                        + inherited_circle_count(
                            n, a, subset_size, stratum, b, 1
                        ) * one[alpha]
                        + inherited_circle_count(
                            n, a, subset_size, stratum, b, 2
                        ) * two[alpha]
                    )
                if gamma >= dimension:
                    inherited += inherited_circle_count(
                        n, a, subset_size, stratum, a, a
                    )
                if inherited < subset_count * lower_bound:
                    return False
        return True

    def passes_conditioned_line_bounds() -> bool:
        """Sum classical connecting-line bounds over every block stratum."""
        for subset_size in range(maximum_block_size + 1, n + 1):
            for stratum in range(a + 1):
                subset_count = (
                    choose(a, stratum)
                    * choose(n - a, subset_size - stratum)
                )
                if subset_count == 0:
                    continue
                multiplicities = {
                    selected: 0
                    for selected in range(3, maximum_block_size + 1)
                }
                for alpha in range(dimension):
                    if available[alpha] == 0:
                        continue
                    b = category_sizes[alpha]
                    zero = available[alpha] - one[alpha] - two[alpha]
                    for selected in range(3, b + 1):
                        multiplicities[selected] += (
                            selected_block_count(
                                n, a, subset_size, stratum, b, 0, selected
                            ) * zero
                            + selected_block_count(
                                n, a, subset_size, stratum, b, 1, selected
                            ) * one[alpha]
                            + selected_block_count(
                                n, a, subset_size, stratum, b, 2, selected
                            ) * two[alpha]
                        )
                if gamma < dimension:
                    for selected in range(3, a + 1):
                        multiplicities[selected] += selected_block_count(
                            n, a, subset_size, stratum, a, a, selected
                        )
                ordinary = (
                    subset_count * choose(subset_size, 2)
                    - sum(
                        (choose(selected, 2) - 1) * number
                        for selected, number in multiplicities.items()
                    )
                )
                if ordinary < (
                    subset_count
                    * ordinary_line_lower_bound(subset_size)
                ):
                    return False
                if ordinary < (
                    3 * subset_count
                    + sum(
                        (selected - 3) * number
                        for selected, number in multiplicities.items()
                        if selected >= 4
                    )
                ):
                    return False
                if multiplicities.get(3, 0) > (
                    subset_count * ORCHARD_UPPER_BOUNDS[subset_size]
                ):
                    return False
                if (subset_size >= 10
                        and maximum_block_size <= subset_size - 2
                        and ordinary + sum(multiplicities.values())
                        < (2 * subset_size - 4) * subset_count):
                    return False
        return True

    def recurse(
        position: int,
        used_two_stratum: int,
        used_one_stratum: int,
        used_zero_stratum: int,
        used_matching: int,
    ) -> bool:
        if position == len(order):
            residual_two = target[0] - used_two_stratum
            residual_one = target[1] - used_one_stratum
            residual_zero = target[2] - used_zero_stratum
            pivot_available = available[pivot]
            if not (
                0 <= residual_two <= pivot_available
                and 0 <= residual_one <= pivot_available - residual_two
                and residual_zero
                == pivot_available - residual_two - residual_one
            ):
                return False
            # The pivot is the 3-circle category.  Its matching coefficient
            # is C(3-2,2)=0, so its two-point intersections consume no
            # outside-pair matching capacity.
            if used_matching > matching_bound:
                return False
            one[pivot] = residual_one
            two[pivot] = residual_two
            if not passes_point_incidence_bounds():
                return False
            if use_deletion_bounds and not passes_deletion_bounds():
                return False
            if (use_subset_inheritance_bounds
                    and not passes_subset_inheritance_bounds()):
                return False
            if (use_conditioned_line_bounds
                    and not passes_conditioned_line_bounds()):
                return False
            rows.append(tuple(
                value for pair in zip(one, two) for value in pair
            ))
            return bool(limit and len(rows) >= limit)

        alpha = order[position]
        b = category_sizes[alpha]
        count = available[alpha]
        maximum_two = (
            0 if gamma < dimension and alpha < dimension else count
        )
        for intersection_two in range(maximum_two + 1):
            contribution_two = (b - 2) * intersection_two
            next_two = used_two_stratum + contribution_two
            if next_two > target[0]:
                break
            matching = choose(b - 2, 2) * intersection_two
            if used_matching + matching > matching_bound:
                break
            for intersection_one in range(count - intersection_two + 1):
                intersection_zero = (
                    count - intersection_one - intersection_two
                )
                next_one = used_one_stratum + (
                    choose(b - 1, 2) * intersection_one
                    + 2 * choose(b - 2, 2) * intersection_two
                )
                if next_one > target[1]:
                    break
                next_zero = used_zero_stratum + (
                    choose(b, 3) * intersection_zero
                    + choose(b - 1, 3) * intersection_one
                    + choose(b - 2, 3) * intersection_two
                )
                if next_zero > target[2]:
                    continue
                one[alpha] = intersection_one
                two[alpha] = intersection_two
                if recurse(
                    position + 1,
                    next_two,
                    next_one,
                    next_zero,
                    used_matching + matching,
                ):
                    return True
        one[alpha] = 0
        two[alpha] = 0
        return False

    recurse(0, 0, 0, 0, 0)
    return tuple(rows)
