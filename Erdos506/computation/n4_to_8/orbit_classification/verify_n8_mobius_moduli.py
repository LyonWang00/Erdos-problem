#!/usr/bin/env python3
"""Exact verification of the complete n=8 generalized-circle Möbius moduli.

The script verifies the rectangular normal form for the realizable generalized-
circle incidence type, its 24 rational relabelling actions, the concurrence
curve giving c(P)=17, and the complete exceptional identifications on that
curve.  All computations are symbolic over Q.
"""

from __future__ import annotations

from itertools import combinations, permutations

import sympy as sp

from family_examples import line_intersection, n8_c17_family, n8_c18_rectangles
from orbit_tools import (
    incidence_summary,
    large_generalized_blocks,
    maximal_incidence,
)
from realize_c18_by_recoloring import circle_equation
from verify_mobius_moduli import (
    cross_ratio,
    exceptional_factors,
    normalized_equal,
)


x, y, lam = sp.symbols("x y lam")
c = 1 - y + x * y
d = x + y - x * y

# The twelve four-point generalized circles in the rectangular labelling.
# The value is the ordered cross ratio in increasing index order.
FOUR_BLOCK_VALUE = {
    (0, 3, 4, 7): x,
    (1, 2, 5, 6): x,
    (0, 1, 2, 3): y,
    (4, 5, 6, 7): y,
    (0, 1, 4, 5): x / c,
    (2, 3, 6, 7): x / c,
    (0, 1, 6, 7): d,
    (2, 3, 4, 5): d,
    (0, 2, 4, 6): x / d,
    (1, 3, 5, 7): x / d,
    (0, 2, 5, 7): c,
    (1, 3, 4, 6): c,
}
FOUR_BLOCKS = frozenset(FOUR_BLOCK_VALUE)


def action_key(action):
    return tuple(sp.sstr(sp.factor(sp.cancel(entry))) for entry in action)


def cr_sphere(values):
    """Cross ratio of four values, with None denoting infinity."""
    a, b, c0, d0 = values
    if a is None:
        return sp.cancel((b - d0) / (b - c0))
    if b is None:
        return sp.cancel((a - c0) / (a - d0))
    if c0 is None:
        return sp.cancel((b - d0) / (a - d0))
    if d0 is None:
        return sp.cancel((a - c0) / (b - c0))
    return sp.cancel((a - c0) * (b - d0) / ((a - d0) * (b - c0)))


CANONICAL = (sp.Integer(0), sp.Integer(1), lam / (lam - 1), None)
ORDER_ACTION = {
    order: sp.factor(cr_sphere([CANONICAL[index] for index in order]))
    for order in permutations(range(4))
}


def ordered_block_value(indices):
    block = tuple(sorted(indices))
    order = tuple(block.index(index) for index in indices)
    return sp.factor(sp.cancel(
        ORDER_ACTION[order].subs(lam, FOUR_BLOCK_VALUE[block])
    ))


def incidence_automorphisms():
    out = []
    for perm in permutations(range(8)):
        image = {
            tuple(sorted(perm[index] for index in block))
            for block in FOUR_BLOCKS
        }
        if image == FOUR_BLOCKS:
            out.append(perm)
    return out


def induced_action(perm):
    return (
        ordered_block_value(tuple(perm[index] for index in (0, 3, 4, 7))),
        ordered_block_value(tuple(perm[index] for index in (0, 1, 2, 3))),
    )


def compose_actions(outer, inner):
    return tuple(sp.factor(sp.cancel(entry.subs(
        {x: inner[0], y: inner[1]}, simultaneous=True
    ))) for entry in outer)


def generated_actions(generators):
    identity = (x, y)
    group = {action_key(identity): identity}
    while True:
        enlarged = dict(group)
        for action in group.values():
            for generator in generators:
                composition = compose_actions(generator, action)
                enlarged[action_key(composition)] = composition
        if enlarged.keys() == group.keys():
            return group
        group = enlarged


def verify_rectangular_normal_form():
    aspect, ratio = sp.symbols("aspect ratio", real=True)
    points = n8_c18_rectangles(aspect, ratio)
    line_triples = {
        tuple(block)
        for four_block in ((0, 3, 4, 7), (1, 2, 5, 6))
        for block in combinations(four_block, 3)
    }
    allowed_four_blocks = set(FOUR_BLOCKS)
    factors = exceptional_factors(
        points, (aspect, ratio), line_triples, allowed_four_blocks
    )
    expected_factors = {
        aspect,
        ratio,
        aspect**2 + 1,
        ratio - 1,
        ratio + 1,
        aspect**2 * (ratio - 1)**2 + (ratio + 1)**2,
        aspect**2 * (ratio + 1)**2 + (ratio - 1)**2,
    }
    assert factors == {sp.factor(value) for value in expected_factors}

    z = [sp.factor(px + sp.I * py) for px, py in points]
    t = sp.symbols("t", real=True)
    substitutions = {
        ratio: (1 + t) / (1 - t),
        aspect**2: y / (1 - y),
    }
    expected = {
        (0, 3, 4, 7): x,
        (1, 2, 5, 6): x,
        (0, 1, 2, 3): y,
        (4, 5, 6, 7): y,
        (0, 1, 4, 5): x / c,
        (2, 3, 6, 7): x / c,
        (0, 1, 6, 7): d,
        (2, 3, 4, 5): d,
        (0, 2, 4, 6): x / d,
        (1, 3, 5, 7): x / d,
        (0, 2, 5, 7): c,
        (1, 3, 4, 6): c,
    }
    for block, target in expected.items():
        value = cross_ratio(z, block)
        value = sp.factor(sp.cancel(value.subs(substitutions)))
        value = value.subs(t**2, x)
        assert sp.factor(value - target) == 0, (block, value, target)

    # Reflection in the x-axis is an anti-Möbius automorphism of every member.
    assert normalized_equal(
        z, z, (1, 0, 3, 2, 5, 4, 7, 6), anti=True
    )


def verify_complete_normalization():
    """Solve the full generalized-circle realization after sending p0 to infinity."""
    a, b, t, k, q = sp.symbols("a b t k q", real=True)
    points = [None] * 8
    points[1] = (sp.Integer(0), sp.Integer(0))
    points[2] = (sp.Integer(1), sp.Integer(0))
    points[4] = (a, b)
    points[5] = (t * a, t * b)
    points[6] = (1 + k * (a - 1), k * b)
    points[7] = line_intersection(
        points[1], points[6], points[2], points[5]
    )
    points[3] = line_intersection(
        points[1], points[2], points[4], points[7]
    )

    circle_blocks = (
        (1, 2, 5, 6), (1, 3, 4, 6), (1, 3, 5, 7),
        (2, 3, 4, 5), (2, 3, 6, 7), (4, 5, 6, 7),
    )
    residuals = []
    for block in circle_blocks:
        determinant = sp.factor(sp.cancel(sp.Matrix([
            [
                points[index][0] ** 2 + points[index][1] ** 2,
                points[index][0], points[index][1], 1,
            ]
            for index in block
        ]).det(method="domain-ge")))
        numerator = sp.factor(determinant.as_numer_denom()[0])
        candidates = [
            factor.subs(b**2, q - a**2)
            for factor, _multiplicity in sp.factor_list(numerator)[1]
            if factor.has(a) and factor.has(k) and factor.has(t)
        ]
        assert len(candidates) == 1
        residuals.append(sp.factor(candidates[0]))

    denominator = 2 * k * t - k - t
    solution = {
        a: t * (k - 1) / denominator,
        q: (k - 1) / denominator,
    }
    assert all(sp.factor(value.subs(solution)) == 0 for value in residuals)
    assert sp.factor(points[3][0] - solution[a]) == 0

    parameter_x = -k * t / (k * t - k - t)
    parameter_y = solution[a]
    inverse = {
        k: x / (x + y - x * y),
        t: x / (1 - y + x * y),
    }
    assert sp.factor(parameter_x.subs(inverse) - x) == 0
    assert sp.factor(parameter_y.subs(inverse) - y) == 0
    squared_height = sp.factor((solution[q] - solution[a] ** 2).subs(inverse))
    assert sp.factor(squared_height - y * (1 - y) / x) == 0


def verify_relabelling_group():
    automorphisms = incidence_automorphisms()
    assert len(automorphisms) == 192
    actual = {
        action_key(induced_action(perm)): induced_action(perm)
        for perm in automorphisms
    }
    generators = (
        (1 / x, y),
        (x, 1 - y),
        ((y - 1) / y, x / (x - 1)),
        (y * (1 - x), x / (1 - y + x * y)),
    )
    generated = generated_actions(generators)
    assert len(actual) == 24
    assert set(actual) == set(generated)
    return actual


def verify_seventeen_circle_locus(actions):
    aspect, ratio = sp.symbols("aspect ratio", real=True)
    points = n8_c18_rectangles(aspect, ratio)
    f_0257 = circle_equation(points, (0, 2, 5, 7))
    f_036 = circle_equation(points, (0, 3, 6))
    f_147 = circle_equation(points, (1, 4, 7))
    coordinate_x, coordinate_y = sp.symbols("x y", real=True)
    common = {
        coordinate_x: 4 * aspect * ratio / (ratio - 1),
        coordinate_y: 2 * ratio / (ratio - 1),
    }
    assert sp.factor((f_0257 - f_036).subs(common)) == 0
    assert sp.factor((f_0257 - f_147).subs(common)) == 0
    concurrence = sp.factor(sp.cancel(f_0257.subs(common)))
    expected = ratio * (
        aspect**2 * ratio**2 + 14 * aspect**2 * ratio
        + aspect**2 + ratio**2 + 2 * ratio + 1
    ) / (ratio - 1)**2
    assert sp.factor(concurrence - expected) == 0

    # In cross-ratio coordinates the concurrence equation is 3xy-3y-1=0.
    t = sp.symbols("t", real=True)
    transformed = sp.factor(sp.cancel(
        expected.subs({
            ratio: (1 + t) / (1 - t),
            aspect**2: y / (1 - y),
        }, simultaneous=True)
    ))
    numerator = sp.factor(transformed.as_numer_denom()[0])
    assert sp.factor(numerator / (1 + t)) == -(3 * t**2 * y - 3 * y - 1)

    # The published rational family parametrizes the whole real locus.
    b, u = sp.symbols("b u", real=True)
    c17_points = n8_c17_family(b)
    z = [sp.factor(px + sp.I * py) for px, py in c17_points]

    # This explicit relabelling identifies the generalized-circle incidence
    # structure of the c=17 family with the rectangular c=18 structure.  Since
    # every triple not contained in a four-block is itself a maximal
    # three-block, checking the four-blocks checks the full structure.
    rectangle_points = n8_c18_rectangles(sp.Integer(2), sp.Integer(3))
    c17_sample = n8_c17_family(sp.Rational(13, 12))
    rectangle_to_c17 = (0, 7, 1, 6, 5, 2, 4, 3)
    rectangle_blocks = large_generalized_blocks(rectangle_points)
    c17_blocks = large_generalized_blocks(c17_sample)
    assert {
        frozenset(rectangle_to_c17[index] for index in block)
        for block in rectangle_blocks
    } == c17_blocks

    # There are no unlisted real incidences anywhere in the asserted domain.
    # The factors b^2+1 and b^2+9 have no real zeros; b=0 is a pole of the
    # coordinate formula and b^2=3 is the only real exceptional-incidence
    # factor.
    c17_lines, c17_circles = maximal_incidence(c17_sample)
    allowed_line_triples = {
        tuple(block)
        for line in c17_lines
        for block in combinations(sorted(line), 3)
    }
    allowed_four_blocks = {
        tuple(block)
        for block in c17_lines | c17_circles
        if len(block) == 4
    }
    assert exceptional_factors(
        c17_points, (b,), allowed_line_triples, allowed_four_blocks
    ) == {b**2 - 3, b**2 + 1, b**2 + 9}
    coordinate_denominators = {
        sp.factor(sp.cancel(coordinate).as_numer_denom()[1])
        for point in c17_points
        for coordinate in point
    }
    assert any(sp.factor(denominator).subs(b, 0) == 0
               for denominator in coordinate_denominators)

    cross_x = cross_ratio(z, (0, 6, 5, 3))
    cross_y = cross_ratio(z, (0, 7, 1, 6))
    assert sp.factor(cross_x + 4 / (b**2 - 3)) == 0
    assert sp.factor(
        cross_y + (b**2 - 3) / (3 * (b**2 + 1))
    ) == 0
    assert sp.factor(3 * cross_x * cross_y - 3 * cross_y - 1) == 0

    # b and -b are directly Möbius equivalent; conjugation also interchanges
    # them.  Consequently direct and extended orbit partitions coincide.
    minus_points = n8_c17_family(-b)
    minus_z = [sp.factor(px + sp.I * py) for px, py in minus_points]
    sign_perm = (1, 0, 3, 2, 5, 4, 7, 6)
    assert normalized_equal(z, minus_z, sign_perm)
    assert normalized_equal(z, minus_z, tuple(range(8)), anti=True)

    x_of_u = 4 / (3 - u)
    y_of_u = (3 - u) / (3 * (1 + u))
    global_maps = set()
    finite_pairs = set()
    curve = 3 * x * y - 3 * y - 1
    for action in actions.values():
        action_x = sp.factor(sp.cancel(action[0].subs(
            {x: x_of_u, y: y_of_u}, simultaneous=True
        )))
        action_y = sp.factor(sp.cancel(action[1].subs(
            {x: x_of_u, y: y_of_u}, simultaneous=True
        )))
        residual = sp.factor(sp.together(curve.subs(
            {x: action_x, y: action_y}, simultaneous=True
        )))
        target_u = sp.factor(3 - 4 / action_x)
        numerator = residual.as_numer_denom()[0]
        if numerator == 0:
            global_maps.add(sp.sstr(target_u))
            continue
        for root in sp.solve(numerator, u):
            if root.is_rational is not True or root <= 0 or root == 3:
                continue
            target = sp.factor(target_u.subs(u, root))
            if target.is_rational is True and target > 0 and target != 3:
                finite_pairs.add((root, target))

    assert global_maps == {"u", "9/u"}
    expected_finite_pairs = {
        (sp.Rational(7), sp.Rational(7)),
        (sp.Rational(15), sp.Rational(15)),
        (sp.Rational(39), sp.Rational(13, 3)),
        (sp.Rational(13, 3), sp.Rational(39)),
        (sp.Rational(3, 5), sp.Rational(3, 5)),
        (sp.Rational(15), sp.Rational(3, 5)),
        (sp.Rational(13, 3), sp.Rational(3, 13)),
        (sp.Rational(27, 13), sp.Rational(3, 13)),
        (sp.Rational(3, 5), sp.Rational(15)),
        (sp.Rational(3, 13), sp.Rational(13, 3)),
        (sp.Rational(9, 7), sp.Rational(7)),
        (sp.Rational(27, 13), sp.Rational(39)),
        (sp.Rational(3, 13), sp.Rational(27, 13)),
        (sp.Rational(7), sp.Rational(9, 7)),
        (sp.Rational(9, 7), sp.Rational(9, 7)),
        (sp.Rational(39), sp.Rational(27, 13)),
    }
    assert finite_pairs == expected_finite_pairs

    # Apart from u -> 9/u, the only enlargement of an orbit is the class
    # {3/13, 27/13, 13/3, 39}.
    special = {
        sp.Rational(3, 13), sp.Rational(27, 13),
        sp.Rational(13, 3), sp.Rational(39),
    }
    for source, target in finite_pairs:
        assert (
            source == target
            or target == sp.factor(9 / source)
            or {source, target} <= special
        )


def verify_exact_representatives():
    representatives = (
        n8_c18_rectangles(sp.Integer(1), sp.Integer(2)),
        [
            (sp.Rational(-39, 10), sp.Rational(-39, 10)),
            (sp.Rational(-131, 34), sp.Rational(-133, 34)),
            (sp.Rational(-133, 34), sp.Rational(-131, 34)),
            (sp.Rational(-23, 6), sp.Rational(-23, 6)),
            (sp.Rational(-47, 12), sp.Rational(-47, 12)),
            (sp.Rational(-77, 20), sp.Rational(-79, 20)),
            (sp.Rational(-79, 20), sp.Rational(-77, 20)),
            (sp.Rational(-15, 4), sp.Rational(-15, 4)),
        ],
        [
            (sp.Rational(13, 120), sp.Rational(13, 40)),
            (sp.Rational(71, 104), sp.Rational(-43, 104)),
            (sp.Rational(199, 104), sp.Rational(85, 104)),
            (sp.Rational(47, 40), sp.Rational(167, 120)),
            (sp.Rational(47, 40), sp.Rational(141, 40)),
            (sp.Rational(35, 8), sp.Rational(59, 24)),
            (sp.Rational(-23, 24), sp.Rational(-23, 8)),
            (sp.Rational(-81, 40), sp.Rational(13, 40)),
        ],
    )
    expected_summaries = (
        {
            "lines": ["0347", "1256"],
            "circles": [
                "035", "036", "056", "124", "127", "147", "247", "356",
                "0123", "0145", "0167", "0246", "0257", "1346", "1357",
                "2345", "2367", "4567",
            ],
        },
        {
            "lines": ["056", "0347"],
            "circles": [
                "035", "036", "124", "127", "147", "247", "356", "0123",
                "0145", "0167", "0246", "0257", "1256", "1346", "1357",
                "2345", "2367", "4567",
            ],
        },
        {
            "lines": ["046", "357"],
            "circles": [
                "024", "026", "135", "137", "157", "246", "0123", "0145",
                "0167", "0257", "0347", "0356", "1247", "1256", "1346",
                "2345", "2367", "4567",
            ],
        },
    )
    for points, expected in zip(representatives, expected_summaries):
        summary = incidence_summary(points)
        assert summary["circle_count"] == 18
        assert summary["lines"] == expected["lines"]
        assert summary["circles"] == expected["circles"]


def main():
    verify_rectangular_normal_form()
    verify_complete_normalization()
    actions = verify_relabelling_group()
    verify_seventeen_circle_locus(actions)
    verify_exact_representatives()
    print("N8_GENERALIZED_MOBIUS_MODULI_PASS")
    print("GENERALIZED_INCIDENCE_AUTOMORPHISMS", 192)
    print("RATIONAL_PARAMETER_ACTIONS", 24)
    print("C17_PARAMETER_DOMAIN", "u>0, u!=3")
    print("C17_GENERIC_IDENTIFICATION", "u~9/u")
    print("C17_SPECIAL_ORBIT", "{3/13,27/13,13/3,39}")


if __name__ == "__main__":
    main()
