"""Symbolically verify the complete n=6 and n=7 parameter spaces and orbit quotients."""

from __future__ import annotations

from itertools import combinations, permutations

import sympy as sp

from family_examples import n6_family, n7_family


def cross_ratio(z, block):
    i, j, k, ell = block
    return sp.factor(
        sp.cancel((z[i] - z[k]) * (z[j] - z[ell]) /
                  ((z[i] - z[ell]) * (z[j] - z[k])))
    )


def compose(first, second):
    return tuple(first[second[i]] for i in range(len(first)))


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    while True:
        enlarged = group | {
            compose(first, second)
            for first in group
            for second in generators
        }
        if enlarged == group:
            return group
        group = enlarged


def normalized_equal(source_z, target_z, perm, anti=False):
    """Exact Möbius/anti-Möbius comparison after normalizing points 0,1,2."""
    if anti:
        source_z = [sp.conjugate(value) for value in source_z]
    za, zb, zc = source_z[0], source_z[1], source_z[2]
    wa, wb, wc = target_z[perm[0]], target_z[perm[1]], target_z[perm[2]]
    for index, value in enumerate(source_z):
        target = target_z[perm[index]]
        source_num = (value - za) * (zb - zc)
        source_den = (value - zc) * (zb - za)
        target_num = (target - wa) * (wb - wc)
        target_den = (target - wc) * (wb - wa)
        if sp.simplify(source_num * target_den - target_num * source_den) != 0:
            return False
    return True


def block_automorphisms(n, blocks):
    block_set = {frozenset(block) for block in blocks}
    return [
        perm for perm in permutations(range(n))
        if {frozenset(perm[i] for i in block) for block in block_set} == block_set
    ]


def coloured_automorphisms(n, lines, circles):
    line_set = {frozenset(block) for block in lines}
    circle_set = {frozenset(block) for block in circles}
    return [
        perm for perm in permutations(range(n))
        if {frozenset(perm[i] for i in block) for block in line_set} == line_set
        and {frozenset(perm[i] for i in block) for block in circle_set} == circle_set
    ]


def exceptional_factors(points, variables, allowed_lines, allowed_circles):
    """Irreducible factors whose vanishing changes the prescribed incidence table."""
    expressions = []
    for i, j in combinations(range(len(points)), 2):
        dx = points[i][0] - points[j][0]
        dy = points[i][1] - points[j][1]
        expressions.append(dx * dx + dy * dy)
    for block in combinations(range(len(points)), 3):
        if block not in allowed_lines:
            expressions.append(sp.Matrix([[*points[i], 1] for i in block]).det(method="domain-ge"))
    for block in combinations(range(len(points)), 4):
        if block not in allowed_circles:
            expressions.append(sp.Matrix([
                [points[i][0] ** 2 + points[i][1] ** 2, *points[i], 1]
                for i in block
            ]).det(method="domain-ge"))
    factors = set()
    for expression in expressions:
        numerator = sp.cancel(expression).as_numer_denom()[0]
        assert numerator != 0
        for factor, _multiplicity in sp.factor_list(numerator)[1]:
            if not any(factor.has(variable) for variable in variables):
                continue
            polynomial = sp.Poly(factor, *variables)
            if polynomial.LC() < 0:
                factor = -factor
            factors.add(sp.factor(factor))
    return factors


def verify_n6():
    a, v, w, x, y = sp.symbols("a v w x y", real=True)
    points = n6_family(a, v)
    lines = {(0, 2, 4), (0, 3, 5), (1, 3, 4)}
    circles = {(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)}

    factors = exceptional_factors(points, (a, v), lines, circles)
    assert factors == {
        a, a - 1, v, a - v ** 2 - 1,
        v ** 2 + 1, a ** 2 - 2 * a + v ** 2 + 1,
    }

    z = [sp.factor(px + sp.I * py) for px, py in points]
    blocks = [(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)]
    lambdas = [cross_ratio(z, block) for block in blocks]
    expected = [
        (a - v ** 2 - 1) / ((a - 1) * (v ** 2 + 1)),
        (a - v ** 2 - 1) / (a - 1),
        (a - 1) * (a - v ** 2 - 1) / ((a - 1) ** 2 + v ** 2),
    ]
    assert all(sp.factor(first - second) == 0 for first, second in zip(lambdas, expected))
    surface = (y - x) / (x * y - 2 * x + 1)
    assert sp.factor(
        expected[2] - surface.subs({
            x: expected[0], y: expected[1],
        }, simultaneous=True)
    ) == 0
    inverse_w = (y - x) / x
    inverse_a = y * (1 - x) / (x * (1 - y))
    assert sp.factor(inverse_w.subs({x: expected[0], y: expected[1]}, simultaneous=True) - v ** 2) == 0
    assert sp.factor(inverse_a.subs({x: expected[0], y: expected[1]}, simultaneous=True) - a) == 0

    automorphisms = block_automorphisms(6, blocks)
    assert len(automorphisms) == 48

    def to_xy(expression):
        return sp.factor(
            expression.subs(v ** 2, w).subs({
                a: y * (1 - x) / (x * (1 - y)),
                w: (y - x) / x,
            }, simultaneous=True)
        )

    compatible = []
    excluded_residual_factors = set()
    induced_actions = set()
    for perm in automorphisms:
        transformed = [
            to_xy(cross_ratio(z, tuple(perm[i] for i in block)))
            for block in blocks
        ]
        residual = sp.factor(sp.together(
            transformed[2] - surface.subs({x: transformed[0], y: transformed[1]}, simultaneous=True)
        ))
        if residual == 0:
            compatible.append(perm)
            induced_actions.add((sp.sstr(transformed[0]), sp.sstr(transformed[1])))
        else:
            numerator = sp.factor(sp.fraction(residual)[0])
            for factor, _multiplicity in sp.factor_list(numerator)[1]:
                if factor.has(x) or factor.has(y):
                    polynomial = sp.Poly(factor, x, y)
                    if polynomial.LC() < 0:
                        factor = -factor
                    excluded_residual_factors.add(sp.factor(factor))

    # The other 36 relabellings can return to the c=8 surface only on its
    # excluded boundary x=0, x=1, y=1, or x=y.
    assert len(compatible) == 12
    assert len(induced_actions) == 6
    assert excluded_residual_factors <= {x, x - 1, y - 1, x - y}

    # Two rational generators give all six induced actions.
    generator_actions = (
        (1 / y, 1 / x),
        (x, (x * y - 2 * x + 1) / (y - x)),
    )

    def action_key(action):
        return tuple(sp.sstr(sp.factor(sp.cancel(entry))) for entry in action)

    action_group = {(sp.sstr(x), sp.sstr(y)): (x, y)}
    while True:
        enlarged = dict(action_group)
        for action in action_group.values():
            for generator in generator_actions:
                composed = tuple(sp.factor(sp.cancel(
                    generator[index].subs({x: action[0], y: action[1]}, simultaneous=True)
                )) for index in range(2))
                enlarged[action_key(composed)] = composed
        if enlarged.keys() == action_group.keys():
            break
        action_group = enlarged
    assert set(action_group) == induced_actions
    assert len(action_group) == 6

    # The generators are realized by exact direct Möbius maps, not merely by
    # cross-ratio coincidences.
    target_first = n6_family((1 + v ** 2) / a, v)
    target_second = n6_family(
        ((a - 1) ** 2 + v ** 2) / (1 + v ** 2 - a),
        -a * v / (1 + v ** 2 - a),
    )
    target_first_z = [sp.factor(px + sp.I * py) for px, py in target_first]
    target_second_z = [sp.factor(px + sp.I * py) for px, py in target_second]
    assert normalized_equal(z, target_first_z, (1, 0, 4, 5, 2, 3))
    assert normalized_equal(z, target_second_z, (3, 2, 1, 0, 4, 5))

    # Every member has a generic anti-Möbius involution swapping the three
    # indicated point pairs, so direct and extended unlabelled orbit partitions agree.
    assert normalized_equal(z, z, (1, 0, 3, 2, 5, 4), anti=True)

    print("N6_PARAMETER_DOMAIN", "a!=0,1; v!=0; a!=1+v^2")
    print("N6_CROSS_RATIO_DOMAIN", "x,y!=0,1; (y-x)/x>0")
    print("N6_ORBIT_ACTION_COUNT", len(induced_actions))


def verify_n7():
    u, v, x, y = sp.symbols("u v x y", real=True)
    points = n7_family(u, v)
    lines = {
        (0, 3, 6), (0, 4, 5), (1, 3, 5),
        (1, 4, 6), (2, 3, 4), (2, 5, 6),
    }
    circles = {
        (0, 1, 3, 4), (0, 1, 5, 6), (0, 2, 3, 5),
        (0, 2, 4, 6), (1, 2, 3, 6), (1, 2, 4, 5),
    }
    factors = exceptional_factors(points, (u, v), lines, circles)
    assert factors == {
        u, u - 1, v, u ** 2 - u + v ** 2,
        u ** 2 + v ** 2, u ** 2 - 2 * u + v ** 2 + 1,
    }

    z = [sp.factor(px + sp.I * py) for px, py in points]
    ordered_circles = sorted(circles)
    actual = [cross_ratio(z, block) for block in ordered_circles]
    delta = u ** 2 - u + v ** 2
    expected_by_block = {
        (0, 1, 3, 4): -v ** 2 / (u * (u - 1)),
        (0, 1, 5, 6): u / (u - 1),
        (0, 2, 3, 5): v ** 2 / (u * delta),
        (0, 2, 4, 6): -u / delta,
        (1, 2, 3, 6): (u - 1) / delta,
        (1, 2, 4, 5): -v ** 2 / ((u - 1) * delta),
    }
    assert all(
        sp.factor(value - expected_by_block[block]) == 0
        for block, value in zip(ordered_circles, actual)
    )

    automorphisms = coloured_automorphisms(7, lines, circles)
    assert len(automorphisms) == 24
    assert set(block_automorphisms(7, circles)) == set(automorphisms)
    generator_first = (0, 2, 1, 3, 5, 4, 6)
    generator_second = (1, 0, 2, 5, 6, 4, 3)
    assert generated_group((generator_first, generator_second)) == set(automorphisms)

    print("N7_PARAMETER_DOMAIN", "u!=0,1; v!=0; u^2-u+v^2!=0")
    print("N7_DIRECT_ORBIT_GROUP_ORDER", len(automorphisms))
    print("N7_EXTENDED_ORBIT_GROUP_ORDER", 2 * len(automorphisms))


def main():
    verify_n6()
    verify_n7()
    print("N6_N7_MOBIUS_PARAMETRIZATION_PASS")


if __name__ == "__main__":
    main()
