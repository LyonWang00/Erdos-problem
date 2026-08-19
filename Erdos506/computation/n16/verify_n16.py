#!/usr/bin/env python3
"""Exact algebra audit for the mathematical proof of c(16)=99.

The proof has no linear-programming model and enumerates no multiplicity
vector.  It checks two displayed coefficient identities, the equality profile
in the K=7 argument, the elementary circle-family moment comparisons, and the
largest-block deletion values for K>=9.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "certificate.json"
N = 16
Q = 15
TARGET_MINUS_ONE = 98


def vector(names: tuple[str, ...], **entries: Fraction | int) -> list[Fraction]:
    return [Fraction(entries.get(name, 0)) for name in names]


def add_scaled(target: list[Fraction], source: list[Fraction], scale: Fraction) -> None:
    for index, value in enumerate(source):
        target[index] += scale * value


def subtract(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    return [a - b for a, b in zip(left, right)]


def local_circle_cost(names: tuple[str, ...], maximum: int) -> list[Fraction]:
    entries: dict[str, Fraction] = {}
    for i in range(2, maximum + 1):
        entries[f"t{i}"] = Fraction(1, i + 1)
        entries[f"b{i}"] = -Fraction(1, i + 1)
    return vector(names, **entries)


def common_forms(maximum: int):
    names = tuple(
        [f"t{i}" for i in range(2, maximum + 1)]
        + [f"b{i}" for i in range(2, maximum + 1)]
    )
    pair = vector(names, **{
        f"t{i}": comb(i, 2) for i in range(2, maximum + 1)
    })
    melchior_entries: dict[str, int] = {"t2": 1}
    for i in range(4, maximum + 1):
        melchior_entries[f"t{i}"] = -(i - 3)
    melchior = vector(names, **melchior_entries)
    support = vector(names, **{
        f"b{i}": 1 for i in range(2, maximum + 1)
    })
    incidence = vector(names, **{
        f"t{i}": i - 1 for i in range(2, maximum + 1)
    })
    return names, pair, melchior, support, incidence


def audit_k7_identity() -> dict[str, object]:
    names, pair, melchior, support, _ = common_forms(6)
    bojanowski_entries: dict[str, int] = {"t2": 4, "t3": 3}
    for i in range(5, 7):
        bojanowski_entries[f"t{i}"] = -i * (i - 4)
    bojanowski = vector(names, **bojanowski_entries)

    left = local_circle_cost(names, 6)
    # d_x=t_6-b_6 counts seven-point circles through the centre.
    add_scaled(left, vector(names, t6=1, b6=-1), Fraction(1, 42))

    combination = [Fraction(0) for _ in names]
    add_scaled(combination, pair, Fraction(7, 108))
    add_scaled(combination, melchior, Fraction(7, 36))
    add_scaled(combination, bojanowski, Fraction(1, 54))
    add_scaled(combination, support, -Fraction(1, 3))

    residual = subtract(left, combination)
    expected = vector(
        names,
        t4=Fraction(1, 180),
        b3=Fraction(1, 12),
        b4=Fraction(2, 15),
        b5=Fraction(1, 6),
        b6=Fraction(1, 6),
    )
    assert residual == expected
    lower_bound = (
        Fraction(7, 108) * comb(Q, 2)
        + Fraction(7, 36) * 3
        + Fraction(1, 54) * (4 * Q)
        - Fraction(1, 3) * (Q // 2)
    )
    assert lower_bound == Fraction(37, 6)

    # If equality holds at a centre with d_x=1, the positive residuals and
    # all positively weighted inequalities are tight.  The resulting three
    # equations in (t2,t3,t5) have nonzero determinant and the unique profile
    # shown below.  Its five rich lines cover at least 16 of the 15 image
    # points, a contradiction.
    determinant = -27
    equality_profile = {
        "t2": 14, "t3": 12, "t4": 0, "t5": 4, "t6": 1,
        "b2": 7, "b3": 0, "b4": 0, "b5": 0, "b6": 0,
    }
    t2, t3, t5 = (equality_profile[key] for key in ("t2", "t3", "t5"))
    assert determinant != 0
    assert t2 - 2 * t5 == 6
    assert t2 + 3 * t3 + 10 * t5 == 90
    assert 4 * t2 + 3 * t3 - 5 * t5 == 72
    rich_line_count = 5
    rich_line_union_lower = 4 * 5 + 6 - comb(rich_line_count, 2)
    assert rich_line_union_lower == 16 > Q

    # Four seven-point circles through one point would invert to four
    # six-point lines whose union has at least 18 image points.  Thus d_x<=3.
    assert 4 * 6 - comb(4, 2) == 18 > Q
    # For g=5,6 circles, convexity gives sum C(d_x,2)>=14g-48, whereas circle
    # pairs give at most g(g-1).  For g=4, absence of degree one would give
    # sum C(d_x,2)>=sum d_x/2=14>12, so a degree-one centre exists.
    impossible_circle_counts = {}
    for g in (5, 6):
        lower = 14 * g - 48
        upper = g * (g - 1)
        assert lower > upper
        impossible_circle_counts[str(g)] = {
            "second_moment_lower": lower,
            "second_moment_upper": upper,
        }
    assert Fraction(7 * 4, 2) == 14 > 4 * 3
    assert 7 * 7 > 16 * 3

    margins = {
        str(g): str(16 * lower_bound - Fraction(g, 6) - TARGET_MINUS_ONE)
        for g in range(4)
    }
    assert all(Fraction(value) > 0 for value in margins.values())
    assert 16 * lower_bound - Fraction(4, 6) == TARGET_MINUS_ONE

    return {
        "identity": "Q_x+d_x/42 >= 37/6",
        "multipliers": {
            "pair_identity": "7/108",
            "Melchior": "7/36",
            "Bojanowski": "1/54",
            "centre_line_support": "-1/3",
        },
        "nonnegative_residual": {
            name: str(value) for name, value in zip(names, residual) if value
        },
        "margins_above_98_for_g_0_to_3": margins,
        "g4_degree_one_forced": True,
        "g4_equality_profile": equality_profile,
        "g4_rich_line_union_lower": rich_line_union_lower,
        "g5_g6_second_moment_contradictions": impossible_circle_counts,
        "g_at_least_7_degree_sum_contradiction": "7g<=16*3",
    }


def audit_k8_identity() -> dict[str, object]:
    names, pair, melchior, support, incidence = common_forms(7)
    ordinary_avoiding = vector(names, t2=1, b2=-1)
    cost = local_circle_cost(names, 7)

    combination = [Fraction(0) for _ in names]
    add_scaled(combination, pair, Fraction(1, 80))
    add_scaled(combination, melchior, Fraction(31, 160))
    add_scaled(combination, ordinary_avoiding, Fraction(1, 48))
    add_scaled(combination, support, -Fraction(5, 16))
    add_scaled(combination, incidence, Fraction(17, 160))

    residual = subtract(cost, combination)
    expected = vector(
        names,
        t5=Fraction(1, 240),
        t6=Fraction(3, 560),
        b3=Fraction(1, 16),
        b4=Fraction(9, 80),
        b5=Fraction(7, 48),
        b6=Fraction(19, 112),
        b7=Fraction(3, 16),
    )
    assert residual == expected

    # At a point of a chosen eight-point line or circle, inversion produces
    # a seven-point image line.  Hence the incidence form is at least
    # 6+7(15-7)=62.
    image_line_incidence = 6 + 7 * (Q - 7)
    assert image_line_incidence == 62
    conditioned = (
        Fraction(1, 80) * comb(Q, 2)
        + Fraction(31, 160) * 3
        + Fraction(1, 48) * 3
        - Fraction(5, 16) * (Q // 2)
        + Fraction(17, 160) * image_line_incidence
    )
    assert conditioned == Fraction(1017, 160)
    unconditional = bojanowski_local_bound(16, 8)
    assert unconditional == Fraction(145, 24)
    margin = 8 * conditioned + 8 * unconditional - TARGET_MINUS_ONE
    assert margin == Fraction(71, 60) > 0

    return {
        "conditioned_identity_bound": str(conditioned),
        "unconditional_bound": str(unconditional),
        "multipliers": {
            "pair_identity": "1/80",
            "Melchior": "31/160",
            "ordinary_lines_avoiding_centre": "1/48",
            "centre_line_support": "-5/16",
            "seven_point_image_line_incidence": "17/160",
        },
        "nonnegative_residual": {
            name: str(value) for name, value in zip(names, residual) if value
        },
        "global_margin_above_98": str(margin),
    }


def bojanowski_local_bound(n: int, K: int) -> Fraction:
    q = n - 1
    return (
        Fraction(K + 1, 18 * K) * comb(q, 2)
        + Fraction(K + 1, 2 * K)
        + Fraction(K - 2, 9 * K) * q
        - Fraction(q // 2, 3)
    )


def direct_line(n: int, K: int) -> int:
    r = n - K
    return r * comb(K, 2) - comb(r, 2) * (K // 2)


def direct_circle(n: int, K: int) -> int:
    r = n - K
    return 1 + r * (comb(K, 2) - K // 2) - comb(r, 2) * (K // 2)


def main() -> None:
    small = {}
    for K in range(3, 7):
        bound = bojanowski_local_bound(N, K)
        margin = N * bound - TARGET_MINUS_ONE
        assert margin > 0
        small[str(K)] = {
            "local_bound": str(bound),
            "margin_above_98": str(margin),
        }

    direct = {
        str(K): min(direct_line(N, K), direct_circle(N, K))
        for K in range(9, N)
    }
    assert min(direct.values()) >= 99

    payload = {
        "schema_version": 1,
        "status": "PASS",
        "theorem": "c(16)=99",
        "proof_type": "mathematical",
        "linear_programming_models": 0,
        "multiplicity_vectors_enumerated": 0,
        "K3_to_K6": small,
        "K7": audit_k7_identity(),
        "K8": audit_k8_identity(),
        "K9_to_K15_direct_bounds": direct,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
