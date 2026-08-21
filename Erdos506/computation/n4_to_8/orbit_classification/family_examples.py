"""Construct exact family members and certify distinct Mobius orbits."""

from __future__ import annotations

import sympy as sp

from orbit_tools import (
    four_block_j_multiset,
    incidence_summary,
    mobius_equivalences,
)


Q = sp.Rational


def line_intersection(p1, p2, p3, p4):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    return (
        sp.factor(((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / den),
        sp.factor(((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / den),
    )


def invert_point(point, center):
    dx = point[0] - center[0]
    dy = point[1] - center[1]
    den = dx * dx + dy * dy
    return (sp.factor(center[0] + dx / den), sp.factor(center[1] + dy / den))


def n6_family(a, v):
    """Similarity normalization of the unique c=8 incidence type."""
    s = sp.factor(a / (1 + v**2))
    t = sp.factor(a * (a - 1) / ((a - 1) ** 2 + v**2))
    p0 = (0, 0)
    p2 = (1, 0)
    p4 = (a, 0)
    p3 = (1, v)
    p5 = (s, s * v)
    p1 = (sp.factor(a + t * (1 - a)), sp.factor(t * v))
    return [p0, p1, p2, p3, p4, p5]


def n7_family(u, v):
    """Two-parameter complete-quadrilateral family of c=11 configurations."""
    p3 = (0, 0)
    p4 = (1, 0)
    p5 = (u, v)
    p6 = (u, sp.factor(u * (1 - u) / v))
    p0 = line_intersection(p3, p6, p4, p5)
    p1 = line_intersection(p3, p5, p4, p6)
    p2 = line_intersection(p3, p4, p5, p6)
    return [p0, p1, p2, p3, p4, p5, p6]


def n8_c17_family(b):
    """One-parameter rational branch containing the published c=17 example."""
    m = sp.factor(3 / (3 - b**2))
    t = sp.factor((b**2 - 3) / (b**2 + 9))
    r = sp.factor(m * (1 + b**2))
    center = (sp.factor((1 - b**2) / (1 + b**2)), sp.factor(2 * b / (1 + b**2)))
    pre = [None] * 7
    pre[0] = (0, 0)
    pre[1] = (1, 0)
    pre[6] = (r, 0)
    pre[2] = (1, b)
    pre[3] = (sp.factor(1 + t * (r - 1)), sp.factor(b * (1 - t)))
    pre[4] = line_intersection(pre[0], pre[3], pre[1], pre[2])
    pre[5] = (m, sp.factor(m * b))
    return [invert_point(point, center) for point in pre] + [center]


def n8_c18_rectangles(aspect, r):
    return [
        (aspect, 1), (aspect, -1), (-aspect, 1), (-aspect, -1),
        (r * aspect, r), (r * aspect, -r), (-r * aspect, r), (-r * aspect, -r),
    ]


def n8_c18_two_squares(r):
    return n8_c18_rectangles(sp.Integer(1), r)


def compare(name, first, second, expected_count):
    first_summary = incidence_summary(first)
    second_summary = incidence_summary(second)
    assert first_summary["circle_count"] == expected_count, first_summary
    assert second_summary["circle_count"] == expected_count, second_summary
    assert sorted(map(len, first_summary["lines"])) == sorted(map(len, second_summary["lines"]))
    result = mobius_equivalences(first, second)
    print(name)
    print(" first", first_summary)
    print(" second", second_summary)
    print(" equivalence", result)
    print(" J_equal", four_block_j_multiset(first) == four_block_j_multiset(second))
    assert not result["mobius"] and not result["anti_mobius"]


def main():
    compare("n6_c8", n6_family(Q(3), Q(1)), n6_family(Q(4), Q(1)), 8)
    compare("n7_c11", n7_family(Q(3, 2), Q(1)), n7_family(Q(5, 3), Q(1)), 11)
    compare("n8_c17", n8_c17_family(Q(13, 12)), n8_c17_family(Q(1)), 17)
    compare("n8_c18", n8_c18_two_squares(Q(2)), n8_c18_two_squares(Q(3)), 18)
    compare(
        "n8_c18_aspect",
        n8_c18_rectangles(Q(1), Q(2)),
        n8_c18_rectangles(Q(2), Q(2)),
        18,
    )


if __name__ == "__main__":
    main()
