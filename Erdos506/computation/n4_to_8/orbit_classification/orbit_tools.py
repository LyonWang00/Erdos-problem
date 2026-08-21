"""Exact incidence enumeration and Mobius/anti-Mobius comparison for n <= 8."""

from __future__ import annotations

from itertools import combinations, permutations

import sympy as sp


def line_det(points):
    return sp.factor(sp.Matrix([[x, y, 1] for x, y in points]).det())


def circle_det(points):
    return sp.factor(
        sp.Matrix([[x * x + y * y, x, y, 1] for x, y in points]).det()
    )


def maximal_incidence(points):
    n = len(points)
    lines = set()
    circles = set()
    for triple in combinations(range(n), 3):
        base = [points[i] for i in triple]
        if line_det(base) == 0:
            block = frozenset(
                i for i in range(n) if line_det([base[0], base[1], points[i]]) == 0
            )
            lines.add(block)
        else:
            block = frozenset(
                i for i in range(n) if circle_det(base + [points[i]]) == 0
            )
            circles.add(block)
    return lines, circles


def block_string(block):
    return "".join(str(i) for i in sorted(block))


def incidence_summary(points):
    lines, circles = maximal_incidence(points)
    return {
        "lines": sorted(map(block_string, lines), key=lambda s: (len(s), s)),
        "circles": sorted(map(block_string, circles), key=lambda s: (len(s), s)),
        "circle_count": len(circles),
    }


def large_generalized_blocks(points):
    lines, circles = maximal_incidence(points)
    return {block for block in lines | circles if len(block) >= 4}


def incidence_isomorphisms(source, target):
    if len(source) != len(target):
        return []
    n = len(source)
    source_blocks = large_generalized_blocks(source)
    target_blocks = large_generalized_blocks(target)
    if sorted(map(len, source_blocks)) != sorted(map(len, target_blocks)):
        return []
    out = []
    for perm in permutations(range(n)):
        image = {frozenset(perm[i] for i in block) for block in source_blocks}
        if image == target_blocks:
            out.append(perm)
    return out


def _normalized_equal(source_z, target_z, perm, anti):
    if anti:
        source_z = [sp.conjugate(z) for z in source_z]
    a, b, c = 0, 1, 2
    za, zb, zc = source_z[a], source_z[b], source_z[c]
    wa, wb, wc = target_z[perm[a]], target_z[perm[b]], target_z[perm[c]]
    for i, zi in enumerate(source_z):
        wi = target_z[perm[i]]
        source_num = (zi - za) * (zb - zc)
        source_den = (zi - zc) * (zb - za)
        target_num = (wi - wa) * (wb - wc)
        target_den = (wi - wc) * (wb - wa)
        if sp.simplify(source_num * target_den - target_num * source_den) != 0:
            return False
    return True


def mobius_equivalences(source, target):
    """Return exact label permutations for direct and orientation-reversing maps."""
    source_z = [sp.simplify(x + sp.I * y) for x, y in source]
    target_z = [sp.simplify(x + sp.I * y) for x, y in target]
    direct = []
    reverse = []
    isomorphisms = incidence_isomorphisms(source, target)
    for perm in isomorphisms:
        if _normalized_equal(source_z, target_z, perm, False):
            direct.append(perm)
        if _normalized_equal(source_z, target_z, perm, True):
            reverse.append(perm)
    return {
        "incidence_isomorphism_count": len(isomorphisms),
        "mobius": direct,
        "anti_mobius": reverse,
    }


def unordered_block_j(points, block):
    i0, i1, i2, i3 = sorted(block)
    z = [sp.simplify(points[i][0] + sp.I * points[i][1]) for i in (i0, i1, i2, i3)]
    lam = sp.cancel((z[0] - z[2]) * (z[1] - z[3]) / ((z[0] - z[3]) * (z[1] - z[2])))
    return sp.factor((1 - lam + lam**2) ** 3 / (lam**2 * (1 - lam) ** 2))


def four_block_j_multiset(points):
    return sorted(
        [unordered_block_j(points, block) for block in large_generalized_blocks(points)],
        key=sp.default_sort_key,
    )
