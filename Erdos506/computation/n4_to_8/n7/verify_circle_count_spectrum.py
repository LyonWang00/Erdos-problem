"""Verify the seven-point extremal incidence type and the gap at 12 circles.

All computations are exact.  The combinatorial part works up to relabelling;
the symbolic identities cover the four five-line q6 orbits not settled by the
two-line lemma in the paper.  The q7 orbit is checked separately by
``verify_q7_exclusion.py``.
"""

from __future__ import annotations

import sympy as sp

from classify_line_extensions import (
    PATTERNS,
    Q7_CORE,
    automorphisms,
    canonical_lines,
    contains_orbit,
    extensions,
)


def classes_at(target):
    answer = {}
    for name, quads in PATTERNS.items():
        autos = automorphisms(quads)
        classes = {}
        for lines, covered in extensions(quads, target):
            key = canonical_lines(lines, autos)
            classes[key] = (len(covered), 35 - len(covered) - 3*len(quads))
        answer[name] = classes
    return answer


at_11 = classes_at(11)
at_12 = classes_at(12)
assert not at_11["q5a"] and not at_11["q5b"]
assert sorted(value for value in at_11["q6"].values() if value[1] == 11) == [(6, 11), (6, 11)]
assert len([1 for value in at_11["q7"].values() if value[1] == 11]) == 2
assert not at_12["q5a"] and not at_12["q5b"]
assert len([1 for value in at_12["q6"].values() if value[1] == 12]) == 5
assert len([1 for value in at_12["q7"].values() if value[1] == 12]) == 1

q6_c11 = {
    lines for lines, (_, circle_count) in at_11["q6"].items()
    if circle_count == 11
}
assert q6_c11 == {
    ((0, 1, 2), (0, 3, 6), (0, 4, 5), (1, 3, 5), (1, 4, 6), (2, 3, 4)),
    ((0, 3, 6), (0, 4, 5), (1, 3, 5), (1, 4, 6), (2, 3, 4), (2, 5, 6)),
}
q7_autos = automorphisms(PATTERNS["q7"])
assert all(
    contains_orbit(lines, Q7_CORE, q7_autos)
    for lines, (_, circle_count) in at_11["q7"].items()
    if circle_count == 11
)


# Orbit q6a.  After a similarity normalization, two points on the circle
# x^2-x+y^2=0 have rational parameters u and v.  The remaining circle
# condition and the coordinate of p2 have the following forms.
u, v = sp.symbols("u v", nonzero=True)
q6a_condition = -2*u*v*(u-v)*(u*v+1)
q6a_p2_parameter = 1/(u*v+1)
assert sp.factor(q6a_condition/(-2*u*v*(u-v))) == u*v+1
assert sp.denom(q6a_p2_parameter) == u*v+1


# Orbit q6b.  The nonzero factors of two circle equations give t=x; the
# remaining equation is then nonzero because y != 0 and x != 1.
x, y, alpha, beta, t = sp.symbols("x y alpha beta t")
R = x*x+y*y
q6b_C = beta*R-2*beta*x+beta+t-1
q6b_F = beta*R-2*beta*x+beta-t+2*x-1
assert sp.expand(q6b_C-q6b_F-2*(t-x)) == 0
S = (x-1)**2+y*y
q6b_A = -2*beta*t*x+2*beta*t+beta*R-beta-t+1
assert sp.factor(q6b_A.subs({t: x, beta: (1-x)/S})-2*y*y*(1-x)/S) == 0


# Orbit q6c.  Three equations imply x=1 and then t=1, which identifies two
# labelled points.
gamma = sp.symbols("gamma")
q6c_D = alpha*R-gamma*(t*t-2*t*x+R)+t*t-2*t*x
subs_c = {alpha: t/R, gamma: t*(t-1)/(t*t-2*t*x+R)}
assert sp.factor(q6c_D.subs(subs_c)) == -2*t*(x-1)
q6c_B = alpha*R-beta*R+2*beta*x-beta-2*x+1
subs_b = {alpha: t/R, x: 1, beta: (1-t)/(R.subs(x, 1)-1)}
assert sp.factor(q6c_B.subs(subs_b)-2*(t-1)) == 0


# Orbit q6d.  In an affine normalization retaining the quadratic part
# x^2+2kxy+h y^2, the circle equations force E=0.  The determinant of the
# nominally missing line is a nonzero factor times exactly E, so this orbit
# actually acquires a sixth line and has eleven rather than twelve circles.
# The same identity excludes the other abstract c=11 extension: that extension
# also contains line 012, and forcing line 256 would put all three diagonal
# points of the complete quadrilateral on line 012.
a, b, h, k = sp.symbols("a b h k")
E = -a*b+a*t-b*t+b+t-1
line_256_numerator = -b*(a-1)*E
assert sp.factor(line_256_numerator/(-b*(a-1))) == E
solution = {h: -b/(a-b+1), k: (1-b)/(a-b+1)}
first_circle = 2*a*b*h-2*a*b*k-a*h-b*h+b+h
second_circle = a*h+b*h-2*b*k+b-h
assert sp.factor(first_circle.subs(solution)) == 0
assert sp.factor(second_circle.subs(solution)) == 0

print("SEVEN-POINT SPECTRUM PASS")
print("c=11 has one realizable incidence orbit; all c=12 incidence orbits are excluded.")
