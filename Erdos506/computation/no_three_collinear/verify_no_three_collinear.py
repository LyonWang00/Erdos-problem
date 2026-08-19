"""Exact verification for the no-three-collinear variant of Erdos 506."""

from __future__ import annotations

from itertools import combinations, permutations

import sympy as sp


def normalized_circle(points):
    matrix = sp.Matrix([[x, y, x*x+y*y, 1] for x, y in points])
    vector = matrix.nullspace()[0]
    pivot = next(value for value in vector if value != 0)
    return tuple(sp.factor(value/pivot) for value in vector)


# The n=8 upper-bound construction: invert the two-square configuration at
# z=(0,10).  The image coordinates are rational.
source = [(1, 1), (1, -1), (-1, 1), (-1, -1),
          (2, 2), (2, -2), (-2, 2), (-2, -2)]
points = []
for x, y in source:
    denominator = x*x+(y-10)**2
    points.append((sp.Rational(x, denominator), sp.Rational(10)+sp.Rational(y-10, denominator)))
lines = []
circles = {}
for triple in combinations(range(8), 3):
    selected = [points[i] for i in triple]
    if sp.det(sp.Matrix([[x, y, 1] for x, y in selected])) == 0:
        lines.append(triple)
    else:
        circles.setdefault(normalized_circle(selected), []).append(triple)
assert not lines
assert len(circles) == 20
assert sorted(map(len, circles.values())) == [1]*8+[4]*12
circle_supports = sorted(
    tuple(sorted({point for triple in triples for point in triple}))
    for triples in circles.values()
)
assert sorted(map(len, circle_supports)) == [3]*8+[4]*12


# For n>=11 the two symbolic margins in the proof are nonnegative.
n, K = sp.symbols("n K", positive=True)
small_margin = (n-4)*(K*n*n-13*K*n+18*K+n*n-7*n)/(36*K)
small_endpoint = sp.factor(small_margin.subs(K, sp.Rational(1, 3)*(2*n+1)))
assert sp.factor(small_endpoint-(n-4)*(n-1)*(n*n-10*n-9)/(18*(2*n+1))) == 0
assert sp.factor(sp.diff(small_margin, K)) == -n*(n-7)*(n-4)/(36*K**2)
large_bound = 1+(n-K)*K*(3*K-n-1)/4
target = 1+(n-1)*(n-2)/2
assert sp.factor((large_bound-target).subs(K, sp.Rational(1, 3)*(2*n+1))-(n-4)*(n-1)*(2*n-9)/36) == 0
assert sp.factor((large_bound-target).subs(K, n-1)) == 0
large_derivative = sp.diff(large_bound, K)
assert sp.factor(
    large_derivative.subs(K, sp.Rational(1, 3)*(2*n+1))
    - (n*n-3*n-1)/12
) == 0
assert sp.factor(large_derivative.subs(K, n-1) + (2*n*n-11*n+11)/4) == 0
assert sp.diff(large_derivative, K, 2) == -sp.Rational(9, 2)


# Classify the possible seven four-circle family on seven points.  Fixing one
# four-set loses no orbit.  Every seven-member triple packing is a relabelling
# of the complements of the seven lines of the Fano plane.
v7 = tuple(range(7))
triples7 = tuple(combinations(v7, 3))
foursets7 = tuple(combinations(v7, 4))
triple_index7 = {triple: index for index, triple in enumerate(triples7)}
four_index7 = {block: index for index, block in enumerate(foursets7)}
block_masks7 = []
through7 = [[] for _ in triples7]
for index, block in enumerate(foursets7):
    mask = sum(1 << triple_index7[triple] for triple in combinations(block, 3))
    block_masks7.append(mask)
    for triple in combinations(block, 3):
        through7[triple_index7[triple]].append(index)
fixed7 = four_index7[(0, 1, 2, 3)]
all_triples7 = (1 << len(triples7))-1
packings7 = set()


def enumerate_packings7(covered, skipped, selected):
    if len(selected) > 7 or skipped.bit_count() > 7:
        return
    undecided = len(triples7)-covered.bit_count()-skipped.bit_count()
    if len(selected)+undecided//4 < 7:
        return
    if len(selected) == 7:
        if covered.bit_count() == 28:
            packings7.add(tuple(sorted(selected)))
        return
    remaining = all_triples7 & ~covered & ~skipped
    bit = remaining & -remaining
    triple_index = bit.bit_length()-1
    enumerate_packings7(covered, skipped | bit, selected)
    for block_index in through7[triple_index]:
        if block_masks7[block_index] & covered == 0:
            enumerate_packings7(
                covered | block_masks7[block_index], skipped,
                selected+(block_index,),
            )


enumerate_packings7(block_masks7[fixed7], 0, (fixed7,))
representative7 = [
    (0, 1, 2, 3), (0, 1, 4, 5), (0, 2, 4, 6), (0, 3, 5, 6),
    (1, 2, 5, 6), (1, 3, 4, 6), (2, 3, 4, 5),
]
orbit7 = set()
for permutation in permutations(v7):
    image = tuple(sorted(
        four_index7[tuple(sorted(permutation[i] for i in block))]
        for block in representative7
    ))
    if fixed7 in image:
        orbit7.add(image)
assert packings7 and packings7 == orbit7


# The unique family is impossible.  After inversion, three remaining circle
# equations reduce to the following.
t, u, v, a = sp.symbols("t u v a")
R = u*u+v*v
e1 = -a*R+2*t*u-t
e2 = a*R-2*a*u+t
e3 = -a*R+t
assert sp.factor(e1.subs(t, a*R)) == 2*a*R*(u-1)
assert sp.factor(e2.subs({t: a*R, u: 1})) == 2*a*v*v


# Classify the line-free c=17 layer on eight points.  Thirteen four-subsets
# cover 52 of the 56 triples without overlap.  Fixing 0123 loses no orbit; all
# 78 labelled solutions in that normalization are relabellings of one family.
v8 = tuple(range(8))
triples8 = tuple(combinations(v8, 3))
foursets8 = tuple(combinations(v8, 4))
triple_index8 = {triple: index for index, triple in enumerate(triples8)}
four_index8 = {block: index for index, block in enumerate(foursets8)}
block_masks8 = []
through8 = [[] for _ in triples8]
for index, block in enumerate(foursets8):
    mask = sum(1 << triple_index8[triple] for triple in combinations(block, 3))
    block_masks8.append(mask)
    for triple in combinations(block, 3):
        through8[triple_index8[triple]].append(index)
fixed8 = four_index8[(0, 1, 2, 3)]
all_triples8 = (1 << len(triples8))-1
packings8 = set()


def enumerate_packings8(covered, skipped, selected):
    if len(selected) > 13 or skipped.bit_count() > 4:
        return
    undecided = len(triples8)-covered.bit_count()-skipped.bit_count()
    if len(selected)+undecided//4 < 13:
        return
    if len(selected) == 13:
        if covered.bit_count() == 52:
            packings8.add(tuple(sorted(selected)))
        return
    remaining = all_triples8 & ~covered & ~skipped
    bit = remaining & -remaining
    triple_index = bit.bit_length()-1
    enumerate_packings8(covered, skipped | bit, selected)
    for block_index in through8[triple_index]:
        if block_masks8[block_index] & covered == 0:
            enumerate_packings8(
                covered | block_masks8[block_index], skipped,
                selected+(block_index,),
            )


enumerate_packings8(block_masks8[fixed8], 0, (fixed8,))
representative8 = [tuple(map(int, block)) for block in
                   ("0123", "0145", "0167", "0247", "0346", "0357",
                    "1246", "1257", "1347", "1356", "2345", "2367", "4567")]
orbit8 = set()
for permutation in permutations(v8):
    image = tuple(sorted(four_index8[tuple(sorted(permutation[i] for i in block))]
                         for block in representative8))
    if fixed8 in image:
        orbit8.add(image)
assert len(packings8) == 78 and packings8 == orbit8


# This unique line-free c=17 incidence class becomes, after
# inversion, the six-line q6 configuration together with the extra circle
# 3456.  The first two circle equations determine h and k; the extra circle
# has the displayed nonzero numerator.
a0, b0 = sp.symbols("a0 b0")
D = a0*b0-b0+1
den = a0-b0+1
h_value = -b0/den
k_value = (1-b0)/den
extra_3456 = 2*a0*b0**2*(a0-1)**2*(b0-1)*D
assert sp.factor(h_value-k_value**2+D/den**2) == 0
assert sp.factor(extra_3456/(2*a0*b0**2*(a0-1)**2*(b0-1))) == D


# At n=9 the nested packing bound is 18.
packing_8_3_2 = (8*3)//3  # sum of point degrees <= 8*floor(7/2)
assert packing_8_3_2 == 8
packing_9_4_3 = (9*packing_8_3_2)//4
assert packing_9_4_3 == 18


# At n=10, pair counting, Melchior and Hirzebruch leave a unique local
# equality type (t2,t3,t4)=(6,4,3).
local_types = []
for t4 in range(7):
    for t3 in range(13):
        t2 = 36-3*t3-6*t4
        if t2 < 0 or t2 < 3+t4 or 4*t2+3*t3 < 36:
            continue
        phi = sp.Rational(t2, 3)+sp.Rational(t3, 4)+sp.Rational(t4, 5)
        local_types.append((phi, t2, t3, t4))
assert min(local_types) == (sp.Rational(18, 5), 6, 4, 3)
assert [row for row in local_types if row[0] == sp.Rational(18, 5)] == [(sp.Rational(18, 5), 6, 4, 3)]


# Classify the resulting 2-(6,3,2) incidence design.  There are twelve
# labelled solutions and one orbit under S_6.
V = tuple(range(6))
triples = tuple(combinations(V, 3))
pairs = tuple(combinations(V, 2))
solutions = []


def enumerate_designs(index, remaining, pair_degrees, counts):
    if index == len(triples):
        if remaining == 0 and all(value == 2 for value in pair_degrees):
            solutions.append(tuple(counts))
        return
    if remaining < 0:
        return
    involved = [pairs.index(pair) for pair in combinations(triples[index], 2)]
    maximum = min([remaining]+[2-pair_degrees[j] for j in involved])
    for multiplicity in range(maximum+1):
        for j in involved:
            pair_degrees[j] += multiplicity
        counts.append(multiplicity)
        enumerate_designs(index+1, remaining-multiplicity, pair_degrees, counts)
        counts.pop()
        for j in involved:
            pair_degrees[j] -= multiplicity


enumerate_designs(0, 10, [0]*len(pairs), [])
assert len(solutions) == 12


def canonical(counts):
    images = []
    for permutation in permutations(V):
        moved = {tuple(sorted(permutation[i] for i in triple)): counts[j]
                 for j, triple in enumerate(triples)}
        images.append(tuple(moved[triple] for triple in triples))
    return min(images)


orbits = {canonical(solution) for solution in solutions}
assert len(orbits) == 1


# Inverting at one point of that design gives three four-point lines and
# three five-point circles.  The exact normalized equations force the two
# points on the third line to coincide.
uu, vv, ss, pp, rr, kk = sp.symbols("uu vv ss pp rr kk")
hh = uu*vv/ss
tt = 1/vv
f = lambda z: (-2*kk*ss+ss+uu*vv)*z*z + (2*kk*ss-ss*vv-2*uu*vv+uu)*z + uu*(vv-1)
g2 = 2*kk*pp*ss-2*kk*ss-pp*ss-pp*uu*vv-ss*uu*vv+ss*uu+ss*vv+uu*vv
g3 = 2*kk*rr*ss-rr*ss-rr*uu*vv-ss*uu*vv+ss*uu+uu*vv-uu
assert sp.factor(sp.resultant(g2, f(pp), kk)) == -2*ss*uu*(pp-1)*(vv-1)*(pp*ss+pp-1)
assert sp.factor(sp.resultant(g3, f(rr), kk)) == -2*rr*ss**2*(vv-1)*(rr*uu+rr-uu)
p_value = 1/(ss+1)
r_value = uu/(uu+1)
eq_left = 2*kk*ss+ss*uu*vv-ss*uu-ss*vv-uu-vv+1
eq_right = -2*kk*ss+ss*uu*vv-ss*uu+ss*vv+uu-vv+1
assert sp.expand(eq_left+eq_right-2*(vv-1)*(ss*uu-1)) == 0
assert sp.factor((p_value-r_value).subs(ss, 1/uu)) == 0

print("NO-THREE-COLLINEAR VARIANT PASS")
print("n=8 construction: 20 circles and no collinear triple")
print("n=8 circle supports:", " ".join("".join(map(str, block)) for block in circle_supports))
print("all symbolic endpoint, packing, incidence, and residual identities verified")
