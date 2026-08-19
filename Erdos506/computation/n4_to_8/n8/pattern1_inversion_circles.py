import sympy as sp


x0, x1, x2, y4, y5, k = sp.symbols("x0 x1 x2 y4 y5 k")


P = {
    0: sp.Matrix([x0, 0, 1]),
    1: sp.Matrix([x1, 0, 1]),
    2: sp.Matrix([x2, 0, 1]),
    4: sp.Matrix([0, y4, 1]),
    5: sp.Matrix([0, y5, 1]),
}
P[6] = sp.Matrix([
    -x0*x2*(y4 - y5),
    -y4*y5*(x0 - x2),
    -x0*y4 + x2*y5,
])
P[7] = sp.Matrix([
    x0*x2*(y4 - y5),
    -y4*y5*(x0 - x2),
    -x0*y5 + x2*y4,
])


def circle_det(block):
    rows = []
    for i in block:
        X, Y, Z = P[i]
        rows.append([X*X + Y*Y + 2*k*X*Y, X*Z, Y*Z, Z*Z])
    return sp.factor(sp.Matrix(rows).det())


circle_blocks = [(0, 2, 4, 5), (0, 1, 4, 6), (1, 2, 5, 6), (0, 1, 5, 7), (1, 2, 4, 7), (0, 2, 6, 7), (4, 5, 6, 7)]
for block in circle_blocks:
    f = circle_det(block)
    print("circle", "".join(map(str, block)), sp.factor(f))

R = x0*x1 - 2*x0*x2 + x1*x2
S = x0*x2 - y4*y5
print("relation_R", R)
circle_0245_reduced = sp.factor(circle_det((0, 2, 4, 5)) / ((x0 - x2)*(y4 - y5)))
assert sp.expand(circle_0245_reduced - S) == 0
print("circle_0245_reduced", circle_0245_reduced)

subs_x1 = sp.solve(sp.Eq(R, 0), x1)[0]
subs_y5 = x0*x2 / y4
print("x1_from_R", subs_x1)
print("y5_from_0245", subs_y5)

reduced_circles = {}
for block in [(0, 1, 4, 6), (1, 2, 5, 6), (0, 1, 5, 7), (1, 2, 4, 7), (0, 2, 6, 7), (4, 5, 6, 7)]:
    f = sp.factor(sp.together(circle_det(block).subs({x1: subs_x1, y5: subs_y5})).as_numer_denom()[0])
    reduced_circles[block] = f
    print("after_R_S", "".join(map(str, block)), f)

A = x0 * x2 + y4**2
B = y4 * (x0 + x2)
common = x0 * x2 - y4**2
for block, factor in [((0, 1, 4, 6), common * (A - k * B)),
                      ((0, 2, 6, 7), common * (k * A - B))]:
    denominator = sp.factor(sp.cancel(reduced_circles[block] / factor).as_numer_denom()[1])
    assert denominator == 1
print("CERTIFICATE PASS: after R=S=0 the two decisive circles give A=kB and kA=B.")
