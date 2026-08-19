import sympy as sp


x0, x1, x2, y4, y5 = sp.symbols("x0 x1 x2 y4 y5")


def line_through_axis_points(x, y):
    # Line through (x,0) and (0,y): X/x + Y/y = 1.
    # Return coefficients A X + B Y + C = 0.
    return sp.Matrix([y, x, -x * y])


def intersect(l1, l2):
    p = l1.cross(l2)
    return sp.Matrix([sp.factor(p[0]), sp.factor(p[1]), sp.factor(p[2])])


L05 = line_through_axis_points(x0, y5)
L24 = line_through_axis_points(x2, y4)
L04 = line_through_axis_points(x0, y4)
L25 = line_through_axis_points(x2, y5)
P6 = intersect(L05, L24)
P7 = intersect(L04, L25)
P1 = sp.Matrix([x1, 0, 1])

det = sp.factor(sp.Matrix.hstack(P1, P6, P7).det())
print("P6", [sp.factor(e) for e in P6])
print("P7", [sp.factor(e) for e in P7])
print("collinearity_1_6_7", det)
reduced = sp.factor(det / (x0 * x2 * y4 * y5 * (x0 - x2) * (y4 - y5)))
R = x0 * x1 - 2 * x0 * x2 + x1 * x2
assert sp.simplify(reduced + R / (x0 * x2)) == 0
print("reduced", reduced)
print("CERTIFICATE PASS: nondegeneracy reduces collinearity 1,6,7 to R=0.")
