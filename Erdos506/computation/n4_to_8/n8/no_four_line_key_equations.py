import sympy as sp


def line(p, q):
    return sp.Matrix([p[0], p[1], 1]).cross(sp.Matrix([q[0], q[1], 1]))


def meet_affine(l, m):
    p = sp.Matrix(l).cross(sp.Matrix(m))
    return (sp.factor(p[0] / p[2]), sp.factor(p[1] / p[2]))


def circle_det(points):
    rows = []
    for x, y in points:
        rows.append([sp.expand(x * x + y * y), x, y, 1])
    return sp.factor(sp.Matrix(rows).det())


def numerator(expr):
    return sp.factor(sp.together(expr).as_numer_denom()[0])


def main():
    a, b, u, s = sp.symbols("a b u s")
    p0 = (a, 0)
    p2 = (b, 0)
    p3 = (0, 0)
    p7 = (u, 1)
    p5 = (sp.expand((1 - s) * a + s * u), s)
    p6 = meet_affine(line(p2, p7), line(p3, p5))
    p4 = meet_affine(line(p0, p6), line(p2, p5))

    e0256 = numerator(circle_det([p0, p2, p5, p6]))
    e0367 = numerator(circle_det([p0, p3, p6, p7]))
    e2357 = numerator(circle_det([p2, p3, p5, p7]))
    f0256 = a**2 * s - a**2 - a * b - 2 * a * s * u + 2 * a * u + s * u**2 + s
    f0367 = -a**2 * s + a**2 + 2 * a * b * s - a * b - 2 * b * s * u + s * u**2 + s
    f2357 = a**2 * s - a**2 + a * b - 2 * a * s * u + s * u**2 + s
    for expression, factor in [(e0256, f0256), (e0367, f0367), (e2357, f2357)]:
        denominator = sp.factor(sp.cancel(expression / factor).as_numer_denom()[1])
        assert denominator == 1
    assert sp.expand(f0256 - f2357 - 2 * a * (u - b)) == 0
    assert sp.expand(
        f0256.subs(u, b) + f0367.subs(u, b) - 2 * s
    ) == 0
    print("p4", [sp.factor(x) for x in p4])
    print("p5", [sp.factor(x) for x in p5])
    print("p6", [sp.factor(x) for x in p6])
    print("E0256", e0256)
    print("E0367", e0367)
    print("E2357", e2357)
    print("F0256", sp.factor(f0256))
    print("F0367", sp.factor(f0367))
    print("F2357", sp.factor(f2357))
    print("F0256_minus_F2357", sp.factor(f0256 - f2357))
    print("F0256_at_u_eq_b", sp.factor(f0256.subs(u, b)))
    print("F0367_at_u_eq_b", sp.factor(f0367.subs(u, b)))
    print("F0256_plus_F0367_after_u_eq_b", sp.factor(f0256.subs(u, b) + f0367.subs(u, b)))
    print("CERTIFICATE PASS: the three circle equations force u=b and then s=0.")


if __name__ == "__main__":
    main()
