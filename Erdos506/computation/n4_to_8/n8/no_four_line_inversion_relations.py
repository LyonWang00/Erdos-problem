import sympy as sp


def line(p, q):
    return sp.Matrix(p).cross(sp.Matrix(q))


def meet(l, m):
    return sp.Matrix(l).cross(sp.Matrix(m))


def det3(a, b, c):
    return sp.Matrix.hstack(sp.Matrix(a), sp.Matrix(b), sp.Matrix(c)).det()


def main():
    a = sp.symbols("a")
    p0 = sp.Matrix([0, 0, 1])
    p2 = sp.Matrix([1, 0, 1])
    p3 = sp.Matrix([a, 0, 1])
    p7 = sp.Matrix([a, 1, 1])

    l057 = line(p0, p7)
    l267 = line(p2, p7)
    l356 = None

    # Write p5 on line 057 and p6 on line 267.  The remaining line condition
    # 3,5,6 determines their parameters up to one scalar.
    s, r = sp.symbols("s r")
    p5 = sp.Matrix([s * a, s, 1])
    p6 = sp.Matrix([1 + r * (a - 1), r, 1])
    condition_356 = sp.factor(det3(p3, p5, p6))
    assert sp.expand(condition_356 - (-a * r + a * s + r * s - s)) == 0
    print("condition_356", condition_356)
    r_solution = sp.solve(sp.Eq(condition_356, 0), r)[0]
    print("r_from_356", sp.factor(r_solution))
    p6 = sp.simplify(p6.subs(r, r_solution))

    p4_from_245 = meet(line(p2, p5), line(p0, p6))
    p4_from_046 = meet(line(p0, p6), line(p2, p5))
    print("p5", [sp.factor(x) for x in p5])
    print("p6", [sp.factor(x) for x in p6])
    print("p4", [sp.factor(x) for x in p4_from_245])
    print("same_p4", [sp.factor(x) for x in p4_from_046])

    col_456 = sp.factor(det3(p4_from_245, p5, p6))
    col_457 = sp.factor(det3(p4_from_245, p5, p7))
    col_467 = sp.factor(det3(p4_from_245, p6, p7))
    col_567 = sp.factor(det3(p5, p6, p7))
    print("col_456", col_456)
    print("col_457", col_457)
    print("col_467", col_467)
    print("col_567", col_567)
    assert sp.simplify(det3(p2, p5, p4_from_245)) == 0
    assert sp.simplify(det3(p0, p6, p4_from_245)) == 0
    print("CERTIFICATE PASS: all derived intersections and prescribed line incidences are exact.")


if __name__ == "__main__":
    main()
