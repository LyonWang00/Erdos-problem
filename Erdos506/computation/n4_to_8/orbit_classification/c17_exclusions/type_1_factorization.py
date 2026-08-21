import sympy as sp


a, b, u, v, t, m, n, s = sp.symbols("a b u v t m n s")


def det_circle(points):
    return sp.factor(sp.Matrix([[x * x + y * y, x, y, 1] for x, y in points]).det())


P = {
    0: (0, 0),
    1: (1, 0),
    2: (a, 0),
    3: (b, 0),
    4: (u, v),
    5: (t * u, t * v),
    6: (m, n),
    7: (s * m, s * n),
}


blocks = "0246 0257 0347 0356 1247 1256 1346 1357 2345 2367 4567".split()


def main():
    for block in blocks:
        ids = [int(ch) for ch in block]
        expr = det_circle([P[i] for i in ids])
        print(block, sp.factor(expr))


if __name__ == "__main__":
    main()
