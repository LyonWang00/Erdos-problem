import sympy as sp


a, b, u, v = sp.symbols("a b u v")
A = u**2 + v**2
n = sp.factor(b * (a - b) * v / (A - b**2))
m = sp.factor((b * (v - n) + n * u) / v)
B = sp.factor(m**2 + n**2)
t = sp.factor(a * b / A)
s = sp.factor(a * b / B)
D = sp.factor(m * v - n * u)


def num(expr):
    return sp.factor(sp.together(expr).as_numer_denom()[0])


eqs = {}
eqs["relation_B"] = num(v * (B - a * b) - n * (A - a * b))
eqs["1247"] = num(
    -a * m * s * v + a * n * s * u - a * n * s + a * v
    + B * s**2 * v - m * s * v - n * s * A + n * s * u
)
eqs["1256"] = num(
    a * m * t * v - a * n * t * u + a * n - a * t * v
    - B * t * v + m * t * v + n * t**2 * A - n * t * u
)
eqs["1346"] = num(b * D + b * n - b * v - B * v + m * v + n * A - n * u)
eqs["1357"] = num(
    -b * m * s * t * v + b * n * s * t * u - b * n * s + b * t * v
    + B * s**2 * t * v - m * s * t * v - n * s * t**2 * A + n * s * t * u
)


def main():
    print("n =", sp.factor(n))
    print("m =", sp.factor(m))
    print("B =", sp.factor(B))
    for name, expr in eqs.items():
        print(name, sp.factor(expr))
        print("factor_list", sp.factor_list(expr))


if __name__ == "__main__":
    main()
