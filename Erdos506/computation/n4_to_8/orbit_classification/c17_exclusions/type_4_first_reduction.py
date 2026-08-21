import sympy as sp


a, b, u, v, B = sp.symbols("a b u v B")
A = u**2 + v**2
W = (u - 1) ** 2 + v**2
C = (a - 1) * (b - 1)
D = sp.factor(v * (a * b - B) / (b - a))
n = sp.factor(b * v * (B - a**2) / (A * (b - a)))
m = sp.factor((D + n * u) / v)
t = sp.factor(C / W)
s = sp.factor(a * b / B)


def num(expr):
    return sp.factor(sp.together(expr).as_numer_denom()[0])


eqB = num(B - m**2 - n**2)
eq0257 = num(
    -a * m * t * v + a * n * t * u - a * n * t + a * n
    + B * s * t * v - n * t**2 * W - 2 * n * t * u + 2 * n * t - n
)
eq0356 = num(
    b * m * t * v - b * n * t * u + b * n * t - b * n
    - B * t * v + n * t**2 * W + 2 * n * t * u - 2 * n * t + n
)
eq1247 = num(
    -a * m * s * v + a * n * s * u - a * n * s + a * v
    + B * s**2 * v - m * s * v - n * s * A + n * s * u
)
eq1256 = num(
    t * (a * m * v - a * n * u + a * n - a * v
         - B * v + m * v + n * t * W + n * u - n)
)
eq1346 = num(b * (m * v - n * u) + b * n - b * v - B * v + m * v + n * A - n * u)
eq1357 = num(
    -t * (-b * m * s * v + b * n * s * u - b * n * s + b * v
          + B * s**2 * v - m * s * v - n * s * t * W - n * s * u + n * s)
)
eq4567 = num(
    -m**3 * s * v + m**2 * n * s * u - m**2 * n * s + m**2 * s * v + m**2 * v
    - m * n**2 * s * v + m * t * W * v - m * v
    + n**3 * s * u - n**3 * s + n**2 * s * v + n**2 * v
    - n * t * (u - 1) * W - n * A + n * u
)


def main():
    for name, expr in [
        ("eqB", eqB),
        ("0257", eq0257),
        ("0356", eq0356),
        ("1247", eq1247),
        ("1256", eq1256),
        ("1346", eq1346),
        ("1357", eq1357),
        ("4567", eq4567),
    ]:
        print(name, sp.factor(expr))
        print("factor_list", sp.factor_list(expr))


if __name__ == "__main__":
    main()
