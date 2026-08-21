import sympy as sp


a, b, A, u = sp.symbols("a b A u")
F1 = (b + 1) * A - 2 * b * u + b**2 * (1 - a)
F2 = (a + 1) * A - 2 * a * u + a**2 * (1 - b)


def main():
    sol = sp.solve([F1, F2], [A, u], dict=True)
    print(sol)
    for s in sol:
        print("A", sp.factor(s[A]))
        print("u", sp.factor(s[u]))
        print("v2", sp.factor(s[A] - s[u] ** 2))


if __name__ == "__main__":
    main()
