import sympy as sp


s, t, g, h = sp.symbols("s t g h")
k = s * t - 2 * s + 1
r = -k / (s - t)
x2 = s * (t - 1) / (t * (s - 1))
y2 = -k / (t * (s - 1))
x3 = t * x2
y3 = t * y2


def N(x, y):
    return x * x + 2 * g * x * y + h * y * y


raw = {
    "1246": N(x2, y2) - x2 - h * y2,
    "1257": N(x2, y2) - s * x2 - h * r * y2,
    "1347": N(x3, y3) - x3 - h * r * y3,
    "1356": N(x3, y3) - s * x3 - h * y3,
    "4567": h * r - s,
}

eqs = {}
print("axis_equations", flush=True)
for name, expr in raw.items():
    num, den = sp.fraction(sp.cancel(expr))
    eqs[name] = sp.factor(num)
    print(name, "num", eqs[name], flush=True)
    print(name, "den", sp.factor(den), flush=True)

sol_gh = sp.solve([eqs["4567"], eqs["1347"] - t * t * eqs["1246"]], [h, g], dict=True)
print("sol_gh_count", len(sol_gh), flush=True)
for sol in sol_gh:
    print("g", sp.factor(sol[g]), flush=True)
    print("h", sp.factor(sol[h]), flush=True)
    for name, expr in eqs.items():
        val = sp.cancel(expr.subs(sol))
        num, den = sp.fraction(val)
        print("reduced", name, "num", sp.factor(num), flush=True)
        print("reduced", name, "den", sp.factor(den), flush=True)
