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


eqs = [
    N(x2, y2) - x2 - h * y2,
    N(x2, y2) - s * x2 - h * r * y2,
    N(x3, y3) - x3 - h * r * y3,
    N(x3, y3) - s * x3 - h * y3,
    h * r - s,
]

branch = {s: t * t + t}
print("branch s=t^2+t", flush=True)
nums = []
for i, expr in enumerate(eqs):
    val = sp.cancel(expr.subs(branch))
    num, den = sp.fraction(val)
    num = sp.factor(num)
    den = sp.factor(den)
    nums.append(num)
    print("eq", i, "num", num, flush=True)
    print("eq", i, "den", den, flush=True)

print("factor_set", flush=True)
for expr in sorted({sp.factor(v) for v in nums}, key=str):
    print(expr, flush=True)

G = sp.groebner(nums, g, h, t, order="lex")
print("groebner_len", len(G.polys), flush=True)
for poly in G.polys:
    print(sp.factor(poly.as_expr()), flush=True)
