import sympy as sp


a, b, u, v, t, s, r, w = sp.symbols("a b u v t s r w")
A = u**2 + v**2
Wa = (u - a) ** 2 + v**2
W1 = (u - 1) ** 2 + v**2

polys = [
    a**2*t - a**2 - 2*a*t*u + a + t*A,
    r**2*v - r*v - s*A*w + 2*s*u*w - s*w - u*w + v*w**2 + w,
    -a**2*t*w + a**2*w - a*r*v + 2*a*t*u*w - a*u*w + r**2*v - t*A*w + v*w**2,
    a + s*A - 2*s*u + s - 1,
    -a**2*s*t**2 + 2*a**2*s*t - a**2*s - a*b*s*t + a*b*s
    + 2*a*s*t**2*u - 2*a*s*t*u + b*s*t - b*t + s**2*t*A
    - s*t**2*A + 2*s*t*u - 2*s*t + t,
    -b*r*v + b*u*w + r**2*v - A*w + v*w**2,
    -a**2*t + a**2 + 2*a*t*u - 2*a*u + s*A - 2*s*u + s - t*A + 2*u - 1,
    -a*r*v + a*u*w + a*v - a*w + r**2*v - r*v - A*w + u*w + v*w**2,
    -a**2*t**2*w + 2*a**2*t*w - a**2*w - a*b*t*w + a*b*w
    + 2*a*t**2*u*w - 2*a*t*u*w - a*t*w + a*w - b*r*t*v
    + b*t*u*w + b*t*v - b*w + r**2*t*v - r*t*v - t**2*A*w + t*u*w + t*v*w**2,
    a*b*s*v - a*b*w - a*r*s*v + a*s*u*w - a*s*w + a*w
    - b*r*s*v + b*s*u*w - b*s*w + b*w + r**2*s*v - s**2*A*w
    + 2*s**2*u*w - s**2*w - 2*s*u*w + s*v*w**2 + 2*s*w - w,
    a**2*r*t*v - a**2*r*v - a**2*t*u*w - a**2*t*v + a**2*t*w + a**2*u*w
    + a**2*v - a**2*w + a*r**2*v - 2*a*r*t*u*v + a*s*u**2*v
    - a*s*u**2*w - 2*a*s*u*v + 2*a*s*u*w + a*s*v**3 - a*s*v**2*w
    + a*s*v - a*s*w + 2*a*t*u**2*w + 2*a*t*u*v - 2*a*t*u*w
    - a*u**2*w - a*v**2*w + a*v*w**2 - a*v + a*w - r**2*v
    - r*s*A*v + 2*r*s*u*v - r*s*v + r*t*A*v + r*v
    + s*u*A*w - 2*s*u**2*w + s*u*w - t*u*A*w - t*u**2*v
    + t*u**2*w - t*v**3 + t*v**2*w + u**2*w - u*w + v**2*w - v*w**2,
]


def main():
    print("polynomial_count", len(polys))
    G = sp.groebner(polys, r, w, s, t, v, u, b, a, order="lex")
    print("basis_len", len(G.polys))
    for g in G.polys[:100]:
        print(sp.factor(g.as_expr()))


if __name__ == "__main__":
    main()
