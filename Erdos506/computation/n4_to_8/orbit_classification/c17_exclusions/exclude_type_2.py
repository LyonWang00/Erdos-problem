import sympy as sp


a, b, u, v, t, s, r, w = sp.symbols("a b u v t s r w")
A = u**2 + v**2

polys = [
    r**2*t*v - r*t*v - s*t**2*A*w + 2*s*t*u*w - s*w - t*u*w + t*v*w**2 + w,
    a*s - a + s**2*t**2*A - 2*s**2*t*u + s**2 - s*t*A + 2*s*t*u - 2*s + 1,
    -a*r*v + a*u*w + r**2*v - t*A*w + v*w**2,
    -b*r*v + b*u*w + r**2*v - A*w + v*w**2,
    b + s*t**2*A - 2*s*t*u + s - 1,
    -a*r*v + a*u*w + a*v - a*w + r**2*v - r*v - A*w + u*w + v*w**2,
    -b*t + b + s*t**2*A - 2*s*t*u + s - t*A + 2*t*u - 1,
    -b*r*t*v + b*t*u*w + b*t*v - b*w + r**2*t*v - r*t*v - t**2*A*w + t*u*w + t*v*w**2,
    t*A - a*b,
    a*b*s*t*v - a*b*w - a*r*s*t*v + a*s*t*u*w - a*s*w + a*w
    - b*r*s*t*v + b*s*t*u*w - b*s*w + b*w + r**2*s*t*v
    - s**2*t**2*A*w + 2*s**2*t*u*w - s**2*w - 2*s*t*u*w + s*t*v*w**2 + 2*s*w - w,
    r**2*v + r*s*t**2*A*v - 2*r*s*t*u*v + r*s*v - r*t*A*v - r*v
    - s*t**2*u*A*w + 2*s*t*u**2*w - s*u*w + t*u*A*w + t*u**2*v
    - t*u**2*w + t*v**3 - t*v**2*w - u**2*w + u*w - v**2*w + v*w**2,
]


def main():
    print("polynomial_count", len(polys))
    G = sp.groebner(polys, r, w, s, t, v, u, b, a, order="lex")
    print("basis_len", len(G.polys))
    for g in G.polys[:80]:
        print(sp.factor(g.as_expr()))


if __name__ == "__main__":
    main()
