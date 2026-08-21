import sympy as sp


a, b, u, v, t, s, r, w = sp.symbols("a b u v t s r w")
A = u**2 + v**2

polys = [
    -b**2*t**2*w + 2*b**2*t*w - b**2*w + 2*b*t**2*u*w - 2*b*t*u*w
    - b*t*w + b*w + r**2*t*v - r*t*v - t**2*A*w + t*u*w + t*v*w**2,
    r**2*s*v - 2*r*s*u*v + 2*r*u*v - r*v + s*A*v - 2*s*v**2*w
    + s*v*w**2 - A*v - u**2*w + u*w,
    -a*r*v + a*u*w + r**2*s*w + r**2*v - 2*r*s*u*w + s*A*w
    - 2*s*v*w**2 + s*w**3 - A*w + v*w**2,
    a*b + b**2*t - b**2 - 2*b*t*u + t*A,
    -b*r*v + b*u*w + r**2*v - A*w + v*w**2,
    b**2*s*t*v - b**2*s*t*w - b**2*s*v + b**2*s*w - b**2*t*v + b**2*v
    - b*r*s*v - 2*b*s*t*u*v + 2*b*s*t*u*w + 2*b*s*u*v - b*s*u*w
    + 2*b*t*u*v - 2*b*u*v + r**2*s**2*v - 2*r*s**2*u*v + 2*r*s*u*v
    + s**2*A*v - 2*s**2*v**2*w + s**2*v*w**2 + s*t*A*v - s*t*A*w
    - 2*s*A*v + 2*s*v**2*w - t*A*v + A*v,
    -a*r*v + a*u*w + a*v - a*w + r**2*v - r*v - A*w + u*w + v*w**2,
    -b*r*v + b*u*w + b*v - b*w + r**2*s*w + r**2*v - 2*r*s*u*w
    - r*v + s*A*w - 2*s*v*w**2 + s*w**3 - u**2*w + u*w - v**2*w + v*w**2,
    a*b*v - a*b*w - a*r*v + a*u*w - b**2*t*w + b**2*w - b*r*v
    + 2*b*t*u*w - b*u*w + r**2*v - t*A*w + v*w**2,
    a*b*v - a*b*w - a*r*v + a*u*w - b*r*v + b*u*w + r**2*s*v
    - 2*r*s*u*v + 2*r*u*v + s*A*v - 2*s*v**2*w + s*v*w**2
    - A*v - u**2*w + v**2*w,
]


def main():
    print("polynomial_count", len(polys))
    G = sp.groebner(polys, r, w, s, t, v, u, b, a, order="lex")
    print("basis_len", len(G.polys))
    for g in G.polys[:120]:
        print(sp.factor(g.as_expr()))


if __name__ == "__main__":
    main()
