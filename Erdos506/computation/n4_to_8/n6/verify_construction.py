import itertools
import sympy as sp


sqrt3 = sp.sqrt(3)

# A six-point configuration over Q(sqrt(3)).
# The points are a scaled version of the coordinates used in the revised proof.
POINTS = [
    (sp.Integer(7), sp.Integer(0)),
    (sp.Integer(1), sp.Integer(0)),
    (sp.Integer(4), 2 * sqrt3),
    (sp.Integer(1), sqrt3 / 2),
    (sp.Integer(1), 4 * sqrt3),
    (sp.Rational(1, 7), 4 * sqrt3 / 7),
]


def det3(rows):
    return sp.factor(sp.Matrix(rows).det())


def line_det(triple):
    return det3([[POINTS[i][0], POINTS[i][1], 1] for i in triple])


def circle_equation(triple):
    """Return (D,E,F) for x^2+y^2+D*x+E*y+F=0, or None if collinear."""
    if line_det(triple) == 0:
        return None
    D, E, F = sp.symbols("D E F")
    equations = []
    for i in triple:
        x, y = POINTS[i]
        equations.append(sp.Eq(x * x + y * y + D * x + E * y + F, 0))
    sol = sp.solve(equations, (D, E, F), dict=True)[0]
    return tuple(sp.factor(sol[v]) for v in (D, E, F))


def points_on_circle(eq):
    D, E, F = eq
    block = []
    for i, (x, y) in enumerate(POINTS):
        if sp.factor(x * x + y * y + D * x + E * y + F) == 0:
            block.append(i)
    return tuple(block)


def fmt(block):
    return "".join(str(i) for i in block)


def main():
    print("points:")
    for i, p in enumerate(POINTS):
        print(f"  p{i}={p}")

    line_blocks = []
    for triple in itertools.combinations(range(6), 3):
        if line_det(triple) == 0:
            line_blocks.append(triple)
    print(f"line_blocks={[fmt(b) for b in line_blocks]}")

    circles = {}
    for triple in itertools.combinations(range(6), 3):
        eq = circle_equation(triple)
        if eq is not None:
            circles.setdefault(eq, set()).add(triple)

    circle_blocks = []
    print(f"circle_count={len(circles)}")
    for eq in sorted(circles, key=lambda item: points_on_circle(item)):
        block = points_on_circle(eq)
        circle_blocks.append(block)
        print(f"  block={fmt(block)} equation={eq}")

    covered = set()
    for block in circle_blocks:
        for triple in itertools.combinations(block, 3):
            if line_det(triple) != 0:
                covered.add(triple)
    noncollinear = {
        triple
        for triple in itertools.combinations(range(6), 3)
        if line_det(triple) != 0
    }
    print(f"noncollinear_triples={len(noncollinear)}")
    print(f"covered_noncollinear_triples={len(covered)}")
    print(f"coverage_complete={covered == noncollinear}")
    print(f"all_on_one_circle={len(max(circle_blocks, key=len)) == 6}")
    print(f"max_circle_size={max(len(b) for b in circle_blocks)}")
    print(f"max_line_size=3")


if __name__ == "__main__":
    main()
