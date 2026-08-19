"""Verify the seven-point construction with exactly eleven circles."""

from itertools import combinations
import sympy as sp


P = (
    (sp.Rational(0), sp.Rational(0)),
    (sp.Rational(1), sp.Rational(0)),
    (sp.Rational(1), sp.Rational(-1)),
    (sp.Rational(0), sp.Rational(1)),
    (sp.Rational(2, 3), sp.Rational(-1, 3)),
    (sp.Rational(2), sp.Rational(-1)),
    (sp.Rational(0), sp.Rational(-1)),
)
k = sp.Rational(1, 2)
h = sp.Rational(1)


def q(x, y):
    return x*x + 2*k*x*y + h*y*y


def normalized_line(triple):
    matrix = sp.Matrix([[P[i][0], P[i][1], 1] for i in triple])
    if matrix.det() != 0:
        return None
    vector = list(matrix.nullspace()[0])
    pivot = next(value for value in vector if value != 0)
    return tuple(sp.factor(value / pivot) for value in vector)


def circle_key(triple):
    A, B, C = sp.symbols("A B C")
    solution = sp.solve(
        [q(*P[i]) + A*P[i][0] + B*P[i][1] + C for i in triple],
        (A, B, C), dict=True,
    )
    assert len(solution) == 1
    return tuple(sp.factor(solution[0][variable]) for variable in (A, B, C))


lines = {}
circles = {}
for triple in combinations(range(7), 3):
    line = normalized_line(triple)
    if line is not None:
        lines.setdefault(line, []).append(triple)
    else:
        circles.setdefault(circle_key(triple), []).append(triple)

expected_lines = {
    (0, 3, 6), (0, 4, 5), (1, 3, 5),
    (1, 4, 6), (2, 3, 4), (2, 5, 6),
}
expected_four_circles = {
    (0, 1, 3, 4), (0, 1, 5, 6), (0, 2, 3, 5),
    (0, 2, 4, 6), (1, 2, 3, 6), (1, 2, 4, 5),
}

assert h - k*k > 0
assert len(set(P)) == 7
assert {triples[0] for triples in lines.values()} == expected_lines
four_circle_blocks = {
    tuple(sorted(set().union(*map(set, triples))))
    for triples in circles.values() if len(triples) == 4
}
assert four_circle_blocks == expected_four_circles
assert sorted(map(len, circles.values())) == [1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4]
assert len(circles) == 11

# The positive-definite form Q is the pullback of the standard Euclidean
# norm under (x,y) -> (x+y/2, sqrt(3)y/2).  Recount once more using the
# standard Euclidean circle determinant, so the certificate does not rely on
# an implicit change-of-metric argument.
sqrt3 = sp.sqrt(3)
euclidean_points = tuple(
    (sp.factor(x + y/2), sp.factor(sqrt3*y/2)) for x, y in P
)
euclidean_circles = {}
euclidean_lines = []
for triple in combinations(range(7), 3):
    matrix = sp.Matrix([
        [euclidean_points[i][0], euclidean_points[i][1], 1]
        for i in triple
    ])
    if sp.simplify(matrix.det()) == 0:
        euclidean_lines.append(triple)
        continue
    A, B, C = sp.symbols("Ae Be Ce")
    solution = sp.solve([
        x*x + y*y + A*x + B*y + C
        for i in triple for x, y in (euclidean_points[i],)
    ], (A, B, C), dict=True)
    assert len(solution) == 1
    key = tuple(sp.simplify(solution[0][z]) for z in (A, B, C))
    euclidean_circles.setdefault(key, []).append(triple)
assert set(euclidean_lines) == expected_lines
assert sorted(map(len, euclidean_circles.values())) == [1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4]
assert len(euclidean_circles) == 11

print("CONSTRUCTION PASS: 7 distinct points, 6 collinear triples,")
print("6 four-point circles, 5 ordinary circles, and c(P)=11.")
print("Euclidean coordinates are (x+y/2, sqrt(3)*y/2) for the listed rational pairs.")
