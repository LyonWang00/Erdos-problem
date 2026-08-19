from fractions import Fraction
from itertools import combinations


def F(a, b=1):
    return Fraction(a, b)


def add(p, q):
    return (p[0] + q[0], p[1] + q[1])


def sub(p, q):
    return (p[0] - q[0], p[1] - q[1])


def mul(c, p):
    return (c * p[0], c * p[1])


def dot(p, q):
    return p[0] * q[0] + p[1] * q[1]


def line_det(points):
    (x1, y1), (x2, y2), (x3, y3) = points
    return x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - y2 * x3


def det4(mat):
    total = F(0)
    for j in range(4):
        submat = [row[:j] + row[j + 1 :] for row in mat[1:]]
        a, b, c = submat[0]
        d, e, f = submat[1]
        g, h, i = submat[2]
        minor = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
        total += ((-1) ** j) * mat[0][j] * minor
    return total


def circle_det(points):
    return det4([[x * x + y * y, x, y, F(1)] for x, y in points])


def line_intersection(p1, p2, p3, p4):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / den
    py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / den
    return (px, py)


def invert_point(p, center):
    d = sub(p, center)
    return add(center, mul(F(1, 1) / dot(d, d), d))


def circle_equation(points):
    # Return D,E,F for x^2+y^2+D x+E y+F=0 through three non-collinear points.
    rows = []
    rhs = []
    for x, y in points[:3]:
        rows.append([x, y, F(1)])
        rhs.append(-(x * x + y * y))
    (a, b, c), (d, e, f), (g, h, i) = rows
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    out = []
    for col in range(3):
        mat = [row[:] for row in rows]
        for r in range(3):
            mat[r][col] = rhs[r]
        (a, b, c), (d, e, f), (g, h, i) = mat
        out.append((a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)) / det)
    return tuple(out)


def fmtq(q):
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def fmtpt(p):
    return f"({fmtq(p[0])}, {fmtq(p[1])})"


def bstr(block):
    return "".join(map(str, block))


def main():
    b = F(13, 12)
    t = F(-263, 1465)
    m = F(432, 263)
    center = (F(-25, 313), F(312, 313))
    r = m * (1 + b * b)

    P = {}
    P[0] = (F(0), F(0))
    P[1] = (F(1), F(0))
    P[6] = (r, F(0))
    P[2] = (F(1), b)
    P[3] = (F(1) + t * (r - 1), b * (1 - t))
    P[4] = line_intersection(P[0], P[3], P[1], P[2])
    P[5] = (m, m * b)

    Q = {7: center}
    for i in range(7):
        Q[i] = invert_point(P[i], center)

    print("parameters")
    print("b", fmtq(b), "t", fmtq(t), "m", fmtq(m), "center", fmtpt(center))
    print("preimage_points")
    for i in range(7):
        print(i, fmtpt(P[i]))
    print("points")
    for i in range(8):
        print(i, fmtpt(Q[i]))

    line_blocks = []
    for comb in combinations(range(8), 3):
        if line_det([Q[i] for i in comb]) == 0:
            line_blocks.append(comb)

    circle_map = {}
    for comb in combinations(range(8), 3):
        if line_det([Q[i] for i in comb]) == 0:
            continue
        eq = circle_equation([Q[i] for i in comb])
        circle_map.setdefault(eq, set()).update(comb)
    circle_blocks = sorted(tuple(sorted(v)) for v in circle_map.values() if len(v) >= 3)

    print("line_blocks", " ".join(bstr(b) for b in line_blocks))
    print("circle_blocks", " ".join(bstr(b) for b in circle_blocks))
    print("circle_count", len(circle_blocks))
    print("circle_equations")
    for block in circle_blocks:
        eq = circle_equation([Q[i] for i in block[:3]])
        print(bstr(block), "D", fmtq(eq[0]), "E", fmtq(eq[1]), "F", fmtq(eq[2]))

    expected_lines = {"012", "013", "023", "046", "123", "357"}
    expected_circles = {
        "0145",
        "0167",
        "024",
        "0257",
        "026",
        "0347",
        "0356",
        "1247",
        "1256",
        "1346",
        "135",
        "137",
        "157",
        "2345",
        "2367",
        "246",
        "4567",
    }
    assert {bstr(b) for b in line_blocks} == expected_lines
    assert {bstr(b) for b in circle_blocks} == expected_circles
    assert len(circle_blocks) == 17
    covered = set()
    for block in line_blocks:
        for triple in combinations(block, 3):
            covered.add(triple)
    for block in circle_blocks:
        for triple in combinations(block, 3):
            if triple in covered:
                raise AssertionError(("duplicate triple", triple, block))
            covered.add(triple)
    assert len(covered) == 56
    print("covered_triples", len(covered))


if __name__ == "__main__":
    main()
