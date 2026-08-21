from fractions import Fraction
from itertools import combinations


def F(a, b=1):
    return Fraction(a, b)


POINTS = {
    0: (F(0), F(0)),
    1: (F(263, 626), F(2178, 4069)),
    2: (F(263, 313), F(4356, 4069)),
    3: (F(789, 626), F(6534, 4069)),
    4: (F(53519, 195938), F(1842342, 1273597)),
    5: (F(184032, 458545), F(7245468, 5961085)),
    6: (F(160557, 917090), F(5527026, 5961085)),
    7: (F(-25, 313), F(312, 313)),
}

EXPECTED_LINES = {"0123", "046", "357"}
EXPECTED_CIRCLES = {
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


def det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def det4(m):
    return sum(
        ((-1) ** j) * m[0][j] * det3([row[:j] + row[j + 1 :] for row in m[1:]])
        for j in range(4)
    )


def line_det(a, b, c):
    return det3([[a[0], a[1], F(1)], [b[0], b[1], F(1)], [c[0], c[1], F(1)]])


def circle_det(ps):
    return det4([[x * x + y * y, x, y, F(1)] for x, y in ps])


def bstr(block):
    return "".join(map(str, block))


def maximal_line_blocks():
    out = set()
    for comb in combinations(range(8), 3):
        if line_det(*(POINTS[i] for i in comb)) != 0:
            continue
        s = set(comb)
        a, b = comb[0], comb[1]
        for j in range(8):
            if line_det(POINTS[a], POINTS[b], POINTS[j]) == 0:
                s.add(j)
        out.add(tuple(sorted(s)))
    return sorted(out, key=lambda z: (len(z), z))


def maximal_circle_blocks():
    out = set()
    for comb in combinations(range(8), 3):
        if line_det(*(POINTS[i] for i in comb)) == 0:
            continue
        base = [POINTS[i] for i in comb]
        s = set(comb)
        for j in range(8):
            if j not in s and circle_det(base + [POINTS[j]]) == 0:
                s.add(j)
        out.add(tuple(sorted(s)))
    return sorted(out, key=lambda z: (len(z), z))


def main():
    assert len(set(POINTS.values())) == 8
    lines = maximal_line_blocks()
    circles = maximal_circle_blocks()
    line_strings = {bstr(b) for b in lines}
    circle_strings = {bstr(b) for b in circles}
    print("line_blocks", " ".join(sorted(line_strings)))
    print("circle_blocks", " ".join(sorted(circle_strings, key=lambda s: (len(s), s))))
    print("circle_count", len(circles))
    assert line_strings == EXPECTED_LINES
    assert circle_strings == EXPECTED_CIRCLES
    assert len(circles) == 17
    covered = {}
    for block in lines:
        for triple in combinations(block, 3):
            covered[triple] = "L" + bstr(block)
    for block in circles:
        for triple in combinations(block, 3):
            if triple in covered:
                raise AssertionError(("triple covered twice", triple, covered[triple], "C" + bstr(block)))
            covered[triple] = "C" + bstr(block)
    assert len(covered) == 56
    print("covered_triples", len(covered))
    print("status OK")


if __name__ == "__main__":
    main()
