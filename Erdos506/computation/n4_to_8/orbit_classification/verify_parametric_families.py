"""Symbolically verify the required incidences in the three parametric families."""

import sympy as sp

from family_examples import n6_family, n7_family, n8_c17_family, n8_c18_rectangles


def line_det(points, block):
    return sp.factor(sp.together(sp.Matrix([[points[i][0], points[i][1], 1] for i in block]).det()))


def circle_det(points, block):
    rows = []
    for i in block:
        x, y = points[i]
        rows.append([x * x + y * y, x, y, 1])
    return sp.factor(sp.together(sp.Matrix(rows).det()))


def check_family(name, points, lines, circles):
    print(name)
    for block in lines:
        value = sp.cancel(line_det(points, block))
        print(" line", "".join(map(str, block)), value)
        assert value == 0
    for block in circles:
        value = sp.cancel(circle_det(points, block))
        print(" circle", "".join(map(str, block)), value)
        assert value == 0


def main():
    a, v = sp.symbols("a v", real=True)
    check_family(
        "n6",
        n6_family(a, v),
        [(0, 2, 4), (0, 3, 5), (1, 3, 4)],
        [(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)],
    )

    u, w = sp.symbols("u w", real=True)
    check_family(
        "n7",
        n7_family(u, w),
        [(0, 3, 6), (0, 4, 5), (1, 3, 5), (1, 4, 6), (2, 3, 4), (2, 5, 6)],
        [(0, 1, 3, 4), (0, 1, 5, 6), (0, 2, 3, 5), (0, 2, 4, 6), (1, 2, 3, 6), (1, 2, 4, 5)],
    )

    b = sp.symbols("b", real=True)
    check_family(
        "n8_c17",
        n8_c17_family(b),
        [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3), (0, 4, 6), (3, 5, 7)],
        [
            (0, 1, 4, 5), (0, 1, 6, 7), (0, 2, 5, 7), (0, 3, 4, 7),
            (0, 3, 5, 6), (1, 2, 4, 7), (1, 2, 5, 6), (1, 3, 4, 6),
            (2, 3, 4, 5), (2, 3, 6, 7), (4, 5, 6, 7),
        ],
    )

    aspect, ratio = sp.symbols("aspect ratio", real=True)
    check_family(
        "n8_c18_rectangles",
        n8_c18_rectangles(aspect, ratio),
        [(0, 3, 4), (0, 3, 7), (0, 4, 7), (3, 4, 7),
         (1, 2, 5), (1, 2, 6), (1, 5, 6), (2, 5, 6)],
        [
            (0, 1, 2, 3), (0, 1, 4, 5), (0, 1, 6, 7), (0, 2, 4, 6),
            (0, 2, 5, 7), (1, 3, 4, 6), (1, 3, 5, 7), (2, 3, 4, 5),
            (2, 3, 6, 7), (4, 5, 6, 7),
        ],
    )
    print("PARAMETRIC INCIDENCE IDENTITIES PASS")


if __name__ == "__main__":
    main()
