"""Find c=18 coloured classes obtained by exact Mobius recolouring of known models."""

from __future__ import annotations

from itertools import combinations
import json
from pathlib import Path

import networkx as nx


BASE = Path(__file__).resolve().parent
import sympy as sp

from family_examples import n8_c17_family, n8_c18_rectangles
from orbit_tools import maximal_incidence


x, y = sp.symbols("x y", real=True)
NODE_MATCH = nx.algorithms.isomorphism.categorical_node_match("kind", None)


def line_equation(points, block):
    p, q = [points[i] for i in sorted(block)[:2]]
    return sp.factor((p[1] - q[1]) * x + (q[0] - p[0]) * y + p[0] * q[1] - q[0] * p[1])


def circle_equation(points, block):
    indices = sorted(block)[:3]
    d, e, f = sp.symbols("d e f")
    solution = sp.solve(
        [points[i][0] ** 2 + points[i][1] ** 2 + d * points[i][0] + e * points[i][1] + f for i in indices],
        (d, e, f), dict=True,
    )[0]
    return sp.factor(x**2 + y**2 + solution[d] * x + solution[e] * y + solution[f])


def object_equations(points):
    lines, circles = maximal_incidence(points)
    objects = []
    for block in sorted(lines, key=lambda b: (len(b), sorted(b))):
        objects.append((block, "L", line_equation(points, block)))
    for block in sorted(circles, key=lambda b: (len(b), sorted(b))):
        objects.append((block, "C", circle_equation(points, block)))
    return objects


def real_intersections(first, second):
    solutions = sp.solve_poly_system([first, second], x, y)
    out = []
    for px, py in solutions or []:
        if sp.simplify(sp.im(px)) == 0 and sp.simplify(sp.im(py)) == 0:
            out.append((sp.simplify(sp.re(px)), sp.simplify(sp.re(py))))
    return out


def same_point(first, second):
    return sp.simplify(first[0] - second[0]) == 0 and sp.simplify(first[1] - second[1]) == 0


def colored_graph(blocks, line_indices):
    graph = nx.Graph()
    for point in range(8):
        graph.add_node(("P", point), kind="P")
    for index, block in enumerate(blocks):
        kind = "L" if index in line_indices else "C"
        node = (kind, index)
        graph.add_node(node, kind=kind + str(len(block)))
        graph.add_edges_from((("P", point), node) for point in block)
    return graph


def record_graph(record):
    graph = nx.Graph()
    for point in range(8):
        graph.add_node(("P", point), kind="P")
    index = 0
    for kind, key in (("L", "lines"), ("C", "circles3"), ("C", "circles_large")):
        for text in record[key]:
            block = tuple(map(int, text))
            node = (kind, index)
            index += 1
            graph.add_node(node, kind=kind + str(len(block)))
            graph.add_edges_from((("P", point), node) for point in block)
    return graph


def match_classes(graph, record_graphs):
    return [
        index for index, target in enumerate(record_graphs, 1)
        if nx.is_isomorphic(graph, target, node_match=NODE_MATCH)
    ]


def analyze(name, points, record_graphs):
    objects = object_equations(points)
    blocks = [entry[0] for entry in objects]
    original_lines = {i for i, entry in enumerate(objects) if entry[1] == "L"}
    seen_centers = []
    realizations = []

    original_graph = colored_graph(blocks, original_lines)
    realizations.append({
        "center": "infinity",
        "line_blocks": ["".join(map(str, sorted(blocks[i]))) for i in sorted(original_lines)],
        "classes": match_classes(original_graph, record_graphs),
    })

    for i, j in combinations(range(len(objects)), 2):
        for center in real_intersections(objects[i][2], objects[j][2]):
            if any(same_point(center, point) for point in points):
                continue
            if any(same_point(center, old) for old in seen_centers):
                continue
            seen_centers.append(center)
            containing = {
                k for k, (_block, _kind, equation) in enumerate(objects)
                if sp.simplify(equation.subs({x: center[0], y: center[1]})) == 0
            }
            if len(containing) != 2:
                continue
            graph = colored_graph(blocks, containing)
            classes = match_classes(graph, record_graphs)
            if classes:
                realizations.append({
                    "center": [sp.sstr(center[0]), sp.sstr(center[1])],
                    "line_blocks": ["".join(map(str, sorted(blocks[k]))) for k in sorted(containing)],
                    "classes": classes,
                })

    print(name, "objects", len(objects), "candidate_centers", len(seen_centers))
    print(name, "classes_seen", sorted({c for rec in realizations for c in rec["classes"]}))
    for record in realizations:
        print(json.dumps(record, sort_keys=True))
    return realizations


def main():
    with (BASE / "c18_backtrack_C4.json").open(encoding="utf-8") as stream:
        records = json.load(stream)["records"]
    record_graphs = [record_graph(record) for record in records]
    payload = {
        "rectangle": analyze("rectangle", n8_c18_rectangles(sp.Rational(1), sp.Rational(2)), record_graphs),
        "c17_family": analyze("c17_family", n8_c17_family(sp.Rational(13, 12)), record_graphs),
    }
    with (BASE / "c18_recoloring_realizations.json").open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2)


if __name__ == "__main__":
    main()
