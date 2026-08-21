"""Construct and verify exact representatives of the three realizable c=18 classes."""

from __future__ import annotations

import json
from pathlib import Path

import networkx as nx
import sympy as sp

from family_examples import invert_point, n8_c17_family, n8_c18_rectangles
from orbit_tools import incidence_summary, maximal_incidence
from realize_c18_by_recoloring import NODE_MATCH, record_graph


BASE = Path(__file__).resolve().parent


def line_intersection(first, second, third, fourth):
    x1, y1 = first
    x2, y2 = second
    x3, y3 = third
    x4, y4 = fourth
    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    return (
        sp.factor(((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denominator),
        sp.factor(((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denominator),
    )


def graph_from_points(points):
    lines, circles = maximal_incidence(points)
    graph = nx.Graph()
    for point in range(8):
        graph.add_node(("P", point), kind="P")
    index = 0
    for kind, blocks in (("L", lines), ("C", circles)):
        for block in sorted(blocks, key=lambda b: (len(b), sorted(b))):
            node = (kind, index)
            index += 1
            graph.add_node(node, kind=kind + str(len(block)))
            graph.add_edges_from((("P", point), node) for point in block)
    return graph


def fmt(value):
    return sp.sstr(sp.factor(value))


def main():
    with (BASE / "c18_backtrack_C4.json").open(encoding="utf-8") as stream:
        records = json.load(stream)["records"]
    targets = [record_graph(record) for record in records]

    rectangle = n8_c18_rectangles(sp.Integer(1), sp.Integer(2))
    class3 = rectangle

    class11_center = (sp.Integer(-4), sp.Integer(-4))
    class11 = [invert_point(point, class11_center) for point in rectangle]

    # b=1 gives a much smaller exact certificate than the published b=13/12
    # representative while retaining the same generalized-circle incidence.
    c17 = n8_c17_family(sp.Integer(1))
    class15_center = line_intersection(c17[0], c17[4], c17[3], c17[5])
    class15 = [invert_point(point, class15_center) for point in c17]

    payload = {}
    for expected, points, center in (
        (3, class3, "infinity"),
        (11, class11, class11_center),
        (15, class15, class15_center),
    ):
        summary = incidence_summary(points)
        assert summary["circle_count"] == 18
        graph = graph_from_points(points)
        matches = [
            index for index, target in enumerate(targets, 1)
            if nx.is_isomorphic(graph, target, node_match=NODE_MATCH)
        ]
        assert matches == [expected], (expected, matches, summary)
        payload[str(expected)] = {
            "source_inversion_center": center if isinstance(center, str) else [fmt(v) for v in center],
            "points": [[fmt(x), fmt(y)] for x, y in points],
            "incidence": summary,
        }
        print("class", expected, "center", center)
        print("points", payload[str(expected)]["points"])
        print("incidence", summary)

    with (BASE / "c18_realizable_representatives.json").open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2)
    print("C18 REALIZABLE REPRESENTATIVES PASS")


if __name__ == "__main__":
    main()
