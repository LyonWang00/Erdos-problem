"""Group the 16 abstract c=18 classes by uncoloured generalized-circle type."""

import json
from pathlib import Path

import networkx as nx

from family_examples import n8_c17_family, n8_c18_rectangles
from orbit_tools import maximal_incidence


BASE = Path(__file__).resolve().parent


NODE_MATCH = nx.algorithms.isomorphism.categorical_node_match("kind", None)


def parse(text):
    return tuple(map(int, text))


def graph_from_blocks(blocks):
    graph = nx.Graph()
    for point in range(8):
        graph.add_node(("P", point), kind="P")
    for index, block in enumerate(blocks):
        node = ("B", index)
        graph.add_node(node, kind="B" + str(len(block)))
        graph.add_edges_from((("P", point), node) for point in block)
    return graph


def graph_from_record(record):
    blocks = [parse(text) for key in ("lines", "circles3", "circles_large") for text in record[key]]
    return graph_from_blocks(blocks)


def graph_from_points(points):
    lines, circles = maximal_incidence(points)
    return graph_from_blocks(sorted(lines | circles, key=lambda block: (len(block), sorted(block))))


def isomorphic(first, second):
    return nx.is_isomorphic(first, second, node_match=NODE_MATCH)


def main():
    with (BASE / "c18_backtrack_C4.json").open(encoding="utf-8") as stream:
        payload = json.load(stream)
    groups = []
    for index, record in enumerate(payload["records"], 1):
        graph = graph_from_record(record)
        for group in groups:
            if isomorphic(graph, group["graph"]):
                group["classes"].append(index)
                break
        else:
            groups.append({"classes": [index], "graph": graph})

    known = {
        "c17_family": graph_from_points(n8_c17_family(__import__("sympy").Rational(13, 12))),
        "rectangle_c18": graph_from_points(n8_c18_rectangles(1, 2)),
    }
    print("uncolored_group_count", len(groups))
    for number, group in enumerate(groups, 1):
        matches = [name for name, graph in known.items() if isomorphic(group["graph"], graph)]
        print("group", number, "classes", group["classes"], "known_matches", matches)


if __name__ == "__main__":
    main()
