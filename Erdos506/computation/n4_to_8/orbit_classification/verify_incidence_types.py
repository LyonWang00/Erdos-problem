"""Independent consistency checks and symbolic exclusions for the c=18 classification."""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import sympy as sp


BASE = Path(__file__).resolve().parent
BRANCHES = ("C7", "L7", "C6", "L6", "C5", "L5", "C4", "L4", "L3_only")
NODE_MATCH_COLOURED = nx.algorithms.isomorphism.categorical_node_match("kind", None)
NODE_MATCH_UNCOLOURED = NODE_MATCH_COLOURED


def parse(text):
    return tuple(map(int, text))


def triples(block):
    return set(combinations(block, 3))


def graph(record, coloured):
    result = nx.Graph()
    for point in range(8):
        result.add_node(("P", point), kind="P")
    index = 0
    for source_kind, key in (("L", "lines"), ("C", "circles3"), ("C", "circles_large")):
        for text in record[key]:
            block = parse(text)
            kind = source_kind if coloured else "B"
            node = (kind, index)
            index += 1
            result.add_node(node, kind=kind + str(len(block)))
            result.add_edges_from((("P", point), node) for point in block)
    return result


def check_record(record):
    lines = [parse(text) for text in record["lines"]]
    circle3 = [parse(text) for text in record["circles3"]]
    circles_large = [parse(text) for text in record["circles_large"]]
    assert len(circle3) + len(circles_large) == 18
    all_blocks = lines + circles_large
    owned = set()
    for block in all_blocks:
        block_triples = triples(block)
        assert not owned.intersection(block_triples)
        owned.update(block_triples)
    assert set(circle3).isdisjoint(owned)
    assert owned | set(circle3) == set(combinations(range(8), 3))
    for first, second in combinations(lines, 2):
        assert len(set(first).intersection(second)) <= 1


def fano_exclusion(records):
    record = records[0]
    blocks4 = [
        parse(text) for key in ("lines", "circles_large")
        for text in record[key] if len(text) == 4
    ]
    through_zero = {
        tuple(point for point in block if point != 0)
        for block in blocks4 if 0 in block
    }
    expected = {
        (1, 2, 3), (1, 4, 5), (1, 6, 7),
        (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 5, 6),
    }
    assert through_zero == expected

    alpha = sp.symbols("alpha", real=True)
    q3 = sp.Matrix([[1, 1, 0]])
    q5 = sp.Matrix([[1, 0, 1]])
    q6 = sp.Matrix([[0, 1, alpha]])
    determinant_356 = sp.factor(sp.Matrix.vstack(q3, q5, q6).det())
    # Lines 167,257,347 force alpha=1; line 356 forces this determinant to vanish.
    assert determinant_356 == -alpha - 1
    assert determinant_356.subs(alpha, 1) == -2
    print("FANO_EXCLUSION", "alpha_from_347=1", "det356_at_alpha1=-2")


def mobius_kantor_exclusion(records):
    record = records[15]
    expected_lines = {"025", "036", "047", "126", "134", "157", "237", "456"}
    assert set(record["lines"]) == expected_lines
    t = sp.symbols("t", real=True)
    p4 = sp.Matrix([[t, 1, t]])
    p5 = sp.Matrix([[1, 0, t]])
    p6 = sp.Matrix([[0, 1, 1]])
    determinant_456 = sp.factor(sp.Matrix.vstack(p4, p5, p6).det())
    assert determinant_456 == -t**2 + t - 1
    assert sp.discriminant(t**2 - t + 1, t) == -3
    print("MOBIUS_KANTOR_EXCLUSION", "det456", determinant_456, "discriminant=-3")


def main():
    payloads = {}
    for branch in BRANCHES:
        with (BASE / f"c18_backtrack_{branch}.json").open(encoding="utf-8") as stream:
            payloads[branch] = json.load(stream)
        assert payloads[branch]["complete"] is True
    assert all(payloads[branch]["orbit_count"] == 0 for branch in BRANCHES if branch != "C4")
    records = payloads["C4"]["records"]
    assert len(records) == 16
    for record in records:
        check_record(record)

    coloured_graphs = [graph(record, True) for record in records]
    for i, j in combinations(range(16), 2):
        assert not nx.is_isomorphic(
            coloured_graphs[i], coloured_graphs[j], node_match=NODE_MATCH_COLOURED
        )

    signatures = Counter()
    for record in records:
        line3 = sum(len(text) == 3 for text in record["lines"])
        line4 = sum(len(text) == 4 for text in record["lines"])
        circle4 = sum(len(text) == 4 for text in record["circles_large"])
        signatures[(line3, line4, circle4)] += 1
    assert signatures == Counter({(1, 1, 11): 8, (2, 0, 12): 4, (0, 2, 10): 3, (8, 0, 10): 1})
    print("SIGNATURES", dict(signatures))

    groups = []
    for index, record in enumerate(records, 1):
        current = graph(record, False)
        for group in groups:
            if nx.is_isomorphic(current, group["graph"], node_match=NODE_MATCH_UNCOLOURED):
                group["classes"].append(index)
                break
        else:
            groups.append({"classes": [index], "graph": current})
    actual_groups = [group["classes"] for group in groups]
    expected_groups = [
        [1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14],
        [3, 11, 15],
        [16],
    ]
    assert actual_groups == expected_groups
    print("UNCOLOURED_GROUPS", actual_groups)

    fano_exclusion(records)
    mobius_kantor_exclusion(records)

    with (BASE / "c18_realizable_representatives.json").open(encoding="utf-8") as stream:
        representatives = json.load(stream)
    assert sorted(map(int, representatives)) == [3, 11, 15]
    print("REALIZABLE_CLASSES", sorted(map(int, representatives)))
    print("C18 CLASSIFICATION VERIFICATION PASS")


if __name__ == "__main__":
    main()
