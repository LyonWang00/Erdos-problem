"""Dependency-free exact backtracking classifier for abstract c(8)=18 types."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
import json
import time

import networkx as nx


POINTS = tuple(range(8))
TRIPLES = tuple(combinations(POINTS, 3))
TRIPLE_INDEX = {triple: i for i, triple in enumerate(TRIPLES)}
PAIRS = tuple(combinations(POINTS, 2))
PAIR_INDEX = {pair: i for i, pair in enumerate(PAIRS)}
BRANCHES = ("C7", "L7", "C6", "L6", "C5", "L5", "C4", "L4", "L3_only")


def bstr(block):
    return "".join(map(str, block))


def bitmask(items, index):
    mask = 0
    for item in items:
        mask |= 1 << index[tuple(item)]
    return mask


@dataclass(frozen=True)
class Candidate:
    kind: str
    block: tuple[int, ...]
    triples: int
    pairs: int
    weight: int


def candidate(kind, block):
    triples = bitmask(combinations(block, 3), TRIPLE_INDEX)
    pairs = bitmask(combinations(block, 2), PAIR_INDEX) if kind == "L" else 0
    weight = len(tuple(combinations(block, 3))) - (kind == "C")
    return Candidate(kind, block, triples, pairs, weight)


def incidence_graph(lines, circles):
    graph = nx.Graph()
    for point in POINTS:
        graph.add_node(("P", point), kind="P")
    for index, block in enumerate(lines):
        node = ("L", index)
        graph.add_node(node, kind="L" + str(len(block)))
        graph.add_edges_from((("P", point), node) for point in block)
    for index, block in enumerate(circles):
        node = ("C", index)
        graph.add_node(node, kind="C" + str(len(block)))
        graph.add_edges_from((("P", point), node) for point in block)
    return graph


NODE_MATCH = nx.algorithms.isomorphism.categorical_node_match("kind", None)


def branch_candidates(branch):
    if branch == "L3_only":
        return None, []
    kind = branch[0]
    size = int(branch[1])
    fixed = candidate(kind, tuple(range(size)))
    out = []
    for k in range(4, size + 1):
        for block in combinations(POINTS, k):
            for this_kind in ("L", "C"):
                if k == size and kind == "L" and this_kind == "C":
                    continue
                item = candidate(this_kind, block)
                if item == fixed:
                    continue
                if item.triples & fixed.triples:
                    continue
                if item.kind == "L" and fixed.kind == "L" and item.pairs & fixed.pairs:
                    continue
                out.append(item)
    # Large weights first reaches the admissible 29..38 window quickly.
    out.sort(key=lambda item: (-item.weight, item.kind, item.block))
    return fixed, out


def available_line3(used_triples, used_line_pairs):
    out = []
    for triple in TRIPLES:
        triple_bit = 1 << TRIPLE_INDEX[triple]
        pair_mask = bitmask(combinations(triple, 2), PAIR_INDEX)
        if not (triple_bit & used_triples) and not (pair_mask & used_line_pairs):
            out.append((triple, triple_bit, pair_mask))
    return out


def enumerate_line3(items, need, start, used_pairs, chosen, emit):
    if need == 0:
        emit(tuple(chosen))
        return
    if len(items) - start < need:
        return
    free_pairs = 28 - used_pairs.bit_count()
    if free_pairs < 3 * need:
        return
    for i in range(start, len(items) - need + 1):
        triple, _triple_bit, pair_mask = items[i]
        if pair_mask & used_pairs:
            continue
        chosen.append(triple)
        enumerate_line3(items, need - 1, i + 1, used_pairs | pair_mask, chosen, emit)
        chosen.pop()


def classify(branch, time_limit, solution_limit):
    fixed, candidates = branch_candidates(branch)
    started = time.time()
    deadline = started + time_limit
    records = []
    graph_buckets = {}
    labelled_solutions = 0
    nodes = 0
    timed_out = False

    if fixed is None:
        initial_selected = []
        initial_triples = 0
        initial_pairs = 0
        initial_weight = 0
    else:
        initial_selected = [fixed]
        initial_triples = fixed.triples
        initial_pairs = fixed.pairs
        initial_weight = fixed.weight

    def add_solution(selected, line3):
        nonlocal labelled_solutions, timed_out
        labelled_solutions += 1
        lines = tuple(sorted(
            [item.block for item in selected if item.kind == "L"] + list(line3),
            key=lambda block: (len(block), block),
        ))
        circles = tuple(sorted(
            [item.block for item in selected if item.kind == "C"],
            key=lambda block: (len(block), block),
        ))
        graph = incidence_graph(lines, circles)
        graph_hash = nx.weisfeiler_lehman_graph_hash(graph, node_attr="kind")
        is_old = any(
            nx.is_isomorphic(graph, old_graph, node_match=NODE_MATCH)
            for old_graph in graph_buckets.get(graph_hash, [])
        )
        if not is_old:
            graph_buckets.setdefault(graph_hash, []).append(graph)
            owned = 0
            for item in selected:
                owned |= item.triples
            for triple in line3:
                owned |= 1 << TRIPLE_INDEX[triple]
            circle3 = [triple for triple in TRIPLES if not (owned >> TRIPLE_INDEX[triple]) & 1]
            record = {
                "lines": [bstr(block) for block in lines],
                "circles3": [bstr(block) for block in circle3],
                "circles_large": [bstr(block) for block in circles],
            }
            records.append(record)
            print(
                branch, "new_class", len(records), "labelled", labelled_solutions,
                "lines", record["lines"],
                "large_circles", record["circles_large"], flush=True,
            )
        if solution_limit and labelled_solutions >= solution_limit:
            timed_out = True

    def visit(selected, start, used_triples, used_line_pairs, weight):
        nonlocal nodes, timed_out
        nodes += 1
        if timed_out:
            return
        if nodes % 10000 == 0 and time.time() >= deadline:
            timed_out = True
            return
        need_line3 = 38 - weight
        if 0 <= need_line3 <= 9:
            items = available_line3(used_triples, used_line_pairs)
            if len(items) >= need_line3 and 28 - used_line_pairs.bit_count() >= 3 * need_line3:
                enumerate_line3(
                    items, need_line3, 0, used_line_pairs, [],
                    lambda line3: add_solution(selected, line3),
                )
                if timed_out:
                    return
        if weight >= 38:
            return
        for i in range(start, len(candidates)):
            item = candidates[i]
            new_weight = weight + item.weight
            if new_weight > 38:
                continue
            if item.triples & used_triples:
                continue
            if item.kind == "L" and item.pairs & used_line_pairs:
                continue
            selected.append(item)
            visit(
                selected,
                i + 1,
                used_triples | item.triples,
                used_line_pairs | item.pairs,
                new_weight,
            )
            selected.pop()
            if timed_out:
                return

    visit(initial_selected, 0, initial_triples, initial_pairs, initial_weight)
    return {
        "branch": branch,
        "complete": not timed_out,
        "elapsed_seconds": time.time() - started,
        "nodes": nodes,
        "labelled_solutions": labelled_solutions,
        "orbit_count": len(records),
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("branch", choices=BRANCHES)
    parser.add_argument("--time-limit", type=float, default=1800)
    parser.add_argument("--solution-limit", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = classify(args.branch, args.time_limit, args.solution_limit)
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
    print("FINAL", result)


if __name__ == "__main__":
    main()
