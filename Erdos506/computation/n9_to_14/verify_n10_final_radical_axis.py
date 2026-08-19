#!/usr/bin/env python3
"""Audited closure of the last corrected ten-point signature.

This verifier has two deliberately small parts.

1.  It enumerates the two symmetry classes of the five K_{2,2} components
    forced by the radical axis and checks the remaining twenty triples.  One
    class has no possible ten-line selection; the other has two complementary
    selections forming one orbit.
2.  For a representative of that orbit it reconstructs a complete projective
    parametrization from the ten three-point lines.  Five selected 3-by-3
    minors of the chord-intersection matrix give an elementary contradiction.

No saved solver verdict or old symbolic certificate is read.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

import sympy as sp

from verify_recursive_inversion_rows import candidate_rows, INFEASIBLE


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "certificates" / "n10_08_final_geometry.json"
RECURSIVE_CERTIFICATE = HERE / "certificates" / "n10_07_recursive.json"
EXPECTED_SIGNATURE = (10, 0, 0, 0, 10, 20, 2, 0, 15)
A = tuple(range(5))
B = tuple(range(5, 10))

# Canonical representative produced by the small radical-axis support model.
LINES3 = (
    (0, 3, 9), (0, 4, 5), (0, 6, 7), (1, 2, 8), (1, 3, 6),
    (1, 7, 9), (2, 4, 7), (2, 5, 6), (3, 5, 8), (4, 8, 9),
)
CIRCLES4 = (
    (0, 1, 5, 8), (0, 1, 6, 9), (0, 2, 5, 7), (0, 2, 8, 9),
    (0, 3, 5, 6), (0, 3, 7, 8), (0, 4, 6, 8), (0, 4, 7, 9),
    (1, 2, 5, 9), (1, 2, 6, 7), (1, 3, 5, 7), (1, 3, 8, 9),
    (1, 4, 5, 6), (1, 4, 7, 8), (2, 3, 6, 8), (2, 3, 7, 9),
    (2, 4, 5, 8), (2, 4, 6, 9), (3, 4, 5, 9), (3, 4, 6, 7),
)


def recursive_partition_certificate():
    """Connect the count-vector reduction to the geometric final case.

    Every realizable point row must have a realizable nine-point child.  The
    recursive certificate leaves four numerical rows, all with five-circle
    degree one and with three-line degree at most three.  Since the signature
    has ten three-lines, their total point incidence is thirty; hence every
    point has three-line degree exactly three.
    """
    document = json.loads(RECURSIVE_CERTIFICATE.read_text(encoding="utf-8"))
    top_records = [
        record for record in document["levels"]["10"]
        if record["status"] not in INFEASIBLE
    ]
    if ([tuple(record["signature"]) for record in top_records]
            != [EXPECTED_SIGNATURE]):
        raise RuntimeError("unexpected recursive ten-point residual")
    child_status = {
        tuple(record["signature"]): record["status"]
        for record in document["levels"]["9"]
    }
    rows, children = candidate_rows(10, EXPECTED_SIGNATURE)
    retained = tuple(
        row for row, child in zip(rows, children)
        if child_status.get(child, "UNKNOWN") not in INFEASIBLE
    )
    if len(retained) != 4:
        raise RuntimeError(("retained point-row count", len(retained)))
    five_circle_degrees = sorted({row[6] for row in retained})
    three_line_degrees = sorted({row[0] for row in retained})
    if five_circle_degrees != [1] or max(three_line_degrees) != 3:
        raise RuntimeError((five_circle_degrees, three_line_degrees))
    balanced_rows = tuple(row for row in retained if row[0] == 3)
    if balanced_rows != ((3, 0, 0, 0, 3, 8, 1, 0),):
        raise RuntimeError(("balanced point row", balanced_rows))
    return {
        "source_file": RECURSIVE_CERTIFICATE.name,
        "source_sha256": sha256(RECURSIVE_CERTIFICATE.read_bytes()).hexdigest(),
        "retained_point_row_count": len(retained),
        "five_circle_degree_set": five_circle_degrees,
        "three_line_degree_set_before_balance": three_line_degrees,
        "forced_balanced_point_row": list(balanced_rows[0]),
        "consequences": [
            "the two five-point circles partition the ten points",
            "every point lies on three selected three-point lines",
            "every point lies on eight selected four-point circles",
        ],
    }


def edge(block, side):
    return tuple(x if side == "A" else x - 5 for x in block
                 if (x < 5) == (side == "A"))


def radical_components(circles4):
    graph = {}
    for block in circles4:
        left = ("A", edge(block, "A"))
        right = ("B", edge(block, "B"))
        graph.setdefault(left, set()).add(right)
        graph.setdefault(right, set()).add(left)
    seen = set()
    answer = []
    for start in sorted(graph):
        if start in seen:
            continue
        stack, component = [start], []
        seen.add(start)
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbour in graph[vertex]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        left = tuple(sorted(v[1] for v in component if v[0] == "A"))
        right = tuple(sorted(v[1] for v in component if v[0] == "B"))
        answer.append((left, right))
    return tuple(sorted(answer))


def image_edge(pair, permutation):
    return tuple(sorted(permutation[x] for x in pair))


def matching_stabilizer(matching):
    target = frozenset(frozenset(part) for part in matching)
    answer = []
    for permutation in permutations(range(5)):
        image = frozenset(
            frozenset(image_edge(edge_value, permutation) for edge_value in part)
            for part in matching
        )
        if image == target:
            answer.append(permutation)
    return tuple(answer)


def component_action(matching, permutation):
    index = {frozenset(part): i for i, part in enumerate(matching)}
    return tuple(index[frozenset(
        image_edge(edge_value, permutation) for edge_value in part
    )] for part in matching)


def compose(first, second):
    return tuple(first[second[i]] for i in range(5))


def inverse(permutation):
    answer = [0] * 5
    for i, value in enumerate(permutation):
        answer[value] = i
    return tuple(answer)


def four_blocks(left_matching, right_matching, sigma):
    return frozenset(
        tuple(sorted(
            left_edge + tuple(x + 5 for x in right_edge)
        ))
        for left_index, left_part in enumerate(left_matching)
        for left_edge in left_part
        for right_edge in right_matching[sigma[left_index]]
    )


def residual_line_designs(blocks):
    triples = tuple(combinations(range(10), 3))
    residual = tuple(
        triple for triple in triples
        if not set(triple) <= set(A)
        and not set(triple) <= set(B)
        and not any(set(triple) <= set(block) for block in blocks)
    )
    designs = []
    for selected in combinations(residual, 10):
        if any(sum(point in triple for triple in selected) != 3
               for point in range(10)):
            continue
        if any(sum(set(pair) <= set(triple) for triple in selected) > 1
               for pair in combinations(range(10), 2)):
            continue
        designs.append(frozenset(selected))
    return residual, tuple(designs)


def transform_blocks(blocks, left_permutation, right_permutation):
    mapping = {x: left_permutation[x] for x in range(5)}
    mapping.update({x + 5: right_permutation[x] + 5 for x in range(5)})
    return frozenset(
        tuple(sorted(mapping[x] for x in block)) for block in blocks
    )


def combinatorial_classification():
    components = radical_components(CIRCLES4)
    if len(components) != 5 or any(
        len(left) != 2 or len(right) != 2 for left, right in components
    ):
        raise RuntimeError("canonical four-circles are not five K_2,2 components")
    left_matching = tuple(left for left, _ in components)
    right_matching = tuple(right for _, right in components)
    if any(set(first) & set(second) for matching in
           (left_matching, right_matching) for first, second in matching):
        raise RuntimeError("a chord pair in one component shares an endpoint")

    # Independently enumerate the perfect matchings of the Petersen graph.
    graph_edges = tuple(combinations(tuple(combinations(range(5), 2)), 2))
    disjoint_pairs = tuple(pair for pair in graph_edges
                           if not set(pair[0]) & set(pair[1]))
    petersen_matchings = set()
    for candidate in combinations(disjoint_pairs, 5):
        used = [edge_value for part in candidate for edge_value in part]
        if len(set(used)) == 10:
            petersen_matchings.add(frozenset(frozenset(x) for x in candidate))

    left_stabilizer = matching_stabilizer(left_matching)
    right_stabilizer = matching_stabilizer(right_matching)
    left_actions = tuple(sorted(set(
        component_action(left_matching, p) for p in left_stabilizer
    )))
    right_actions = tuple(sorted(set(
        component_action(right_matching, p) for p in right_stabilizer
    )))

    remaining = set(permutations(range(5)))
    identity = tuple(range(5))
    representatives = []
    double_orbits = []
    while remaining:
        representative = identity if identity in remaining else min(remaining)
        orbit = {
            compose(right, compose(representative, inverse(left)))
            for left in left_actions for right in right_actions
        }
        representatives.append(representative)
        double_orbits.append(orbit)
        remaining.difference_update(orbit)

    records = []
    for sigma in representatives:
        blocks = four_blocks(left_matching, right_matching, sigma)
        residual, designs = residual_line_designs(blocks)
        automorphisms = []
        for left in left_stabilizer:
            for right in right_stabilizer:
                if transform_blocks(blocks, left, right) == blocks:
                    automorphisms.append((left, right))
        unseen = set(designs)
        orbits = []
        while unseen:
            design = min(unseen, key=lambda value: tuple(sorted(value)))
            orbit = {
                transform_blocks(design, left, right)
                for left, right in automorphisms
            } & set(designs)
            unseen.difference_update(orbit)
            orbits.append(orbit)
        records.append({
            "component_bijection": list(sigma),
            "four_circle_count": len(blocks),
            "residual_triple_count": len(residual),
            "line_design_count": len(designs),
            "line_design_orbit_count": len(orbits),
            "line_design_orbit_sizes": sorted(map(len, orbits)),
            "contains_canonical_representative": (
                blocks == frozenset(CIRCLES4)
                and frozenset(LINES3) in set(designs)
            ),
        })

    expected = {
        "perfect_matching_count": 6,
        "left_stabilizer_size": 20,
        "right_stabilizer_size": 20,
        "double_orbit_count": 2,
        "double_orbit_sizes": [20, 100],
        "line_design_counts": [2, 0],
        "line_design_orbit_counts": [1, 0],
    }
    actual = {
        "perfect_matching_count": len(petersen_matchings),
        "left_stabilizer_size": len(left_stabilizer),
        "right_stabilizer_size": len(right_stabilizer),
        "double_orbit_count": len(representatives),
        "double_orbit_sizes": sorted(map(len, double_orbits)),
        "line_design_counts": [record["line_design_count"] for record in records],
        "line_design_orbit_counts": [
            record["line_design_orbit_count"] for record in records
        ],
    }
    if actual != expected or not records[0]["contains_canonical_representative"]:
        raise RuntimeError({"expected": expected, "actual": actual, "records": records})
    return {**actual, "records": records}


def primitive_row(points, block):
    left = [x for x in block if x < 5]
    right = [x for x in block if x >= 5]
    raw = points[left[0]].cross(points[left[1]]).cross(
        points[right[0]].cross(points[right[1]])
    )
    common = sp.factor(sp.gcd_list(list(raw)))
    if common not in (0, 1, -1):
        raw = raw.applyfunc(lambda value: sp.factor(value / common))
    return raw.T, common


def algebraic_certificate():
    a, b, c, d, e, f = sp.symbols("a b c d e f")
    initial = {
        0: sp.Matrix((1, 0, 0)),
        1: sp.Matrix((0, 1, 0)),
        2: sp.Matrix((0, 0, 1)),
        3: sp.Matrix((1, 1, 1)),
        9: sp.Matrix((a + 1, a, a)),
        8: sp.Matrix((0, 1, b)),
        6: sp.Matrix((c, c + 1, c)),
        7: sp.Matrix((d * (a + 1), d * a + 1, d * a)),
        4: sp.Matrix((e * d * (a + 1), e * (d * a + 1), e * d * a + 1)),
        5: sp.Matrix((f * c, f * (c + 1), f * c + 1)),
    }
    line_equations = tuple(sorted(set(
        sp.factor(sp.Matrix.hstack(*(initial[x] for x in line)).det())
        for line in LINES3
        if sp.factor(sp.Matrix.hstack(*(initial[x] for x in line)).det()) != 0
    ), key=str))
    expected_line_equations = {
        -c + d * a,
        f * b - 1,
        (a + 1) * (e * b - 1),
        (e * f * c - e * f * d * a + e * d * a + e - f * c - f),
    }
    if set(line_equations) != expected_line_equations:
        raise RuntimeError((line_equations, expected_line_equations))

    # c=da and f=1/b.  If eb != 1, the last two equations force a=-1
    # and da=-1; then point 7 equals point 2.  Thus every nondegenerate
    # realization belongs to the branch e=f=1/b.
    points = {
        0: sp.Matrix((1, 0, 0)),
        1: sp.Matrix((0, 1, 0)),
        2: sp.Matrix((0, 0, 1)),
        3: sp.Matrix((1, 1, 1)),
        9: sp.Matrix((a + 1, a, a)),
        8: sp.Matrix((0, 1, b)),
        6: sp.Matrix((d * a, d * a + 1, d * a)),
        7: sp.Matrix((d * (a + 1), d * a + 1, d * a)),
        4: sp.Matrix((d * (a + 1), d * a + 1, d * a + b)),
        5: sp.Matrix((d * a, d * a + 1, d * a + b)),
    }
    if any(sp.factor(sp.Matrix.hstack(*(points[x] for x in line)).det()) != 0
           for line in LINES3):
        raise RuntimeError("specialized parametrization misses a line equation")

    rows, common_factors = zip(*(
        primitive_row(points, block) for block in CIRCLES4
    ))
    selected = {
        (0, 4, 14): a * d + b**2,
        (0, 14, 15): -a * b + a * d - b,
        (0, 4, 15): a * b - a * d - a - 1,
        (0, 15, 16): -a * d + b * d - b - d,
        (0, 11, 15): b * (a * b + 1),
    }
    verified_minors = []
    for indices, expected in selected.items():
        determinant = sp.factor(sp.Matrix.vstack(*(rows[i] for i in indices)).det())
        if sp.factor(determinant - expected) != 0:
            raise RuntimeError((indices, determinant, expected))
        verified_minors.append({
            "row_indices": list(indices),
            "circle_blocks": [list(CIRCLES4[i]) for i in indices],
            "polynomial": str(expected),
        })

    # Only rows 0,4,11,14,15,16 occur above.  Their removed common factors
    # are monomials in a,b,d, all nonzero in a valid realization.
    selected_row_factors = {
        index: sp.factor(common_factors[index])
        for index in (0, 4, 11, 14, 15, 16)
    }
    if selected_row_factors != {
        0: a * d, 4: b, 11: 1, 14: 1, 15: 1, 16: a * d,
    }:
        raise RuntimeError(selected_row_factors)

    p1, p2, p3, p4, p5 = (selected[key] for key in selected)
    relation_a = sp.factor(-(p2 + p3))  # a+b+1
    relation_b = sp.factor((p5 / b).subs(a, -b - 1))
    relation_c = sp.factor(
        (b + 1) * p4.subs({a: -b - 1, d: b**2 / (b + 1)})
    )
    univariate_first = sp.Poly(-relation_b, b)       # b^2+b-1
    univariate_second = sp.Poly(relation_c / b, b)  # (b-1)(2b+1)
    gcd_value = sp.gcd(univariate_first, univariate_second)
    if (relation_a != a + b + 1
            or univariate_first.as_expr() != b**2 + b - 1
            or sp.factor(univariate_second.as_expr()) != (b - 1) * (2 * b + 1)
            or gcd_value.as_expr() != 1):
        raise RuntimeError({
            "relation_a": relation_a,
            "first": univariate_first.as_expr(),
            "second": univariate_second.as_expr(),
            "gcd": gcd_value.as_expr(),
        })

    return {
        "projective_frame": [0, 1, 2, 3],
        "line_parameter_count_before_relations": 6,
        "line_equations": [str(value) for value in line_equations],
        "nondegenerate_branch_parameters": ["a", "b", "d"],
        "required_nonzero_parameters": ["a", "b", "d"],
        "selected_row_common_factors": {
            str(index): str(value) for index, value in selected_row_factors.items()
        },
        "selected_minors": verified_minors,
        "univariate_consequences": [
            str(univariate_first.as_expr()),
            str(sp.factor(univariate_second.as_expr())),
        ],
        "univariate_gcd": str(gcd_value.as_expr()),
    }


def main():
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    recursive_partition = recursive_partition_certificate()
    classification = combinatorial_classification()
    certificate = algebraic_certificate()
    document = {
        "schema_version": 1,
        "status": "PASS",
        "claim": (
            "The final corrected n=10 signature is impossible: after the "
            "radical-axis classification there is one abstract line-system "
            "orbit, and five chord-intersection minors contradict its complete "
            "projective parametrization."
        ),
        "signature": list(EXPECTED_SIGNATURE),
        "recursive_partition_certificate": recursive_partition,
        "classification": classification,
        "algebraic_certificate": certificate,
    }
    payload = json.dumps(document, sort_keys=True).encode()
    document["payload_sha256"] = sha256(payload).hexdigest()
    OUTPUT.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": document["status"],
        "double_orbit_count": classification["double_orbit_count"],
        "line_design_counts": classification["line_design_counts"],
        "selected_minor_count": len(certificate["selected_minors"]),
        "univariate_gcd": certificate["univariate_gcd"],
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
