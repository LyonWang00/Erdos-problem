import argparse
from itertools import combinations, permutations


PTS = tuple(range(7))
TRIPLES = tuple(combinations(PTS, 3))
FOURS = tuple(combinations(PTS, 4))
LINE_BLOCKS = tuple(combinations(PTS, 3)) + FOURS

PATTERNS = {
    "q5a": [(0, 1, 3, 5), (0, 1, 4, 6), (0, 2, 3, 6), (0, 2, 4, 5), (1, 2, 3, 4)],
    "q5b": [(0, 1, 5, 6), (0, 2, 4, 6), (0, 3, 4, 5), (1, 2, 3, 6), (1, 2, 4, 5)],
    "q6": [(0, 1, 3, 4), (0, 1, 5, 6), (0, 2, 3, 5), (0, 2, 4, 6), (1, 2, 3, 6), (1, 2, 4, 5)],
    "q7": [(0, 1, 2, 3), (0, 1, 4, 5), (0, 2, 4, 6), (0, 3, 5, 6), (1, 2, 5, 6), (1, 3, 4, 6), (2, 3, 4, 5)],
}

TARGET = 10
Q7_CORE = ((0, 1, 6), (0, 2, 5))


def tris(block):
    return set(combinations(block, 3))


def pblock(block, p):
    return tuple(sorted(p[i] for i in block))


def automorphisms(quads):
    target = set(tuple(sorted(q)) for q in quads)
    out = []
    for p in permutations(PTS):
        if set(pblock(q, p) for q in quads) == target:
            out.append(p)
    return out


def canonical_lines(lines, autos):
    best = None
    for p in autos:
        key = tuple(sorted(pblock(b, p) for b in lines))
        if best is None or key < best:
            best = key
    return best


def extensions(quads, target):
    qtris = set().union(*(tris(q) for q in quads))
    allowed = []
    for b in LINE_BLOCKS:
        if tris(b) & qtris:
            continue
        if any(len(set(b) & set(q)) > 2 for q in quads):
            continue
        allowed.append(b)
    # c(P)=35-ell-3q.  To test a putative c(P)<=TARGET we therefore
    # need at least 35-TARGET-3q collinear triples.
    need = max(0, 35 - target - 3 * len(quads))
    out = []

    def rec(i, fam, covered):
        if i == len(allowed):
            if len(covered) >= need:
                out.append((tuple(fam), frozenset(covered)))
            return
        possible = set(covered)
        for b in allowed[i:]:
            possible |= tris(b)
        if len(possible) < need:
            return
        rec(i + 1, fam, covered)
        b = allowed[i]
        if all(len(set(b) & set(a)) <= 1 for a in fam):
            fam.append(b)
            rec(i + 1, fam, covered | tris(b))
            fam.pop()

    rec(0, [], frozenset())
    return out


def contains_orbit(lines, core, autos):
    line_set = set(lines)
    return any(set(pblock(b, p) for b in core) <= line_set for p in autos)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=int, default=TARGET)
    args = parser.parse_args()
    summaries = {}
    for name, quads in PATTERNS.items():
        autos = automorphisms(quads)
        raw = extensions(quads, args.target)
        classes = {}
        for lines, covered in raw:
            key = canonical_lines(lines, autos)
            classes.setdefault(key, (key, len(covered), 35 - len(covered) - 3 * len(quads)))
        stats = {}
        for _, ell, c in classes.values():
            stats[(ell, c)] = stats.get((ell, c), 0) + 1
        print("pattern", name)
        print(" automorphisms", len(autos))
        print(" raw_extensions", len(raw))
        print(" classes", len(classes))
        print(" stats ell c class_count")
        for (ell, c), count in sorted(stats.items(), key=lambda x: (x[0][1], x[0][0])):
            print(" ", ell, c, count)
        for idx, (_, (lines, ell, c)) in enumerate(sorted(classes.items(), key=lambda kv: (kv[1][2], kv[1][1], kv[0])), 1):
            print(f" #{idx}: ell={ell} c={c}")
            print("  lines", " ".join("".join(map(str, b)) for b in lines) or "-")
        print()
        summaries[name] = (autos, classes)

    # Both five-circle orbits have no line extension capable of reaching ten
    # circles.  The six-circle orbit has a unique such class; its seven
    # three-point lines cover all 21 pairs and hence have no ordinary line,
    # which is impossible over the real plane by Sylvester--Gallai.  Every
    # surviving seven-circle extension contains the same two-line core, up to
    # an automorphism of the circle family.
    if args.target != 10:
        return
    assert len(summaries["q5a"][1]) == 0
    assert len(summaries["q5b"][1]) == 0
    q6_classes = list(summaries["q6"][1].values())
    assert len(q6_classes) == 1
    q6_lines, q6_ell, q6_c = q6_classes[0]
    assert (q6_ell, q6_c) == (7, 10)
    assert set().union(*(set(combinations(line, 2)) for line in q6_lines)) == set(combinations(PTS, 2))
    q7_autos, q7_classes = summaries["q7"]
    assert q7_classes
    assert all(contains_orbit(lines, Q7_CORE, q7_autos)
               for lines, _, _ in q7_classes.values())
    print("COVERAGE PASS: q5 has no candidate; q6 reduces to the real-nonrepresentable")
    print("seven-line system; every q7 candidate contains lines 016 and 025.")


if __name__ == "__main__":
    main()
