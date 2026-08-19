"""Classify abstract six-point incidence types with exactly eight circles."""

from itertools import combinations, permutations


V = tuple(range(6))
TRIPLES = tuple(combinations(V, 3))
QUADS = tuple(combinations(V, 4))
TINDEX = {block: index for index, block in enumerate(TRIPLES)}
QMASK = {
    block: sum(1 << TINDEX[triple] for triple in combinations(block, 3))
    for block in QUADS
}
ALL = (1 << len(TRIPLES)) - 1


def canonical(lines, circles4):
    images = []
    for permutation in permutations(V):
        moved_lines = tuple(sorted(
            tuple(sorted(permutation[i] for i in block)) for block in lines
        ))
        moved_circles = tuple(sorted(
            tuple(sorted(permutation[i] for i in block)) for block in circles4
        ))
        images.append((moved_lines, moved_circles))
    return min(images)


def triple_packings(allowed, size, incompatible_with=()):
    allowed = tuple(allowed)
    answer = []

    def visit(start, selected):
        if len(selected) == size:
            answer.append(tuple(selected))
            return
        for index in range(start, len(allowed)):
            block = allowed[index]
            if any(len(set(block) & set(old)) > 1 for old in selected):
                continue
            if any(len(set(block) & set(old)) > 1 for old in incompatible_with):
                continue
            visit(index + 1, selected + [block])

    visit(0, [])
    return answer


types = set()
labelled = 0
for q4 in range(4):
    for circles4 in combinations(QUADS, q4):
        masks = [QMASK[block] for block in circles4]
        if any(masks[i] & masks[j] for i in range(q4) for j in range(i)):
            continue
        circle_mask = sum(masks)
        line_triples_needed = 12 - 3*q4
        for line4 in (None,) + QUADS:
            line4s = () if line4 is None else (line4,)
            if line4 is not None and QMASK[line4] & circle_mask:
                continue
            l3 = line_triples_needed - 4*len(line4s)
            if l3 < 0:
                continue
            allowed = [
                triple for triple in TRIPLES
                if not (1 << TINDEX[triple]) & circle_mask
                and (line4 is None or not (set(triple) <= set(line4)))
            ]
            for lines3 in triple_packings(allowed, l3, line4s):
                line_mask = sum(1 << TINDEX[block] for block in lines3)
                if line4 is not None:
                    line_mask |= QMASK[line4]
                if line_mask & circle_mask:
                    continue
                remaining = ALL & ~(line_mask | circle_mask)
                if remaining.bit_count() != 8-q4:
                    continue
                lines = tuple(sorted(line4s + lines3))
                labelled += 1
                types.add(canonical(lines, circles4))

assert labelled == 120
assert len(types) == 1
print("SIX-POINT EQUALITY CLASSIFICATION PASS")
print(f"labelled={labelled}")
print(f"orbits={len(types)}")
for index, (lines, circles4) in enumerate(sorted(types), 1):
    covered = sum(QMASK[block] for block in circles4)
    for line in lines:
        covered |= QMASK[line] if len(line) == 4 else 1 << TINDEX[line]
    circles3 = [block for block in TRIPLES if not covered & (1 << TINDEX[block])]
    print(index, "lines", lines, "circles4", circles4, "circles3", circles3)
