from itertools import combinations, permutations


PTS = tuple(range(7))
FOURS = tuple(combinations(PTS, 4))
INDEX = {b: i for i, b in enumerate(FOURS)}
BY_INDEX = {i: b for b, i in INDEX.items()}
PERMS = tuple(permutations(PTS))


def perm_block(block, p):
    return tuple(sorted(p[i] for i in block))


PERM_BIT = tuple(
    tuple(1 << INDEX[perm_block(b, p)] for b in FOURS)
    for p in PERMS
)


def canonical(fam):
    ids = tuple(INDEX[b] for b in fam)
    best = None
    for pi in range(len(PERMS)):
        mask = 0
        for i in ids:
            mask |= PERM_BIT[pi][i]
        if best is None or mask < best:
            best = mask
    return best


def bits(mask):
    i = 0
    out = []
    while mask:
        if mask & 1:
            out.append(BY_INDEX[i])
        mask >>= 1
        i += 1
    return tuple(out)


def main():
    raw = []

    def rec(i, fam):
        if i == len(FOURS):
            if len(fam) >= 5:
                raw.append(tuple(fam))
            return
        if len(fam) + (len(FOURS) - i) < 5:
            return
        rec(i + 1, fam)
        b = FOURS[i]
        bs = set(b)
        if all(len(bs & set(a)) <= 2 for a in fam):
            fam.append(b)
            rec(i + 1, fam)
            fam.pop()

    rec(0, [])
    classes = {}
    for fam in raw:
        classes.setdefault(canonical(fam), fam)

    raw_stats = {}
    class_stats = {}
    for fam in raw:
        raw_stats[len(fam)] = raw_stats.get(len(fam), 0) + 1
    for mask in classes:
        class_stats[mask.bit_count()] = class_stats.get(mask.bit_count(), 0) + 1

    print("raw_stats q count")
    for q in sorted(raw_stats):
        print(q, raw_stats[q])
    print("class_stats q count")
    for q in sorted(class_stats):
        print(q, class_stats[q])
    print("classes")
    for idx, mask in enumerate(sorted(classes, key=lambda m: (m.bit_count(), m)), 1):
        print(f"#{idx}: q={mask.bit_count()}")
        print(" ", " ".join("".join(map(str, b)) for b in bits(mask)))


if __name__ == "__main__":
    main()
