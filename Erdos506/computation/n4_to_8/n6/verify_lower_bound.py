from itertools import combinations

A = range(4)
pairs = [tuple(c) for c in combinations(A, 2)]

def disjoint_pair_sets(k):
    out = []
    for S in combinations(pairs, k):
        used = []
        ok = True
        for p in S:
            if set(p) & set(used):
                ok = False
                break
            used.extend(p)
        if ok:
            out.append(tuple(sorted(S)))
    return out

matchings = disjoint_pair_sets(2)

low_parameter_rows = []
for tau4 in range(3):
    for tau5 in range(3):
        for s in range(2):
            for b in range(3):
                # The extra four-point circles through both exterior points
                # pair points of the base circle by a partial matching.
                if b > 2:
                    continue
                N = 20 - tau4 - tau5 - s
                c = N - 3 - 3 * b
                if c <= 7:
                    low_parameter_rows.append((tau4, tau5, s, b, N, c))

print("parameter rows satisfying only tau<=2, s<=1, b<=2 and c<=7")
for row in low_parameter_rows:
    print("tau4=%d tau5=%d s=%d b=%d noncollinear_triples=%d circle_count=%d" % row)

print()
print("after the geometric consequence b=2 => s=0:")
for row in low_parameter_rows:
    tau4, tau5, s, b, N, c = row
    if b == 2 and s == 0:
        print("tau4=%d tau5=%d s=%d b=%d noncollinear_triples=%d circle_count=%d" % row)

print()
print("perfect matchings of the four points on the fixed circle:")
for i, m in enumerate(matchings):
    print(i, m)

print()
print("If b=2, the two four-point circles use one perfect matching B.")
print("A point with tau=2 must be one of the two other perfect matchings,")
print("because lying on a pair from B would make that four-point circle degenerate.")
for B in matchings:
    alternatives = [M for M in matchings if M != B]
    print("B=%s -> tau=2 alternatives %s" % (B, alternatives))
