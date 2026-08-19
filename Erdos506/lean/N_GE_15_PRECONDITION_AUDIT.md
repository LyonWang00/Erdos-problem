# Preconditions in the mathematical proof for `n >= 15`

## Uniform branch, `n >= 17`

For inversion at `x`, the image contains `q=n-1` noncollinear points and every
connecting line has multiplicity at most `K-1`.  Therefore
`3(K-1)<2(n-1)` is exactly sufficient for the form of Bojanowski's inequality
used in the proof.  Melchior applies because the image is a non-pencil real
projective arrangement.  Lines through the inversion centre have disjoint
image-point supports and each uses at least two image points, giving the
bound by `floor(q/2)`.

In the complementary branch the original proof weakened the available
condition.  The Lean proof instead keeps the exact consequence
`3r <= n-1`, where `r=n-K`, and proves the coarse circle bound directly by a
nonnegative factorization.  The case `r=1` (`K=n-1`) remains separate and uses
the exact largest-circle formula, so no rounding is lost there.

## The sixteen-point layer

For `K<=8`, the maximum image-line multiplicity is at most seven, strictly
below `2q/3=10`, so Bojanowski applies.  In the `K=7` branch, four six-point
image lines would have union at least `4*6-C(4,2)=18>15`; hence at most three
seven-point circles pass through a point.  All incidence-moment estimates use
only this bound and the fact that two circles meet in at most two points.

In the equality case with four seven-point circles, the moment inequality
forces a point of circle-degree one.  The coefficient identity then forces
the unique local profile `(14,12,0,4,1)`.  Its four five-point lines and one
six-point line have union at least `4*5+6-C(5,2)=16>15`, a contradiction.

For `K=8`, the special seven-point image line exists only at the eight centres
lying on the chosen largest line or circle.  The remaining eight centres use
the unconditional local bound, whose multiplicity hypothesis still holds.

## The fifteen-point layer

After the largest-block reduction, `K<=7`, so the fourteen-line dual
arrangement has vertex multiplicity at most six and satisfies Bojanowski's
hypothesis.  Because the number of lines is even, the product of their
homogeneous equations descends to the projective plane and two-colours its
faces.  The two colour-class excesses are nonnegative and congruent modulo
three; consequently their sum, the Melchior defect, cannot equal one.

When the defect is zero, Cuntz's complete fourteen-line classification is
used.  The near-pencil and the row with a sevenfold vertex are excluded by
`K<=7`.  The three remaining multiplicity vectors are the only external
finite input in this layer.  The global pair, triple, and local-incidence
identities count maximal blocks, so every pair/triple occurs exactly once in
the appropriate line or circle block.  Finally, Csima--Sawyer applies to the
original noncollinear fifteen-point set and gives `l2>=7`.
