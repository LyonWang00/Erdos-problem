#!/usr/bin/env python3
"""Exact algebra audit for the mathematical proof covering every n >= 17.

The proof combines the point-pair identity, Melchior's inequality, and
Bojanowski's inequality after inversion about an arbitrary point.  It has no
enumeration of n, K, line-multiplicity vectors, or geometric configurations.
Only the endpoint n=17 is checked separately after the general parity
calculation.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "certificate.json"


def main() -> None:
    n, K, i, u = sp.symbols("n K i u", integer=True, positive=True)

    # Positive multipliers of the point-pair identity, Melchior's
    # inequality, and Bojanowski's inequality, respectively.
    alpha = (K + 1) / (18 * K)
    beta = (K + 1) / (6 * K)
    gamma = (K - 2) / (36 * K)

    # The resulting coefficient is exactly 1/(i+1) for i=2,3.  For i=4
    # and 5<=i<=K-1, the nonnegative slack factors as follows.
    slack_2 = sp.factor(
        sp.Rational(1, 3) - (alpha + beta + 4 * gamma)
    )
    slack_3 = sp.factor(
        sp.Rational(1, 4) - (3 * alpha + 3 * gamma)
    )
    slack_4 = sp.factor(
        sp.Rational(1, 5) - (6 * alpha - beta)
    )
    slack_i = sp.factor(
        1 / (i + 1)
        - (
            alpha * i * (i - 1) / 2
            - beta * (i - 3)
            - gamma * i * (i - 4)
        )
    )
    assert slack_2 == 0
    assert slack_3 == 0
    assert slack_4 == (K - 5) / (30 * K)
    expected_slack_i = (i - 3) * (i - 2) * (K - i - 1) / (
        12 * K * (i + 1)
    )
    assert sp.factor(slack_i - expected_slack_i) == 0

    # If q=n-1 and e=floor(q/2), the local contribution at every inversion
    # centre is at least the following quantity.
    q, e = sp.symbols("q e", integer=True, nonnegative=True)
    local = sp.factor(
        alpha * q * (q - 1) / 2
        + 3 * beta
        + 4 * q * gamma
        - e / 3
    )
    expected_local = -(
        K * (12 * e - (q + 1) ** 2 - (q + 1) - 16)
        - (q + 1) ** 2 + 11 * (q + 1) - 28
    ) / (36 * K)
    assert sp.factor(local - expected_local) == 0

    # The margin n*local-(F(n)-1) decreases with K.  In the branch where
    # Bojanowski applies, K <= (2n+1)/3, so the continuous right endpoint is
    # a valid lower estimate.  Separate formulas remove the parity in e.
    qn = n - 1
    choose_q_2 = qn * (qn - 1) / 2
    local_n = local.subs(q, qn)
    even_e = (n - 2) / 2
    odd_e = (n - 1) / 2
    even_target_minus_one = choose_q_2 - even_e
    odd_target_minus_one = choose_q_2 - odd_e
    even_margin = sp.factor(
        n * local_n.subs(e, even_e) - even_target_minus_one
    )
    odd_margin = sp.factor(
        n * local_n.subs(e, odd_e) - odd_target_minus_one
    )
    expected_even = (
        K * n**3 - 23 * K * n**2 + 100 * K * n - 72 * K
        + n**3 - 11 * n**2 + 28 * n
    ) / (36 * K)
    expected_odd = (
        K * n**3 - 23 * K * n**2 + 94 * K * n - 54 * K
        + n**3 - 11 * n**2 + 28 * n
    ) / (36 * K)
    assert sp.factor(even_margin - expected_even) == 0
    assert sp.factor(odd_margin - expected_odd) == 0
    assert sp.factor(sp.diff(even_margin, K)) == (
        -n * (n - 7) * (n - 4) / (36 * K**2)
    )
    assert sp.factor(sp.diff(odd_margin, K)) == (
        -n * (n - 7) * (n - 4) / (36 * K**2)
    )

    endpoint = (2 * n + 1) / 3
    even_endpoint = sp.factor(even_margin.subs(K, endpoint))
    odd_endpoint = sp.factor(odd_margin.subs(K, endpoint))
    even_polynomial = n**4 - 21 * n**3 + 72 * n**2 + 20 * n - 36
    odd_polynomial = n**4 - 21 * n**3 + 66 * n**2 + 35 * n - 27
    assert sp.factor(
        even_endpoint - even_polynomial / (18 * (2 * n + 1))
    ) == 0
    assert sp.factor(
        odd_endpoint - odd_polynomial / (18 * (2 * n + 1))
    ) == 0

    even_coefficients = [
        int(value)
        for value in sp.Poly(sp.expand(even_polynomial.subs(n, 18 + 2 * u)), u)
        .all_coeffs()
    ]
    odd_coefficients = [
        int(value)
        for value in sp.Poly(sp.expand(odd_polynomial.subs(n, 19 + 2 * u)), u)
        .all_coeffs()
    ]
    assert even_coefficients == [16, 408, 3528, 11056, 6156]
    assert odd_coefficients == [16, 440, 4140, 14472, 10746]
    assert all(value > 0 for value in even_coefficients + odd_coefficients)

    # At n=17 the integer bound is K<=11, not the continuous endpoint 35/3.
    n17_margin = sp.factor(odd_margin.subs({n: 17, K: 11}))
    assert n17_margin == sp.Rational(10, 33)

    # The complementary branch has 3K>=n+12.  The largest-line/circle
    # estimate is bounded below by B(n,K).  Its derivative is positive at
    # the left endpoint for n>=17, and the cubic can therefore have no
    # interior minimum.  Both endpoints exceed the target.
    k_real = sp.symbols("k_real", real=True)
    B = 1 + (n - k_real) * k_real * (3 * k_real - n - 3) / 4
    r_real = n - k_real
    # Replacing floor(K/2) by K/2 in the subtracted terms of the two
    # largest-block estimates gives these lower relaxations.  The circle
    # relaxation is B itself; the line relaxation is still stronger.
    line_relaxed = r_real * k_real * (3 * k_real - n - 1) / 4
    circle_relaxed = 1 + r_real * k_real * (3 * k_real - n - 3) / 4
    assert sp.factor(
        line_relaxed - B - (r_real * k_real / 2 - 1)
    ) == 0
    assert sp.factor(circle_relaxed - B) == 0
    upper_target = (n**2 - 4 * n + 6) / 2
    left = (n + 12) / 3
    derivative = sp.diff(B, k_real)
    derivative_left = sp.factor(derivative.subs(k_real, left))
    assert sp.factor(
        derivative_left - (2 * n**2 + 21 * n - 360) / 12
    ) == 0
    assert sp.factor(derivative_left.subs(n, 17)) > 0
    assert sp.Poly(derivative, k_real).LC() < 0
    left_margin = sp.factor(B.subs(k_real, left) - upper_target)
    right_margin = sp.factor(B.subs(k_real, n - 2) - upper_target)
    assert sp.factor(left_margin - (5 * n - 38)) == 0
    assert sp.factor(right_margin - (n - 7) * (n - 2) / 2) == 0

    payload = {
        "schema_version": 1,
        "status": "PASS",
        "theorem": "c(n)=F(n) for every n>=17",
        "parameter_or_profile_enumeration": 0,
        "local_coefficient_slacks": {
            "i=2": str(slack_2),
            "i=3": str(slack_3),
            "i=4": str(slack_4),
            "5<=i<=K-1": str(slack_i),
        },
        "small_maximum_block_branch": {
            "condition": "3(K-1)<2(n-1)",
            "n17_margin": str(n17_margin),
            "even_n_ge_18_positive_coefficients": even_coefficients,
            "odd_n_ge_19_positive_coefficients": odd_coefficients,
        },
        "large_maximum_block_branch": {
            "condition": "3K>=n+12",
            "largest_line_relaxation_minus_B": str(
                sp.factor(line_relaxed - B)),
            "largest_circle_relaxation_minus_B": str(
                sp.factor(circle_relaxed - B)),
            "derivative_at_left_endpoint": str(derivative_left),
            "left_endpoint_margin": str(left_margin),
            "right_endpoint_margin": str(right_margin),
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
