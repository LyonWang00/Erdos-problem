#!/usr/bin/env python3
"""Exact arithmetic audit of the short projected-arrangement proof for n=15.

The proof has no orbit, support, coordinate, or solver enumeration.  Its only
finite external input is Cuntz's complete list of stretchable simplicial
arrangements on fourteen lines.
"""

from __future__ import annotations

import json
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "certificate.json"
REFERENCE_SHA256 = "5a765f81d1ee1a527e7ef3b5b4bfa02b354c52903f80ae11d3fbfda90ece3eed"

# The type with one sevenfold point is excluded by K<=7, because a local
# multiplicity r corresponds to a global block of size r+1.
CUNTZ_FOURTEEN_LINE_TYPES = {
    "A(14,1)": (7, 21, 0, 0, 0, 1),
    "A(14,2)": (11, 12, 4, 2, 0, 0),
    "A(14,3)": (9, 16, 4, 1, 0, 0),
    "A(14,4)": (10, 14, 4, 0, 1, 0),
}


def local_residual(row: tuple[int, ...]) -> int:
    t2, t3, t4, t5, _t6, *_ = row
    return 136 * t2 + 93 * t3 + 60 * t4 + 30 * t5 - 2902


def melchior_defect(row: tuple[int, ...]) -> int:
    return row[0] - 3 - sum((r - 3) * row[r - 2]
                            for r in range(4, len(row) + 2))


def main() -> None:
    allowed_simplicial = {
        label: row for label, row in CUNTZ_FOURTEEN_LINE_TYPES.items()
        if row[5] == 0
    }
    simplicial_checks = []
    for label, row in allowed_simplicial.items():
        assert sum(comb(r, 2) * row[r - 2]
                   for r in range(2, 8)) == comb(14, 2)
        assert melchior_defect(row) == 0
        residual = local_residual(row)
        assert residual in {0, 10, 80}
        simplicial_checks.append({
            "label": label,
            "type_t2_through_t7": list(row),
            "residual": residual,
        })
    assert sorted(item["residual"] for item in simplicial_checks) == [0, 10, 80]

    # For positive Melchior defect delta, the even-line face colouring rules
    # out delta=1.  Bojanowski gives
    # 3t4+9t5+18t6 <= 44+3delta.  Therefore, for delta>=2,
    # R >= 234+105delta-(25/3)(44+3delta)
    #   = (240delta-398)/3 > 0.
    assert 3 * (234 + 105 * 2) - 25 * (44 + 3 * 2) == 82
    # The affine numerator is increasing in the defect, so checking its first
    # allowed value proves the estimate for every defect >= 2.
    assert 240 * 2 - 398 == 82
    assert 240 > 0

    # Audit the coefficient identity
    # 420C-35270 = 140(l2-7)+420l4+980l5+1680l6+2520l7
    #               + sum_x R(t(x)).
    sizes = range(3, 8)
    line_pair_coefficients = {k: comb(k, 2) for k in sizes}
    triple_coefficients = {k: comb(k, 3) for k in sizes}
    local_incidence_coefficients = {
        3: 136 * 3,
        4: 93 * 4,
        5: 60 * 5,
        6: 30 * 6,
        7: 0,
    }
    explicit_line_slack = {3: 0, 4: 420, 5: 980,
                           6: 1680, 7: 2520}
    for k in sizes:
        line_coefficient = (
            -140 * line_pair_coefficients[k]
            + explicit_line_slack[k]
            + local_incidence_coefficients[k]
        )
        circle_coefficient = local_incidence_coefficients[k]
        assert line_coefficient == -12 * triple_coefficients[k]
        assert circle_coefficient == 420 - 12 * triple_coefficients[k]
    assert (140 * comb(15, 2) - 140 * 7
            - 15 * 2902 - 12 * comb(15, 3)) == -35270

    # If C=84, the left side after subtracting the LP bound is 10.  Every
    # positive line slack is at least 140.  The local analysis above shows
    # that the only residuals below 28 are 0 and 10, belonging respectively
    # to A(14,4) and A(14,2).  Hence fourteen centres would have the first
    # type and one the second; their t2 sum 151 is not divisible by 3, while
    # it must equal three times the number of global three-point blocks.
    assert 420 * 84 - 35270 == 10
    assert 14 * 10 + 11 == 151
    assert 151 % 3 == 1

    payload = {
        "schema_version": 1,
        "status": "PASS",
        "n": 15,
        "branch": "maximum block size at most seven",
        "solver_or_orbit_enumeration": False,
        "reference_pdf_sha256": REFERENCE_SHA256,
        "local_inequality": (
            "136*t2+93*t3+60*t4+30*t5-2902 >= 0"
        ),
        "positive_defect_minimum_residual": 28,
        "simplicial_type_checks": simplicial_checks,
        "global_identity": (
            "420*C-35270=140*(l2-7)+420*l4+980*l5+1680*l6+"
            "2520*l7+sum_x R_x"
        ),
        "C84_exclusion": "unique residual 10 gives sum(t2)=151 mod 3",
        "conclusion": "circle_count_at_least_85",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
