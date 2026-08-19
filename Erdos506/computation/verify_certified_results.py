#!/usr/bin/env python3
"""Run the computations used for the main theorem and its no-three-collinear corollary.

The run checks all mathematical identities, exact constructions, and finite
classifiers.  Every subprocess must terminate with exit code zero.
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

from runtime_dependencies import activate_ortools

try:
    import sympy
except ModuleNotFoundError as exc:
    raise SystemExit(
        "SymPy is required for the exact algebraic checks. Install it in the "
        "Python environment used to run this entry point."
    ) from exc


HERE = Path(__file__).resolve().parent
PYTHON = sys.executable


SHORT_STEPS = [
    "n_ge_17/verify_n_ge_17.py",
    "n16/verify_n16.py",
    "n9_to_14/verify_largest_block_reduction.py",
    "n15/verify_n15.py",
    "n4_to_8/n6/verify_construction.py",
    "n4_to_8/n8/verify_construction.py",
    "n4_to_8/n6/verify_lower_bound.py",
    "n4_to_8/n6/classify_extremal_types.py",
    "n4_to_8/n7/classify_quad_packings.py",
    "n4_to_8/n7/classify_line_extensions.py",
    "n4_to_8/n7/verify_circle_count_spectrum.py",
    "n4_to_8/n7/verify_q7_exclusion.py",
    "n4_to_8/n7/verify_extremal_configuration.py",
    "n4_to_8/n8/classify_boundary_layers.py",
    "n4_to_8/n8/pattern1_inversion_relation.py",
    "n4_to_8/n8/pattern1_inversion_circles.py",
    "n4_to_8/n8/no_four_line_inversion_relations.py",
    "n4_to_8/n8/no_four_line_key_equations.py",
    "no_three_collinear/verify_no_three_collinear.py",
]


def run(relative: str, *arguments: str) -> None:
    path = HERE / relative
    if not path.is_file():
        raise FileNotFoundError(path)
    print(f"\n===== {relative} =====", flush=True)
    subprocess.run(
        [PYTHON, str(path), *arguments],
        cwd=HERE.parent,
        check=True,
    )


def main() -> None:
    print(f"Python {sys.version.split()[0]}: {sys.executable}")
    print(f"SymPy {sympy.__version__}")
    print(f"OR-Tools {activate_ortools()}: checked by computation/runtime_dependencies.py")
    for step in SHORT_STEPS:
        run(step)
    print("\nALL MAIN THEOREM AND NO-THREE-COLLINEAR CHECKS PASSED")


if __name__ == "__main__":
    main()
