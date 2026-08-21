# Exceptional configurations for n = 6, 7, 8

This directory contains the exact checks used in the paper's classification of
the exceptional configurations and in its description of their Möbius moduli.
Run commands from this directory with the versions in
`../../requirements.txt`.

## Direct verification

- `verify_incidence_types.py` checks the sixteen abstract line--circle
  incidence types at `n = 8, c(P) = 18`, their three isomorphism classes after
  lines and circles are both regarded as generalized circles, the Fano and
  Möbius--Kantor real-realizability obstructions, and the three realizable
  line--circle incidence types.
- `verify_parametric_families.py` verifies all required line and circle
  determinant identities for the displayed parametric families.
- `verify_mobius_moduli.py` verifies the complete parameter domains and the
  finite quotient actions stated for `n = 6, 7`.
- `verify_n8_mobius_moduli.py` verifies the complete generalized-circle
  Möbius moduli at `n = 8, c(P) = 17, 18`, including the rectangular normal
  form, all 24 rational relabelling actions, the exact representatives and
  incidence lists, and the exceptional identifications on the
  seventeen-circle locus.
- `family_examples.py` enumerates all incidences of exact rational sample
  members and proves that the selected pairs are neither Möbius nor
  anti-Möbius equivalent.

The usual run order is

```text
python verify_incidence_types.py
python verify_parametric_families.py
python verify_mobius_moduli.py
python verify_n8_mobius_moduli.py
python family_examples.py
```

The files `classify_c18_backtrack.py`, `group_c18_uncolored.py`,
`realize_c18_by_recoloring.py`, and `construct_c18_representatives.py`
regenerate the `c18_*.json` records. They are not needed when checking the
released records. The word `uncolored` in the legacy filename means precisely
that lines and circles are both regarded as generalized circles; it is not
used as mathematical terminology in the paper.

## The seventeen-circle layer

The subdirectory `c17_exclusions` contains the seven abstract incidence
records, the exact rational verifier for the unique realized record, and the
symbolic exclusions of the other six records. First run

```text
python c17_exclusions/verify_target_records.py
python c17_exclusions/verify_realized_type.py
```

The scripts beginning with `type_1_`, `type_4_`, and `type_7_` give the hand
eliminations for types 1, 4, and 7. The three `exclude_type_*.py` files compute
the Gröbner-basis certificates for types 2, 5, and 6. Every computation is
exact over the integers, rationals, or a symbolic polynomial ring; no
numerical tolerance is used.
