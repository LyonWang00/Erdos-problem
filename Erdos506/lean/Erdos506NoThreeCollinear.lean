import Mathlib

/-!
# Arithmetic-combinatorial formalization of the no-three-collinear variant

`circleMinimumNoThree n` in the final theorem denotes the geometric minimum
called `c_nc(n)` in the manuscript.  Euclidean incidence statements and the
outputs of exhaustive finite classifiers enter as explicit hypotheses.  All
arithmetic, congruence, interval, and case-closing deductions are checked here.
There are no proof placeholders or added axioms.
-/

namespace Erdos506NoThree

def benchmark (n : ℕ) : ℕ := 1 + Nat.choose (n - 1) 2

def target (n : ℕ) : ℕ := if n = 8 then 20 else benchmark n

theorem benchmark_small_values :
    benchmark 4 = 4 ∧ benchmark 5 = 7 ∧ benchmark 6 = 11 ∧
    benchmark 7 = 16 ∧ benchmark 8 = 22 ∧ benchmark 9 = 29 ∧
    benchmark 10 = 37 := by
  native_decide

/- If all eight-point circles have size three or four, triple ownership gives
`56 = c3 + 4*c4` and `C = c3+c4`.  Thus no total strictly between 17 and 20
is possible. -/
theorem eight_point_gap
    (C c3 c4 : ℕ)
    (hcount : C = c3 + c4)
    (htriples : 56 = c3 + 4 * c4)
    (hlower : 17 ≤ C)
    (hupper : C < 20) :
    C = 17 := by
  omega

/- The elementary degree count behind the 3-(9,4,1) packing bound. -/
theorem nine_point_four_circle_bound
    (c4 : ℕ) (hincidence : 4 * c4 ≤ 9 * 8) :
    c4 ≤ 18 := by
  omega

theorem nine_point_circle_lower_bound
    (C c3 c4 : ℕ)
    (hcount : C = c3 + c4)
    (htriples : 84 = c3 + 4 * c4)
    (hpacking : c4 ≤ 18) :
    30 ≤ C := by
  omega

/- The six-point four-circle packing calculation. -/
theorem six_point_circle_lower_bound
    (C c3 c4 : ℕ)
    (hcount : C = c3 + c4)
    (htriples : 20 = c3 + 4 * c4)
    (hpacking : c4 ≤ 3) :
    11 ≤ C := by
  omega

/- The local n=10 inequality after pair counting, Melchior, Hirzebruch, and
the divisibility condition have been translated to integer coefficients. -/
theorem ten_point_local_bound
    (t2 t3 t4 : ℕ)
    (hpairs : t2 + 3 * t3 + 6 * t4 = 36)
    (hmelchior : 3 + t4 ≤ t2)
    (hhirzebruch : 36 ≤ 4 * t2 + 3 * t3)
    (hdivisible : t2 % 3 = 0) :
    72 ≤ 20 * t2 + 15 * t3 + 12 * t4 := by
  omega

theorem ten_point_local_equality
    (t2 t3 t4 : ℕ)
    (hpairs : t2 + 3 * t3 + 6 * t4 = 36)
    (hhirzebruch : 36 ≤ 4 * t2 + 3 * t3)
    (hdivisible : t2 % 3 = 0)
    (hequality : 20 * t2 + 15 * t3 + 12 * t4 = 72) :
    t2 = 6 ∧ t3 = 4 ∧ t4 = 3 := by
  omega

/- Positivity of the two endpoint factors in the n>=11 argument. -/
theorem small_endpoint_factor_positive
    (n : ℤ) (hn : 11 ≤ n) :
    0 < (n - 4) * (n - 1) * (n ^ 2 - 10 * n - 9) := by
  have h1 : 0 < n - 4 := by omega
  have h2 : 0 < n - 1 := by omega
  have h3 : 0 < n ^ 2 - 10 * n - 9 := by
    nlinarith [sq_nonneg (n - 11)]
  positivity

theorem large_endpoint_factor_nonnegative
    (n : ℤ) (hn : 4 ≤ n) :
    0 ≤ (n - 4) * (n - 1) * (2 * n - 9) := by
  by_cases hn4 : n = 4
  · subst n
    norm_num
  · have h1 : 0 < n - 4 := by omega
    have h2 : 0 < n - 1 := by omega
    have h3 : 0 < 2 * n - 9 := by omega
    positivity

/- End-to-end closure.  The hypotheses have the same semantic roles as in the
manuscript: explicit constructions give the upper bounds; the mathematical
n>=11 argument gives `hLargeLower`; exact finite verifications give the six
small lower bounds and the exclusion of the eight-point 17-circle class. -/
theorem no_three_collinear_result
    (circleMinimumNoThree : ℕ → ℕ)
    (hGeneralUpper : ∀ n, 4 ≤ n → circleMinimumNoThree n ≤ benchmark n)
    (hEightUpper : circleMinimumNoThree 8 ≤ 20)
    (hLargeLower : ∀ n, 11 ≤ n → benchmark n ≤ circleMinimumNoThree n)
    (hFourLower : 4 ≤ circleMinimumNoThree 4)
    (hFiveLower : 7 ≤ circleMinimumNoThree 5)
    (hSixLower : 11 ≤ circleMinimumNoThree 6)
    (hSevenLower : 16 ≤ circleMinimumNoThree 7)
    (hEightMainLower : 17 ≤ circleMinimumNoThree 8)
    (hEightTripleData : circleMinimumNoThree 8 < 20 →
      ∃ c3 c4, circleMinimumNoThree 8 = c3 + c4 ∧ 56 = c3 + 4 * c4)
    (hEightNoSeventeen : circleMinimumNoThree 8 ≠ 17)
    (hNineLower : 29 ≤ circleMinimumNoThree 9)
    (hTenLower : 37 ≤ circleMinimumNoThree 10) :
    ∀ n, 4 ≤ n → circleMinimumNoThree n = target n := by
  intro n hn
  by_cases h8 : n = 8
  · subst n
    have hnotlt : ¬circleMinimumNoThree 8 < 20 := by
      intro hlt
      obtain ⟨c3, c4, hcount, htriples⟩ := hEightTripleData hlt
      exact hEightNoSeventeen
        (eight_point_gap (circleMinimumNoThree 8) c3 c4
          hcount htriples hEightMainLower hlt)
    have hlower : 20 ≤ circleMinimumNoThree 8 := by omega
    change circleMinimumNoThree 8 = 20
    omega
  · have hupper := hGeneralUpper n hn
    have htarget : target n = benchmark n := by simp [target, h8]
    rw [htarget]
    by_cases hlarge : 11 ≤ n
    · exact Nat.le_antisymm hupper (hLargeLower n hlarge)
    · interval_cases n <;>
        norm_num [benchmark, Nat.choose] at hupper ⊢ <;> omega

end Erdos506NoThree
