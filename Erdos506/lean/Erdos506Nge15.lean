import Mathlib

/-!
# Arithmetic formalization for the `n ≥ 15` proof of Erdős problem 506

This file formalizes the exact algebraic and integral consequences used in
the manuscript.  The real-projective line-arrangement theorems (Melchior,
Bojanowski, Csima--Sawyer, and Cuntz's classification) enter the mathematical
proof as named hypotheses; none of them is silently reintroduced as an axiom
here.  There are no proof placeholders.
-/

namespace Erdos506

open scoped BigOperators

def benchmark (n : ℕ) : ℕ :=
  1 + Nat.choose (n - 1) 2 - (n - 1) / 2

theorem benchmark_even (m : ℕ) (hm : 1 ≤ m) :
    benchmark (2 * m) = 2 * (m - 1) ^ 2 + 1 := by
  rw [benchmark, Nat.choose_two_right]
  have hsub : 2 * m - 1 = 2 * (m - 1) + 1 := by omega
  rw [hsub]
  simp only [Nat.add_sub_cancel]
  have hprod : (2 * (m - 1) + 1) * (2 * (m - 1)) =
      2 * ((2 * (m - 1) + 1) * (m - 1)) := by ring
  rw [hprod]
  have hexpand : (2 * (m - 1) + 1) * (m - 1) =
      2 * (m - 1) ^ 2 + (m - 1) := by ring
  rw [hexpand]
  omega

theorem benchmark_odd (m : ℕ) (hm : 1 ≤ m) :
    benchmark (2 * m + 1) = 2 * m * (m - 1) + 1 := by
  rw [benchmark, Nat.choose_two_right]
  simp only [Nat.add_sub_cancel]
  have hsub : 2 * m - 1 = 2 * (m - 1) + 1 := by omega
  rw [hsub]
  have hprod : 2 * m * (2 * (m - 1) + 1) =
      2 * (m * (2 * (m - 1) + 1)) := by ring
  rw [hprod]
  have hexpand : m * (2 * (m - 1) + 1) =
      2 * m * (m - 1) + m := by ring
  rw [hexpand]
  omega

theorem integer_strict_margin_closes
    (C F : ℕ) (hmargin : (F : ℚ) - 1 < C) : F ≤ C := by
  cases F with
  | zero => simp
  | succ F =>
      by_contra h
      have hnat : C ≤ F := by omega
      have hcast : (C : ℚ) ≤ (F : ℚ) := by exact_mod_cast hnat
      norm_num at hmargin
      linarith

section LocalCoefficientCombination

/- The coefficient calculation in the Bojanowski--Melchior combination. -/
theorem coefficient_difference_ge_five
    (K i : ℚ) (hK : K ≠ 0) (hi : i + 1 ≠ 0) :
    1 / (i + 1) -
        ((K + 1) / (18 * K) * (i * (i - 1) / 2) +
          (K + 1) / (6 * K) * (-(i - 3)) +
          (K - 2) / (36 * K) * (-(i * (i - 4)))) =
      (i - 3) * (i - 2) * (K - i - 1) / (12 * K * (i + 1)) := by
  field_simp [hK, hi]
  ring

theorem coefficient_difference_four (K : ℚ) (hK : K ≠ 0) :
    1 / 5 -
        ((K + 1) / (18 * K) * 6 +
          (K + 1) / (6 * K) * (-1)) =
      (K - 5) / (30 * K) := by
  field_simp
  ring

end LocalCoefficientCombination

section LargeN

def localBojanowskiBound (n K floorHalf : ℚ) : ℚ :=
  ((K + 1) / (18 * K)) * ((n - 1) * (n - 2) / 2) +
    (K + 1) / (2 * K) + ((K - 2) / (9 * K)) * (n - 1) -
    floorHalf / 3

def evenMargin (n K : ℚ) : ℚ :=
  (K * n ^ 3 - 23 * K * n ^ 2 + 100 * K * n - 72 * K +
      n ^ 3 - 11 * n ^ 2 + 28 * n) / (36 * K)

def oddMargin (n K : ℚ) : ℚ :=
  (K * n ^ 3 - 23 * K * n ^ 2 + 94 * K * n - 54 * K +
      n ^ 3 - 11 * n ^ 2 + 28 * n) / (36 * K)

/- These identities check the two parity expansions after summing the local
bound over all inversion centres. -/
theorem even_margin_identity (n K : ℚ) (hK : K ≠ 0) :
    n * localBojanowskiBound n K ((n - 2) / 2) -
      ((n ^ 2 - 4 * n + 6) / 2 - 1) =
      evenMargin n K := by
  simp only [evenMargin, localBojanowskiBound]
  field_simp
  ring

theorem odd_margin_identity (n K : ℚ) (hK : K ≠ 0) :
    n * localBojanowskiBound n K ((n - 1) / 2) -
      ((n ^ 2 - 4 * n + 5) / 2 - 1) =
      oddMargin n K := by
  simp only [oddMargin, localBojanowskiBound]
  field_simp
  ring

theorem nonnegative_div_antitone
    (D K U : ℚ) (hD : 0 ≤ D) (hK : 0 < K) (hU : 0 < U)
    (hKU : K ≤ U) : D / U ≤ D / K := by
  exact (div_le_div_iff₀ hU hK).2 (mul_le_mul_of_nonneg_left hKU hD)

theorem even_endpoint_identity (n : ℚ) (hn : 2 * n + 1 ≠ 0) :
    evenMargin n ((2 * n + 1) / 3) =
      (n ^ 4 - 21 * n ^ 3 + 72 * n ^ 2 + 20 * n - 36) /
        (18 * (2 * n + 1)) := by
  simp only [evenMargin]
  field_simp
  ring

theorem odd_endpoint_identity (n : ℚ) (hn : 2 * n + 1 ≠ 0) :
    oddMargin n ((2 * n + 1) / 3) =
      (n ^ 4 - 21 * n ^ 3 + 66 * n ^ 2 + 35 * n - 27) /
        (18 * (2 * n + 1)) := by
  simp only [oddMargin]
  field_simp
  ring

theorem even_endpoint_numerator_positive (u : ℚ) (hu : 0 ≤ u) :
    0 < (18 + 2 * u) ^ 4 - 21 * (18 + 2 * u) ^ 3 +
      72 * (18 + 2 * u) ^ 2 + 20 * (18 + 2 * u) - 36 := by
  have h : 0 < 16 * u ^ 4 + 408 * u ^ 3 + 3528 * u ^ 2 +
      11056 * u + 6156 := by positivity
  nlinarith

theorem odd_endpoint_numerator_positive (u : ℚ) (hu : 0 ≤ u) :
    0 < (19 + 2 * u) ^ 4 - 21 * (19 + 2 * u) ^ 3 +
      66 * (19 + 2 * u) ^ 2 + 35 * (19 + 2 * u) - 27 := by
  have h : 0 < 16 * u ^ 4 + 440 * u ^ 3 + 4140 * u ^ 2 +
      14472 * u + 10746 := by positivity
  nlinarith

/- A calculus-free version of the monotonicity step: the margin is an affine
term plus `n(n-7)(n-4)/(36K)`. -/
theorem even_margin_antitone
    (n K U : ℚ) (hn : 7 ≤ n) (hK : 0 < K) (hU : 0 < U)
    (hKU : K ≤ U) : evenMargin n U ≤ evenMargin n K := by
  have hn0 : 0 ≤ n := by linarith
  have hn7 : 0 ≤ n - 7 := by linarith
  have hn4 : 0 ≤ n - 4 := by linarith
  have hD : 0 ≤ n * (n - 7) * (n - 4) :=
    mul_nonneg (mul_nonneg hn0 hn7) hn4
  have hratio := nonnegative_div_antitone
    (n * (n - 7) * (n - 4)) K U hD hK hU hKU
  simp only [evenMargin] at *
  field_simp at hratio ⊢
  nlinarith

theorem odd_margin_antitone
    (n K U : ℚ) (hn : 7 ≤ n) (hK : 0 < K) (hU : 0 < U)
    (hKU : K ≤ U) : oddMargin n U ≤ oddMargin n K := by
  have hn0 : 0 ≤ n := by linarith
  have hn7 : 0 ≤ n - 7 := by linarith
  have hn4 : 0 ≤ n - 4 := by linarith
  have hD : 0 ≤ n * (n - 7) * (n - 4) :=
    mul_nonneg (mul_nonneg hn0 hn7) hn4
  have hratio := nonnegative_div_antitone
    (n * (n - 7) * (n - 4)) K U hD hK hU hKU
  simp only [oddMargin] at *
  field_simp at hratio ⊢
  nlinarith

theorem even_margin_positive
    (u K : ℚ) (hu : 0 ≤ u) (hK : 0 < K)
    (hRange : 3 * K ≤ 2 * (18 + 2 * u) + 1) :
    0 < evenMargin (18 + 2 * u) K := by
  let n : ℚ := 18 + 2 * u
  let U : ℚ := (2 * n + 1) / 3
  have hn : 7 ≤ n := by dsimp [n]; linarith
  have hU : 0 < U := by dsimp [U, n]; positivity
  have hKU : K ≤ U := by dsimp [U, n] at *; linarith
  have hmono := even_margin_antitone n K U hn hK hU hKU
  have hnum := even_endpoint_numerator_positive u hu
  have hden : 0 < 18 * (2 * n + 1) := by dsimp [n]; positivity
  have hendpoint : 0 < evenMargin n U := by
    rw [even_endpoint_identity n (by positivity)]
    exact div_pos hnum hden
  exact lt_of_lt_of_le hendpoint hmono

theorem odd_margin_positive
    (u K : ℚ) (hu : 0 ≤ u) (hK : 0 < K)
    (hRange : 3 * K ≤ 2 * (19 + 2 * u) + 1) :
    0 < oddMargin (19 + 2 * u) K := by
  let n : ℚ := 19 + 2 * u
  let U : ℚ := (2 * n + 1) / 3
  have hn : 7 ≤ n := by dsimp [n]; linarith
  have hU : 0 < U := by dsimp [U, n]; positivity
  have hKU : K ≤ U := by dsimp [U, n] at *; linarith
  have hmono := odd_margin_antitone n K U hn hK hU hKU
  have hnum := odd_endpoint_numerator_positive u hu
  have hden : 0 < 18 * (2 * n + 1) := by dsimp [n]; positivity
  have hendpoint : 0 < oddMargin n U := by
    rw [odd_endpoint_identity n (by positivity)]
    exact div_pos hnum hden
  exact lt_of_lt_of_le hendpoint hmono

theorem n17_margin_positive (K : ℚ) (hK : 0 < K) (hKmax : K ≤ 11) :
    0 < oddMargin 17 K := by
  have hmono := odd_margin_antitone (17 : ℚ) K 11 (by norm_num) hK
    (by norm_num) hKmax
  have hvalue : oddMargin 17 11 = 10 / 33 := by
    norm_num [oddMargin]
  rw [hvalue] at hmono
  linarith

def coarseCircleBound (n r : ℚ) : ℚ :=
  1 + r * (n - r) * (2 * n - 3 * r - 3) / 4

def quadraticBenchmarkUpper (n : ℚ) : ℚ :=
  (n ^ 2 - 4 * n + 6) / 2

/- This replaces the manuscript's derivative argument in the large-block
branch.  It uses the stronger inequality actually available there,
`3r ≤ n-1`, and an explicit factorization around `r=2`. -/
theorem coarse_circle_bound_ge
    (n r : ℚ) (hn : 17 ≤ n) (hr : 2 ≤ r) (hrange : 3 * r ≤ n - 1) :
    quadraticBenchmarkUpper n ≤ coarseCircleBound n r := by
  let Q : ℚ := 2 * n ^ 2 - 5 * n * r - 13 * n +
    3 * r ^ 2 + 9 * r + 18
  have hA1 : 0 ≤ n - 2 := by linarith
  have hA2 : 0 ≤ 2 * n - 23 := by linarith
  have hA : 0 ≤ (n - 2) * (2 * n - 23) := mul_nonneg hA1 hA2
  have hB1 : 0 ≤ 4 * n - 3 * r - 8 := by linarith
  have hB2 : 0 ≤ n - 3 * r - 1 := by linarith
  have hB : 0 ≤ (4 * n - 3 * r - 8) * (n - 3 * r - 1) :=
    mul_nonneg hB1 hB2
  have hQ : 0 ≤ Q := by
    dsimp [Q]
    nlinarith [hA, hB]
  have hprod : 0 ≤ (r - 2) * Q := mul_nonneg (by linarith) hQ
  dsimp [quadraticBenchmarkUpper, coarseCircleBound]
  dsimp [Q] at hprod
  nlinarith [mul_nonneg (show 0 ≤ n - 7 by linarith)
    (show 0 ≤ n - 2 by linarith)]

def largestLineBound (n K : ℕ) : ℕ :=
  let r := n - K
  r * Nat.choose K 2 - Nat.choose r 2 * (K / 2)

def largestCircleBound (n K : ℕ) : ℕ :=
  let r := n - K
  1 + r * (Nat.choose K 2 - K / 2) - Nat.choose r 2 * (K / 2)

theorem n16_large_block_table (K : ℕ) (hmin : 9 ≤ K) (hmax : K ≤ 15) :
    99 ≤ largestLineBound 16 K ∧ 99 ≤ largestCircleBound 16 K := by
  interval_cases K <;> native_decide

theorem n15_large_block_table (K : ℕ) (hmin : 8 ≤ K) (hmax : K ≤ 14) :
    85 ≤ largestLineBound 15 K ∧ 85 ≤ largestCircleBound 15 K := by
  interval_cases K <;> native_decide

theorem n16_small_block_table (K : ℕ) (hmin : 3 ≤ K) (hmax : K ≤ 6) :
    98 < 16 * localBojanowskiBound 16 K 7 := by
  interval_cases K <;> norm_num [localBojanowskiBound]

theorem n16_k8_unconditional_value :
    localBojanowskiBound 16 8 7 = 145 / 24 := by
  norm_num [localBojanowskiBound]

/- End-to-end arithmetic closures for the small-maximum-block branch.  The
hypothesis `hcount` is exactly the summed local incidence bound; the remaining
steps, including passage from a strict rational margin to an integral circle
count, are checked here. -/
theorem even_small_block_closes
    (u C : ℕ) (K : ℚ) (hK : 0 < K)
    (hRange : 3 * K ≤ 2 * (18 + 2 * (u : ℚ)) + 1)
    (hcount :
      ((benchmark (18 + 2 * u) : ℕ) : ℚ) - 1 +
          evenMargin (18 + 2 * (u : ℚ)) K ≤ C) :
    benchmark (18 + 2 * u) ≤ C := by
  have hmargin := even_margin_positive (u : ℚ) K (by positivity) hK hRange
  apply integer_strict_margin_closes
  exact lt_of_lt_of_le (by linarith) hcount

theorem odd_small_block_closes
    (u C : ℕ) (K : ℚ) (hK : 0 < K)
    (hRange : 3 * K ≤ 2 * (19 + 2 * (u : ℚ)) + 1)
    (hcount :
      ((benchmark (19 + 2 * u) : ℕ) : ℚ) - 1 +
          oddMargin (19 + 2 * (u : ℚ)) K ≤ C) :
    benchmark (19 + 2 * u) ≤ C := by
  have hmargin := odd_margin_positive (u : ℚ) K (by positivity) hK hRange
  apply integer_strict_margin_closes
  exact lt_of_lt_of_le (by linarith) hcount

theorem n17_small_block_closes
    (C : ℕ) (K : ℚ) (hK : 0 < K) (hKmax : K ≤ 11)
    (hcount : ((benchmark 17 : ℕ) : ℚ) - 1 + oddMargin 17 K ≤ C) :
    benchmark 17 ≤ C := by
  have hmargin := n17_margin_positive K hK hKmax
  apply integer_strict_margin_closes
  exact lt_of_lt_of_le (by linarith) hcount

theorem even_large_block_closes
    (m r C : ℕ) (hm : 9 ≤ m) (hr : 2 ≤ r)
    (hrange : 3 * r ≤ 2 * m - 1)
    (hcount : coarseCircleBound (2 * m) r ≤ C) :
    benchmark (2 * m) ≤ C := by
  have hn : (17 : ℚ) ≤ 2 * (m : ℚ) := by exact_mod_cast (by omega : 17 ≤ 2 * m)
  have hrr : (2 : ℚ) ≤ r := by exact_mod_cast hr
  have hrange' : 3 * r + 1 ≤ 2 * m := by omega
  have hrangeq' : (3 : ℚ) * r + 1 ≤ 2 * m := by exact_mod_cast hrange'
  have hrangeq : (3 : ℚ) * r ≤ 2 * m - 1 := by
    linarith
  have hcoarse := coarse_circle_bound_ge (2 * (m : ℚ)) (r : ℚ)
    hn hrr hrangeq
  have hcastsub : ((m - 1 : ℕ) : ℚ) = (m : ℚ) - 1 := by
    rw [Nat.cast_sub (by omega)]
    norm_num
  have hbenchmark :
      ((benchmark (2 * m) : ℕ) : ℚ) = quadraticBenchmarkUpper (2 * m) := by
    rw [benchmark_even m (by omega)]
    simp only [quadraticBenchmarkUpper]
    push_cast
    rw [hcastsub]
    ring
  have hfinal : ((benchmark (2 * m) : ℕ) : ℚ) ≤ C := by
    rw [hbenchmark]
    exact hcoarse.trans hcount
  exact_mod_cast hfinal

theorem odd_large_block_closes
    (m r C : ℕ) (hm : 8 ≤ m) (hr : 2 ≤ r)
    (hrange : 3 * r ≤ 2 * m)
    (hcount : coarseCircleBound (2 * m + 1) r ≤ C) :
    benchmark (2 * m + 1) ≤ C := by
  have hn : (17 : ℚ) ≤ 2 * (m : ℚ) + 1 := by
    exact_mod_cast (by omega : 17 ≤ 2 * m + 1)
  have hrr : (2 : ℚ) ≤ r := by exact_mod_cast hr
  have hrangeq : (3 : ℚ) * r ≤ 2 * m := by exact_mod_cast hrange
  have hcoarse := coarse_circle_bound_ge (2 * (m : ℚ) + 1) (r : ℚ)
    hn hrr (by norm_num; exact hrangeq)
  have hcastsub : ((m - 1 : ℕ) : ℚ) = (m : ℚ) - 1 := by
    rw [Nat.cast_sub (by omega)]
    norm_num
  have hbenchmark :
      ((benchmark (2 * m + 1) : ℕ) : ℚ) ≤
        quadraticBenchmarkUpper (2 * m + 1) := by
    rw [benchmark_odd m (by omega)]
    simp only [quadraticBenchmarkUpper]
    push_cast
    rw [hcastsub]
    nlinarith
  have hfinal : ((benchmark (2 * m + 1) : ℕ) : ℚ) ≤ C :=
    hbenchmark.trans (hcoarse.trans hcount)
  exact_mod_cast hfinal

end LargeN

section N16

def qLocal6 (t2 t3 t4 t5 t6 b2 b3 b4 b5 b6 : ℚ) : ℚ :=
  (t2 - b2) / 3 + (t3 - b3) / 4 + (t4 - b4) / 5 +
    (t5 - b5) / 6 + (t6 - b6) / 7

theorem n16_k7_coefficient_identity
    (t2 t3 t4 t5 t6 b2 b3 b4 b5 b6 : ℚ) :
    qLocal6 t2 t3 t4 t5 t6 b2 b3 b4 b5 b6 + (t6 - b6) / 42 =
      (7 / 108) * (t2 + 3 * t3 + 6 * t4 + 10 * t5 + 15 * t6) +
      (7 / 36) * (t2 - t4 - 2 * t5 - 3 * t6) +
      (1 / 54) * (4 * t2 + 3 * t3 - 5 * t5 - 12 * t6) -
      (1 / 3) * (b2 + b3 + b4 + b5 + b6) +
      (t4 / 180 + b3 / 12 + 2 * b4 / 15 + b5 / 6 + b6 / 6) := by
  simp only [qLocal6]
  ring

theorem n16_k7_local_lower_bound
    (T M B S R Q d : ℚ)
    (hidentity : Q + d / 42 =
      (7 / 108) * T + (7 / 36) * M + (1 / 54) * B - (1 / 3) * S + R)
    (hT : T = 105) (hM : 3 ≤ M) (hB : 60 ≤ B) (hS : S ≤ 7)
    (hR : 0 ≤ R) : 37 / 6 ≤ Q + d / 42 := by
  linarith

theorem n16_k7_equality_profile
    (t2 t3 t4 t5 t6 b3 b4 b5 b6 : ℚ)
    (ht4 : 0 ≤ t4) (hb3 : 0 ≤ b3) (hb4 : 0 ≤ b4)
    (hb5 : 0 ≤ b5) (hb6 : 0 ≤ b6)
    (hT : t2 + 3 * t3 + 6 * t4 + 10 * t5 + 15 * t6 = 105)
    (hM : t2 - t4 - 2 * t5 - 3 * t6 = 3)
    (hB : 4 * t2 + 3 * t3 - 5 * t5 - 12 * t6 = 60)
    (hresidual : t4 / 180 + b3 / 12 + 2 * b4 / 15 +
      b5 / 6 + b6 / 6 = 0)
    (hd : t6 - b6 = 1) :
    t2 = 14 ∧ t3 = 12 ∧ t4 = 0 ∧ t5 = 4 ∧ t6 = 1 ∧
      b3 = 0 ∧ b4 = 0 ∧ b5 = 0 ∧ b6 = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> linarith

theorem n16_seven_circle_family_at_most_four
    (g : ℕ) (hincidence : 7 * g ≤ 48)
    (hmoment : 14 * g ≤ g * (g - 1) + 48) : g ≤ 4 := by
  have hg : g ≤ 6 := by omega
  by_contra h
  have hcases : g = 5 ∨ g = 6 := by omega
  rcases hcases with rfl | rfl <;> norm_num at hmoment

theorem n16_union_arithmetic :
    4 * 6 - Nat.choose 4 2 > 15 ∧
      4 * 5 + 6 - Nat.choose 5 2 > 15 := by
  native_decide

theorem n16_four_circles_force_degree_one
    (d : Fin 16 → ℕ) (hdegree : ∀ x, d x ≤ 3)
    (hincidence : ∑ x, d x = 28)
    (hmoment : ∑ x, Nat.choose (d x) 2 ≤ 12) :
    ∃ x, d x = 1 := by
  by_contra h
  push Not at h
  have hpoint : ∀ x, d x ≤ 2 * Nat.choose (d x) 2 := by
    intro x
    have hx := hdegree x
    have hxne := h x
    have hcases : d x = 0 ∨ d x = 2 ∨ d x = 3 := by omega
    rcases hcases with hd | hd | hd <;> simp [hd]
  have hsumle : (∑ x, d x) ≤ ∑ x, 2 * Nat.choose (d x) 2 :=
    Finset.sum_le_sum fun x _ => hpoint x
  have htwice : (∑ x, 2 * Nat.choose (d x) 2) =
      2 * ∑ x, Nat.choose (d x) 2 := by
    rw [Finset.mul_sum]
  omega

def qLocal7
    (t2 t3 t4 t5 t6 t7 b2 b3 b4 b5 b6 b7 : ℚ) : ℚ :=
  (t2 - b2) / 3 + (t3 - b3) / 4 + (t4 - b4) / 5 +
    (t5 - b5) / 6 + (t6 - b6) / 7 + (t7 - b7) / 8

theorem n16_k8_coefficient_identity
    (t2 t3 t4 t5 t6 t7 b2 b3 b4 b5 b6 b7 : ℚ) :
    qLocal7 t2 t3 t4 t5 t6 t7 b2 b3 b4 b5 b6 b7 =
      (1 / 80) * (t2 + 3 * t3 + 6 * t4 + 10 * t5 + 15 * t6 + 21 * t7) +
      (31 / 160) * (t2 - t4 - 2 * t5 - 3 * t6 - 4 * t7) +
      (1 / 48) * (t2 - b2) -
      (5 / 16) * (b2 + b3 + b4 + b5 + b6 + b7) +
      (17 / 160) * (t2 + 2 * t3 + 3 * t4 + 4 * t5 + 5 * t6 + 6 * t7) +
      (t5 / 240 + 3 * t6 / 560 + b3 / 16 + 9 * b4 / 80 +
        7 * b5 / 48 + 19 * b6 / 112 + 3 * b7 / 16) := by
  simp only [qLocal7]
  ring

theorem n16_k8_local_lower_bound
    (T M O S I R Q : ℚ)
    (hidentity : Q = (1 / 80) * T + (31 / 160) * M + (1 / 48) * O -
      (5 / 16) * S + (17 / 160) * I + R)
    (hT : T = 105) (hM : 3 ≤ M) (hO : 3 ≤ O) (hS : S ≤ 7)
    (hI : 62 ≤ I) (hR : 0 ≤ R) : 1017 / 160 ≤ Q := by
  linarith

theorem n16_final_margin :
    0 < 8 * (1017 / 160 : ℚ) + 8 * (145 / 24 : ℚ) - 98 := by
  norm_num

end N16

section N15

theorem even_coloring_excess_not_one
    (delta epsPlus epsMinus : ℕ)
    (hsum : delta = epsPlus + epsMinus)
    (hcongr : epsPlus % 3 = epsMinus % 3) : delta ≠ 1 := by
  omega

theorem n15_pair_delta_elimination
    (t2 t3 t4 t5 t6 delta : ℤ)
    (hpairs : t2 + 3 * t3 + 6 * t4 + 10 * t5 + 15 * t6 = 91)
    (hdelta : t2 = 3 + t4 + 2 * t5 + 3 * t6 + delta) :
    3 * t3 + 7 * t4 + 12 * t5 + 18 * t6 = 88 - delta := by
  linarith

theorem n15_rich_bound
    (t2 t3 t4 t5 t6 delta : ℤ)
    (helim : 3 * t3 + 7 * t4 + 12 * t5 + 18 * t6 = 88 - delta)
    (hdelta : t2 = 3 + t4 + 2 * t5 + 3 * t6 + delta)
    (hbojanowski : 56 ≤ 4 * t2 + 3 * t3 - 5 * t5 - 12 * t6) :
    3 * t4 + 9 * t5 + 18 * t6 ≤ 44 + 3 * delta := by
  linarith

theorem n15_residual_identity
    (t2 t3 t4 t5 t6 delta : ℤ)
    (helim : 3 * t3 + 7 * t4 + 12 * t5 + 18 * t6 = 88 - delta)
    (hdelta : t2 = 3 + t4 + 2 * t5 + 3 * t6 + delta) :
    136 * t2 + 93 * t3 + 60 * t4 + 30 * t5 - 2902 =
      234 + 105 * delta - 21 * t4 - 70 * t5 - 150 * t6 := by
  linarith

theorem n15_positive_defect_residual
    (t4 t5 t6 delta R : ℤ)
    (ht4 : 0 ≤ t4) (ht5 : 0 ≤ t5)
    (hdelta : 2 ≤ delta)
    (hrich : 3 * t4 + 9 * t5 + 18 * t6 ≤ 44 + 3 * delta)
    (hR : R = 234 + 105 * delta - 21 * t4 - 70 * t5 - 150 * t6) :
    28 ≤ R := by
  have hweight :
      21 * t4 + 70 * t5 + 150 * t6 ≤
        25 * (t4 + 3 * t5 + 6 * t6) := by
    linarith
  omega

def n15SimplicialType (t2 t3 t4 t5 t6 : ℤ) : Prop :=
  (t2 = 11 ∧ t3 = 12 ∧ t4 = 4 ∧ t5 = 2 ∧ t6 = 0) ∨
  (t2 = 9 ∧ t3 = 16 ∧ t4 = 4 ∧ t5 = 1 ∧ t6 = 0) ∨
  (t2 = 10 ∧ t3 = 14 ∧ t4 = 4 ∧ t5 = 0 ∧ t6 = 1)

theorem n15_simplicial_residual_gap
    (t2 t3 t4 t5 t6 : ℤ) (htype : n15SimplicialType t2 t3 t4 t5 t6) :
    let R := 136 * t2 + 93 * t3 + 60 * t4 + 30 * t5 - 2902
    0 ≤ R ∧ (0 < R → R = 10 ∨ 28 ≤ R) := by
  rcases htype with h | h | h <;> rcases h with ⟨rfl, rfl, rfl, rfl, rfl⟩ <;>
    norm_num

theorem n15_local_residual_gap_from_inputs
    (t2 t3 t4 t5 t6 delta R : ℤ)
    (ht4 : 0 ≤ t4) (ht5 : 0 ≤ t5) (hdelta0 : 0 ≤ delta)
    (hdelta1 : delta ≠ 1)
    (hpairs : t2 + 3 * t3 + 6 * t4 + 10 * t5 + 15 * t6 = 91)
    (hdelta : t2 = 3 + t4 + 2 * t5 + 3 * t6 + delta)
    (hbojanowski : 56 ≤ 4 * t2 + 3 * t3 - 5 * t5 - 12 * t6)
    (hR : R = 136 * t2 + 93 * t3 + 60 * t4 + 30 * t5 - 2902)
    (hclassification : delta = 0 → n15SimplicialType t2 t3 t4 t5 t6) :
    0 ≤ R ∧ (0 < R → R = 10 ∨ 28 ≤ R) := by
  have helim := n15_pair_delta_elimination t2 t3 t4 t5 t6 delta hpairs hdelta
  by_cases hz : delta = 0
  · have hs := n15_simplicial_residual_gap t2 t3 t4 t5 t6
      (hclassification hz)
    dsimp only at hs
    rw [hR]
    exact hs
  · have hd2 : 2 ≤ delta := by omega
    have hrich := n15_rich_bound t2 t3 t4 t5 t6 delta
      helim hdelta hbojanowski
    have hRid := n15_residual_identity t2 t3 t4 t5 t6 delta
      helim hdelta
    have hRform : R = 234 + 105 * delta - 21 * t4 - 70 * t5 - 150 * t6 := by
      linarith
    have hlower := n15_positive_defect_residual t4 t5 t6 delta R
      ht4 ht5 hd2 hrich hRform
    exact ⟨by omega, fun _ => Or.inr hlower⟩

theorem n15_global_identity_from_counts
    (C l2 l3 l4 l5 l6 l7 c3 c4 c5 c6 c7 Rsum : ℤ)
    (hC : C = c3 + c4 + c5 + c6 + c7)
    (hpairs : l2 + 3 * l3 + 6 * l4 + 10 * l5 + 15 * l6 + 21 * l7 = 105)
    (htriples : l3 + c3 + 4 * (l4 + c4) + 10 * (l5 + c5) +
      20 * (l6 + c6) + 35 * (l7 + c7) = 455)
    (hRsum : Rsum = 408 * (l3 + c3) + 372 * (l4 + c4) +
      300 * (l5 + c5) + 180 * (l6 + c6) - 15 * 2902) :
    420 * C - 35270 = 140 * (l2 - 7) + 420 * l4 +
      980 * l5 + 1680 * l6 + 2520 * l7 + Rsum := by
  linarith

theorem n15_global_identity_gives_84
    (C l2 l4 l5 l6 l7 Rsum : ℤ)
    (hl2 : 7 ≤ l2) (hl4 : 0 ≤ l4) (hl5 : 0 ≤ l5)
    (hl6 : 0 ≤ l6) (hl7 : 0 ≤ l7) (hR : 0 ≤ Rsum)
    (hid : 420 * C - 35270 = 140 * (l2 - 7) + 420 * l4 +
      980 * l5 + 1680 * l6 + 2520 * l7 + Rsum) :
    84 ≤ C := by
  omega

theorem n15_equality_forces_slack
    (l2 l4 l5 l6 l7 Rsum : ℤ)
    (hl2 : 7 ≤ l2) (hl4 : 0 ≤ l4) (hl5 : 0 ≤ l5)
    (hl6 : 0 ≤ l6) (hl7 : 0 ≤ l7) (hR : 0 ≤ Rsum)
    (hid : 10 = 140 * (l2 - 7) + 420 * l4 +
      980 * l5 + 1680 * l6 + 2520 * l7 + Rsum) :
    l2 = 7 ∧ l4 = 0 ∧ l5 = 0 ∧ l6 = 0 ∧ l7 = 0 ∧ Rsum = 10 := by
  omega

theorem n15_final_congruence_obstruction
    (e : Fin 15 → ℕ) (hsum : ∑ x, e x = 1)
    (hdiv : 3 ∣ ∑ x, (10 + e x)) : False := by
  have hcount : ∑ _x : Fin 15, 10 = 150 := by norm_num
  have htotal : ∑ x, (10 + e x) = 151 := by
    rw [Finset.sum_add_distrib, hcount, hsum]
  rw [htotal] at hdiv
  norm_num at hdiv

end N15

end Erdos506
