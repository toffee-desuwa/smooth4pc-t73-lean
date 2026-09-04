import Mathlib

namespace Smooth4PC

/-!
# Filtered cubic naturality

An `h`-adic formal series is modelled here by its coefficient sequence
`ℕ → _`.  This is the filtered model: no completion and no `PowerSeries`
structure is used, so every statement below is a finite identity between the
degree-`n` coefficients of a product of three series.

* a vector series `u : ℕ → E` stands for `u(h) = ∑ (u k) h ^ k`;
* an operator series `A : ℕ → (E →ₗ[ℚ] E)` stands for `A(h) = ∑ (A j) h ^ j`;
* a row series `ell : ℕ → (E →ₗ[ℚ] ℚ)` stands for `ℓ(h) = ∑ (ell i) h ^ i`.

Two facts used by the T73 endpoint-transport argument are recorded.

1. If the operator series starts in degree three then the degree-three
   coefficient of `ℓ(h) A(h) u(h)` is `ell 0 (A 3 (u 0))`, and the degree-zero,
   degree-one and degree-two coefficients vanish.
2. Conjugating the operator series, the vector series and the row series by one
   and the same `h`-independent linear automorphism leaves every coefficient of
   `ℓ(h) A(h) u(h)` unchanged.
-/

noncomputable section

variable {E : Type*} [AddCommGroup E] [Module ℚ E]

/-- Degree-`n` coefficient of the scalar series `ℓ(h) A(h) u(h)`.

The Cauchy product is written as a nested finite sum: the outer index `i` runs
over `Finset.range (n + 1)`, the inner index `j` runs over
`Finset.range (n + 1 - i)`, and the remaining index is forced to be
`k = n - i - j`.  The pairs `(i, j)` in that nested range are exactly the pairs
with `i + j ≤ n`, so the sum ranges exactly over the triples `(i, j, k)` of
natural numbers with `i + j + k = n`. -/
def pairingCoeff (ell : ℕ → (E →ₗ[ℚ] ℚ)) (A : ℕ → (E →ₗ[ℚ] E)) (u : ℕ → E)
    (n : ℕ) : ℚ :=
  ∑ i ∈ Finset.range (n + 1), ∑ j ∈ Finset.range (n + 1 - i),
    ell i (A j (u (n - i - j)))

/-- The operator series carries no term below degree three, i.e.
`A(h) ∈ h ^ 3 * End(E)`. -/
def StartsAtThree (A : ℕ → (E →ₗ[ℚ] E)) : Prop :=
  A 0 = 0 ∧ A 1 = 0 ∧ A 2 = 0

/-- Every coefficient of `ℓ(h) A(h) u(h)` below degree three vanishes when the
operator series starts in degree three: each surviving summand carries an
operator index `j ≤ n < 3`. -/
theorem pairingCoeff_lt_three_eq_zero
    (ell : ℕ → (E →ₗ[ℚ] ℚ)) (A : ℕ → (E →ₗ[ℚ] E)) (u : ℕ → E)
    (hA : StartsAtThree A) (n : ℕ) (hn : n < 3) :
    pairingCoeff ell A u n = 0 := by
  obtain ⟨h0, h1, h2⟩ := hA
  interval_cases n <;>
    simp [pairingCoeff, Finset.sum_range_succ, h0, h1, h2]

/-- The divided cubic of the pairing.

If `A(h) = ρ(W) - I` lies in `h ^ 3 * End(E)`, if `u(h) = u₀ + O(h)` is the cup
and `ℓ(h) = ℓ₀ + O(h)` is the cap, then the degree-three coefficient
`[h ^ 3] ℓ(h) (ρ(W) - I) u(h)` equals `ℓ₀ A₃ u₀`.  In particular it does not
depend on the positive-order corrections of the cup, of the cap, or of the
basis change: those corrections can only pair against `A 0`, `A 1` or `A 2`,
which vanish. -/
theorem pairingCoeff_three_of_startsAtThree
    (ell : ℕ → (E →ₗ[ℚ] ℚ)) (A : ℕ → (E →ₗ[ℚ] E)) (u : ℕ → E)
    (hA : StartsAtThree A) :
    pairingCoeff ell A u 3 = ell 0 (A 3 (u 0)) := by
  obtain ⟨h0, h1, h2⟩ := hA
  simp [pairingCoeff, Finset.sum_range_succ, h0, h1, h2]

/-- Conjugation of an operator series by an `h`-independent automorphism:
`A(h) ↦ P A(h) P⁻¹`. -/
def transportOperator (P : E ≃ₗ[ℚ] E) (A : ℕ → (E →ₗ[ℚ] E)) :
    ℕ → (E →ₗ[ℚ] E) :=
  fun j => (P.toLinearMap.comp (A j)).comp P.symm.toLinearMap

/-- Transport of a vector series by an `h`-independent automorphism:
`u(h) ↦ P u(h)`. -/
def transportVector (P : E ≃ₗ[ℚ] E) (u : ℕ → E) : ℕ → E :=
  fun k => P (u k)

/-- Transport of a row series by an `h`-independent automorphism:
`ℓ(h) ↦ ℓ(h) P⁻¹`. -/
def transportRow (P : E ≃ₗ[ℚ] E) (ell : ℕ → (E →ₗ[ℚ] ℚ)) :
    ℕ → (E →ₗ[ℚ] ℚ) :=
  fun i => (ell i).comp P.symm.toLinearMap

/-- Simultaneous transport of the row, the operator and the vector by one and
the same `h`-independent automorphism leaves every coefficient of the pairing
unchanged, term by term: the inner `P⁻¹ P` and `P P⁻¹` cancel inside each
summand `ell i (A j (u k))`. -/
theorem pairingCoeff_transport
    (P : E ≃ₗ[ℚ] E) (ell : ℕ → (E →ₗ[ℚ] ℚ)) (A : ℕ → (E →ₗ[ℚ] E))
    (u : ℕ → E) (n : ℕ) :
    pairingCoeff (transportRow P ell) (transportOperator P A)
        (transportVector P u) n = pairingCoeff ell A u n := by
  simp [pairingCoeff, transportRow, transportOperator, transportVector]

/-- Conjugation by an `h`-independent automorphism preserves the property of
starting in degree three. -/
theorem startsAtThree_transport
    (P : E ≃ₗ[ℚ] E) (A : ℕ → (E →ₗ[ℚ] E)) (hA : StartsAtThree A) :
    StartsAtThree (transportOperator P A) := by
  obtain ⟨h0, h1, h2⟩ := hA
  refine ⟨?_, ?_, ?_⟩ <;> simp [transportOperator, h0, h1, h2]

/-- The divided cubic is invariant under a simultaneous change of basis: the
transported data still starts in degree three, and its degree-three pairing
coefficient is the untransported `ℓ₀ A₃ u₀`. -/
theorem cubic_invariant_under_simultaneous_transport
    (P : E ≃ₗ[ℚ] E) (ell : ℕ → (E →ₗ[ℚ] ℚ)) (A : ℕ → (E →ₗ[ℚ] E))
    (u : ℕ → E) (hA : StartsAtThree A) :
    pairingCoeff (transportRow P ell) (transportOperator P A)
        (transportVector P u) 3 = ell 0 (A 3 (u 0)) := by
  rw [pairingCoeff_transport P ell A u 3]
  exact pairingCoeff_three_of_startsAtThree ell A u hA

/-- Source-side transport by a series `Q(h) = ∑ (Q m) h ^ m` acting on the
vector series: `(Q u)(h) = Q(h) u(h)`, whose degree-`k` coefficient is
`∑_{m ≤ k} Q m (u (k - m))`. -/
def transportVectorSeries (Q : ℕ → (E →ₗ[ℚ] E)) (u : ℕ → E) : ℕ → E :=
  fun k => ∑ m ∈ Finset.range (k + 1), Q m (u (k - m))

/-- Degree-zero coefficient of a source-side series transport whose degree-zero
term is the identity. -/
theorem transportVectorSeries_zero_of_id
    (Q : ℕ → (E →ₗ[ℚ] E)) (u : ℕ → E) (hQ : Q 0 = LinearMap.id) :
    transportVectorSeries Q u 0 = u 0 := by
  simp [transportVectorSeries, hQ]

/-- Mirror of `cubicComposite_right_identity_transport`: a source-side
transport by a series that is the identity in degree zero does not alter the
divided cubic of the pairing. -/
theorem pairingCoeff_three_transportVectorSeries
    (ell : ℕ → (E →ₗ[ℚ] ℚ)) (A : ℕ → (E →ₗ[ℚ] E)) (u : ℕ → E)
    (Q : ℕ → (E →ₗ[ℚ] E)) (hA : StartsAtThree A) (hQ : Q 0 = LinearMap.id) :
    pairingCoeff ell A (transportVectorSeries Q u) 3 = ell 0 (A 3 (u 0)) := by
  rw [pairingCoeff_three_of_startsAtThree ell A (transportVectorSeries Q u) hA,
    transportVectorSeries_zero_of_id Q u hQ]

end

end Smooth4PC
