import Mathlib

namespace Smooth4PC

/-!
# Reynolds averaging over placements and a cocone of extended rows

This module records the finite algebra behind a statewise "Reynolds cocone"
argument for a cabled system indexed by owner states.  A state is a function
`r : Fin n → ℕ` (the intended application has `n = 5` owners) and states carry
the product order, so `r ≤ r'` means `r i ≤ r' i` for every owner `i`.

For each state `r` the data consists of a rational vector space `M r`, a finite
type `Placement r` of placements, and one row `row r ω : M r →ₗ[ℚ] ℚ` per
placement.  The Reynolds average `reynolds S r` is the arithmetic mean of the
rows over the placements available at `r`.

Three groups of statements are proved.

1. `reynolds_comp_of_row_perm`: if precomposition with `β` permutes the rows at
   a fixed state, then the Reynolds average is unchanged by `β`.
2. `reynolds_comp_dotted` and `reynolds_comp_undotted`: if a pair addition
   pulls the rows at `r'` back to the rows at `r` along a placement map all of
   whose fibres have one and the same positive size `k`, then it pulls the
   Reynolds average at `r'` back to the Reynolds average at `r`; a pair
   addition annihilated by every row at `r'` is annihilated by the average.
3. `extendedRow_eq_comp_of_le`, `extendedRow_eq_of_threshold` and
   `extendedRow_comp_psi`: rows given above a threshold state `s0` and
   descending along pair additions extend to every state, the extension does
   not depend on which common upper bound is used, and the descent relation
   then holds at every state, low states included.

## What is hypothesis and what is proved

Everything geometric enters as a hypothesis of a statement or as a field of a
structure, and none of it is proved here.  In particular this module neither
proves nor asserts:

* that the actual MWW `β` action permutes the placements of a state (this is
  the hypothesis `hbeta` of `reynolds_comp_of_row_perm`);
* that an actual pair addition forgets one distinguished pair, so that the
  induced placement map has fibres all of one and the same size (this is the
  hypothesis `hfib` of `card_placement_of_uniform_fibres`, together with `h1`
  and `h0` relating the rows);
* that the actual pair additions commute with one another, or that they
  assemble into a family of path composites that is functorial and that is
  built from the pair additions themselves (these are the fields `Psi_self`
  and `Psi_succ` of `CableCocone`; the one-step law `CableCocone.Psi_step`,
  the composition law `CableCocone.Psi_comp` and the two-step law
  `CableCocone.Psi_two_steps` are theorems derived from those two, not further
  hypotheses);
* that the actual rows above a threshold satisfy the descent relation (the
  field `descent` of `CableCocone`).

Only the algebraic consequences of those inputs are established below.
-/

noncomputable section

universe u

/-- One pair addition raises the state coordinatewise: `r ≤ r + e i`. -/
theorem le_add_single {n : ℕ} (r : Fin n → ℕ) (i : Fin n) :
    r ≤ r + Pi.single i 1 :=
  Pi.le_def.mpr fun j => by
    rw [Pi.add_apply]
    exact Nat.le_add_right _ _

/-- A cabled system over `n` owners.

`M r` is the cabled summand at the state `r`, `Placement r` is the finite set of
placements available at `r`, and `row r ω` is the row `D_{r,ω}` cut out by the
placement `ω`.  No compatibility between different states is required here;
that is the business of `CableCocone`. -/
structure CableSystem (n : ℕ) where
  /-- The cabled summand at each state. -/
  M : (Fin n → ℕ) → Type u
  [instAddCommGroup : ∀ r, AddCommGroup (M r)]
  [instModule : ∀ r, Module ℚ (M r)]
  /-- The finite set of placements available at each state. -/
  Placement : (Fin n → ℕ) → Type u
  [instFintype : ∀ r, Fintype (Placement r)]
  [instDecidableEq : ∀ r, DecidableEq (Placement r)]
  /-- The row cut out by one placement at one state. -/
  row : ∀ r, Placement r → (M r →ₗ[ℚ] ℚ)

attribute [instance] CableSystem.instAddCommGroup CableSystem.instModule
  CableSystem.instFintype CableSystem.instDecidableEq

variable {n : ℕ}

/-- The Reynolds average of the rows at a state: the arithmetic mean of
`row r ω` over all placements `ω` at `r`. -/
def reynolds (S : CableSystem n) (r : Fin n → ℕ) : S.M r →ₗ[ℚ] ℚ :=
  (Fintype.card (S.Placement r) : ℚ)⁻¹ • ∑ ω : S.Placement r, S.row r ω

/-- Pointwise description of the Reynolds average. -/
theorem reynolds_apply (S : CableSystem n) (r : Fin n → ℕ) (x : S.M r) :
    reynolds S r x
      = (Fintype.card (S.Placement r) : ℚ)⁻¹
          * ∑ ω : S.Placement r, S.row r ω x := by
  simp only [reynolds, LinearMap.smul_apply, LinearMap.sum_apply, smul_eq_mul]

/-! ### Theorem 1: invariance of the average under a placement-permuting map -/

/-- **Beta invariance by placement transport.**  If precomposition with the
linear map `β` carries the row at each placement `ω` to the row at `σ ω`, for
one permutation `σ` of the placements of the state `r`, then `β` leaves the
Reynolds average at `r` unchanged.

That the actual MWW `β` action does permute placements is a hypothesis
(`hbeta`), not a conclusion of this module. -/
theorem reynolds_comp_of_row_perm (S : CableSystem n) (r : Fin n → ℕ)
    (β : S.M r →ₗ[ℚ] S.M r) (σ : Equiv.Perm (S.Placement r))
    (hbeta : ∀ ω : S.Placement r, (S.row r ω).comp β = S.row r (σ ω)) :
    (reynolds S r).comp β = reynolds S r := by
  ext x
  have hpoint : ∀ ω : S.Placement r, S.row r ω (β x) = S.row r (σ ω) x :=
    fun ω => LinearMap.congr_fun (hbeta ω) x
  have hsum : ∑ ω : S.Placement r, S.row r ω (β x)
      = ∑ ω : S.Placement r, S.row r ω x := by
    calc ∑ ω : S.Placement r, S.row r ω (β x)
        = ∑ ω : S.Placement r, S.row r (σ ω) x :=
          Finset.sum_congr rfl fun ω _ => hpoint ω
      _ = ∑ ω : S.Placement r, S.row r ω x :=
          Fintype.sum_equiv σ (fun ω => S.row r (σ ω) x)
            (fun ω => S.row r ω x) (fun ω => rfl)
  rw [LinearMap.comp_apply, reynolds_apply, reynolds_apply, hsum]

/-! ### Theorem 2: pair additions and fibre counting -/

/-- If every fibre of a placement map `π : Placement r' → Placement r` has the
same size `k`, then the placement count at `r'` is `k` times the placement count
at `r`.

That an actual pair addition induces such a placement map is a hypothesis. -/
theorem card_placement_of_uniform_fibres (S : CableSystem n) (r r' : Fin n → ℕ)
    (π : S.Placement r' → S.Placement r) (k : ℕ)
    (hfib : ∀ ω : S.Placement r,
      (Finset.univ.filter (fun ω' => π ω' = ω)).card = k) :
    Fintype.card (S.Placement r') = Fintype.card (S.Placement r) * k := by
  have hfw := Finset.card_eq_sum_card_fiberwise
    (f := π) (s := (Finset.univ : Finset (S.Placement r')))
    (t := (Finset.univ : Finset (S.Placement r))) (fun ω' _ => Finset.mem_univ _)
  simpa [hfib, Finset.card_univ] using hfw

/-- **Pair addition, dotted case.**  Let `ψ1 : M r →ₗ M r'` pull the row at each
placement `ω'` of `r'` back to the row at `π ω'`, where every fibre of
`π : Placement r' → Placement r` has one and the same positive size `k`.  Then
`ψ1` pulls the Reynolds average at `r'` back to the Reynolds average at `r`:
the factor `k` produced by regrouping the sum over the fibres of `π` cancels
against the ratio of the two normalisations. -/
theorem reynolds_comp_dotted (S : CableSystem n) (r r' : Fin n → ℕ)
    (ψ1 : S.M r →ₗ[ℚ] S.M r') (π : S.Placement r' → S.Placement r) (k : ℕ)
    (hk : 0 < k)
    (hfib : ∀ ω : S.Placement r,
      (Finset.univ.filter (fun ω' => π ω' = ω)).card = k)
    (h1 : ∀ ω' : S.Placement r', (S.row r' ω').comp ψ1 = S.row r (π ω')) :
    (reynolds S r').comp ψ1 = reynolds S r := by
  have hcard := card_placement_of_uniform_fibres S r r' π k hfib
  have hk' : (k : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr hk.ne'
  have harith : ∀ N A : ℚ, (N * (k : ℚ))⁻¹ * ((k : ℚ) * A) = N⁻¹ * A := by
    intro N A
    rw [mul_inv]
    calc N⁻¹ * (k : ℚ)⁻¹ * ((k : ℚ) * A)
        = N⁻¹ * ((k : ℚ)⁻¹ * (k : ℚ)) * A := by ring
      _ = N⁻¹ * A := by rw [inv_mul_cancel₀ hk', mul_one]
  ext x
  have hstep : ∀ ω' : S.Placement r', S.row r' ω' (ψ1 x) = S.row r (π ω') x :=
    fun ω' => LinearMap.congr_fun (h1 ω') x
  have hinner : ∀ ω : S.Placement r,
      ∑ ω' ∈ Finset.univ.filter (fun ω' => π ω' = ω), S.row r (π ω') x
        = (k : ℚ) * S.row r ω x := by
    intro ω
    calc ∑ ω' ∈ Finset.univ.filter (fun ω' => π ω' = ω), S.row r (π ω') x
        = ∑ _ω' ∈ Finset.univ.filter (fun ω' => π ω' = ω), S.row r ω x :=
          Finset.sum_congr rfl fun ω' hω' => by
            rw [(Finset.mem_filter.mp hω').2]
      _ = (Finset.univ.filter (fun ω' => π ω' = ω)).card • S.row r ω x :=
          Finset.sum_const _
      _ = (k : ℚ) * S.row r ω x := by rw [hfib ω, nsmul_eq_mul]
  have hsum : ∑ ω' : S.Placement r', S.row r' ω' (ψ1 x)
      = (k : ℚ) * ∑ ω : S.Placement r, S.row r ω x := by
    calc ∑ ω' : S.Placement r', S.row r' ω' (ψ1 x)
        = ∑ ω' : S.Placement r', S.row r (π ω') x :=
          Finset.sum_congr rfl fun ω' _ => hstep ω'
      _ = ∑ ω : S.Placement r,
            ∑ ω' ∈ Finset.univ.filter (fun ω' => π ω' = ω), S.row r (π ω') x :=
          (Finset.sum_fiberwise_of_maps_to
            (fun ω' _ => Finset.mem_univ (π ω')) _).symm
      _ = ∑ ω : S.Placement r, (k : ℚ) * S.row r ω x :=
          Finset.sum_congr rfl fun ω _ => hinner ω
      _ = (k : ℚ) * ∑ ω : S.Placement r, S.row r ω x := by rw [Finset.mul_sum]
  rw [LinearMap.comp_apply, reynolds_apply, reynolds_apply, hsum, hcard]
  push_cast
  exact harith _ _

/-- **Pair addition, undotted case.**  A map annihilated by every row at `r'` is
annihilated by the Reynolds average at `r'`. -/
theorem reynolds_comp_undotted (S : CableSystem n) (r r' : Fin n → ℕ)
    (ψ0 : S.M r →ₗ[ℚ] S.M r')
    (h0 : ∀ ω' : S.Placement r', (S.row r' ω').comp ψ0 = 0) :
    (reynolds S r').comp ψ0 = 0 := by
  ext x
  have hzero : ∀ ω' : S.Placement r', S.row r' ω' (ψ0 x) = 0 := by
    intro ω'
    have hval := LinearMap.congr_fun (h0 ω') x
    simpa using hval
  rw [LinearMap.comp_apply, reynolds_apply, LinearMap.zero_apply,
    Finset.sum_eq_zero (fun ω' _ => hzero ω'), mul_zero]

/-! ### Theorem 3: low states and independence of the chosen upper bound -/

/-- A cabled system together with once-dotted pair additions `psi r i` and a
family of path composites `Psi r r'` for `r ≤ r'`, plus rows `D` descending
along the pair additions above a threshold state `s0`.

This is the **second** version of the third theorem offered in the
specification: the path composite `Psi` is supplied as a field satisfying the
functoriality laws `Psi r r = id`, `Psi r' r'' ∘ Psi r r' = Psi r r''` and
`Psi r (r + e i) = psi r i`, rather than being built by hand from `psi` along a
chosen coordinate order.  Commutation of the pair additions is exactly what the
existence of such a `Psi` encodes, and it is an input here, not a theorem.

The only laws assumed of `Psi` are `Psi_self` (the empty path is the identity)
and `Psi_succ` (lengthening a path on top applies one more pair addition).
Together they pin every path composite down to an iterate of `psi`: the
one-step law `Psi_step`, the general composition law `Psi_comp` and the
two-step law `Psi_two_steps` are proved from them below.

`D` is given at every state, but only its values at states above `s0` are ever
constrained or used: `descent` speaks only about states `r ≥ s0`, and
`extendedRow` evaluates `D` only at states of the form `r ⊔ s0 ≥ s0`. -/
structure CableCocone (n : ℕ) (S : CableSystem n) where
  /-- The once-dotted pair addition raising the `i`-th owner state by one. -/
  psi : ∀ (r : Fin n → ℕ) (i : Fin n), S.M r →ₗ[ℚ] S.M (r + Pi.single i 1)
  /-- The composite along a path of pair additions from `r` up to `r'`. -/
  Psi : ∀ (r r' : Fin n → ℕ), r ≤ r' → (S.M r →ₗ[ℚ] S.M r')
  /-- The empty path is the identity. -/
  Psi_self : ∀ r : Fin n → ℕ, Psi r r le_rfl = LinearMap.id
  /-- Lengthening a path by one pair addition on top applies that pair
  addition after the shorter path.  This pins every path composite down to an
  iterate of `psi`; the one-step law, the composition law and the two-step law
  are all derived from it below rather than assumed. -/
  Psi_succ : ∀ (r r' : Fin n → ℕ) (h : r ≤ r') (i : Fin n),
    Psi r (r' + Pi.single i 1) (h.trans (le_add_single r' i))
      = (psi r' i).comp (Psi r r' h)
  /-- The threshold state above which the rows are given. -/
  s0 : Fin n → ℕ
  /-- The rows; only their values at states above `s0` are used. -/
  D : ∀ r : Fin n → ℕ, S.M r →ₗ[ℚ] ℚ
  /-- Descent of the rows along pair additions, above the threshold. -/
  descent : ∀ (r : Fin n → ℕ), s0 ≤ r → ∀ i : Fin n,
    (D (r + Pi.single i 1)).comp (psi r i) = D r

variable {S : CableSystem n}

/-- A one-step path is the corresponding pair addition.  Derived from
`Psi_succ` at the empty path together with `Psi_self`, so it is not assumed. -/
theorem CableCocone.Psi_step (C : CableCocone n S) (r : Fin n → ℕ) (i : Fin n) :
    C.Psi r (r + Pi.single i 1) (le_add_single r i) = C.psi r i := by
  have hstep := C.Psi_succ r r le_rfl i
  rw [C.Psi_self r, LinearMap.comp_id] at hstep
  exact hstep

/-- Composition along a path whose upper half is empty. -/
theorem CableCocone.Psi_comp_of_eq (C : CableCocone n S) (r r2 r3 : Fin n → ℕ)
    (h : r ≤ r2) (h2 : r2 ≤ r3) (heq : r2 = r3) :
    (C.Psi r2 r3 h2).comp (C.Psi r r2 h) = C.Psi r r3 (h.trans h2) := by
  subst heq
  rw [show C.Psi _ _ h2 = LinearMap.id from C.Psi_self _, LinearMap.id_comp]

/-- Induction carrier for `CableCocone.Psi_comp`: the composition law holds once
the upper half of the path is at most `N` pair additions long. -/
theorem CableCocone.Psi_comp_aux (C : CableCocone n S) :
    ∀ (N : ℕ) (r r2 r3 : Fin n → ℕ) (h : r ≤ r2) (h2 : r2 ≤ r3),
      (∑ m, (r3 m - r2 m)) ≤ N →
      (C.Psi r2 r3 h2).comp (C.Psi r r2 h) = C.Psi r r3 (h.trans h2) := by
  intro N
  induction N with
  | zero =>
      intro r r2 r3 h h2 hN
      refine C.Psi_comp_of_eq r r2 r3 h h2 (le_antisymm h2 ?_)
      refine Pi.le_def.mpr fun m => ?_
      have hle : r3 m - r2 m ≤ ∑ m2, (r3 m2 - r2 m2) :=
        Finset.single_le_sum (f := fun m2 => r3 m2 - r2 m2)
          (fun _ _ => Nat.zero_le _) (Finset.mem_univ m)
      omega
  | succ N ih =>
      intro r r2 r3 h h2 hN
      by_cases heq : r2 = r3
      · exact C.Psi_comp_of_eq r r2 r3 h h2 heq
      obtain ⟨j, hj⟩ := Function.ne_iff.mp heq
      have h2m := Pi.le_def.mp h2
      have hlt : r2 j < r3 j := lt_of_le_of_ne (h2m j) hj
      obtain ⟨w, hw⟩ : ∃ w : Fin n → ℕ, w + Pi.single j 1 = r3 := by
        refine ⟨fun m => r3 m - (Pi.single j 1 : Fin n → ℕ) m, ?_⟩
        funext m
        simp only [Pi.add_apply]
        by_cases hm : m = j
        · subst hm
          rw [Pi.single_eq_same]
          omega
        · rw [Pi.single_eq_of_ne hm]
          omega
      subst hw
      have hjval : (w + Pi.single j 1 : Fin n → ℕ) j = w j + 1 := by
        rw [Pi.add_apply, Pi.single_eq_same]
      have hrw : r2 ≤ w := by
        refine Pi.le_def.mpr fun m => ?_
        have hm := h2m m
        rw [Pi.add_apply] at hm
        by_cases hmj : m = j
        · subst hmj
          rw [hjval] at hlt
          omega
        · rw [Pi.single_eq_of_ne hmj] at hm
          omega
      have hrwj : r2 j ≤ w j := Pi.le_def.mp hrw j
      have hdist : (∑ m : Fin n, (w m - r2 m)) ≤ N := by
        have e1 : (w j + 1 - r2 j)
              + ∑ m ∈ Finset.univ.erase j,
                  ((w + Pi.single j 1 : Fin n → ℕ) m - r2 m)
            = ∑ m : Fin n, ((w + Pi.single j 1 : Fin n → ℕ) m - r2 m) := by
          rw [← hjval]
          exact Finset.add_sum_erase Finset.univ
            (fun m => (w + Pi.single j 1 : Fin n → ℕ) m - r2 m) (Finset.mem_univ j)
        have e2 : (w j - r2 j) + ∑ m ∈ Finset.univ.erase j, (w m - r2 m)
            = ∑ m : Fin n, (w m - r2 m) :=
          Finset.add_sum_erase Finset.univ (fun m => w m - r2 m) (Finset.mem_univ j)
        have e3 : ∑ m ∈ Finset.univ.erase j,
              ((w + Pi.single j 1 : Fin n → ℕ) m - r2 m)
            = ∑ m ∈ Finset.univ.erase j, (w m - r2 m) := by
          refine Finset.sum_congr rfl fun m hm => ?_
          rw [Pi.add_apply, Pi.single_eq_of_ne (Finset.ne_of_mem_erase hm),
            Nat.add_zero]
        omega
      have hIH := ih r r2 w h hrw hdist
      have hs1 := C.Psi_succ r2 w hrw j
      have hs2 := C.Psi_succ r w (h.trans hrw) j
      ext x
      have p1 := LinearMap.congr_fun hs1 (C.Psi r r2 h x)
      have p2 := LinearMap.congr_fun hs2 x
      have p3 := LinearMap.congr_fun hIH x
      simp only [LinearMap.comp_apply] at p1 p2 p3
      show C.Psi r2 (w + Pi.single j 1) (hrw.trans (le_add_single w j))
            (C.Psi r r2 h x)
          = C.Psi r (w + Pi.single j 1)
              ((h.trans hrw).trans (le_add_single w j)) x
      rw [p1, p3, ← p2]

/-- **Path composites compose.**  Derived from `Psi_self` and `Psi_succ` by
induction on the number of pair additions in the upper half of the path, so the
composition law is not an independent hypothesis. -/
theorem CableCocone.Psi_comp (C : CableCocone n S) (r r2 r3 : Fin n → ℕ)
    (h : r ≤ r2) (h2 : r2 ≤ r3) :
    (C.Psi r2 r3 h2).comp (C.Psi r r2 h) = C.Psi r r3 (h.trans h2) :=
  C.Psi_comp_aux (∑ m, (r3 m - r2 m)) r r2 r3 h h2 le_rfl

/-- **Two pair additions in a row.**  The length-two path composite is the
composite of the two pair additions along it, in order. -/
theorem CableCocone.Psi_two_steps (C : CableCocone n S) (r : Fin n → ℕ)
    (i j : Fin n) :
    C.Psi r (r + Pi.single i 1 + Pi.single j 1)
        ((le_add_single r i).trans (le_add_single (r + Pi.single i 1) j))
      = (C.psi (r + Pi.single i 1) j).comp (C.psi r i) := by
  rw [C.Psi_succ r (r + Pi.single i 1) (le_add_single r i) j, C.Psi_step r i]

/-- Descent along the empty path. -/
theorem row_comp_Psi_of_eq (C : CableCocone n S) (a b : Fin n → ℕ) (hab : a ≤ b)
    (heq : a = b) : (C.D b).comp (C.Psi a b hab) = C.D a := by
  subst heq
  rw [show C.Psi _ _ hab = LinearMap.id from C.Psi_self _, LinearMap.comp_id]

/-- **Descent along an arbitrary increasing path.**  Above the threshold the
rows descend along one pair addition by hypothesis; by induction on the total
number of pair additions separating two states they descend along every
increasing path. -/
theorem row_descent_along_path (C : CableCocone n S) :
    ∀ (N : ℕ) (a b : Fin n → ℕ) (hab : a ≤ b),
      (∑ j, (b j - a j)) ≤ N → C.s0 ≤ a →
      (C.D b).comp (C.Psi a b hab) = C.D a := by
  intro N
  induction N with
  | zero =>
      intro a b hab hN _
      refine row_comp_Psi_of_eq C a b hab (le_antisymm hab ?_)
      refine Pi.le_def.mpr fun j => ?_
      have hle : b j - a j ≤ ∑ j' , (b j' - a j') :=
        Finset.single_le_sum (f := fun j' => b j' - a j')
          (fun _ _ => Nat.zero_le _) (Finset.mem_univ j)
      omega
  | succ N ih =>
      intro a b hab hN hs
      by_cases heq : a = b
      · exact row_comp_Psi_of_eq C a b hab heq
      obtain ⟨i, hi⟩ := Function.ne_iff.mp heq
      have hab' := Pi.le_def.mp hab
      have hlt : a i < b i := lt_of_le_of_ne (hab' i) hi
      have hstep : a + Pi.single i 1 ≤ b := by
        refine Pi.le_def.mpr fun j => ?_
        rw [Pi.add_apply]
        by_cases hj : j = i
        · subst hj
          rw [Pi.single_eq_same]
          omega
        · rw [Pi.single_eq_of_ne hj]
          have hbound := hab' j
          omega
      have hdist :
          (∑ j : Fin n, (b j - (a + Pi.single i 1 : Fin n → ℕ) j)) ≤ N := by
        have e1 : (b i - a i)
              + ∑ j ∈ Finset.univ.erase i, (b j - a j)
            = ∑ j : Fin n, (b j - a j) :=
          Finset.add_sum_erase Finset.univ (fun j => b j - a j)
            (Finset.mem_univ i)
        have e2 : (b i - (a + Pi.single i 1 : Fin n → ℕ) i)
              + ∑ j ∈ Finset.univ.erase i,
                  (b j - (a + Pi.single i 1 : Fin n → ℕ) j)
            = ∑ j : Fin n, (b j - (a + Pi.single i 1 : Fin n → ℕ) j) :=
          Finset.add_sum_erase Finset.univ
            (fun j => b j - (a + Pi.single i 1 : Fin n → ℕ) j)
            (Finset.mem_univ i)
        have e3 : ∑ j ∈ Finset.univ.erase i,
              (b j - (a + Pi.single i 1 : Fin n → ℕ) j)
            = ∑ j ∈ Finset.univ.erase i, (b j - a j) := by
          refine Finset.sum_congr rfl fun j hj => ?_
          rw [Pi.add_apply, Pi.single_eq_of_ne (Finset.ne_of_mem_erase hj),
            Nat.add_zero]
        have e4 : (a + Pi.single i 1 : Fin n → ℕ) i = a i + 1 := by
          rw [Pi.add_apply, Pi.single_eq_same]
        rw [e4] at e2
        omega
      have hs' : C.s0 ≤ a + Pi.single i 1 := hs.trans (le_add_single a i)
      have hrec := ih (a + Pi.single i 1) b hstep hdist hs'
      ext x
      have hc := LinearMap.congr_fun
        (C.Psi_comp a (a + Pi.single i 1) b (le_add_single a i) hstep) x
      have hr := LinearMap.congr_fun hrec
        (C.Psi a (a + Pi.single i 1) (le_add_single a i) x)
      have hd := LinearMap.congr_fun (C.descent a hs i) x
      have hp := LinearMap.congr_fun (C.Psi_step a i) x
      simp only [LinearMap.comp_apply] at hc hr hd ⊢
      calc C.D b (C.Psi a b hab x)
          = C.D b (C.Psi (a + Pi.single i 1) b hstep
              (C.Psi a (a + Pi.single i 1) (le_add_single a i) x)) := by rw [hc]
        _ = C.D (a + Pi.single i 1)
              (C.Psi a (a + Pi.single i 1) (le_add_single a i) x) := hr
        _ = C.D (a + Pi.single i 1) (C.psi a i x) := by rw [hp]
        _ = C.D a x := hd

/-- The extended row at an arbitrary state: push the state up to `r ⊔ s0`, which
lies above the threshold, and pull the given row there back along the path
composite. -/
def extendedRow (C : CableCocone n S) (r : Fin n → ℕ) : S.M r →ₗ[ℚ] ℚ :=
  (C.D (r ⊔ C.s0)).comp (C.Psi r (r ⊔ C.s0) le_sup_left)

/-- **(b) Independence of the chosen common upper bound.**  For every state `r`
and every `r'` that lies above both `r` and the threshold, the extended row at
`r` is the given row at `r'` pulled back along the path composite. -/
theorem extendedRow_eq_comp_of_le (C : CableCocone n S) (r r' : Fin n → ℕ)
    (hr : r ≤ r') (hs : C.s0 ≤ r') :
    extendedRow C r = (C.D r').comp (C.Psi r r' hr) := by
  have hsup : r ⊔ C.s0 ≤ r' := sup_le hr hs
  have hrsup : r ≤ r ⊔ C.s0 := le_sup_left
  have hdesc : (C.D r').comp (C.Psi (r ⊔ C.s0) r' hsup) = C.D (r ⊔ C.s0) :=
    row_descent_along_path C (∑ j, (r' j - (r ⊔ C.s0) j)) (r ⊔ C.s0) r' hsup
      le_rfl le_sup_right
  ext x
  have h1 := LinearMap.congr_fun hdesc (C.Psi r (r ⊔ C.s0) hrsup x)
  have h2 := LinearMap.congr_fun (C.Psi_comp r (r ⊔ C.s0) r' hrsup hsup) x
  simp only [LinearMap.comp_apply] at h1 h2
  show C.D (r ⊔ C.s0) (C.Psi r (r ⊔ C.s0) hrsup x) = C.D r' (C.Psi r r' hr x)
  rw [← h1, h2]

/-- **(a) The extension agrees with the given rows above the threshold.** -/
theorem extendedRow_eq_of_threshold (C : CableCocone n S) (r : Fin n → ℕ)
    (hs : C.s0 ≤ r) : extendedRow C r = C.D r := by
  rw [extendedRow_eq_comp_of_le C r r le_rfl hs,
    show C.Psi r r le_rfl = LinearMap.id from C.Psi_self r, LinearMap.comp_id]

/-- **(c) The descent relation holds at every state.**  Pulling the extended row
at `r + e i` back along the pair addition `psi r i` returns the extended row at
`r`, with no assumption relating `r` to the threshold: the low states are
covered as well. -/
theorem extendedRow_comp_psi (C : CableCocone n S) (r : Fin n → ℕ) (i : Fin n) :
    (extendedRow C (r + Pi.single i 1)).comp (C.psi r i) = extendedRow C r := by
  have hA : extendedRow C (r + Pi.single i 1)
      = (C.D ((r + Pi.single i 1) ⊔ C.s0)).comp
          (C.Psi (r + Pi.single i 1) ((r + Pi.single i 1) ⊔ C.s0) le_sup_left) :=
    extendedRow_eq_comp_of_le C (r + Pi.single i 1)
      ((r + Pi.single i 1) ⊔ C.s0) le_sup_left le_sup_right
  have hB : extendedRow C r
      = (C.D ((r + Pi.single i 1) ⊔ C.s0)).comp
          (C.Psi r ((r + Pi.single i 1) ⊔ C.s0)
            ((le_add_single r i).trans le_sup_left)) :=
    extendedRow_eq_comp_of_le C r ((r + Pi.single i 1) ⊔ C.s0)
      ((le_add_single r i).trans le_sup_left) le_sup_right
  rw [hA, hB, ← C.Psi_step r i]
  ext x
  have hcomp := LinearMap.congr_fun
    (C.Psi_comp r (r + Pi.single i 1) ((r + Pi.single i 1) ⊔ C.s0)
      (le_add_single r i) le_sup_left) x
  simp only [LinearMap.comp_apply] at hcomp ⊢
  rw [hcomp]

end

end Smooth4PC
