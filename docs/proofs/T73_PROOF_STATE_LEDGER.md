# T73 proof state ledger

**Current authority:** the conditional paper together with
`docs/proofs/T73_GEOMETRIC_PREMISE_STATUS_20260903.md`.  Computational P0
is discharged for the Johnson replacement.  C is discharged for the
replacement collar.  S is discharged for the Johnson replacement reversed
1-handle picture.  The MWW four-handle layer is discharged.  The
identification `X_J ≅ Sigma_A^0` is discharged as the constructed CS
handle picture.  Lean ExternalGeometry remains uninhabited.  Entries
below describing unavailable historical PD/TH objects are recovery history,
not current premises.  Frozen cubic is `2624`.

**Updated:** 2026-09-03
**Purpose:** recovery source of truth for the mathematical chain. Read this
before rerunning any long sphere computation.

## Current carrying class

```text
class:        v_T = eta_R[T_1], T_1 = B_act^{-1} U_(0,5)
not class:    xi
not class:    mixed Z
base value:   D_3(v_T) = 2624
final degree: (homological, quantum) = (0,494)
scope:        h^3 associated graded / ordinary divided functional
```

The outer detector contains one factor `rho_h(W)-I`. Applying it to `xi`
would introduce a second factor and starts in order six. The mixed-Z branch
has an actual raw lift but is contained in the full action-closed two-cup
`psi0` ideal and is retired as the carrying class.

## Layer audit

| layer | current status | current evidence |
|---|---|---|
| E1 paired-annulus Hattori coefficient | PROVED | [`ACTUAL_PD_CABLE_UNIT_CERT.json`](../../evidence/public_geometry/ACTUAL_PD_CABLE_UNIT_CERT.json), SHA `7F3D3618D6A790A9B60EE8085B647AC2AB742E1BC9C15841F1BEF015034217B5` |
| E2 selected diagonal class | PROVED | main proof Sections 2--3 |
| E3 coefficient q-trace and raw-state binding | PROVED | `MR_REES_QUANTUM_TRACE_SOURCE_LEDGER.md`, SHA `99913F0AF70CF5FF650492C2740F98B7053FEBBB58212DA32B506D863836EA84`; `RAW_STATE_BINDING_RESULT.md`, SHA `BCBFBC9F8A5350D323EC2876352EC25C6690E5BB6AD4F2C264441BA4D73CDEF6` |
| E5/E6 one-/two-handle beta/psi quotient | PROVED at divided cubic | `T73_EVIDENCE_ONE_CUP_E5_E6.md`, source SHA `6A63978D734EFF30B8F6C2E6F9800F9338B3499BA1CA4D6AF4880716468FE865` |
| G1 actual-Gompf/DIAGRAM framed lift | OPEN | historical DIAGRAM lift unavailable; Johnson P0 is the replacement presentation |
| E8 three chosen sphere rows | PROVED at divided cubic | `T73_EVIDENCE_E8_CHOSEN_SPHERES.md`, source SHA `C33AF1FCAE8F0056A75D1841E1151786C2EBB756AD8A768C4E3E6704274E2B43` |

The E5/E6 proof uses the one-cup through-86 cell, the ordinary divided
functional, beta's `O(h)` defect against `W-I=O(h^3)`, and split-injectivity of
the dotted direct system after quotienting the undotted image. It does not use
the retired owner-only mixed-Z model.

The first G1 PASS was premature because it treated `Ae_i` as a complete
attaching circle.  The later slogan used uniqueness of regular neighborhoods
to place a detector ball after the word identities `t empty` and `x to z`.
That uniqueness argument is not used.  The Johnson replacement now has a
certified discrete Voronoi Heegaard pair, two local product cancellation
movies, and a reconstructed 44-strand braid recovering the public word.

## Sphere evidence chronology

The early `OPEN` files are historical diagnostics, not the latest state.

| local time, 2026-08-31 | artifact | status |
|---|---|---|
| 02:01--02:02 | `THREE_SPHERE_MAP_LEVEL_FINAL_VERDICT.md`; `ACTUAL_THREE_SPHERE_H0_FINAL_DECISION.md` | EARLY OPEN; superseded by later chosen-sphere artifacts |
| 12:10--12:16 | TH1 chosen decision and hardened hostile replay | PASS |
| 12:13--12:18 | TH2 chosen decision and hardened hostile replay | PASS |
| 12:35 | THXY full macro successor and re-audit | PASS |
| 12:48 | `FULL_CHAIN_INDEPENDENT_FINAL.md` | later combined audit; its old `xi` wording is not consumed |
| 20:09 | global chain certificates | corroborating combined records; current upstream class is replaced by `v_T` proof |

### Do not rerun before checking these identities

```text
TH1 geometry
EE620E6B085A5F9E1C73CFDD1AD04FC0682CEC74DA3DBF8AFE70DD19C038E3A0
  evidence/public_geometry/TH1_EXHAUSTIVE_ALL_ROW_GEOMETRY.json

TH2 geometry
4D1B627C0343A1C464319704EAFADCA127902A2E4E90CF1C283004359B7ADC24
  evidence/public_geometry/TH2_EXHAUSTIVE_ALL_ROW_GEOMETRY.json

THXY geometry
EABF67C0D1AE309F0281297710A283B8ED5A11D88C8E810837932313F4C53227
  evidence/public_geometry/THXY_FULL_MACRO_P3FREE_HJ_CERT.json
```

The three chosen spheres lie in one actual `W2`, in mutually disjoint movie
sectors, and their class matrix has determinant `+1`. HJ's theorem relates the
chosen system to the historical attaching system by isotopy, permutation and
3--3 handle slides. A 3--3 handle slide gives an equivalent upper cobordism
relative to incoming `W2`; it is not an arbitrary boundary-diffeomorphism
substitution.

For each chosen sphere, every noninvertible event is new--new and the old
one-cup block is an identity cylinder. Actual core counits therefore give the
whole-old-source divided rows

\[
\Lambda_3\Sigma_j(0)=0,
\qquad
\Lambda_3(\Sigma_j(1)-\operatorname{Id})=0,
\qquad j=1,2,3.
\]

## Deliberately unclaimed

- historical attaching-sphere movies were not recovered;
- no full-`q` series cocone is claimed;
- no conclusion should use `D_3(xi)` or the retired mixed-Z all-level model;
- the earlier end-to-end `HOLE` verdict was tied to the incorrect short G1
  bridge; after replacement by the AR product-ribbon proof, a fresh full-chain
  independent review is still required.

## Erratum — 2026-09-02 endpoint-coordinate unification

The historical base value above records the 2 September public draft and is
not the current value.  That draft paired the THXY-indexed vector
`e_0-e_5` with the collar-indexed covector `e_87^*-e_2^*`.  The cabled word
uses the collar convention.  The unique derived position table now maps the
same physical cup to collar positions `(2,87)`, so the current pair is

```text
u   = e_2 - e_87
ell = e_87^* - e_2^*
```

The exact replay gives

```text
ell (rho(W)-I) u, epsilon^0..6 = [0,0,0,-328,14596,-410246,9595271]
ell (rho(W)-I) u, h^0..6       = [0,0,0,2624,221728,11760112,520583560]
D_3(eta_R[T_1])                = 2624
D_3(xi)                        = 0
```

The nonvanishing status is unchanged, but the number remains conditional on
the frozen collar word and therefore on A1.  The position-table SHA-256 is
`119C7E9E74AE6C820DA72A84CDFD5D445D81E6C3AACCC209C25D37C323961508`;
the public receipt SHA-256 is
`765124A7625A5CA06BF4BCCD24B7FF24DEFD20DE9F12BE9F31F0AE0B6C5EF907`.
