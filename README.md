# Smooth4PC T73 — conditional skein-lasagna obstruction

[中文说明](README.zh-CN.md)

This repository supports a **conditional** skein-lasagna obstruction for the
trace-73 Cappell--Shaneson homotopy 4-sphere associated with
\[
A=\begin{pmatrix}0&269&1240\\0&41&189\\1&0&32\end{pmatrix}
\]
(Iwaki standard form \(X_{41,189,73}\)).

**No counterexample to the smooth four-dimensional Poincaré conjecture is
claimed.** The controlling manuscript is
[`paper/spc4-t73-candidate/main.tex`](paper/spc4-t73-candidate/main.tex)
(*A conditional skein-lasagna obstruction for a trace-73 Cappell--Shaneson
sphere*).

## What is proved, and what is open

For an explicit **Johnson-generator** handle presentation the paper formulates
the geometric inputs **P0, C, S, and P3** needed for a skein-lasagna comparison
at quantum degree \(494\), together with the identification
\(X_J\cong\Sigma_A^0\). **P0, C, and S are open**: the committed certificates
check finite combinatorial models of these inputs (control strands, local
model movies, model spheres), not the actual Cappell--Shaneson geometry. The
paper's status table (Section 3) is the claim boundary.

An exact finite calculation gives the nonzero integer \(2624\) as the divided
cubic of a frozen Burau computation in the collar endpoint convention;
identifying it with the MWW divided cubic of an actual geometric class is part
of the open input C. An Artin--Magnus certificate and the pure-braid
Andreadakis theorem establish a third-order property of the public braid word.

A Lean development formalizes the **abstract quotient argument**: given
interface data assembling the MWW quotients and four-handle transport
(`ExternalGeometry`), a nonzero degree-\(494\) class would obstruct
diffeomorphism with \(S^4\). Those geometric interfaces are **not** constructed
in Lean.

| Layer | Status |
| --- | --- |
| Finite algebra (\(2624\), degree \(494\), \(\det A=\det(A-I)=1\)) | Checked in Lean |
| Abstract conditional implication | Checked in Lean |
| Johnson P0 / C / S (actual geometry) | **Open** (finite-model certificates only) |
| P3 computed part (E12: empty-link degree 494 on \(S^4\)) | Checked; E11/E13 depend on P0 |
| Lean inhabitant of `ExternalGeometry` | **Open** |

The exact Lean boundary is
[`Smooth4PC/T73External.lean`](Smooth4PC/T73External.lean). Premises audit:

```text
python3 scripts/audit_t73_premises.py --check
```

Expected summary: `P0/C/S/P3=PASS`, `OVERALL=OPEN`, `COUNTEREXAMPLE=False`.
The `PASS` values are certificate-internal (the finite models replay); the
paper's status table records P0, C, and S as open.

**Erratum (2 September 2026).** An earlier draft mixed two endpoint index
tables and reported \(-59072\). With both objects in the collar table used by
the braid word, the exact value is \(+2624\) (still nonzero).

## PDFs for review

- English: [`output/pdf/spc4-t73-candidate.pdf`](output/pdf/spc4-t73-candidate.pdf)
- Chinese: [`output/pdf/spc4-t73-candidate-zh.pdf`](output/pdf/spc4-t73-candidate-zh.pdf)

Paper sources and build notes:
[`paper/spc4-t73-candidate/README.md`](paper/spc4-t73-candidate/README.md).

## Start here

1. Read the paper abstract and §3 (precise statements) in
   [`paper/spc4-t73-candidate/main.tex`](paper/spc4-t73-candidate/main.tex)
   or the English PDF above.
2. Follow [`REPRODUCING.md`](REPRODUCING.md) to build from a fresh clone, audit
   reported axioms, and recompute the detector.
3. Replay the finite detector and Johnson P0/C/S/P3 certificates below
   (or the fuller checklist in `REPRODUCING.md`).
4. For the independent-review boundary map, see
   [`docs/INDEPENDENT_REVIEW.md`](docs/INDEPENDENT_REVIEW.md).

## Replay calculation scripts

Requires **Python 3.10+** from the repository root. On Windows, `python` is
fine wherever `python3` appears. Each `--check` regenerates the certificate
in memory and compares it to the committed JSON under `audit/`.

### Finite detector (\(D_3=2624\))

```text
python -I -B scripts/recompute_t73_delta3.py --check
```

Expect `DELTA3_ETA_T1=2624`, `DELTA3_XI=0`, and `VERIFY=PASS`. The cup
vector `u` and cap row `ell` are no longer hand-written: they are derived from
the single endpoint authority `data/T73_ENDPOINT_CONVENTION.json` (physical
endpoint ids, orientations, geometric and public orders, pivotal
coefficients) by the endpoint transport program, which also proves
`W_public = P W_geometric P^-1` letter by letter along the actual cabled word:

```text
python -B scripts/build_t73_endpoint_transport.py --check
python -B scripts/verify_t73_endpoint_transport.py
```

Expect `ENDPOINT_TRANSPORT=PASS`, `NO_UNRESOLVED_SIGNS=PASS`, `DELTA3=2624`,
and the coordinate controls `-59072` / `-2496` reproduced as illegal partial
transports. The filtered cubic lemma behind this (`[h^3] ell A u = ell_0 A_3 u_0`
when `A` starts in order three, and simultaneous-conjugation invariance) is
formalized in `Smooth4PC/FilteredCubicNaturality.lean`
(`python -B tests/test_t73_filtered_cubic_naturality.py -v`). The statewise
Reynolds-average algebra of the two-handle cocone (beta invariance by placement
transport, pair-addition relations by fibre counting, extension to all cable
states below the selected one along iterated once-dotted additions,
independence of the common upper bound) is formalized for every state in
`Smooth4PC/ReynoldsCableCocone.lean`
(`python -B tests/test_t73_reynolds_cable_cocone.py -v`); its geometric inputs
are hypotheses there and belong to the open input C.

### Johnson P0 → C → S → P3 (recommended order)

```text
# P0a: dual-block regular neighbourhoods, explicit collapses, Regina cross-check
python -B scripts/verify_elementary_collapse.py --check
python -B scripts/verify_t73_handlebody_bridge_regina.py --check   # needs the regina module

# P0 (~1–2 min): AR bridge, cancellations, geometric braid
python -B scripts/certify_t73_p0_johnson.py --check

# Optional explicit P0 reconstruction input (not committed by default)
python -B scripts/build_t73_p0_reconstruction_input.py --write
python -B scripts/reconstruct_t73_p0.py audit/t73_p0_reconstruction_input.json

# C: product rectangles, comparison supports, assembled witness
python -B scripts/certify_t73_c1_cut_link.py --check
python -B scripts/certify_t73_c2_comparison.py --check
python -B scripts/generate_t73_c_comparison_witness.py --check

# S: reversed belt spheres and relative-move ledger
python -B scripts/certify_t73_s_standard_spheres.py --check
python -B scripts/certify_t73_s_relative_moves.py --check

# P3: four-handle picture, standard-S^4 degree 494, CS identification
python -B scripts/certify_t73_p3_four_handle.py --check
python -B scripts/certify_t73_e12_s4.py --check
python -B scripts/certify_t73_e13_close.py --check
python -B scripts/certify_t73_e13_identification.py --check

# Premise summary (must remain OVERALL=OPEN / no counterexample)
python -B scripts/audit_t73_premises.py --check
python -B scripts/check_t73_claim_boundary.py
```

| Script | Role | Expect |
| --- | --- | --- |
| `certify_t73_p0_johnson.py` | P0 Johnson replacement | `T73_P0_JOHNSON_CERTIFICATE=PASS` |
| `reconstruct_t73_p0.py` | Strict PL collar vs public word | `P0_RECONSTRUCTION=PASS`, `B44_LENGTH=11340` |
| `certify_t73_c1_cut_link.py` | 44 ribbons + 227 leftover \(z\) | `RECTANGLES=44`, `LEFTOVER_Z_CIRCLES=227` |
| `certify_t73_c2_comparison.py` | Disjoint C2 supports / \(H\) movies | `T73_C2_COMPARISON=PASS` |
| `generate_t73_c_comparison_witness.py` | C ledger bound to P0/C1/C2 | `C_STATUS=PASS` |
| `certify_t73_s_*.py` | S sphere system + moves | `PASS` |
| `certify_t73_p3_four_handle.py` | \(X_J\) four-handle layer | `E11`/`E12` PASS; `E13=PARTIAL` by design |
| `certify_t73_e12_s4.py` | Empty-link degree \(494\) on standard \(S^4\) | `S4_DEGREE_494_ZERO=True` |
| `certify_t73_e13_*.py` | \(X_J\cong\Sigma_A^0\) pipeline | `IDENTIFIED_WITH_SIGMA=True` |
| `audit_t73_premises.py` | Aggregate status | `P0/C/S/P3=PASS`, `COUNTEREXAMPLE=False` |
| `check_t73_claim_boundary.py` | Paper keeps conditional claim | `T73_CLAIM_BOUNDARY=OPEN_GEOMETRY` |

**Notes.**

- Certificates are SHA-chained: C binds to P0, S to P0+C, P3 to P0+C+S. Changing an upstream certificate without regenerating dependents will fail `--check`.
- `certify_t73_p3_four_handle.py` reporting `E13=PARTIAL` is expected; full \(\Sigma_A^0\) identification is in the `e13_*` scripts.
- C1/C2 check **combinatorial / PL model** obligations used by the paper; they do not re-prove MWW/BPW categorical comparison by themselves.
- Lean compile (`tests/test_t73_minimal_formalization.py`) is separate and slower (~5–10 min); see [`REPRODUCING.md`](REPRODUCING.md).

## Reproducibility contract

- Lean toolchain: `leanprover/lean4:v4.32.1`
- mathlib revision: `520045ab14e26149ee970e2e617ca04b09bde5d6`
- Python: 3.10 or later
- expected axiom reports: `38`
- allowed reported axioms: `propext`, `Classical.choice`, `Quot.sound`
- `sorryAx`: absent
- expected detector value: `2624`

The committed `lake-manifest.json` pins Lean dependencies. Build products and
local dependency copies are not part of the source contract.

## Scope

This is a public verification package for a **conditional** obstruction, not a
claim of peer acceptance and not a claimed counterexample. The most useful
adverse review attacks the Johnson geometric bindings and the remaining Lean
`ExternalGeometry` assembly, not the already checked integer arithmetic.

## Why this is being released on GitHub first

GitHub is the first public release channel for access, speed and
reproducibility—not a substitute for scholarly review. Existing arXiv history
is in computer science; mathematics-category endorsement may be unavailable.
The full argument, Lean sources, certificates and replay instructions are
inspectable here. After substantive scrutiny, a conventional preprint
submission is intended.

## License

The repository is released under the [MIT License](LICENSE).
