# Reproducing the T73 audit

These commands are intended for a fresh clone with no pre-existing `.lake/`
directory or local `deps/` tree.

## 1. Prerequisites

- Git
- [elan](https://github.com/leanprover/elan), which installs the Lean version
  named by `lean-toolchain`
- Python 3.10 or later

On Windows, if Git fails with `SEC_E_NO_CREDENTIALS`, set the OpenSSL backend
for the current shell before fetching dependencies:

```powershell
$env:GIT_CONFIG_COUNT = '1'
$env:GIT_CONFIG_KEY_0 = 'http.sslBackend'
$env:GIT_CONFIG_VALUE_0 = 'openssl'
```

## 2. Clone and materialize the pinned dependencies

```text
git clone <repository-url> smooth4pc-t73-lean
cd smooth4pc-t73-lean
lake update
git diff --exit-code -- lake-manifest.json
lake exe cache get
```

`git diff --exit-code` must return exit code `0`; otherwise the dependency
lockfile does not match the published source.

## 3. Compile the complete audited chain

```text
python -B tests/test_t73_minimal_formalization.py -v
lake lean T73Audit.lean
```

Expected results:

- the Python suite reports `Ran 2 tests` and `OK`;
- `T73Audit.lean` exits `0`;
- exactly 38 `#print axioms` reports appear;
- every report is a subset of
  `propext`, `Classical.choice`, `Quot.sound`;
- the output contains no `sorryAx`.

The Python test enforces those conditions, scans the audited Lean sources for
forbidden proof escapes and builds every project module into a fresh temporary
olean root.

## 4. Recompute the detector independently

```text
python -I -B scripts/recompute_t73_delta3.py --check
```

The script reconstructs the registered point-push word from primitive crossing
rows, checks the SHA-pinned B88 position-to-passage table, rebuilds the
two-cable action and performs exact truncated-polynomial arithmetic. Its input
does not contain the expected cubic as a result field.
The bearing output is:

```text
POSITION_TABLE_SHA256=119C7E9E74AE6C820DA72A84CDFD5D445D81E6C3AACCC209C25D37C323961508
ELL_RHOW_MINUS_I_U_EPS=[0,0,0,-328,14596,-410246,9595271]
ELL_RHOW_MINUS_I_SQUARED_U_EPS=[0,0,0,0,0,0,-102729600]
DELTA3_ETA_T1=2624
DELTA3_XI=0
VERIFY=PASS
```

The cup and cap terms are derived from the endpoint authority
`data/T73_ENDPOINT_CONVENTION.json`; regenerate and verify the transport with

```text
python -B scripts/build_t73_endpoint_transport.py --check
python -B scripts/verify_t73_endpoint_transport.py
python -B tests/test_t73_endpoint_transport.py -v
```

The bearing lines are `ENDPOINT_TRANSPORT=PASS`, `NO_UNRESOLVED_SIGNS=PASS`
and `DELTA3=2624`.

## 5. Replay Johnson P0 / C / S / P3 certificates

Run from the repository root with Python 3.10+. On Windows, `python` may be
used in place of `python3`. Each `--check` regenerates the certificate and
compares it to the committed file under `audit/`. Prefer this order (SHA
chain: P0 → C → S → P3):

```text
python -B scripts/certify_t73_p0_johnson.py --check

# Optional: materialize and verify the strict P0 reconstruction input
python -B scripts/build_t73_p0_reconstruction_input.py --write
python -B scripts/reconstruct_t73_p0.py audit/t73_p0_reconstruction_input.json

python -B scripts/certify_t73_c1_cut_link.py --check
python -B scripts/certify_t73_c2_comparison.py --check
python -B scripts/generate_t73_c_comparison_witness.py --check

python -B scripts/certify_t73_s_standard_spheres.py --check
python -B scripts/certify_t73_s_relative_moves.py --check

python -B scripts/certify_t73_p3_four_handle.py --check
python -B scripts/certify_t73_e12_s4.py --check
python -B scripts/certify_t73_e13_close.py --check
python -B scripts/certify_t73_e13_identification.py --check

python -B scripts/audit_t73_premises.py --check
python -B scripts/check_t73_claim_boundary.py
```

Expected highlights:

- P0: `T73_P0_JOHNSON_CERTIFICATE=PASS` (full geometric braid; ~1–2 minutes)
- optional reconstruct: `P0_RECONSTRUCTION=PASS`, `B44_LENGTH=11340`
- C1: `RECTANGLES=44`, `LEFTOVER_Z_CIRCLES=227`
- C2 / C witness: `PASS`
- S: `PASS` on spheres and relative moves
- P3 four-handle: `E11`/`E12` PASS and `E13=PARTIAL` **by design**
  (identification is completed by the `e13_*` scripts)
- E12: `S4_DEGREE_494_ZERO=True`, `ABOUT_STANDARD_S4=True`
- E13: `IDENTIFIED_WITH_SIGMA=True`
- audit: `P0/C/S/P3=PASS` (certificate-internal finite models), `PAPER_STATUS=P0:OPEN,C:OPEN,S:OPEN,P3_E11:OPEN,P3_E12:PROVED,P3_E13:OPEN` (paper claim boundary), `OVERALL=OPEN`, `COUNTEREXAMPLE=False`
- claim boundary: `T73_CLAIM_BOUNDARY=OPEN_GEOMETRY`

A short copy of this checklist also appears in the root `README.md` /
`README.zh-CN.md`.

## 6. Verify the public geometry evidence

Before that, verify the public geometry evidence and recompute the
2,126,291-crossing global-descending certificate:

```text
python -I -B scripts/verify_public_geometry_evidence.py
```

The command must end with `GLOBAL_DESCENDING=PASS` and `VERIFY=PASS`.

## 7. Check the convention freeze

The legality criteria were committed before the detailed convention search:

```text
git merge-base --is-ancestor cf4a990 9d75dcd
```

The command must exit `0`. The two records are:

- `docs/proofs/QSTAR_R7_LEGALITY_CRITERIA_PRESEARCH_20260901.md`
- `docs/proofs/QSTAR_R7_LEGALITY_ADJUDICATION_20260901.md`

## 8. Interpret the result correctly

Compilation verifies the finite algebra and the implication from the
`ExternalGeometry` and Cappell--Shaneson interfaces to the final conclusion.
It does not turn those interfaces into kernel-checked differential topology.
Their published sources, exact scope and candidate-specific evidence are
listed in `docs/INDEPENDENT_REVIEW.md`.

For the erratum compilation and complete axiom output, see
`docs/proofs/QSTAR_R8_LEAN_COMPILE_RECEIPT_20260903.md` and
`docs/proofs/QSTAR_R8_T73AUDIT_RAW_20260903.txt`.  The earlier fresh-clone
record remains at `docs/proofs/PUBLIC_RELEASE_CLEAN_REPLAY_20260901.md`.

## 9. Build the review paper

From `paper/spc4-t73-candidate` in WSL or Linux:

```text
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
mkdir -p ../../output/pdf
cp main.pdf ../../output/pdf/spc4-t73-candidate.pdf
```

See `paper/spc4-t73-candidate/README.md` for the Chinese edition build.
