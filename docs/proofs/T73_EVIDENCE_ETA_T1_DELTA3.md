# Single-input decision: `delta_3(eta_R[T1])`

Date: 2026-08-31
Erratum: 2026-09-02

## Binary verdict after endpoint-index unification

```text
delta_3(eta_R[T1]) = 2624 != 0
delta_3(xi)         =    0

ENDPOINT-INDEX ERRATUM: NONZERO SURVIVES
```

The functional is defined once:

\[
\Delta_h(x)=\ell(\rho_h(W)-I)\Phi_h\operatorname{Sh}(x),
\qquad
\delta_3(x)=[h^3]\Delta_h(x).
\]

## 1. Erratum history

`ERRATUM_HISTORY_BEGIN` The 2 September public draft reported `-59072`.
That value mixed the THXY endpoint labels used for `u` with the collar/Burau
labels used for `ell`.  The Artin word is generated in the collar coordinates,
so the physical cup must first be translated by
`data/B88_POSITION_TO_PASSAGE_TABLE.json`.  `ERRATUM_HISTORY_END`

The table binds the physical cup

```text
r_xy negative passage c_r_xy_neg:0000  -> B88 index 2;
m_2  positive passage c_m_2_pos:0310   -> B88 index 87.
```

Thus, in the single coordinate convention used by the word and cap,

\[
u=e_2-e_{87},
\qquad
\ell=e_{87}^*-e_2^*.
\]

The source/target identification remains a candidate-specific geometric
premise.  This erratum fixes only the finite endpoint coordinates.

## 2. The two inputs remain distinct

For the fixed coefficient, let `eta_R[T]` be the MWW coefficient-trace class
of the Hattori identity at object `T`.  With

\[
B_\Omega=WF_\Omega,
\quad T_0=F_\Omega^{-1}U_1,
\quad T_1=F_\Omega^{-1}W^{-1}U_1,
\]

we have `B_Omega T0 = W U1` and `B_Omega T1 = U1`.  The corrected endpoint
binding is

\[
\Phi_h\operatorname{Sh}(\eta_R[T_1])=u+O(h).
\]

For

\[
\xi=\eta_R[T_0]-s_{\mathrm{inv}}\eta_R[T_1],
\]

the same convention gives

\[
\Phi_h\operatorname{Sh}(\xi)=(\rho_h(W)-I)u.
\]

Hence the detector applies one copy of `rho_h(W)-I` to `eta_R[T1]` and two
copies to `xi`.

## 3. Exact calculation

`scripts/recompute_t73_delta3.py` reconstructs the 45,360-letter B88 word
from primitive crossing rows and evaluates the unreduced Burau action over
`Z[epsilon]/(epsilon^7)`.  With the corrected `u` and the unchanged `ell` it
computes

```text
ell (rho(W)-I) u, epsilon degrees 0..6:
  [0, 0, 0, -328, 14596, -410246, 9595271]

ell (rho(W)-I)^2 u, epsilon degrees 0..6:
  [0, 0, 0, 0, 0, 0, -102729600]
```

Because `epsilon=-2h+3h^2-4h^3+...`,

\[
[h^3]\ell(\rho_h(W)-I)u=(-2)^3(-328)=2624,
\]

whereas the square starts in order six.  Therefore

\[
\boxed{\delta_3(\eta_R[T_1])=2624\ne0},
\qquad
\boxed{\delta_3(\xi)=0}.
\]

The complete `h` profiles through degree six are

```text
eta: [0, 0, 0, 2624, 221728, 11760112, 520583560]
xi:  [0, 0, 0, 0, 0, 0, -6574694400]
```

## 4. Reproduction

```powershell
python -I -B scripts/recompute_t73_delta3.py --check --write-receipt
```

The immutable outputs are:

```text
B88_POSITION_TO_PASSAGE_TABLE SHA256:
  119C7E9E74AE6C820DA72A84CDFD5D445D81E6C3AACCC209C25D37C323961508

T73_DELTA3_PUBLIC_RECEIPT.json SHA256:
  765124A7625A5CA06BF4BCCD24B7FF24DEFD20DE9F12BE9F31F0AE0B6C5EF907
```

This decision settles the endpoint-coordinate arithmetic only.  It does not
discharge the candidate-specific Hattori, two-handle, sphere-map, or closed
four-manifold assumptions.

## 5. Endpoint authority (2026-09-04)

The cup and cap terms are no longer hand-written. They are derived from the
single endpoint authority `data/T73_ENDPOINT_CONVENTION.json` by
`scripts/build_t73_endpoint_transport.py`, which records for every one of the
88 endpoints its physical identity, owner, orientation, geometric order (frozen
MWW cut slot), public order (collar position), THXY index and pivotal
coefficient, verifies `W_public = P W_geometric P^-1` letter by letter along
the actual 45360-letter cabled word, and reproduces the withdrawn values as
illegal partial transports (THXY cup with collar cap: -59072; THXY cup and cap
with the untransported word: -2496). The public input now records only the name
and SHA-256 of the authority file; the receipt schema is v2 and its SHA-256 is
the value quoted in section 4.
