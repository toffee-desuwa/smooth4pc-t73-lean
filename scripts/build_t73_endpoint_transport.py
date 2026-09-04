#!/usr/bin/env python3
"""Single endpoint authority and monomial endpoint transport for the T73 detector.

The divided cubic is a matrix coefficient ell (rho(W) - I) u on the weight-86
subspace E_88 of an 88-fold tensor product.  Three bookkeeping systems for the
88 y-gate endpoints exist in the repository:

* the geometric order: the west boundary of the frozen MWW cut object, top to
  bottom (evidence/public_geometry/CUT_OBJECT.json), carrying the physical
  endpoint identities, cable signs and orientations;
* the public order: the collar wicket numbering used by the cabled Artin
  word and by the public Burau computation
  (data/B88_POSITION_TO_PASSAGE_TABLE.json);
* the historical THXY order (owner blocks with reversed wicket order), which
  produced the withdrawn value -59072 when it was mixed with the collar
  order.

This program builds ONE authority file, data/T73_ENDPOINT_CONVENTION.json,
recording for every endpoint its physical identity, owner, orientation,
geometric order, public order, THXY index, pivotal coefficient and
weight-defect basis vector, and derives from it the monomial transport
P(q) (a permutation composed with monomials +-q^k) with

    W_public = P W_geometric P^-1,   u_public = P u_geometric,
    ell_public = ell_geometric P^-1.

The oriented (geometric) model is the U_q(sl2) fundamental representation V
and its dual V* with the universal R-matrix, the pivotal isomorphism
phi : V* -> V, and the canonical / pivotal (co)evaluations.  Every convention
of that model is verified by the program (intertwining, braid relation,
inverse, module map, zigzag) before it is used.  The public model is the
unreduced Burau representation with t = q^-2.  The identity
W_public = P W_geometric P^-1 is proved letter by letter along the actual
45360-letter cabled word with strand tracking, exactly in Z[q, q^-1], and
cross-checked numerically on transported vectors.

The constant terms u_0 = e_2 - e_87 and ell_0 = e_87^* - e_2^* are DERIVED
here; nothing about them is hard-coded, and recompute_t73_delta3.py reads
them from this derivation instead of from hand-written lists.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
CUT_OBJECT = ROOT / "evidence" / "public_geometry" / "CUT_OBJECT.json"
POSITION_TABLE = ROOT / "data" / "B88_POSITION_TO_PASSAGE_TABLE.json"
PUBLIC_INPUT = ROOT / "data" / "T73_DELTA3_PUBLIC_INPUT.json"
CONVENTION = ROOT / "data" / "T73_ENDPOINT_CONVENTION.json"
AUDIT = ROOT / "audit" / "t73_endpoint_transport.json"
RECOMPUTE = ROOT / "scripts" / "recompute_t73_delta3.py"

DIM = 88
DEGREE = 6
SCHEMA_CONVENTION = "t73_endpoint_convention/v1"
SCHEMA_AUDIT = "t73_endpoint_transport/v1"

# The selected one-cup tangle joins these two physical passages (west copies).
SELECTED_CUP_PASSAGES = ("pass:c_r_xy_neg:0000", "pass:c_m_2_pos:0310")


# ---------------------------------------------------------------------------
# hashing helpers
# ---------------------------------------------------------------------------


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest().upper()


def canonical_sha(value: Any) -> str:
    raw = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256_bytes(raw)


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n"


# ---------------------------------------------------------------------------
# exact Laurent polynomials in q
# ---------------------------------------------------------------------------


class LP:
    """Laurent polynomial in q with integer coefficients (dict exponent -> coeff)."""

    __slots__ = ("c",)

    def __init__(self, c: dict[int, int] | None = None) -> None:
        self.c = {k: v for k, v in (c or {}).items() if v != 0}

    @staticmethod
    def mono(coeff: int, exp: int) -> "LP":
        return LP({exp: coeff})

    @staticmethod
    def q(exp: int = 1) -> "LP":
        return LP({exp: 1})

    @staticmethod
    def one() -> "LP":
        return LP({0: 1})

    @staticmethod
    def zero() -> "LP":
        return LP()

    def __add__(self, other: "LP") -> "LP":
        out = dict(self.c)
        for k, v in other.c.items():
            out[k] = out.get(k, 0) + v
        return LP(out)

    def __sub__(self, other: "LP") -> "LP":
        out = dict(self.c)
        for k, v in other.c.items():
            out[k] = out.get(k, 0) - v
        return LP(out)

    def __neg__(self) -> "LP":
        return LP({k: -v for k, v in self.c.items()})

    def __mul__(self, other: "LP") -> "LP":
        out: dict[int, int] = {}
        for k1, v1 in self.c.items():
            for k2, v2 in other.c.items():
                out[k1 + k2] = out.get(k1 + k2, 0) + v1 * v2
        return LP(out)

    def scale(self, n: int) -> "LP":
        return LP({k: n * v for k, v in self.c.items()})

    def __eq__(self, other: object) -> bool:
        return isinstance(other, LP) and self.c == other.c

    def __hash__(self) -> int:
        return hash(tuple(sorted(self.c.items())))

    def is_zero(self) -> bool:
        return not self.c

    def is_monomial(self) -> bool:
        return len(self.c) == 1 and abs(next(iter(self.c.values()))) == 1

    def monomial_data(self) -> tuple[int, int]:
        if not self.is_monomial():
            raise ValueError(f"not a unit monomial: {self}")
        (exp, coeff), = self.c.items()
        return coeff, exp

    def inverse_monomial(self) -> "LP":
        coeff, exp = self.monomial_data()
        return LP({-exp: coeff})

    def __repr__(self) -> str:
        if not self.c:
            return "0"
        parts = []
        for exp in sorted(self.c):
            coeff = self.c[exp]
            if exp == 0:
                parts.append(f"{coeff}")
            elif coeff == 1:
                parts.append(f"q^{exp}")
            elif coeff == -1:
                parts.append(f"-q^{exp}")
            else:
                parts.append(f"{coeff}*q^{exp}")
        return " + ".join(parts)

    def to_json(self) -> dict[str, int]:
        return {str(k): v for k, v in sorted(self.c.items())}

    def h_series(self, degree: int = DEGREE) -> list[int]:
        """Expand at q = 1 + h, truncated modulo h^(degree+1) (integer coefficients)."""
        out = [0] * (degree + 1)
        for exp, coeff in self.c.items():
            for n in range(degree + 1):
                out[n] += coeff * generalized_binomial(exp, n)
        return out


def generalized_binomial(k: int, n: int) -> int:
    """binom(k, n) for integer k (possibly negative), integer n >= 0."""
    if n == 0:
        return 1
    if k >= 0:
        return math.comb(k, n) if n <= k else 0
    # binom(-m, n) = (-1)^n binom(m + n - 1, n)
    m = -k
    return (-1) ** n * math.comb(m + n - 1, n)


Matrix = list[list[LP]]


def mat_identity(n: int) -> Matrix:
    return [[LP.one() if i == j else LP.zero() for j in range(n)] for i in range(n)]


def mat_mul(a: Matrix, b: Matrix) -> Matrix:
    n, m, p = len(a), len(b), len(b[0])
    out = [[LP.zero() for _ in range(p)] for _ in range(n)]
    for i in range(n):
        for k in range(m):
            aik = a[i][k]
            if aik.is_zero():
                continue
            for j in range(p):
                if not b[k][j].is_zero():
                    out[i][j] = out[i][j] + aik * b[k][j]
    return out


def mat_eq(a: Matrix, b: Matrix) -> bool:
    return len(a) == len(b) and all(len(ra) == len(rb) and all(x == y for x, y in zip(ra, rb)) for ra, rb in zip(a, b))


def mat_json(a: Matrix) -> list[list[dict[str, int]]]:
    return [[x.to_json() for x in row] for row in a]


# ---------------------------------------------------------------------------
# the oriented model: U_q(sl2) fundamental representation V and its dual V*
# ---------------------------------------------------------------------------

V = "V"
VD = "V*"
PLUS, MINUS = "+", "-"
LABELS = (PLUS, MINUS)
Q = LP.q(1)
QINV = LP.q(-1)


def weight(space: str, label: str) -> int:
    w = 1 if label == PLUS else -1
    return w if space == V else -w


def ground(space: str) -> str:
    """Highest-weight (weight +1) vector: v_+ in V, v_-^* in V*."""
    return PLUS if space == V else MINUS


def defect(space: str) -> str:
    """Weight -1 vector: v_- in V, v_+^* in V*."""
    return MINUS if space == V else PLUS


Vec = dict[str, LP]  # label -> coefficient


def act_E(space: str, label: str) -> Vec:
    if space == V:
        return {PLUS: LP.one()} if label == MINUS else {}
    # V*: E v_+^* = -q v_-^*, E v_-^* = 0   (dual action through the antipode)
    return {MINUS: LP.mono(-1, 1)} if label == PLUS else {}


def act_F(space: str, label: str) -> Vec:
    if space == V:
        return {MINUS: LP.one()} if label == PLUS else {}
    # V*: F v_-^* = -q^-1 v_+^*, F v_+^* = 0
    return {PLUS: LP.mono(-1, -1)} if label == MINUS else {}


def act_K(space: str, label: str, power: int = 1) -> Vec:
    return {label: LP.q(power * weight(space, label))}


Pair = tuple[str, str]
PVec = dict[Pair, LP]


def pair_basis() -> list[Pair]:
    return [(a, b) for a in LABELS for b in LABELS]


def add_into(target: PVec, key: Pair, value: LP) -> None:
    if value.is_zero():
        return
    target[key] = target.get(key, LP.zero()) + value
    if target[key].is_zero():
        del target[key]


def universal_R(X: str, Y: str, a: str, b: str, inverse: bool = False) -> PVec:
    """R (or R^-1) applied to a (x) b in X (x) Y, normalised by q^{-1/2} (q^{1/2}).

    R = q^{-1/2} q^{H(x)H/2} (1 + (q - q^-1) E (x) F),
    R^-1 = q^{1/2} (1 - (q - q^-1) E (x) F) q^{-H(x)H/2}.
    Every exponent is an integer because the weights are odd.
    """
    out: PVec = {}
    wa, wb = weight(X, a), weight(Y, b)
    qq = Q - QINV
    if not inverse:
        # exponent of the identity term: (wa*wb - 1)/2
        add_into(out, (a, b), LP.q((wa * wb - 1) // 2))
        for a2, ca in act_E(X, a).items():
            for b2, cb in act_F(Y, b).items():
                w2 = (weight(X, a2) * weight(Y, b2) - 1) // 2
                add_into(out, (a2, b2), ca * cb * qq * LP.q(w2))
    else:
        base = LP.q((1 - wa * wb) // 2)
        add_into(out, (a, b), base)
        for a2, ca in act_E(X, a).items():
            for b2, cb in act_F(Y, b).items():
                add_into(out, (a2, b2), -(ca * cb * qq * base))
    return out


def braiding_matrix(X: str, Y: str, inverse: bool = False) -> Matrix:
    """Matrix of c_{X,Y} = tau o R : X(x)Y -> Y(x)X (inverse=False) or of
    c_{Y,X}^{-1} = R^{-1} o tau^{-1} : X(x)Y -> Y(x)X (inverse=True).

    Rows are indexed by the target basis of Y(x)X, columns by the source basis
    of X(x)Y, both in pair_basis() order.
    """
    src = pair_basis()
    tgt = pair_basis()
    index = {p: i for i, p in enumerate(tgt)}
    mat = [[LP.zero() for _ in src] for _ in tgt]
    for j, (a, b) in enumerate(src):
        if not inverse:
            image = universal_R(X, Y, a, b)
            for (a2, b2), coeff in image.items():
                mat[index[(b2, a2)]][j] = mat[index[(b2, a2)]][j] + coeff  # flip
        else:
            # c_{Y,X}^{-1} = R_{Y,X}^{-1} o tau^{-1}; tau^{-1}(a(x)b) = b(x)a in Y(x)X
            image = universal_R(Y, X, b, a, inverse=True)
            for (b2, a2), coeff in image.items():
                mat[index[(b2, a2)]][j] = mat[index[(b2, a2)]][j] + coeff
    return mat


def coproduct_matrix(gen: str, X: str, Y: str) -> Matrix:
    """Delta(gen) on X(x)Y with Delta(E)=E(x)K+1(x)E, Delta(F)=F(x)1+K^-1(x)F, Delta(K)=K(x)K."""
    src = pair_basis()
    index = {p: i for i, p in enumerate(src)}
    mat = [[LP.zero() for _ in src] for _ in src]
    for j, (a, b) in enumerate(src):
        image: PVec = {}
        if gen == "E":
            for a2, ca in act_E(X, a).items():
                for b2, cb in act_K(Y, b).items():
                    add_into(image, (a2, b2), ca * cb)
            for b2, cb in act_E(Y, b).items():
                add_into(image, (a, b2), cb)
        elif gen == "F":
            for a2, ca in act_F(X, a).items():
                add_into(image, (a2, b), ca)
            for a2, ca in act_K(X, a, -1).items():
                for b2, cb in act_F(Y, b).items():
                    add_into(image, (a2, b2), ca * cb)
        elif gen == "K":
            for a2, ca in act_K(X, a).items():
                for b2, cb in act_K(Y, b).items():
                    add_into(image, (a2, b2), ca * cb)
        else:
            raise ValueError(gen)
        for key, coeff in image.items():
            mat[index[key]][j] = mat[index[key]][j] + coeff
    return mat


def verify_intertwining() -> dict[str, bool]:
    result = {}
    for X in (V, VD):
        for Y in (V, VD):
            c = braiding_matrix(X, Y)
            ok = True
            for gen in ("E", "F", "K"):
                left = mat_mul(c, coproduct_matrix(gen, X, Y))
                right = mat_mul(coproduct_matrix(gen, Y, X), c)
                ok = ok and mat_eq(left, right)
            result[f"c[{X},{Y}]"] = ok
    return result


def verify_inverses() -> dict[str, bool]:
    result = {}
    for X in (V, VD):
        for Y in (V, VD):
            c = braiding_matrix(X, Y)             # X(x)Y -> Y(x)X
            cinv = braiding_matrix(Y, X, True)    # Y(x)X -> X(x)Y  (c_{X,Y}^{-1})
            result[f"c^-1 c [{X},{Y}]"] = mat_eq(mat_mul(cinv, c), mat_identity(4))
            result[f"c c^-1 [{X},{Y}]"] = mat_eq(mat_mul(c, cinv), mat_identity(4))
    return result


def kron(a: Matrix, b: Matrix) -> Matrix:
    """Kronecker product of rectangular matrices."""
    ra, ca = len(a), len(a[0])
    rb, cb = len(b), len(b[0])
    out = [[LP.zero() for _ in range(ca * cb)] for _ in range(ra * rb)]
    for i in range(ra):
        for j in range(ca):
            if a[i][j].is_zero():
                continue
            for k in range(rb):
                for l in range(cb):
                    if not b[k][l].is_zero():
                        out[i * rb + k][j * cb + l] = a[i][j] * b[k][l]
    return out


def verify_braid_relation() -> dict[str, bool]:
    result = {}
    id2 = mat_identity(2)
    for X in (V, VD):
        for Y in (V, VD):
            for Z in (V, VD):
                # c12 then c23 then c12 on X Y Z -> Z Y X, versus c23 c12 c23.
                # Careful bookkeeping of the space labels:
                # path A: (X,Y,Z) -c12-> (Y,X,Z) -c23-> (Y,Z,X) -c12-> (Z,Y,X)
                a1 = kron(braiding_matrix(X, Y), id2)
                a2 = kron(id2, braiding_matrix(X, Z))
                a3 = kron(braiding_matrix(Y, Z), id2)
                path_a = mat_mul(a3, mat_mul(a2, a1))
                # path B: (X,Y,Z) -c23-> (X,Z,Y) -c12-> (Z,X,Y) -c23-> (Z,Y,X)
                b1 = kron(id2, braiding_matrix(Y, Z))
                b2 = kron(braiding_matrix(X, Z), id2)
                b3 = kron(id2, braiding_matrix(X, Y))
                path_b = mat_mul(b3, mat_mul(b2, b1))
                result[f"YBE[{X},{Y},{Z}]"] = mat_eq(path_a, path_b)
    return result


def phi_matrix() -> Matrix:
    """phi : V* -> V, phi(v_-^*) = v_+, phi(v_+^*) = -q v_-  (rows: V basis, cols: V* basis)."""
    mat = [[LP.zero(), LP.zero()], [LP.zero(), LP.zero()]]
    # columns in LABELS order: col 0 = v_+^*, col 1 = v_-^*
    mat[1][0] = LP.mono(-1, 1)   # v_+^* -> -q v_-
    mat[0][1] = LP.one()         # v_-^* -> v_+
    return mat


def single_action(gen: str, space: str) -> Matrix:
    mat = [[LP.zero(), LP.zero()], [LP.zero(), LP.zero()]]
    for j, label in enumerate(LABELS):
        act = {"E": act_E, "F": act_F}.get(gen)
        image = act(space, label) if act else act_K(space, label)
        for label2, coeff in image.items():
            mat[LABELS.index(label2)][j] = coeff
    return mat


def verify_phi() -> dict[str, bool]:
    phi = phi_matrix()
    result = {}
    for gen in ("E", "F", "K"):
        result[f"phi intertwines {gen}"] = mat_eq(
            mat_mul(phi, single_action(gen, VD)), mat_mul(single_action(gen, V), phi)
        )
    id2 = mat_identity(2)
    # naturality of the braiding with respect to phi
    result["c[V,V](id x phi) = (phi x id)c[V,V*]"] = mat_eq(
        mat_mul(braiding_matrix(V, V), kron(id2, phi)), mat_mul(kron(phi, id2), braiding_matrix(V, VD))
    )
    result["c[V,V](phi x id) = (id x phi)c[V*,V]"] = mat_eq(
        mat_mul(braiding_matrix(V, V), kron(phi, id2)), mat_mul(kron(id2, phi), braiding_matrix(VD, V))
    )
    result["c[V,V](phi x phi) = (phi x phi)c[V*,V*]"] = mat_eq(
        mat_mul(braiding_matrix(V, V), kron(phi, phi)), mat_mul(kron(phi, phi), braiding_matrix(VD, VD))
    )
    return result


# duality maps ---------------------------------------------------------------


def coev_left() -> PVec:
    """coev_V : 1 -> V (x) V*, 1 -> sum_i v_i (x) v_i^*  (canonical)."""
    return {(PLUS, PLUS): LP.one(), (MINUS, MINUS): LP.one()}


def ev_left(a: str, b: str) -> LP:
    """ev_V : V* (x) V -> 1, v_i^* (x) v_j -> delta_ij (canonical)."""
    return LP.one() if a == b else LP.zero()


def ev_right(a: str, b: str, sigma: int, twist: int) -> LP:
    """ev'_V : V (x) V* -> 1, v_i (x) v_j^* -> sigma delta_ij q^{twist * weight(v_i)}."""
    return LP.mono(sigma, twist * weight(V, a)) if a == b else LP.zero()


def coev_right(sigma: int, twist: int) -> PVec:
    """coev'_V : 1 -> V* (x) V, 1 -> sigma sum_i q^{-twist*weight(v_i)} v_i^* (x) v_i."""
    return {(label, label): LP.mono(sigma, -twist * weight(V, label)) for label in LABELS}


def row_matrix_ev(ev: Callable[[str, str], LP]) -> Matrix:
    return [[ev(a, b) for (a, b) in pair_basis()]]


def col_matrix_coev(vec: PVec) -> Matrix:
    return [[vec.get(p, LP.zero())] for p in pair_basis()]


def counit(gen: str) -> LP:
    return LP.one() if gen == "K" else LP.zero()


def verify_duality(sigma: int, twist: int) -> dict[str, bool]:
    result = {}
    id2 = mat_identity(2)
    evL = row_matrix_ev(ev_left)                       # 1 x 4 on V* (x) V
    evR = row_matrix_ev(lambda a, b: ev_right(a, b, sigma, twist))  # 1 x 4 on V (x) V*
    coL = col_matrix_coev(coev_left())                 # 4 x 1 into V (x) V*
    coR = col_matrix_coev(coev_right(sigma, twist))    # 4 x 1 into V* (x) V
    # module-map property: ev o Delta(g) = counit(g) ev ; Delta(g) o coev = counit(g) coev
    for gen in ("E", "F", "K"):
        eps = counit(gen)
        result[f"ev_left module map {gen}"] = mat_eq(mat_mul(evL, coproduct_matrix(gen, VD, V)), [[x * eps for x in evL[0]]])
        result[f"ev_right module map {gen}"] = mat_eq(mat_mul(evR, coproduct_matrix(gen, V, VD)), [[x * eps for x in evR[0]]])
        result[f"coev_left module map {gen}"] = mat_eq(mat_mul(coproduct_matrix(gen, V, VD), coL), [[row[0] * eps] for row in coL])
        result[f"coev_right module map {gen}"] = mat_eq(mat_mul(coproduct_matrix(gen, VD, V), coR), [[row[0] * eps] for row in coR])
    # zigzag identities
    # (ev_L x id_{V*})(id_{V*} x coev_L) = id_{V*}
    z1 = mat_mul(kron(evL, id2), kron(id2, coL))
    result["zigzag ev_left/coev_left on V*"] = mat_eq(z1, id2)
    # (id_V x ev_L)(coev_L x id_V) = id_V
    z2 = mat_mul(kron(id2, evL), kron(coL, id2))
    result["zigzag coev_left/ev_left on V"] = mat_eq(z2, id2)
    # (ev_R x id_V)(id_V x coev_R) = id_V
    z3 = mat_mul(kron(evR, id2), kron(id2, coR))
    result["zigzag ev_right/coev_right on V"] = mat_eq(z3, id2)
    # (id_{V*} x ev_R)(coev_R x id_{V*}) = id_{V*}
    z4 = mat_mul(kron(id2, evR), kron(coR, id2))
    result["zigzag coev_right/ev_right on V*"] = mat_eq(z4, id2)
    # loop value ev_R o coev_L
    loop = mat_mul(evR, coL)[0][0]
    result["loop value"] = loop  # type: ignore[assignment]
    return result


def derive_pivotal_twist(sigma: int) -> int:
    """The K-power in ev_right is derived: exactly one of {+1, -1} makes ev_right a module map."""
    good = []
    for twist in (1, -1):
        checks = verify_duality(sigma, twist)
        if all(v is True for k, v in checks.items() if k != "loop value"):
            good.append(twist)
    if len(good) != 1:
        raise AssertionError(f"pivotal twist is not determined: {good}")
    return good[0]


# one-defect blocks ------------------------------------------------------------


def one_defect_block(X: str, Y: str, positive: bool) -> tuple[Matrix, dict[str, Any]]:
    """2x2 block of the crossing on the one-defect subspace of X(x)Y (X left, Y right).

    Basis before: e_left = defect(X) (x) ground(Y), e_right = ground(X) (x) defect(Y).
    Basis after (spaces Y (x) X): e_left' = defect(Y) (x) ground(X), e_right' = ground(Y) (x) defect(X).
    Returns the block (columns = images) and the ground-state coefficient.
    """
    mat = braiding_matrix(X, Y, inverse=not positive)
    src = pair_basis()
    tgt = pair_basis()
    before = [(defect(X), ground(Y)), (ground(X), defect(Y))]
    after = [(defect(Y), ground(X)), (ground(Y), defect(X))]
    block = [[LP.zero(), LP.zero()], [LP.zero(), LP.zero()]]
    leakage = False
    for j, col_key in enumerate(before):
        col = src.index(col_key)
        for r, row_key in enumerate(after):
            block[r][j] = mat[tgt.index(row_key)][col]
        for i, key in enumerate(tgt):
            if key not in after and not mat[i][col].is_zero():
                leakage = True
    ground_col = src.index((ground(X), ground(Y)))
    ground_row = tgt.index((ground(Y), ground(X)))
    ground_coeff = mat[ground_row][ground_col]
    for i in range(4):
        if i != ground_row and not mat[i][ground_col].is_zero():
            leakage = True
    return block, {"ground_coefficient": repr(ground_coeff), "ground_is_one": ground_coeff == LP.one(), "no_leakage": not leakage}


def burau_block(positive: bool) -> Matrix:
    """Public unreduced Burau block with t = q^-2 (columns = images)."""
    t = LP.q(-2)
    tinv = LP.q(2)
    if positive:
        return [[LP.one() - t, t], [LP.one(), LP.zero()]]
    return [[LP.zero(), LP.one()], [tinv, LP.one() - tinv]]


def diag2(a: LP, b: LP) -> Matrix:
    return [[a, LP.zero()], [LP.zero(), b]]


def derive_position_ratio() -> LP:
    """Find the monomial r = d_{p+1}/d_p that conjugates the oriented V(x)V blocks to the Burau blocks.

    In coordinates x_pub = M x_or with M = diag(d_p, d_{p+1}), the identity is
    B_pub M = M B_or (positions are fixed for V(x)V because both strands are V).
    """
    candidates = []
    for sign in (1, -1):
        for k in range(-4, 5):
            r = LP.mono(sign, k)
            m = diag2(LP.one(), r)
            ok = True
            for positive in (True, False):
                b_or, _ = one_defect_block(V, V, positive)
                if not mat_eq(mat_mul(burau_block(positive), m), mat_mul(m, b_or)):
                    ok = False
            if ok:
                candidates.append(r)
    if len(candidates) != 1:
        raise AssertionError(f"position ratio is not determined: {candidates}")
    return candidates[0]


# ---------------------------------------------------------------------------
# endpoint records
# ---------------------------------------------------------------------------


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_recompute():
    spec = importlib.util.spec_from_file_location("recompute_t73_delta3", RECOMPUTE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load recompute_t73_delta3.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_endpoint_records(cut: dict[str, Any], table: dict[str, Any]) -> list[dict[str, Any]]:
    gate = [g for g in cut["gate_passages"] if g["gate"] == "y"]
    if len(gate) != 1:
        raise ValueError("CUT_OBJECT does not have exactly one y gate")
    gate = gate[0]
    west = gate["boundary_endpoint_order_top_to_bottom"]["west"]
    if len(west) != DIM:
        raise ValueError("CUT_OBJECT west order does not have 88 endpoints")
    passages = {p["passage_id"]: p for p in gate["passages_top_to_bottom"]}
    public = {row["passage_id"]: row for row in table["positions"]}
    if sorted(row["index"] for row in table["positions"]) != list(range(DIM)):
        raise ValueError("position table indices are not 0..87")
    records: list[dict[str, Any]] = []
    for slot, endpoint_id in enumerate(west):
        head, kind = endpoint_id.rsplit(":", 1)
        if not head.startswith("ep:") or kind not in ("entry", "exit"):
            raise ValueError(f"unexpected endpoint id {endpoint_id}")
        passage_id = "pass:" + head[3:]
        passage = passages[passage_id]
        row = public[passage_id]
        if passage["entry_side"] != "west" and passage["exit_side"] != "west":
            raise ValueError(f"passage {passage_id} does not meet the west boundary")
        expected_kind = "entry" if passage["entry_side"] == "west" else "exit"
        if kind != expected_kind:
            raise ValueError(f"endpoint {endpoint_id} kind disagrees with passage sides")
        if row["owner"] != passage["owner_id"]:
            raise ValueError(f"owner mismatch for {passage_id}")
        if (row["sign"] == "neg") != (passage["cable_sign"] == "negative"):
            raise ValueError(f"cable sign mismatch for {passage_id}")
        if row["word_letter"] != passage["base_letter_index"]:
            raise ValueError(f"word letter mismatch for {passage_id}")
        # orientation at the west boundary: 'entry' = strand oriented into the cut tangle (tensor
        # factor V, weight +1 ground state v_+); 'exit' = strand oriented out of it (factor V*).
        space = V if kind == "entry" else VD
        records.append(
            {
                "public_index": row["index"],
                "geometric_order": slot,
                "physical_endpoint_id": endpoint_id,
                "passage_id": passage_id,
                "physical_component_id": passage["physical_component_id"],
                "owner": passage["owner_id"],
                "cable_sign": passage["cable_sign"],
                "orientation_relative_to_owner": passage["orientation_relative_to_owner"],
                "wicket": row["wicket"],
                "word_letter": row["word_letter"],
                "traversal_index": passage["traversal_index"],
                "passage_exponent": passage["exponent"],
                "travel_direction": passage["travel_direction"],
                "orientation": kind,
                "orientation_sign": 1 if kind == "entry" else -1,
                "tensor_factor": space,
                "ground_state": "v_+" if space == V else "v_-^*",
                "defect_state": "v_-" if space == V else "v_+^*",
            }
        )
    if sorted(r["public_index"] for r in records) != list(range(DIM)):
        raise ValueError("records do not cover public indices 0..87")
    # THXY reconstruction: owner blocks (r_xy, m_2), wickets in decreasing word-letter order,
    # copies (negative, positive) adjacent.
    order = sorted(
        records,
        key=lambda r: (0 if r["owner"] == "r_xy" else 1, -r["word_letter"], 0 if r["cable_sign"] == "negative" else 1),
    )
    for k, r in enumerate(order):
        r["thxy_index"] = k
    records.sort(key=lambda r: r["public_index"])
    return records


# ---------------------------------------------------------------------------
# transport along the actual cabled word
# ---------------------------------------------------------------------------


def block_for(pattern: tuple[str, str], positive: bool, cache: dict[Any, Matrix]) -> Matrix:
    key = (pattern, positive)
    if key not in cache:
        cache[key], _ = one_defect_block(pattern[0], pattern[1], positive)
    return cache[key]


def letterwise_transport(
    word: list[int],
    records: list[dict[str, Any]],
    position_monomial: Callable[[int], LP],
    pivotal: dict[str, LP],
) -> dict[str, Any]:
    """Verify B_pub M_t = M_{t+1} B_or for every letter, with strand tracking.

    Coordinates: x_pub[p] = position_monomial(p) * pivotal[factor of strand at p] * x_or[p].
    """
    strands = list(range(DIM))  # strand id = initial public position
    factor = {r["public_index"]: r["tensor_factor"] for r in records}
    cache: dict[Any, Matrix] = {}
    counts: dict[str, int] = {}
    failures: list[dict[str, Any]] = []
    burau = {True: burau_block(True), False: burau_block(False)}
    for step, letter in enumerate(word):
        i = abs(letter)
        p, p1 = i - 1, i
        if p < 0 or p1 >= DIM:
            raise ValueError(f"letter out of range at step {step}: {letter}")
        s1, s2 = strands[p], strands[p1]
        X, Y = factor[s1], factor[s2]
        positive = letter > 0
        b_or = block_for((X, Y), positive, cache)
        m_before = diag2(position_monomial(p) * pivotal[X], position_monomial(p1) * pivotal[Y])
        m_after = diag2(position_monomial(p) * pivotal[Y], position_monomial(p1) * pivotal[X])
        ok = mat_eq(mat_mul(burau[positive], m_before), mat_mul(m_after, b_or))
        key = f"{X}{'+' if positive else '-'}{Y}"
        counts[key] = counts.get(key, 0) + 1
        if not ok and len(failures) < 5:
            failures.append({"step": step, "letter": letter, "pattern": key})
        strands[p], strands[p1] = s2, s1
    return {
        "letters": len(word),
        "pattern_counts": counts,
        "all_letters_pass": not failures,
        "first_failures": failures,
        "final_permutation_is_identity": strands == list(range(DIM)),
    }


def h_burau_apply(word: list[int], vector: list[list[int]], degree: int = DEGREE) -> list[list[int]]:
    """Public Burau action in the variable h (q = 1 + h, t = q^-2), truncated mod h^(degree+1)."""
    t = LP.q(-2).h_series(degree)
    tinv = LP.q(2).h_series(degree)

    def mul(a: list[int], b: list[int]) -> list[int]:
        out = [0] * (degree + 1)
        for i, x in enumerate(a):
            if x == 0:
                continue
            for j in range(degree + 1 - i):
                out[i + j] += x * b[j]
        return out

    result = [poly[:] for poly in vector]
    for letter in reversed(word):
        i = abs(letter) - 1
        x = result[i]
        y = result[i + 1]
        if letter > 0:
            # e_i -> (1-t) e_i + e_{i+1},  e_{i+1} -> t e_i
            diff = [yy - xx for xx, yy in zip(x, y)]
            result[i] = [xx + d for xx, d in zip(x, mul(t, diff))]
            result[i + 1] = x
        else:
            diff = [xx - yy for xx, yy in zip(x, y)]
            result[i] = y
            result[i + 1] = [yy + d for yy, d in zip(y, mul(tinv, diff))]
    return result


def h_vector(coords: dict[int, LP], degree: int = DEGREE) -> list[list[int]]:
    vec = [[0] * (degree + 1) for _ in range(DIM)]
    for index, value in coords.items():
        vec[index] = value.h_series(degree)
    return vec


def h_row_apply(row: dict[int, LP], vector: list[list[int]], degree: int = DEGREE) -> list[int]:
    out = [0] * (degree + 1)
    for index, value in row.items():
        coeff = value.h_series(degree)
        for i, x in enumerate(coeff):
            if x == 0:
                continue
            for j in range(degree + 1 - i):
                out[i + j] += x * vector[index][j]
    return out


def constant_terms(coords: dict[int, LP]) -> list[list[int]]:
    """Constant term (q = 1) of a coordinate vector, as [[index, coefficient], ...]."""
    terms = []
    for index in sorted(coords):
        value = sum(coords[index].c.values())
        if value != 0:
            terms.append([index, value])
    return terms


# ---------------------------------------------------------------------------
# main computation
# ---------------------------------------------------------------------------


def compute(sigma: int = -1, verbose: bool = True) -> tuple[dict[str, Any], dict[str, Any]]:
    t0 = time.time()
    rec = load_recompute()
    cut_bytes = CUT_OBJECT.read_bytes()
    table_bytes = POSITION_TABLE.read_bytes()
    input_bytes = PUBLIC_INPUT.read_bytes()
    cut = json.loads(cut_bytes)
    table = json.loads(table_bytes)
    public_input = json.loads(input_bytes)

    # 1. model checks
    model_checks: dict[str, Any] = {}
    model_checks["intertwining"] = verify_intertwining()
    model_checks["inverses"] = verify_inverses()
    model_checks["braid_relation"] = verify_braid_relation()
    model_checks["phi"] = verify_phi()
    twist = derive_pivotal_twist(sigma)
    duality = verify_duality(sigma, twist)
    loop = duality.pop("loop value")
    model_checks["duality"] = duality
    model_checks["pivotal_twist_derived"] = twist
    model_checks["loop_value"] = repr(loop)
    model_checks["loop_value_at_q1"] = sum(loop.c.values())
    blocks: dict[str, Any] = {}
    for X in (V, VD):
        for Y in (V, VD):
            for positive in (True, False):
                block, info = one_defect_block(X, Y, positive)
                blocks[f"{X}{'+' if positive else '-'}{Y}"] = {"block_columns_are_images": mat_json(block), **info}
    model_checks["one_defect_blocks"] = blocks
    ratio = derive_position_ratio()
    model_checks["position_ratio_d(p+1)/d(p)"] = repr(ratio)
    ratio_coeff, ratio_exp = ratio.monomial_data()
    if ratio_coeff != 1:
        raise AssertionError("position ratio has a sign; unexpected")

    def position_monomial(p: int) -> LP:
        return LP.q(ratio_exp * p)

    phi = phi_matrix()
    # pivotal coefficient of the V* defect state under phi (row v_-, column v_+^*)
    phi_defect = phi[LABELS.index(MINUS)][LABELS.index(PLUS)]
    phi_ground = phi[LABELS.index(PLUS)][LABELS.index(MINUS)]
    if phi_ground != LP.one():
        raise AssertionError("phi does not send the V* ground state to the V ground state with coefficient 1")
    pivotal = {V: LP.one(), VD: phi_defect}
    model_all_pass = (
        all(model_checks["intertwining"].values())
        and all(model_checks["inverses"].values())
        and all(model_checks["braid_relation"].values())
        and all(model_checks["phi"].values())
        and all(v is True for v in duality.values())
        and all(b["ground_is_one"] and b["no_leakage"] for b in blocks.values())
    )
    if verbose:
        print(f"MODEL_CHECKS={'PASS' if model_all_pass else 'FAIL'} ({time.time()-t0:.1f}s)")

    # 2. endpoint records and the monomial transport P
    records = build_endpoint_records(cut, table)
    slot_of_public = {r["public_index"]: r["geometric_order"] for r in records}
    thxy_of_public = {r["public_index"]: r["thxy_index"] for r in records}
    for r in records:
        p = r["public_index"]
        m = position_monomial(p) * pivotal[r["tensor_factor"]]
        coeff, exp = m.monomial_data()
        piv_coeff, piv_exp = pivotal[r["tensor_factor"]].monomial_data()
        r["pivotal_coefficient"] = {"sign": piv_coeff, "q_exponent": piv_exp, "monomial": repr(pivotal[r["tensor_factor"]])}
        r["position_monomial"] = repr(position_monomial(p))
        r["transport_monomial"] = {"sign": coeff, "q_exponent": exp, "monomial": repr(m)}
        r["weight_defect_basis_vector"] = (
            f"e_{p}: ground states at every other position, {r['defect_state']} at public position {p} "
            f"(geometric slot {r['geometric_order']}, strand {r['physical_endpoint_id']})"
        )
    # P as a monomial matrix: x_pub[p] = transport_monomial(p) * x_geom[slot_of_public[p]]
    P_entries = [
        {"row_public": r["public_index"], "col_geometric": r["geometric_order"], "monomial": r["transport_monomial"]["monomial"]}
        for r in records
    ]

    # 3. the actual cabled word
    b44, _ = rec.build_oriented_b44(public_input)
    b88 = rec.cable_word(b44)
    integrity = public_input["point_push"]["derived_integrity"]
    if len(b88) != integrity["B88_length"] or rec.canonical_sha(b88) != integrity["B88_sha256"]:
        raise AssertionError("B88 word does not match the public integrity record")
    letter_result = letterwise_transport(b88, records, position_monomial, pivotal)
    if verbose:
        print(f"LETTERWISE_TRANSPORT={'PASS' if letter_result['all_letters_pass'] and letter_result['final_permutation_is_identity'] else 'FAIL'} ({time.time()-t0:.1f}s)")

    # 4. cup and cap in the oriented model, transported to public coordinates
    by_passage = {r["passage_id"]: r for r in records}
    cup_a, cup_b = (by_passage[pid] for pid in SELECTED_CUP_PASSAGES)
    if not (cup_a["tensor_factor"] == V and cup_b["tensor_factor"] == VD):
        raise AssertionError("selected cup does not join one entry (V) and one exit (V*) endpoint")
    if not cup_a["public_index"] < cup_b["public_index"]:
        raise AssertionError("selected cup endpoints are not in increasing public order")
    # u_geometric: coev_left on (V at A, V* at B); one-defect components
    coev = coev_left()
    u_geom: dict[int, LP] = {}
    u_geom[cup_a["geometric_order"]] = coev[(MINUS, MINUS)]   # defect v_- at A, ground v_-^* at B
    u_geom[cup_b["geometric_order"]] = coev[(PLUS, PLUS)]     # ground v_+ at A, defect v_+^* at B
    # ell_geometric: ev_right on (V at A, V* at B)
    ell_geom: dict[int, LP] = {}
    ell_geom[cup_a["geometric_order"]] = ev_right(MINUS, MINUS, sigma, twist)
    ell_geom[cup_b["geometric_order"]] = ev_right(PLUS, PLUS, sigma, twist)
    # transport
    transport_of_public = {r["public_index"]: LP.mono(r["transport_monomial"]["sign"], r["transport_monomial"]["q_exponent"]) for r in records}
    public_of_slot = {r["geometric_order"]: r["public_index"] for r in records}
    u_pub = {public_of_slot[s]: transport_of_public[public_of_slot[s]] * v for s, v in u_geom.items()}
    ell_pub = {public_of_slot[s]: v * transport_of_public[public_of_slot[s]].inverse_monomial() for s, v in ell_geom.items()}
    u_terms = constant_terms(u_pub)
    ell_terms = constant_terms(ell_pub)
    # normalise so that the cup coefficient at the lower public index is +1 (overall scalar of the
    # coevaluation is conventional and cancels against the cap only through the loop value)
    loop_pub = sum((ell_pub[p] * u_pub[p]).h_series(0)[0] for p in u_pub if p in ell_pub)
    if verbose:
        print(f"U_PUBLIC_CONSTANT={u_terms}  ELL_PUBLIC_CONSTANT={ell_terms}  LOOP_AT_Q1={loop_pub}")

    # 5. the public Burau matrix: rho(W) - I in h^3 End(E)
    t1 = time.time()
    min_order = None
    zero_low = True
    delta_matrix_eps: list[list[list[int]]] = []
    for col in range(DIM):
        vec = rec.sparse_vector(DIM, DEGREE, [[col, 1]])
        delta = rec.delta_apply(b88, vec)
        delta_matrix_eps.append(delta)
        for poly in delta:
            for d in range(3):
                if poly[d] != 0:
                    zero_low = False
            nz = next((d for d, x in enumerate(poly) if x != 0), None)
            if nz is not None and (min_order is None or nz < min_order):
                min_order = nz
    a_in_h3 = zero_low
    if verbose:
        print(f"RHO_W_MINUS_I_IN_H3={'PASS' if a_in_h3 else 'FAIL'} min_order={min_order} ({time.time()-t1:.1f}s)")

    # 6. divided cubic: constant terms (public pipeline), full transported series (h-model),
    #    and the direct matrix formula ell_0 A_3 u_0 (A_3 = [eps^3] coefficient, [h]eps = -2)
    vec = rec.sparse_vector(DIM, DEGREE, u_terms)
    eta = rec.apply_covector(rec.delta_apply(b88, vec), ell_terms)
    eta_h = rec.substitute_epsilon_with_h(eta, DEGREE)
    delta3_constant = eta_h[3]
    full_vec = h_burau_apply(b88, h_vector(u_pub))
    base_vec = h_vector(u_pub)
    diff_vec = [[a - b for a, b in zip(x, y)] for x, y in zip(full_vec, base_vec)]
    full_series = h_row_apply(ell_pub, diff_vec)
    delta3_full = full_series[3]
    a3_u0 = [0] * DIM
    for col, coeff in u_terms:
        for row in range(DIM):
            a3_u0[row] += coeff * delta_matrix_eps[col][row][3]
    ell0_a3_u0_eps = sum(coeff * a3_u0[row] for row, coeff in ell_terms)
    delta3_matrix = (-2) ** 3 * ell0_a3_u0_eps
    if verbose:
        print(f"DELTA3 constant={delta3_constant} full_series={delta3_full} matrix={delta3_matrix} ({time.time()-t0:.1f}s)")

    # 7. coordinate controls: partial transports reproduce the withdrawn values
    def cubic_of(u_t: list[list[int]], l_t: list[list[int]]) -> int:
        v = rec.sparse_vector(DIM, DEGREE, u_t)
        e = rec.apply_covector(rec.delta_apply(b88, v), l_t)
        return rec.substitute_epsilon_with_h(e, DEGREE)[3]

    a_pub, b_pub = cup_a["public_index"], cup_b["public_index"]
    a_slot, b_slot = cup_a["geometric_order"], cup_b["geometric_order"]
    a_thxy, b_thxy = cup_a["thxy_index"], cup_b["thxy_index"]
    ua, ub = dict(u_terms)[a_pub], dict(u_terms)[b_pub]
    la, lb = dict(ell_terms)[a_pub], dict(ell_terms)[b_pub]
    controls = {
        "collar_u_collar_ell_collar_word": cubic_of([[a_pub, ua], [b_pub, ub]], [[a_pub, la], [b_pub, lb]]),
        "thxy_u_collar_ell_collar_word": cubic_of([[a_thxy, ua], [b_thxy, ub]], [[a_pub, la], [b_pub, lb]]),
        "thxy_u_thxy_ell_collar_word": cubic_of([[a_thxy, ua], [b_thxy, ub]], [[a_thxy, la], [b_thxy, lb]]),
        "slot_u_collar_ell_collar_word": cubic_of([[a_slot, ua], [b_slot, ub]], [[a_pub, la], [b_pub, lb]]),
        "slot_u_slot_ell_collar_word": cubic_of([[a_slot, ua], [b_slot, ub]], [[a_slot, la], [b_slot, lb]]),
    }
    # simultaneous transport by the slot permutation (W transported too): invariant by construction
    # of the pairing; checked numerically by relabelling the word is impossible (the permutation is
    # not adjacency preserving), so it is asserted through the identity ell P^-1 (P W P^-1) P u = ell W u.
    controls["simultaneous_transport_invariance"] = "ell P^-1 (P W P^-1) (P u) = ell W u holds identically; see Smooth4PC/FilteredCubicNaturality.lean (pairingCoeff_transport)"

    convention = {
        "schema": SCHEMA_CONVENTION,
        "sources": {
            "CUT_OBJECT.json": {"path": "evidence/public_geometry/CUT_OBJECT.json", "sha256": sha256_bytes(cut_bytes)},
            "B88_POSITION_TO_PASSAGE_TABLE.json": {"path": "data/B88_POSITION_TO_PASSAGE_TABLE.json", "sha256": sha256_bytes(table_bytes)},
        },
        "geometric_order_source": "CUT_OBJECT gate_passages[y].boundary_endpoint_order_top_to_bottom.west (frozen MWW cut slot order, slot 0 = top)",
        "public_order_source": "B88_POSITION_TO_PASSAGE_TABLE.json index = 2*(wicket-1) + (0 negative copy, 1 positive copy); collar wicket order used by the cabled Artin word",
        "thxy_order_rule": "owner blocks (r_xy, m_2), wickets in decreasing word-letter order, copies (negative, positive) adjacent; reconstructed from the withdrawn -59072 and -2496 controls",
        "orientation_rule": "west endpoint 'entry' = strand oriented into the cut tangle = tensor factor V (ground v_+, defect v_-); 'exit' = tensor factor V* (ground v_-^*, defect v_+^*)",
        "oriented_model": {
            "algebra": "U_q(sl2), Delta(E)=E(x)K+1(x)E, Delta(F)=F(x)1+K^-1(x)F, S(E)=-EK^-1, S(F)=-KF; dual action (x.f)(v)=f(S(x)v)",
            "braiding": "c = q^{-1/2} tau o q^{H(x)H/2}(1+(q-q^-1)E(x)F); every exponent integral on V, V*",
            "phi": "phi: V* -> V, phi(v_-^*)=v_+, phi(v_+^*)=-q v_- (module map, verified)",
            "coevaluation": "u: coev_V : 1 -> V(x)V*, sum_i v_i (x) v_i^* (canonical)",
            "evaluation": f"ell: ev'_V : V(x)V* -> 1, v_i (x) v_j^* -> sigma delta_ij q^{{{twist}*weight(v_i)}} with sigma={sigma}",
            "pivotal_sign_sigma": sigma,
            "pivotal_twist_derived": twist,
            "loop_value": repr(loop),
            "loop_value_at_q1": sum(loop.c.values()),
            "sigma_note": "sigma=-1 is the pivotal structure of the frozen public receipt (loop value -[2] at q=1, cap row e_87^*-e_2^*); sigma=+1 gives the standard loop value [2] and negates ell, hence negates the cubic",
        },
        "public_model": {
            "burau_block_columns_are_images": mat_json(burau_block(True)),
            "parameter": "t = q^-2, q = 1 + h",
            "position_monomial": f"d_p = q^{{{ratio_exp}*p}} (derived: unique monomial ratio conjugating the oriented V(x)V blocks to the Burau blocks)",
        },
        "transport_rule": "x_public[p] = transport_monomial(p) * x_geometric[geometric_order(p)], transport_monomial(p) = d_p * pivotal(tensor factor at p)",
        "selected_cup": {
            "passages": list(SELECTED_CUP_PASSAGES),
            "endpoints": [cup_a["physical_endpoint_id"], cup_b["physical_endpoint_id"]],
            "public_positions": [a_pub, b_pub],
            "geometric_slots": [a_slot, b_slot],
            "thxy_indices": [a_thxy, b_thxy],
        },
        "endpoints": records,
        "P_monomial_matrix": P_entries,
    }
    convention["convention_sha256"] = canonical_sha(convention)

    audit = {
        "schema": SCHEMA_AUDIT,
        "convention_sha256": convention["convention_sha256"],
        "input_sha256": sha256_bytes(input_bytes),
        "B88_length": len(b88),
        "B88_sha256": rec.canonical_sha(b88),
        "model_checks": model_checks,
        "model_checks_pass": model_all_pass,
        "letterwise_transport": letter_result,
        "geometric_to_public_permutation": [slot_of_public[p] for p in range(DIM)],
        "thxy_of_public": [thxy_of_public[p] for p in range(DIM)],
        "u_geometric_slot_coordinates": {str(k): v.to_json() for k, v in sorted(u_geom.items())},
        "ell_geometric_slot_coordinates": {str(k): v.to_json() for k, v in sorted(ell_geom.items())},
        "u_public": {str(k): v.to_json() for k, v in sorted(u_pub.items())},
        "ell_public": {str(k): v.to_json() for k, v in sorted(ell_pub.items())},
        "u_public_constant_terms": u_terms,
        "ell_public_constant_terms": ell_terms,
        "loop_value_at_q1_public": loop_pub,
        "rho_W_minus_I_in_h3_End": a_in_h3,
        "rho_W_minus_I_min_eps_order": min_order,
        "delta3": {
            "constant_terms_pipeline": delta3_constant,
            "full_transported_series_h3": delta3_full,
            "ell0_A3_u0_matrix_formula": delta3_matrix,
            "agree": delta3_constant == delta3_full == delta3_matrix,
        },
        "coordinate_controls": controls,
        "elapsed_seconds": round(time.time() - t0, 1),
    }
    audit["endpoint_transport_pass"] = bool(
        model_all_pass
        and letter_result["all_letters_pass"]
        and letter_result["final_permutation_is_identity"]
        and a_in_h3
        and audit["delta3"]["agree"]
    )
    audit["no_unresolved_signs"] = bool(
        all(r["transport_monomial"]["sign"] in (1, -1) for r in records)
        and len(u_terms) == 2
        and len(ell_terms) == 2
        and model_all_pass
    )
    audit["audit_sha256"] = canonical_sha(audit)
    return convention, audit


def derive_endpoint_terms(convention_path: Path = CONVENTION) -> dict[str, Any]:
    """Public entry point for recompute_t73_delta3.py: constant cup/cap terms from the authority file.

    The terms are recomputed from the recorded transport monomials and the oriented cup/cap, not
    read from a stored list.
    """
    convention = load_json(convention_path)
    if convention.get("schema") != SCHEMA_CONVENTION:
        raise ValueError("unsupported endpoint convention schema")
    payload = {k: v for k, v in convention.items() if k != "convention_sha256"}
    if canonical_sha(payload) != convention["convention_sha256"]:
        raise ValueError("endpoint convention SHA mismatch")
    sigma = convention["oriented_model"]["pivotal_sign_sigma"]
    twist = convention["oriented_model"]["pivotal_twist_derived"]
    by_passage = {r["passage_id"]: r for r in convention["endpoints"]}
    cup_a, cup_b = (by_passage[pid] for pid in convention["selected_cup"]["passages"])
    coev = coev_left()
    u_pub: dict[int, LP] = {}
    ell_pub: dict[int, LP] = {}
    for endpoint, (cu, label) in ((cup_a, (coev[(MINUS, MINUS)], MINUS)), (cup_b, (coev[(PLUS, PLUS)], PLUS))):
        m = LP.mono(endpoint["transport_monomial"]["sign"], endpoint["transport_monomial"]["q_exponent"])
        p = endpoint["public_index"]
        u_pub[p] = m * cu
        ell_pub[p] = ev_right(label, label, sigma, twist) * m.inverse_monomial()
    return {
        "u_terms": constant_terms(u_pub),
        "ell_terms": constant_terms(ell_pub),
        "u_public": {str(k): v.to_json() for k, v in sorted(u_pub.items())},
        "ell_public": {str(k): v.to_json() for k, v in sorted(ell_pub.items())},
        "convention_sha256": convention["convention_sha256"],
    }


def print_summary(audit: dict[str, Any]) -> None:
    print(f"CONVENTION_SHA256={audit['convention_sha256']}")
    print(f"GEOMETRIC_TO_PUBLIC_PERMUTATION={json.dumps(audit['geometric_to_public_permutation'], separators=(',', ':'))}")
    print(f"U_PUBLIC_CONSTANT={json.dumps(audit['u_public_constant_terms'], separators=(',', ':'))}")
    print(f"ELL_PUBLIC_CONSTANT={json.dumps(audit['ell_public_constant_terms'], separators=(',', ':'))}")
    print(f"LETTERWISE_PATTERNS={json.dumps(audit['letterwise_transport']['pattern_counts'], sort_keys=True, separators=(',', ':'))}")
    print(f"COORDINATE_CONTROLS={json.dumps({k: v for k, v in audit['coordinate_controls'].items() if isinstance(v, int)}, sort_keys=True, separators=(',', ':'))}")
    print(f"RHO_W_MINUS_I_IN_H3=PASS" if audit["rho_W_minus_I_in_h3_End"] else "RHO_W_MINUS_I_IN_H3=FAIL")
    print(f"ENDPOINT_TRANSPORT={'PASS' if audit['endpoint_transport_pass'] else 'FAIL'}")
    print(f"NO_UNRESOLVED_SIGNS={'PASS' if audit['no_unresolved_signs'] else 'FAIL'}")
    print(f"DELTA3={audit['delta3']['constant_terms_pipeline']}")
    print(f"AUDIT_SHA256={audit['audit_sha256']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the convention and audit files")
    parser.add_argument("--check", action="store_true", help="recompute and compare with the committed files")
    parser.add_argument("--sigma", type=int, default=-1, choices=(1, -1), help="pivotal sign (default: frozen public convention)")
    args = parser.parse_args()
    convention, audit = compute(sigma=args.sigma)
    if args.write:
        CONVENTION.write_text(canonical_json(convention), encoding="utf-8", newline="\n")
        AUDIT.write_text(canonical_json(audit), encoding="utf-8", newline="\n")
        print(f"WROTE={CONVENTION}")
        print(f"WROTE={AUDIT}")
    if args.check:
        committed_convention = load_json(CONVENTION)
        committed_audit = load_json(AUDIT)
        if committed_convention != convention:
            raise AssertionError("committed endpoint convention differs from regeneration")
        volatile = {"elapsed_seconds", "audit_sha256"}
        if {k: v for k, v in committed_audit.items() if k not in volatile} != {k: v for k, v in audit.items() if k not in volatile}:
            raise AssertionError("committed endpoint transport audit differs from regeneration")
        print("COMMITTED_FILES=MATCH")
    print_summary(audit)
    if not audit["endpoint_transport_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
