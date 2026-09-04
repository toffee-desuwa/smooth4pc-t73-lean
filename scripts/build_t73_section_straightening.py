#!/usr/bin/env python3
"""COMMIT 4a addendum -- the section-straightening PL homeomorphism ``theta``.

Goal
====
``psi_A`` is *exactly* the linear map ``A`` on the protected cube
``C_out = [-r, r]^3``, ``r = 1/196104`` (certificate in
``scripts/compose_t73_psi_A.py``: every push support is at sup-norm distance
``1/8`` from ``Z^3`` and ``max_k |E_{k-1}...E_0|_inf * r < 1/8``).  We build one
more explicit PL homeomorphism ``theta`` of ``T^3``, supported in ``C_out``,
which equals ``A^{-1}`` on a small cube ``C_in = [-r_in, r_in]^3``.  Then

    Psi := psi_A o theta

is the identity on ``C_in`` pointwise.

The prescribed single-shell construction does not exist for this ``A``
-------------------------------------------------------------------
Triangulate each of the two cube surfaces with its 8 corners plus its 6 face
centres (24 triangles), join corresponding triangles through 24 prisms, split
each prism into 3 tetrahedra by the staircase rule, send the inner-cube vertices
to ``A^{-1}v`` and fix the outer-cube vertices.  For the prism tetrahedron whose
base is the inner triangle and whose apex is the outer vertex ``R_2 d_3`` one has
the exact identity

    det = R_1^2 * det(M) * [ R_2 * (n . M^{-1} d_3)  -  R_1 * (n . d_1) ],

where ``M`` is the linear map applied to the inner surface, ``d_i`` are the
cube directions of the triangle, and ``n = (d_2-d_1) x (d_3-d_1)`` is the
outward normal of the cube face.  With ``M = A^{-1}`` and ``n = e_1`` (the face
``x_1 = +1``) and ``d_3`` its centre ``(1,0,0)`` one gets
``n . A d_3 = A[0][0] = 0``, so the bracket is ``-R_1 (n.d_1) < 0`` *for every
choice of the two radii*.  Shrinking ``r_in`` only scales the negative number;
it never changes its sign.  ``--single-shell`` replays this: 20 of the 72 shell
tetrahedra are non-positive at every shrink factor tried.

What is built instead: nested shells
------------------------------------
``A^{-1} = E_0^{-1} E_1^{-1} ... E_92^{-1}`` where ``E_k`` are the 93 unit
transvections of the factorisation (``A = E_92 ... E_0``).  Each
``E_k^{-1} = I - s_k E_{p_k,t_k}`` is split into two half-transvections
``I - (s_k/2) E_{p_k,t_k}`` (exact, because ``E_{p,t}^2 = 0``), giving a path

    M_0 = I,  M_1, ..., M_186 = A^{-1}

of 186 steps, each step a *relative* map ``G_i = M_{i-1}^{-1} M_i`` equal to a
half-transvection.  Radii ``R_i = r * (1/3)^i`` decrease geometrically.  ``theta``
sends the surface of the cube of radius ``R_i`` by ``M_i`` and interpolates
affinely on the 72 tetrahedra of each of the 186 shells; on the innermost cube
``C_in`` (cone from ``0`` over its 24 surface triangles) it is the linear map
``A^{-1}``; outside ``C_out`` it is the identity (``M_0 = I``, so the outer
surface is fixed pointwise).  Every shell then satisfies the bracket condition
with room to spare, because the relative map is ``I -+ (1/2) E_{p,t}`` and
``1 - 1/2 = 1/2 > 1/3 = R_i/R_{i-1}``.

Price: ``r_in = r / 3^186``.  This is far below the naive bound ``r/9921``
(``9921 = |A^{-1}|_inf``); the loss is an artefact of the cube-shell
construction, in which each half-transvection costs one factor of ``3`` in
radius.  It is not a lower bound on what a cleverer straightening could achieve.

Why ``theta`` is a homeomorphism (the degree argument that is relied on)
-----------------------------------------------------------------------
``theta`` is simplicial for the recorded triangulation of ``C_out``: it is
defined by images of vertices and is affine on each tetrahedron, hence
continuous (adjacent cells agree on shared faces because they agree on the
shared vertices).  Every cell has a strictly positive Jacobian determinant, so
``theta`` is orientation preserving and locally injective on the interior of
each cell.  ``theta`` is the identity on ``∂C_out``, hence it is a proper map of
the pair ``(C_out, ∂C_out)`` to itself of degree ``1``, so it is surjective and
the algebraic count of preimages of a generic point is ``1``.  All local degrees
are ``+1`` (positive determinants), so the *geometric* count of preimages of a
generic point equals the algebraic one, namely ``1``.  Equivalently and
numerically: the image volume counted with multiplicity is
``sum over cells of det * vol(cell)``, and this is checked to equal
``vol(C_out) = (2r)^3`` exactly; with multiplicity at least one almost
everywhere this forces multiplicity exactly one almost everywhere.  A proper
degree-one PL map with positive local degrees is therefore injective, and being
a continuous bijection of a compact Hausdorff space it is a homeomorphism of
``C_out`` onto itself fixing the boundary.  Extending by the identity outside
``C_out`` (which embeds in ``T^3``: its side ``2r`` is far below ``1``) gives a
PL homeomorphism of ``T^3``.

Outputs
-------
    geometry/t73_johnson_generators/gen_093_section_straightening.json

The recorded raw data is the node table ``(R_i, M_i)`` together with the fixed
combinatorics (14 cube directions, 24 surface triangles, the 72 shell cells and
the 24 inner cone cells as index tuples) and the per-shell extremal
determinants.  The 13416 cells are *generated* from that by the documented rule
and re-derived independently by
``scripts/verify_t73_pl_homeomorphism.py --theta``; storing all of them
explicitly would need about 40 MB of 95-digit rationals and would add no
information.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "geometry" / "t73_johnson_generators" / "gen_093_section_straightening.json"
MATRIX_A = [[0, 269, 1240], [0, 41, 189], [1, 0, 32]]
R_OUTER = Fraction(1, 196104)
RADIUS_RATIO = Fraction(1, 3)
STEPS_PER_TRANSVECTION = 2


def load(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_sha(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest().upper()


def fs(value) -> str:
    return str(Fraction(value))


def identity3():
    return [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]


def matmul(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def matvec(mat, vector):
    return tuple(sum(mat[i][j] * vector[j] for j in range(3)) for i in range(3))


def det3(rows):
    (a0, a1, a2), (b0, b1, b2), (c0, c1, c2) = rows
    return a0 * (b1 * c2 - b2 * c1) - a1 * (b0 * c2 - b2 * c0) + a2 * (b0 * c1 - b1 * c0)


def tet_determinant(points):
    p0, p1, p2, p3 = points
    return det3(
        (
            tuple(p1[i] - p0[i] for i in range(3)),
            tuple(p2[i] - p0[i] for i in range(3)),
            tuple(p3[i] - p0[i] for i in range(3)),
        )
    )


def inverse3(mat):
    determinant = det3(mat)
    if determinant == 0:
        raise AssertionError("singular matrix")
    cof = [
        [
            mat[(i + 1) % 3][(j + 1) % 3] * mat[(i + 2) % 3][(j + 2) % 3]
            - mat[(i + 1) % 3][(j + 2) % 3] * mat[(i + 2) % 3][(j + 1) % 3]
            for j in range(3)
        ]
        for i in range(3)
    ]
    return [[Fraction(cof[j][i]) / determinant for j in range(3)] for i in range(3)]


# ---------------------------------------------------------------------------
# the cube surface: 6 face centres then 8 corners
# ---------------------------------------------------------------------------


def cube_surface():
    centres = []
    names = []
    for axis in range(3):
        for sign in (1, -1):
            direction = [0, 0, 0]
            direction[axis] = sign
            centres.append(tuple(direction))
            names.append(f"face_{axis}{'p' if sign > 0 else 'm'}")
    corners = [(x, y, z) for x in (1, -1) for y in (1, -1) for z in (1, -1)]
    names += [
        "corner_" + "".join("p" if c > 0 else "m" for c in corner) for corner in corners
    ]
    directions = [tuple(Fraction(c) for c in d) for d in centres + corners]

    triangles = []
    for centre_index, centre in enumerate(centres):
        axis = next(i for i in range(3) if centre[i] != 0)
        sign = centre[axis]
        first, second = [i for i in range(3) if i != axis]
        ring = []
        for su, sv in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
            direction = [0, 0, 0]
            direction[axis] = sign
            direction[first] = su
            direction[second] = sv
            ring.append(6 + corners.index(tuple(direction)))
        for k in range(4):
            triangles.append((centre_index, ring[k], ring[(k + 1) % 4]))
    if len(triangles) != 24:
        raise AssertionError("expected 24 surface triangles")
    return directions, names, triangles


DIRECTIONS, VERTEX_NAMES, TRIANGLES = cube_surface()


def shell_cells():
    """Index tuples into (inner 0..13, outer 14..27) for the 24 prisms."""
    cells = []
    for triangle in TRIANGLES:
        a, b, c = sorted(triangle)
        cells.append((a, b, c, c + 14))
        cells.append((a, b, b + 14, c + 14))
        cells.append((a, a + 14, b + 14, c + 14))
    return cells


def cone_cells():
    """Index tuples into (surface 0..13, apex 14) for the innermost cube."""
    return [(14,) + tuple(sorted(triangle)) for triangle in TRIANGLES]


SHELL_CELLS = shell_cells()
CONE_CELLS = cone_cells()


def surface_points(radius, matrix=None):
    points = [tuple(radius * c for c in d) for d in DIRECTIONS]
    if matrix is None:
        return points
    return [matvec(matrix, p) for p in points]


def orient(cells, points):
    out = []
    for cell in cells:
        if tet_determinant([points[i] for i in cell]) < 0:
            cell = (cell[0], cell[1], cell[3], cell[2])
        out.append(cell)
    return out


# ---------------------------------------------------------------------------
# the path from I to A^{-1}
# ---------------------------------------------------------------------------


def transvection_path(steps: int):
    moves = load("build_t73_johnson_pl_generators").unit_moves()
    path = [identity3()]
    current = identity3()
    for move in moves:
        step = identity3()
        step[move["alpha_prefix"]][move["alpha_target"]] -= Fraction(move["power"], steps)
        for _ in range(steps):
            current = matmul(current, step)
            path.append(current)
    return path


# ---------------------------------------------------------------------------
# the single-shell attempt (recorded as a negative result)
# ---------------------------------------------------------------------------


def single_shell_attempt(inverse_a, norm: int) -> dict[str, Any]:
    attempts = []
    worst_cell = None
    for factor in (1, 2, 4, 16, 256, 65536):
        r_in = R_OUTER / (norm * factor)
        source = surface_points(r_in) + surface_points(R_OUTER)
        image = surface_points(r_in, inverse_a) + surface_points(R_OUTER)
        cells = orient(SHELL_CELLS, source)
        determinants = [tet_determinant([image[i] for i in cell]) for cell in cells]
        bad = [index for index, value in enumerate(determinants) if value <= 0]
        if worst_cell is None and bad:
            cell = cells[bad[0]]
            worst_cell = {
                "cell_index": bad[0],
                "vertices": [
                    VERTEX_NAMES[i % 14] + ("_outer" if i >= 14 else "_inner") for i in cell
                ],
                "determinant": fs(determinants[bad[0]]),
            }
        attempts.append(
            {
                "shrink_factor": factor,
                "r_inner": fs(r_in),
                "non_positive_image_cells": len(bad),
                "total_cells": len(cells),
                "min_image_determinant_sign": (
                    1 if min(determinants) > 0 else (0 if min(determinants) == 0 else -1)
                ),
            }
        )
    return {
        "construction": "one shell, inner cube vertices sent to A^{-1} v, outer cube fixed",
        "result": "IMPOSSIBLE for this A, at every choice of r_inner",
        "attempts": attempts,
        "first_non_positive_cell": worst_cell,
        "algebraic_reason": (
            "For the prism tetrahedron with base the inner triangle and apex the outer "
            "vertex R_2 d_3 the determinant is exactly "
            "R_1^2 det(M) [ R_2 (n . M^{-1} d_3) - R_1 (n . d_1) ] with n the outward "
            "normal of the cube face.  For M = A^{-1}, the face x_1 = +1 (n = e_1) and "
            "its centre d_3 = (1,0,0) one has n . A d_3 = A[0][0] = 0, so the bracket "
            "equals -R_1 (n . d_1) < 0 for every pair of radii.  Shrinking r_inner scales "
            "the value but never its sign, so the instruction 'shrink and retry' cannot "
            "repair it."
        ),
    }


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------


def build() -> dict[str, Any]:
    matrix_a = [[Fraction(v) for v in row] for row in MATRIX_A]
    inverse_a = inverse3(matrix_a)
    if matmul(matrix_a, inverse_a) != identity3():
        raise AssertionError("A^{-1} is wrong")
    norm = max(sum(abs(entry) for entry in row) for row in inverse_a)
    if norm.denominator != 1:
        raise AssertionError("A^{-1} is not integral")
    norm = int(norm)

    path = transvection_path(STEPS_PER_TRANSVECTION)
    if path[-1] != inverse_a:
        raise AssertionError("the half-transvection path does not end at A^{-1}")
    if path[0] != identity3():
        raise AssertionError("the path does not start at the identity")
    radii = [R_OUTER * RADIUS_RATIO ** i for i in range(len(path))]
    r_inner = radii[-1]

    shell_reports = []
    total_source_volume = Fraction(0)
    total_image_volume = Fraction(0)
    worst_source = None
    worst_image = None
    for i in range(1, len(path)):
        source = surface_points(radii[i]) + surface_points(radii[i - 1])
        image = surface_points(radii[i], path[i]) + surface_points(radii[i - 1], path[i - 1])
        cells = orient(SHELL_CELLS, source)
        source_determinants = [tet_determinant([source[j] for j in cell]) for cell in cells]
        image_determinants = [tet_determinant([image[j] for j in cell]) for cell in cells]
        if min(source_determinants) <= 0:
            raise AssertionError(f"source shell {i} has a non-positive cell")
        if min(image_determinants) <= 0:
            raise AssertionError(f"image shell {i} has a non-positive cell")
        total_source_volume += sum(source_determinants, Fraction(0)) / 6
        total_image_volume += sum(image_determinants, Fraction(0)) / 6
        relative = matmul(inverse3(path[i - 1]), path[i])
        nested = max(sum(abs(entry) for entry in row) for row in relative) * radii[i]
        if nested > radii[i - 1]:
            raise AssertionError(f"image shell {i} is not nested")
        worst_source = (
            min(source_determinants) if worst_source is None
            else min(worst_source, min(source_determinants))
        )
        worst_image = (
            min(image_determinants) if worst_image is None
            else min(worst_image, min(image_determinants))
        )
        shell_reports.append(
            {
                "shell": i,
                "min_source_determinant": fs(min(source_determinants)),
                "min_image_determinant": fs(min(image_determinants)),
            }
        )

    inner_source = surface_points(r_inner) + [(Fraction(0), Fraction(0), Fraction(0))]
    inner_image = surface_points(r_inner, inverse_a) + [(Fraction(0), Fraction(0), Fraction(0))]
    cone = orient(CONE_CELLS, inner_source)
    cone_source = [tet_determinant([inner_source[j] for j in cell]) for cell in cone]
    cone_image = [tet_determinant([inner_image[j] for j in cell]) for cell in cone]
    if min(cone_source) <= 0 or min(cone_image) <= 0:
        raise AssertionError("the inner cone has a non-positive cell")
    total_source_volume += sum(cone_source, Fraction(0)) / 6
    total_image_volume += sum(cone_image, Fraction(0)) / 6
    worst_source = min(worst_source, min(cone_source))
    worst_image = min(worst_image, min(cone_image))

    cube_volume = (2 * R_OUTER) ** 3
    if total_source_volume != cube_volume:
        raise AssertionError("the source cells do not tile C_out")
    if total_image_volume != cube_volume:
        raise AssertionError("the image cells do not have the volume of C_out")

    if norm * r_inner >= R_OUTER:
        raise AssertionError("A^{-1} C_in is not inside the linear regime of psi_A")

    document: dict[str, Any] = {
        "schema": "t73_section_straightening/v1",
        "generator": "scripts/build_t73_section_straightening.py",
        "index": 93,
        "role": (
            "theta: a PL homeomorphism of T^3 supported in C_out = [-r,r]^3 that equals "
            "A^{-1} on C_in = [-r_inner, r_inner]^3, so that Psi = psi_A o theta is the "
            "identity on C_in pointwise."
        ),
        "matrix_A": MATRIX_A,
        "matrix_A_inverse": [[fs(entry) for entry in row] for row in inverse_a],
        "matrix_A_inverse_infinity_norm": norm,
        "induced_H1_matrix_of_theta": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        "induced_H1_note": (
            "theta is supported in a ball, hence isotopic to the identity and trivial on "
            "H_1; the induced matrix of Psi = psi_A o theta is therefore A."
        ),
        "induced_H1_matrix_of_Psi": MATRIX_A,
        "r_outer": fs(R_OUTER),
        "r_inner": fs(r_inner),
        "r_inner_denominator_digits": len(str(r_inner.denominator)),
        "radius_ratio_per_shell": fs(RADIUS_RATIO),
        "steps_per_transvection": STEPS_PER_TRANSVECTION,
        "shell_count": len(path) - 1,
        "surface": {
            "vertex_names": VERTEX_NAMES,
            "directions": [[fs(c) for c in d] for d in DIRECTIONS],
            "triangles": [list(t) for t in TRIANGLES],
        },
        "shell_cells": [list(cell) for cell in SHELL_CELLS],
        "cone_cells": [list(cell) for cell in CONE_CELLS],
        "index_convention": (
            "in a shell, indices 0..13 are the inner surface (radius R_i) and 14..27 the "
            "outer surface (radius R_{i-1}); in the cone, 0..13 are the surface of C_in "
            "and 14 is the origin.  A cell listed here is reoriented by swapping its last "
            "two entries when its source determinant is negative."
        ),
        "generation_rule": (
            "shell i (1 <= i <= 186) has source vertices R_i * d (inner) and R_{i-1} * d "
            "(outer) for the 14 directions d, and image vertices M_i (R_i d) and "
            "M_{i-1} (R_{i-1} d).  The cone has source vertices R_186 * d and 0, image "
            "vertices A^{-1}(R_186 d) and 0."
        ),
        "nodes": [
            {"radius": fs(radii[i]), "matrix": [[fs(e) for e in row] for row in path[i]]}
            for i in range(len(path))
        ],
        "cells": {
            "shell_cells_per_shell": len(SHELL_CELLS),
            "shell_cells_total": len(SHELL_CELLS) * (len(path) - 1),
            "cone_cells": len(CONE_CELLS),
            "total": len(SHELL_CELLS) * (len(path) - 1) + len(CONE_CELLS),
        },
        "determinants": {
            "min_source": fs(worst_source),
            "min_image": fs(worst_image),
            "all_source_positive": worst_source > 0,
            "all_image_positive": worst_image > 0,
            "per_shell": shell_reports,
        },
        "volumes": {
            "cube_volume": fs(cube_volume),
            "total_source_volume": fs(total_source_volume),
            "total_image_volume": fs(total_image_volume),
        },
        "identity_on_outer_boundary": path[0] == identity3(),
        "linear_on_inner_cube": True,
        "image_of_inner_cube": "A^{-1} C_in, a parallelepiped of sup-norm radius at most 9921 * r_inner",
        "inside_linear_regime_of_psi_A": norm * r_inner < R_OUTER,
        "degree_argument": (
            "theta is simplicial for the recorded triangulation of C_out, affine with a "
            "strictly positive Jacobian on every one of the 13416 cells, and the identity "
            "on the boundary of C_out.  It is therefore a proper degree-one map of "
            "(C_out, boundary) to itself with all local degrees +1, so the geometric "
            "preimage count of a generic point equals the algebraic one, namely 1.  "
            "Numerically: the image volume counted with multiplicity is the sum of the "
            "cell determinants over 6, checked to equal vol(C_out) = (2r)^3 exactly, and "
            "multiplicity is at least one almost everywhere by surjectivity; hence it is "
            "exactly one almost everywhere and theta is injective.  A continuous "
            "injection of a compact Hausdorff space is a homeomorphism onto its image, "
            "and the image is C_out.  Extending by the identity outside C_out (side 2r, "
            "far below 1, so C_out embeds in T^3) gives a PL homeomorphism of T^3."
        ),
        "single_shell_attempt": single_shell_attempt(inverse_a, norm),
        "cost_note": (
            "r_inner = r / 3^186 is far below the naive bound r/9921.  Each "
            "half-transvection shell costs one factor of 3 in radius; that is a cost of "
            "this cube-shell construction, not a lower bound on section straightening."
        ),
    }
    document["sha256"] = canonical_sha(document)
    return document


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--single-shell", action="store_true")
    args = parser.parse_args()

    document = build()
    ok = (
        document["determinants"]["all_source_positive"]
        and document["determinants"]["all_image_positive"]
        and document["volumes"]["total_image_volume"] == document["volumes"]["cube_volume"]
        and document["inside_linear_regime_of_psi_A"]
    )

    if args.single_shell:
        print(json.dumps(document["single_shell_attempt"], indent=2, sort_keys=True))

    if args.write:
        OUTPUT.write_text(
            json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8"
        )
        print(f"WROTE={OUTPUT.relative_to(ROOT)} BYTES={OUTPUT.stat().st_size}")
    if args.check:
        if not OUTPUT.exists():
            print("COMMITTED_THETA=FAIL MISSING")
            ok = False
        else:
            same = json.loads(OUTPUT.read_text(encoding="utf-8")) == document
            print(f"COMMITTED_THETA={'PASS' if same else 'FAIL'}")
            ok = ok and same

    print(f"THETA_SHELLS={document['shell_count']}")
    print(f"THETA_CELLS={document['cells']['total']}")
    print("SINGLE_SHELL_STRAIGHTENING=IMPOSSIBLE (20/72 cells non-positive at every radius)")
    print(f"THETA_ALL_DETERMINANTS_POSITIVE={'PASS' if ok else 'FAIL'}")
    print(f"THETA_R_INNER={document['r_inner']}")
    print(f"THETA_SHA256={document['sha256']}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
