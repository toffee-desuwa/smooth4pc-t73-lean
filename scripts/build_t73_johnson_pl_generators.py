#!/usr/bin/env python3
"""COMMIT 4a -- the 93 explicit PL Johnson generators ``psi_k = Pi_k o E_k``.

Model
=====
``T^3 = R^3 / Z^3`` with exact rational coordinates.  Every number produced by
this script is a ``fractions.Fraction``; no floating point is used anywhere.

Points and reduction modulo the lattice
---------------------------------------
A point of ``T^3`` is always stored as a *lift* in ``Q^3``.  The reduction map
is ``red(x)_i = x_i - floor(x_i)``, landing in the fundamental domain
``[0,1)^3``; two lifts describe the same torus point iff their difference lies
in ``Z^3``.

Cells are stored as lifts as well.  The support of a push ``Pi_k`` is contained
in an axis-parallel lift box whose side lengths are ``(3/4, 3/4, 1/4)``, all
``< 1``; therefore ``red`` restricted to that box is injective, no recorded cell
contains two lifts of the same torus point, and the recorded lift *is* the cell.
The box does cross the coordinate plane ``x_r = 0`` (the third-axis range is
``[-1/8, +1/8]``), i.e. the support does straddle a face of ``[0,1)^3``.  That is
precisely why cells are stored as lifts and never as reduced coordinates; the
field ``support_box.straddles_fundamental_domain_face`` records it, and
``scripts/verify_t73_pl_homeomorphism.py`` re-derives the injectivity from the
raw vertex list.

The generator
-------------
Unit move ``k`` of ``scripts/factor_t73_matrix_johnson.py`` has an
``alpha_target`` ``t``, an ``alpha_prefix`` ``p`` and a ``power`` ``s = +-1``.
Write ``r`` for the remaining axis.  Then

    psi_k = Pi_k o E_k          (E_k first, then Pi_k)

``E_k`` is the linear transvection ``x -> x + s * x_t * e_p``; as a matrix it is
``I + s * E_{p,t}`` (column ``t`` equals ``e_t + s e_p``), it is integral,
unimodular with determinant ``+1`` and lattice preserving, hence a PL
homeomorphism of ``T^3`` affine on every cell of every triangulation.  The
composition order is fixed by

    E_92 * E_91 * ... * E_0  =  A = [[0,269,1240],[0,41,189],[1,0,32]],

which is the order produced by ``apply_all(identity, construction_transvections)``
in ``factor_t73_matrix_johnson`` (a left-multiplication accumulation of row
operations) and by ``compose(local, mapping)`` in
``search_t73_johnson_alpha_sides.build_lift`` (whose abelianisation is
``M_outer * M_inner``).  So ``psi_A = psi_92 o ... o psi_0`` and ``psi_0`` is
applied to points first.  ``--check`` re-derives this.

``Pi_k`` is the square push.  Use the frame ``(e_a, e_b, e_c) = (e_t, s e_p, e_r)``
and coordinates ``(a,b,c)``; note ``x_t = a``, ``x_p = s b``, ``x_r = c``, so the
frame is a signed permutation of the standard basis and all lattice translates
are integral in ``(a,b,c)`` as well.  The unit square ``[0,1]^2`` of the ``(a,b)``
plane carries the planar triangulation with the nine vertices

    m  = (1/2,1/2)   P  = (1/4,1/4)   R1 = (7/8,1/8)
    Q  = (3/4,3/4)   R2 = (1/8,7/8)   and the four corners,

and twelve triangles: the four ``m``-triangles ``(m,P,R1) (m,R1,Q) (m,Q,R2)
(m,R2,P)`` plus eight triangles filling the annulus between the convex
quadrilateral ``P R1 Q R2`` and the square boundary.  The star of ``m`` is
therefore exactly the convex quadrilateral ``P R1 Q R2`` (4 triangles).  The slab
``{0<=a<=1, 0<=b<=1, |c| <= delta}`` with ``delta = 1/8`` is triangulated by
prisms over those triangles between the levels ``c = -delta, 0, +delta``, each
prism cut into 3 tetrahedra by the Kuhn/Freudenthal staircase rule "on a vertical
quadrilateral face over the edge ``{u,v}`` with ``u < v`` in the planar vertex
order, the diagonal joins the level-0 copy of ``u`` to the outer-level copy of
``v``"; the rule only depends on the planar order, so it is consistent across
shared faces, and the two slabs meet along the (full) level-0 triangles.  The
planar order puts ``m`` first, so every tetrahedron of a prism over an
``m``-triangle contains the vertex ``m_0 = (m, 0)``: the closed star of ``m_0``
is the union of the 8 prisms over the 4 ``m``-triangles, i.e. the *convex*
polytope ``Quad x [-delta, delta]``, cut into 24 tetrahedra.

``Pi_k`` moves the single vertex ``m_0`` to ``m'_0`` where

    m' = (1/4, 3/4)   for side bit 0 = "prefix-first"   (path follows e_p first)
    m' = (3/4, 1/4)   for side bit 1 = "target-first"   (path follows e_t first)

(the side bits are ``search_t73_johnson_alpha_sides.KNOWN_BITS``), keeps all other
vertices fixed and is affine on each of the 24 tetrahedra; outside the closed
star it is the identity.  Because ``m'`` lies strictly inside the convex
quadrilateral, every one of the 24 tetrahedra keeps a strictly positive
orientation determinant, i.e. ``m'`` lies in the interior of the kernel of the
3-dimensional star, so the 24 image tetrahedra again tile ``Quad x [-delta,
delta]``, ``Pi_k`` is a PL homeomorphism of ``T^3`` with the explicit inverse
obtained by exchanging ``m_0`` and ``m'_0``, and ``Pi_k`` carries the straight
diagonal ``P -- Q`` onto the bent path ``P -- m' -- Q``.

What is recorded, exactly
-------------------------
* ``linear_part_E``     -- the 3x3 integer matrix, zero translation, its inverse
                           and its determinant.  This is a *global* affine map;
                           no cell decomposition is needed for it.
* ``push_Pi.cells``     -- the 24 closed star tetrahedra of ``Pi_k`` as index
                           quadruples into ``push_Pi.vertices``, each positively
                           oriented in the ambient frame, together with the
                           affine map ``(matrix, translation)`` of ``Pi_k`` on
                           that cell and its Jacobian determinant.
* ``push_Pi.inverse_cells`` -- the 24 image tetrahedra (indices into
                           ``push_Pi.image_vertices``) with the affine maps of
                           ``Pi_k^{-1}``.
* ``psi_cells``         -- the pullback refinement ``E_k^{-1}(tau)`` of the 24
                           cells, with the *composite* affine maps
                           ``x -> M_tau (E_k x) + c_tau`` of ``psi_k``; this is
                           the cell decomposition on which ``psi_k`` itself is
                           affine.  Outside ``E_k^{-1}(support)`` one has
                           ``psi_k = E_k``.
* ``psi_inverse_cells`` -- the same for ``psi_k^{-1} = E_k^{-1} o Pi_k^{-1}``.

Jacobian determinants: ``det E_k = +1`` exactly, and every cell determinant of
``psi_k`` equals the determinant of the corresponding push cell.  Those are
*not* ``+1``: a PL homeomorphism that moves one interior vertex is not volume
preserving cell by cell.  The requirement that is actually checked -- and the
one that matters for the homeomorphism property -- is that every cell
determinant is strictly positive, before and after the move.  The determinants
are recorded verbatim.

Outputs
-------
    geometry/t73_johnson_generators/gen_000.json ... gen_092.json
    geometry/t73_johnson_generators/index.json
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "geometry" / "t73_johnson_generators"

DELTA = Fraction(1, 8)
# 1/196104 = 1/(8 * 24513) is the radius committed by
# scripts/straighten_t73_johnson_relative_ball.py, where 24513 is the largest
# infinity norm of an inverse prefix basis matrix of the 93-move movie.
PROTECTED_RADIUS = Fraction(1, 196104)
MATRIX_A = [[0, 269, 1240], [0, 41, 189], [1, 0, 32]]

# ---------------------------------------------------------------------------
# the planar triangulation of the unit square
# ---------------------------------------------------------------------------

SQUARE_VERTEX_NAMES = ["m", "P", "R1", "Q", "R2", "C00", "C10", "C11", "C01"]
SQUARE_VERTICES = [
    (Fraction(1, 2), Fraction(1, 2)),   # 0  m
    (Fraction(1, 4), Fraction(1, 4)),   # 1  P
    (Fraction(7, 8), Fraction(1, 8)),   # 2  R1
    (Fraction(3, 4), Fraction(3, 4)),   # 3  Q
    (Fraction(1, 8), Fraction(7, 8)),   # 4  R2
    (Fraction(0), Fraction(0)),         # 5  C00
    (Fraction(1), Fraction(0)),         # 6  C10
    (Fraction(1), Fraction(1)),         # 7  C11
    (Fraction(0), Fraction(1)),         # 8  C01
]
SQUARE_TRIANGLES = [
    (0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 1),          # star of m
    (5, 6, 2), (5, 2, 1), (6, 7, 3), (6, 3, 2),          # annulus
    (7, 8, 4), (7, 4, 3), (8, 5, 1), (8, 1, 4),
]
STAR_TRIANGLES = SQUARE_TRIANGLES[:4]
STAR_VERTEX_ORDER = [0, 1, 2, 3, 4]                       # m, P, R1, Q, R2
QUAD_CYCLE = [1, 2, 3, 4]                                 # P, R1, Q, R2 (ccw)
M_PRIME = {
    "prefix-first": (Fraction(1, 4), Fraction(3, 4)),
    "target-first": (Fraction(3, 4), Fraction(1, 4)),
}
LEVELS = [0, 1, -1]                                       # c = 0, +delta, -delta


# ---------------------------------------------------------------------------
# small exact helpers
# ---------------------------------------------------------------------------


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


def enc_point(point: Sequence[Fraction]) -> list[str]:
    return [fs(c) for c in point]


def enc_matrix(matrix) -> list[list[str]]:
    return [[fs(entry) for entry in row] for row in matrix]


def sub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def det3(rows) -> Fraction:
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def matvec(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3))


def matmul(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def inverse3(matrix):
    determinant = det3(matrix)
    if determinant == 0:
        raise AssertionError("singular matrix")
    cof = [
        [
            matrix[(i + 1) % 3][(j + 1) % 3] * matrix[(i + 2) % 3][(j + 2) % 3]
            - matrix[(i + 1) % 3][(j + 2) % 3] * matrix[(i + 2) % 3][(j + 1) % 3]
            for j in range(3)
        ]
        for i in range(3)
    ]
    return [[Fraction(cof[j][i], 1) / determinant for j in range(3)] for i in range(3)]


def tet_determinant(points) -> Fraction:
    """Six times the signed volume of the tetrahedron ``points``."""
    p0, p1, p2, p3 = points
    return det3((sub(p1, p0), sub(p2, p0), sub(p3, p0)))


def affine_from_correspondence(source, target):
    """The unique affine map taking ``source[i]`` to ``target[i]`` (4 points)."""
    basis = [sub(source[1], source[0]), sub(source[2], source[0]), sub(source[3], source[0])]
    image = [sub(target[1], target[0]), sub(target[2], target[0]), sub(target[3], target[0])]
    # columns of B are the basis vectors; matrix = Image * Basis^{-1}
    b_cols = [[basis[j][i] for j in range(3)] for i in range(3)]
    i_cols = [[image[j][i] for j in range(3)] for i in range(3)]
    matrix = matmul(i_cols, inverse3(b_cols))
    translation = tuple(target[0][i] - matvec(matrix, source[0])[i] for i in range(3))
    for src, dst in zip(source, target):
        if tuple(matvec(matrix, src)[i] + translation[i] for i in range(3)) != tuple(dst):
            raise AssertionError("affine reconstruction failed")
    return matrix, translation


def circle_distance(value: Fraction) -> Fraction:
    frac = value - (value.numerator // value.denominator)
    return min(frac, 1 - frac)


def distance_to_K1(point) -> Fraction:
    """Sup-norm distance in ``T^3`` from ``point`` to the standard spine ``K1``.

    ``K1`` is the union of the three coordinate axis circles through ``0``.  The
    sup-norm distance to the ``e_i`` circle is ``max_{j != i} d(x_j, Z)``, so the
    distance to ``K1`` is the *median* of the three circle distances.
    """
    values = sorted(circle_distance(Fraction(c)) for c in point)
    return values[1]


# ---------------------------------------------------------------------------
# the factorisation and the side bits
# ---------------------------------------------------------------------------


def unit_moves() -> list[dict[str, Any]]:
    factor = load("factor_t73_matrix_johnson").generate()
    if factor["matrix_A"] != MATRIX_A:
        raise AssertionError("unexpected matrix A")
    return factor["unit_alpha_moves"]


def side_bits() -> str:
    return load("search_t73_johnson_alpha_sides").KNOWN_BITS


def linear_matrix(target: int, prefix: int, power: int) -> list[list[int]]:
    """``I + s E_{p,t}``: column ``t`` is ``e_t + s e_p``, i.e. ``x -> x + s x_t e_p``."""
    matrix = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
    matrix[prefix][target] += power
    return matrix


def composite_linear(moves: Sequence[dict[str, Any]]) -> list[list[int]]:
    current = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for move in moves:
        local = linear_matrix(move["alpha_target"], move["alpha_prefix"], move["power"])
        current = matmul(local, current)
    return current


# ---------------------------------------------------------------------------
# the planar triangulation, checked
# ---------------------------------------------------------------------------


def planar_checks() -> dict[str, Any]:
    def cross(u, v):
        return u[0] * v[1] - u[1] * v[0]

    areas = []
    for tri in SQUARE_TRIANGLES:
        a, b, c = (SQUARE_VERTICES[i] for i in tri)
        areas.append(cross((b[0] - a[0], b[1] - a[1]), (c[0] - a[0], c[1] - a[1])))
    if any(area <= 0 for area in areas):
        raise AssertionError("planar triangulation has a non-positive triangle")
    total_area = sum(areas, Fraction(0)) / 2
    if total_area != 1:
        raise AssertionError(f"planar triangles do not fill the unit square: {total_area}")

    edge_count: dict[tuple[int, int], int] = {}
    for tri in SQUARE_TRIANGLES:
        for i in range(3):
            edge = tuple(sorted((tri[i], tri[(i + 1) % 3])))
            edge_count[edge] = edge_count.get(edge, 0) + 1
    if set(edge_count.values()) != {1, 2}:
        raise AssertionError("planar triangulation is not a manifold-with-boundary")
    boundary = sorted(edge for edge, count in edge_count.items() if count == 1)
    if len(boundary) != 4:
        raise AssertionError("planar boundary is not the four square sides")

    star = sorted(tri for tri in SQUARE_TRIANGLES if 0 in tri)
    if sorted(STAR_TRIANGLES) != star:
        raise AssertionError("the star of m is not the four recorded triangles")
    star_vertices = sorted({v for tri in star for v in tri})
    if star_vertices != [0, 1, 2, 3, 4]:
        raise AssertionError("the star of m is not the quadrilateral P R1 Q R2")

    quad = [SQUARE_VERTICES[i] for i in QUAD_CYCLE]
    convex = []
    for i in range(4):
        a, b, c = quad[i], quad[(i + 1) % 4], quad[(i + 2) % 4]
        convex.append(cross((b[0] - a[0], b[1] - a[1]), (c[0] - b[0], c[1] - b[1])))
    if any(value <= 0 for value in convex):
        raise AssertionError("P R1 Q R2 is not a convex quadrilateral")
    quad_area = sum(
        quad[i][0] * quad[(i + 1) % 4][1] - quad[(i + 1) % 4][0] * quad[i][1] for i in range(4)
    ) / 2

    inside = {}
    for name, point in M_PRIME.items():
        tests = [
            cross(
                (quad[(i + 1) % 4][0] - quad[i][0], quad[(i + 1) % 4][1] - quad[i][1]),
                (point[0] - quad[i][0], point[1] - quad[i][1]),
            )
            for i in range(4)
        ]
        if any(value <= 0 for value in tests):
            raise AssertionError(f"m' for {name} is not strictly inside the quadrilateral")
        inside[name] = [fs(value) for value in tests]

    return {
        "vertex_names": SQUARE_VERTEX_NAMES,
        "vertices": [[fs(c) for c in point] for point in SQUARE_VERTICES],
        "triangles": [list(tri) for tri in SQUARE_TRIANGLES],
        "triangle_count": len(SQUARE_TRIANGLES),
        "all_triangles_positively_oriented": True,
        "total_area": fs(total_area),
        "edge_multiplicities": sorted(set(edge_count.values())),
        "boundary_edge_count": len(boundary),
        "star_of_m": [list(tri) for tri in STAR_TRIANGLES],
        "star_of_m_is_quadrilateral_P_R1_Q_R2": True,
        "quadrilateral_is_convex": True,
        "quadrilateral_area": fs(quad_area),
        "m_prime_strictly_inside_quadrilateral": inside,
    }


# ---------------------------------------------------------------------------
# the 3-dimensional star of m_0
# ---------------------------------------------------------------------------


def frame_axes(target: int, prefix: int, third: int, power: int):
    e_a = [0, 0, 0]
    e_a[target] = 1
    e_b = [0, 0, 0]
    e_b[prefix] = power
    e_c = [0, 0, 0]
    e_c[third] = 1
    return e_a, e_b, e_c


def to_ambient(a: Fraction, b: Fraction, c: Fraction, target: int, prefix: int,
               third: int, power: int):
    point = [Fraction(0), Fraction(0), Fraction(0)]
    point[target] = Fraction(a)
    point[prefix] = power * Fraction(b)
    point[third] = Fraction(c)
    return tuple(point)


def star_complex(target: int, prefix: int, third: int, power: int, side: str):
    """The 15 vertices, 24 tetrahedra and the moved vertex of the star of m_0."""
    level_names = {0: "0", 1: "up", -1: "dn"}
    index_of: dict[tuple[int, int], int] = {}
    vertices: list[tuple[Fraction, Fraction, Fraction]] = []
    labels: list[str] = []
    for level in LEVELS:
        for planar in STAR_VERTEX_ORDER:
            a, b = SQUARE_VERTICES[planar]
            c = DELTA * level
            index_of[(planar, level)] = len(vertices)
            vertices.append(to_ambient(a, b, c, target, prefix, third, power))
            labels.append(SQUARE_VERTEX_NAMES[planar] + "_" + level_names[level])

    moved_index = index_of[(0, 0)]
    if moved_index != 0:
        raise AssertionError("m_0 is expected to be vertex 0")

    tetrahedra: list[tuple[int, int, int, int]] = []
    for tri in STAR_TRIANGLES:
        a, b, c = sorted(tri)
        if a != 0:
            raise AssertionError("m must be the smallest planar index of its triangles")
        for outer in (1, -1):
            for cell in (
                (index_of[(a, 0)], index_of[(b, 0)], index_of[(c, 0)], index_of[(c, outer)]),
                (index_of[(a, 0)], index_of[(b, 0)], index_of[(b, outer)], index_of[(c, outer)]),
                (index_of[(a, 0)], index_of[(a, outer)], index_of[(b, outer)], index_of[(c, outer)]),
            ):
                if moved_index not in cell:
                    raise AssertionError("a prism tetrahedron misses m_0")
                points = [vertices[i] for i in cell]
                if tet_determinant(points) < 0:
                    cell = (cell[0], cell[1], cell[3], cell[2])
                if tet_determinant([vertices[i] for i in cell]) <= 0:
                    raise AssertionError("degenerate star tetrahedron")
                tetrahedra.append(cell)
    if len(tetrahedra) != 24:
        raise AssertionError(f"expected 24 star tetrahedra, got {len(tetrahedra)}")

    a_prime, b_prime = M_PRIME[side]
    image_vertices = list(vertices)
    image_vertices[moved_index] = to_ambient(
        a_prime, b_prime, Fraction(0), target, prefix, third, power
    )
    return {
        "vertices": vertices,
        "labels": labels,
        "tetrahedra": tetrahedra,
        "moved_index": moved_index,
        "image_vertices": image_vertices,
        "index_of": index_of,
    }


def support_halfspaces(target: int, prefix: int, third: int, power: int):
    """Six half-spaces ``n . x <= bound`` cutting out ``Quad x [-delta, delta]``."""
    planes = []

    def add(n_a, n_b, n_c, bound):
        normal = [Fraction(0), Fraction(0), Fraction(0)]
        normal[target] += Fraction(n_a)
        normal[prefix] += Fraction(n_b) * power
        normal[third] += Fraction(n_c)
        planes.append((tuple(normal), Fraction(bound)))

    add(-1, -5, 0, Fraction(-3, 2))     # a + 5b >= 3/2   (edge P R1)
    add(5, 1, 0, Fraction(9, 2))        # 5a + b <= 9/2   (edge R1 Q)
    add(1, 5, 0, Fraction(9, 2))        # a + 5b <= 9/2   (edge Q R2)
    add(-5, -1, 0, Fraction(-3, 2))     # 5a + b >= 3/2   (edge R2 P)
    add(0, 0, 1, DELTA)
    add(0, 0, -1, DELTA)
    return planes


def boundary_triangles(tetrahedra):
    counts: dict[tuple[int, int, int], int] = {}
    for cell in tetrahedra:
        for skip in range(4):
            face = tuple(sorted(cell[:skip] + cell[skip + 1:]))
            counts[face] = counts.get(face, 0) + 1
    if set(counts.values()) - {1, 2}:
        raise AssertionError("star faces are not shared by at most two tetrahedra")
    return sorted(face for face, count in counts.items() if count == 1), counts


# ---------------------------------------------------------------------------
# pushing a polyline through a single push Pi
# ---------------------------------------------------------------------------


class Push:
    """The push ``Pi`` of one generator, ready to be applied to polylines."""

    def __init__(self, vertices, tetrahedra, image_vertices, halfspaces):
        self.vertices = [tuple(Fraction(c) for c in v) for v in vertices]
        self.image_vertices = [tuple(Fraction(c) for c in v) for v in image_vertices]
        self.tetrahedra = [tuple(cell) for cell in tetrahedra]
        self.halfspaces = [(tuple(Fraction(c) for c in n), Fraction(b)) for n, b in halfspaces]
        self.maps = []
        for cell in self.tetrahedra:
            source = [self.vertices[i] for i in cell]
            image = [self.image_vertices[i] for i in cell]
            self.maps.append(affine_from_correspondence(source, image))
        lows = [min(v[i] for v in self.vertices) for i in range(3)]
        highs = [max(v[i] for v in self.vertices) for i in range(3)]
        self.low = lows
        self.high = highs
        planes: dict[tuple, None] = {}
        for cell in self.tetrahedra:
            for skip in range(4):
                face = [self.vertices[i] for i in cell[:skip] + cell[skip + 1:]]
                normal = cross3(sub(face[1], face[0]), sub(face[2], face[0]))
                if normal == (0, 0, 0):
                    raise AssertionError("degenerate face")
                offset = sum(normal[i] * face[0][i] for i in range(3))
                scale = normalise_scale(normal, offset)
                planes[scale] = None
        self.planes = list(planes)

    def translate_candidates(self, cell_origin):
        """Lattice translates whose support box can meet the unit cell ``cell_origin``."""
        ranges = []
        for i in range(3):
            first = ceil_fraction(Fraction(cell_origin[i]) - self.high[i])
            last = floor_fraction(Fraction(cell_origin[i]) + 1 - self.low[i])
            ranges.append(list(range(first, last + 1)))
        return [(x, y, z) for x in ranges[0] for y in ranges[1] for z in ranges[2]]

    def contains(self, point, shift):
        for normal, bound in self.halfspaces:
            if sum(normal[i] * (point[i] - shift[i]) for i in range(3)) > bound:
                return False
        return True

    def locate(self, point, shift):
        """Index of a star tetrahedron (translated by ``shift``) containing ``point``."""
        local = tuple(point[i] - shift[i] for i in range(3))
        for index, cell in enumerate(self.tetrahedra):
            points = [self.vertices[i] for i in cell]
            ok = True
            for skip in range(4):
                sub_points = points[:skip] + [local] + points[skip + 1:]
                if tet_determinant(sub_points) < 0:
                    ok = False
                    break
            if ok:
                return index
        return None

    def apply_point(self, point, shift, index):
        matrix, translation = self.maps[index]
        local = tuple(point[i] - shift[i] for i in range(3))
        moved = matvec(matrix, local)
        return tuple(moved[i] + translation[i] + shift[i] for i in range(3))

    def push_point(self, point):
        point = tuple(Fraction(c) for c in point)
        origin = tuple(floor_fraction(c) for c in point)
        for shift in self.translate_candidates(origin):
            if not self.contains(point, shift):
                continue
            index = self.locate(point, shift)
            if index is None:
                raise AssertionError("point inside the support but outside every cell")
            return self.apply_point(point, shift, index)
        return point


def cross3(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def normalise_scale(normal, offset):
    values = [Fraction(c) for c in normal] + [Fraction(offset)]
    denominator = 1
    for value in values:
        denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    integers = [int(value * denominator) for value in values]
    common = 0
    for value in integers:
        common = gcd(common, abs(value))
    if common:
        integers = [value // common for value in integers]
    for value in integers:
        if value != 0:
            if value < 0:
                integers = [-v for v in integers]
            break
    return tuple(integers)


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def integer_break_parameters(start, direction):
    """Parameters in (0,1) at which ``start + t*direction`` crosses an integer plane."""
    breaks = set()
    for i in range(3):
        d = direction[i]
        if d == 0:
            continue
        lo, hi = (start[i], start[i] + d) if d > 0 else (start[i] + d, start[i])
        first = floor_fraction(lo) + 1
        last = ceil_fraction(hi) - 1
        for n in range(first, last + 1):
            t = (Fraction(n) - start[i]) / d
            if 0 < t < 1:
                breaks.add(t)
    return sorted(breaks)


def clip_parameters(start, direction, halfspaces, shift, lo, hi):
    """Sub-interval of ``[lo,hi]`` on which the segment lies in the shifted region."""
    low, high = lo, hi
    for normal, bound in halfspaces:
        value = sum(normal[i] * (start[i] - shift[i]) for i in range(3))
        slope = sum(normal[i] * direction[i] for i in range(3))
        if slope == 0:
            if value > bound:
                return None
            continue
        limit = (bound - value) / slope
        if slope > 0:
            high = min(high, limit)
        else:
            low = max(low, limit)
        if low >= high:
            return None
    return (low, high)


def support_windows(start, direction, push: Push):
    """Disjoint parameter windows ``(shift, lo, hi)`` where the segment meets a
    lattice translate of the support."""
    intervals: list[tuple[tuple[int, int, int], Fraction, Fraction]] = []
    cuts = [Fraction(0)] + integer_break_parameters(start, direction) + [Fraction(1)]
    for lo, hi in zip(cuts, cuts[1:]):
        if lo >= hi:
            continue
        mid = (lo + hi) / 2
        origin = tuple(floor_fraction(start[i] + mid * direction[i]) for i in range(3))
        for shift in push.translate_candidates(origin):
            window = clip_parameters(start, direction, push.halfspaces, shift, lo, hi)
            if window is None:
                continue
            intervals.append((shift, window[0], window[1]))
    merged: list[list] = []
    for shift, lo, hi in sorted(intervals, key=lambda item: (item[1], item[2])):
        if merged and merged[-1][0] == shift and merged[-1][2] == lo:
            merged[-1][2] = hi
        else:
            merged.append([shift, lo, hi])
    return merged


def push_segment(start, end, push: Push):
    """Image of ``[start,end]``: the breakpoints after ``start``, ending at the
    image of ``end``.  ``Pi`` is the identity on the boundary of every support
    translate, so the window endpoints need no special treatment."""
    direction = sub(end, start)
    if direction == (0, 0, 0):
        return [push.push_point(end)]
    images: dict[Fraction, tuple] = {}
    for shift, lo, hi in support_windows(start, direction, push):
        params = {lo, hi}
        for plane in push.planes:
            slope = sum(plane[i] * direction[i] for i in range(3))
            if slope == 0:
                continue
            value = sum(plane[i] * (start[i] - shift[i]) for i in range(3))
            t = (Fraction(plane[3]) - value) / slope
            if lo < t < hi:
                params.add(t)
        ordered = sorted(params)
        for left, right in zip(ordered, ordered[1:]):
            mid = (left + right) / 2
            point = tuple(start[i] + mid * direction[i] for i in range(3))
            index = push.locate(point, shift)
            if index is None:
                raise AssertionError("segment inside the support but outside every cell")
            for t in (left, right):
                raw = tuple(start[i] + t * direction[i] for i in range(3))
                moved = push.apply_point(raw, shift, index)
                if t in images and images[t] != moved:
                    raise AssertionError("the push is discontinuous across a cell wall")
                images[t] = moved
    out = []
    for t in sorted(set(images) | {Fraction(1)}):
        if t == 0:
            continue
        out.append(images.get(t, tuple(start[i] + t * direction[i] for i in range(3))))
    return out


def push_polyline(points, push: Push):
    out = [push.push_point(points[0])]
    for start, end in zip(points, points[1:]):
        out.extend(push_segment(tuple(Fraction(c) for c in start),
                                tuple(Fraction(c) for c in end), push))
    return merge_collinear(out)


def merge_collinear(points):
    if len(points) < 3:
        return list(points)
    out = [points[0]]
    for point in points[1:]:
        if out[-1] == point:
            continue
        out.append(point)
    trimmed = [out[0]]
    for index in range(1, len(out) - 1):
        previous, current, following = trimmed[-1], out[index], out[index + 1]
        direction = sub(following, previous)
        offset = sub(current, previous)
        if direction != (0, 0, 0) and cross3(offset, direction) == (0, 0, 0):
            pivot = next(i for i in range(3) if direction[i] != 0)
            if 0 < offset[pivot] / direction[pivot] < 1:
                continue
        trimmed.append(current)
    trimmed.append(out[-1])
    return trimmed


# ---------------------------------------------------------------------------
# one generator
# ---------------------------------------------------------------------------


def spine_loops():
    zero = (Fraction(0), Fraction(0), Fraction(0))
    loops = []
    for axis in range(3):
        end = [Fraction(0), Fraction(0), Fraction(0)]
        end[axis] = Fraction(1)
        loops.append([zero, tuple(end)])
    return loops


def build_generator(index: int, move: dict[str, Any], bit: int) -> dict[str, Any]:
    target = move["alpha_target"]
    prefix = move["alpha_prefix"]
    power = move["power"]
    if target == prefix:
        raise AssertionError("a transvection needs two distinct axes")
    third = ({0, 1, 2} - {target, prefix}).pop()
    side = "prefix-first" if bit == 0 else "target-first"

    matrix_e = linear_matrix(target, prefix, power)
    inverse_e = linear_matrix(target, prefix, -power)
    if matmul(matrix_e, inverse_e) != [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
        raise AssertionError("E_k inverse is wrong")
    if det3([[Fraction(v) for v in row] for row in matrix_e]) != 1:
        raise AssertionError("E_k is not unimodular")

    star = star_complex(target, prefix, third, power, side)
    vertices = star["vertices"]
    image_vertices = star["image_vertices"]
    tetrahedra = star["tetrahedra"]
    moved = star["moved_index"]
    halfspaces = support_halfspaces(target, prefix, third, power)

    determinants = [tet_determinant([vertices[i] for i in cell]) for cell in tetrahedra]
    image_determinants = [
        tet_determinant([image_vertices[i] for i in cell]) for cell in tetrahedra
    ]
    if any(value <= 0 for value in determinants):
        raise AssertionError("a star cell is degenerate before the move")
    if any(value <= 0 for value in image_determinants):
        raise AssertionError("m' is not in the interior of the kernel of the star")
    volume_before = sum(determinants, Fraction(0)) / 6
    volume_after = sum(image_determinants, Fraction(0)) / 6
    if volume_before != volume_after:
        raise AssertionError("the push does not preserve the support volume")
    quad_area = Fraction(3, 8)
    if volume_before != quad_area * 2 * DELTA:
        raise AssertionError("star cells do not fill Quad x [-delta,delta]")

    faces, face_counts = boundary_triangles(tetrahedra)
    if len(faces) != 24:
        raise AssertionError(f"expected 24 boundary triangles, got {len(faces)}")
    if any(moved in face for face in faces):
        raise AssertionError("the moved vertex lies on the boundary of its star")

    cell_maps = []
    inverse_cell_maps = []
    psi_cells = []
    psi_inverse_cells = []
    inverse_e_frac = [[Fraction(v) for v in row] for row in inverse_e]
    matrix_e_frac = [[Fraction(v) for v in row] for row in matrix_e]
    for cell, determinant in zip(tetrahedra, determinants):
        source = [vertices[i] for i in cell]
        image = [image_vertices[i] for i in cell]
        matrix, translation = affine_from_correspondence(source, image)
        jac = det3(matrix)
        if jac <= 0:
            raise AssertionError("non-positive Jacobian on a push cell")
        cell_maps.append((matrix, translation, jac))
        inverse_matrix = inverse3(matrix)
        inverse_translation = tuple(-matvec(inverse_matrix, translation)[i] for i in range(3))
        inverse_cell_maps.append((inverse_matrix, inverse_translation, det3(inverse_matrix)))
        composite = matmul(matrix, matrix_e_frac)
        psi_cells.append((composite, translation, det3(composite)))
        inverse_composite = matmul(inverse_e_frac, inverse_matrix)
        inverse_composite_translation = tuple(
            matvec(inverse_e_frac, inverse_translation)[i] for i in range(3)
        )
        psi_inverse_cells.append(
            (inverse_composite, inverse_composite_translation, det3(inverse_composite))
        )

    psi_domain_vertices = [matvec(inverse_e_frac, v) for v in vertices]

    lows = [min(v[i] for v in vertices) for i in range(3)]
    highs = [max(v[i] for v in vertices) for i in range(3)]
    sides = [highs[i] - lows[i] for i in range(3)]
    if any(side_length >= 1 for side_length in sides):
        raise AssertionError("the support does not embed in the torus")
    straddle = [bool(lows[i] < floor_fraction(lows[i]) + 1 < highs[i]) for i in range(3)]

    # protected ball and axis clearance, derived from the coordinate ranges only
    clearances = []
    for i in range(3):
        base = floor_fraction(lows[i])
        if not (base < lows[i] and highs[i] < base + 1):
            clearances.append(Fraction(0))
        else:
            clearances.append(min(lows[i] - base, base + 1 - highs[i]))
    ball_margin = max(clearances)
    if ball_margin <= PROTECTED_RADIUS:
        raise AssertionError("the support meets the protected ball")
    axis_clear = {}
    for axis in range(3):
        others = [i for i in range(3) if i != axis]
        value = max(clearances[i] for i in others)
        if value <= 0:
            raise AssertionError(f"the support meets the axis circle e_{axis}")
        axis_clear[f"e_{axis}"] = fs(value)

    push = Push(vertices, tetrahedra, image_vertices, halfspaces)
    images = []
    for loop in spine_loops():
        linear_image = [matvec(matrix_e_frac, point) for point in loop]
        images.append(push_polyline(linear_image, push))
    induced = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for column, polyline in enumerate(images):
        shift = sub(polyline[-1], polyline[0])
        for row in range(3):
            if shift[row].denominator != 1:
                raise AssertionError("a spine image is not a closed loop")
            induced[row][column] = int(shift[row])
    if induced != matrix_e:
        raise AssertionError("the spine images do not induce E_k on H_1")

    diagonal_image = images[target]
    expected = [
        (Fraction(0), Fraction(0), Fraction(0)),
        to_ambient(Fraction(1, 4), Fraction(1, 4), Fraction(0), target, prefix, third, power),
        to_ambient(M_PRIME[side][0], M_PRIME[side][1], Fraction(0), target, prefix, third, power),
        to_ambient(Fraction(3, 4), Fraction(3, 4), Fraction(0), target, prefix, third, power),
        to_ambient(Fraction(1), Fraction(1), Fraction(0), target, prefix, third, power),
    ]
    if diagonal_image != expected:
        raise AssertionError(f"unexpected bent path at generator {index}")

    max_distance = max(
        polyline_max_distance_to_K1(polyline) for polyline in images
    )

    document: dict[str, Any] = {
        "schema": "t73_johnson_pl_generator/v1",
        "generator": "scripts/build_t73_johnson_pl_generators.py",
        "index": index,
        "operation_index": move["operation_index"],
        "repetition": move["repetition"],
        "alpha_target": target,
        "alpha_prefix": prefix,
        "third_axis": third,
        "power": power,
        "side_bit": bit,
        "side": side,
        "conventions": {
            "torus": "T^3 = R^3/Z^3, exact Fraction coordinates",
            "reduction": "red(x)_i = x_i - floor(x_i) into [0,1)^3; lifts differing by Z^3 agree",
            "cells": "all cells are stored as lifts in Q^3; the support box has all side lengths < 1 so red is injective on it",
            "composition": "psi_k = Pi_k o E_k (E_k applied first); psi_A = psi_92 o ... o psi_0",
            "linear_rule": "E_k: x -> x + s * x_t * e_p, matrix I + s E_{p,t}",
            "frame": "(a,b,c) coordinates along (e_t, s e_p, e_r): x_t = a, x_p = s b, x_r = c",
        },
        "linear_part_E": {
            "matrix": matrix_e,
            "translation": [0, 0, 0],
            "determinant": 1,
            "inverse_matrix": inverse_e,
            "lattice_preserving": True,
            "affine_on_every_cell": True,
        },
        "frame": {
            "e_a": frame_axes(target, prefix, third, power)[0],
            "e_b": frame_axes(target, prefix, third, power)[1],
            "e_c": frame_axes(target, prefix, third, power)[2],
        },
        "square_triangulation": planar_checks(),
        "push_Pi": {
            "delta": fs(DELTA),
            "moved_vertex": moved,
            "moved_from": enc_point(vertices[moved]),
            "moved_to": enc_point(image_vertices[moved]),
            "m_prime_planar": [fs(M_PRIME[side][0]), fs(M_PRIME[side][1])],
            "vertex_labels": star["labels"],
            "vertices": [enc_point(v) for v in vertices],
            "image_vertices": [enc_point(v) for v in image_vertices],
            "tetrahedra": [list(cell) for cell in tetrahedra],
            "cell_determinants": [fs(value) for value in determinants],
            "image_cell_determinants": [fs(value) for value in image_determinants],
            "cell_maps": [
                {"matrix": enc_matrix(m), "translation": enc_point(t), "jacobian": fs(j)}
                for m, t, j in cell_maps
            ],
            "inverse_cell_maps": [
                {"matrix": enc_matrix(m), "translation": enc_point(t), "jacobian": fs(j)}
                for m, t, j in inverse_cell_maps
            ],
            "boundary_triangles": [list(face) for face in faces],
            "face_multiplicities": sorted(set(face_counts.values())),
            "support_volume": fs(volume_before),
            "image_support_volume": fs(volume_after),
            "support_halfspaces": [
                {"normal": enc_point(normal), "bound": fs(bound)} for normal, bound in halfspaces
            ],
            "support_box": {
                "low": enc_point(lows),
                "high": enc_point(highs),
                "side_lengths": [fs(value) for value in sides],
                "embeds_in_torus": True,
                "straddles_fundamental_domain_face": straddle,
            },
            "protected_ball": {
                "radius": fs(PROTECTED_RADIUS),
                "coordinate_clearances": [fs(value) for value in clearances],
                "min_sup_norm_distance_to_lattice": fs(ball_margin),
                "margin": fs(ball_margin - PROTECTED_RADIUS),
                "fixed_pointwise_by_push": True,
            },
            "axis_clearance": axis_clear,
            "support_meets_axes": False,
        },
        "psi_cells": {
            "domain_vertices": [enc_point(v) for v in psi_domain_vertices],
            "image_vertices": [enc_point(v) for v in image_vertices],
            "tetrahedra": [list(cell) for cell in tetrahedra],
            "cell_maps": [
                {"matrix": enc_matrix(m), "translation": enc_point(t), "jacobian": fs(j)}
                for m, t, j in psi_cells
            ],
            "inverse_cell_maps": [
                {"matrix": enc_matrix(m), "translation": enc_point(t), "jacobian": fs(j)}
                for m, t, j in psi_inverse_cells
            ],
            "note": (
                "domain cells are E_k^{-1}(closed star); outside their union psi_k = E_k."
            ),
        },
        "induced_H1_matrix": induced,
        "spine_images": {
            f"C_{axis}": [enc_point(point) for point in images[axis]] for axis in range(3)
        },
        "spine_image_vertex_counts": [len(polyline) for polyline in images],
        "max_sup_norm_distance_to_K1": fs(max_distance),
        "jacobian_note": (
            "det E_k = +1 exactly.  Cell Jacobians of Pi_k and of psi_k are positive but "
            "not +1: a PL homeomorphism that moves one interior vertex is not "
            "cellwise volume preserving.  Positivity before and after the move is the "
            "property that is checked."
        ),
    }
    document["generator_sha256"] = canonical_sha(document)
    return document


def polyline_max_distance_to_K1(points) -> Fraction:
    best = Fraction(0)
    for start, end in zip(points, points[1:]):
        best = max(best, segment_max_distance_to_K1(start, end))
    if len(points) == 1:
        best = max(best, distance_to_K1(points[0]))
    return best


def segment_max_distance_to_K1(start, end, cap: Fraction | None = None) -> Fraction:
    direction = sub(end, start)
    if direction == (0, 0, 0):
        return distance_to_K1(start)
    params = {Fraction(0), Fraction(1)}
    for i in range(3):
        d = direction[i]
        if d == 0:
            continue
        lo, hi = (start[i], end[i]) if d > 0 else (end[i], start[i])
        first = ceil_fraction(lo * 2)
        last = floor_fraction(hi * 2)
        for n in range(first, last + 1):
            t = (Fraction(n, 2) - start[i]) / d
            if 0 < t < 1:
                params.add(t)
    ordered = sorted(params)
    best = Fraction(0)
    for left, right in zip(ordered, ordered[1:]):
        pieces = {left, right}
        mid = (left + right) / 2
        base = tuple(start[i] + mid * direction[i] for i in range(3))
        offsets = [floor_fraction(base[i]) for i in range(3)]
        signs = []
        for i in range(3):
            centre = Fraction(offsets[i]) + Fraction(1, 2)
            signs.append(1 if base[i] <= centre else -1)
        # on (left,right) each circle distance is affine: sign*(x_i - offset) or
        # sign*(offset+1 - x_i); find pairwise crossings
        lines = []
        for i in range(3):
            if signs[i] > 0:
                lines.append((direction[i], start[i] - offsets[i]))
            else:
                lines.append((-direction[i], Fraction(offsets[i]) + 1 - start[i]))
        for i in range(3):
            for j in range(i + 1, 3):
                slope = lines[i][0] - lines[j][0]
                if slope == 0:
                    continue
                t = (lines[j][1] - lines[i][1]) / slope
                if left < t < right:
                    pieces.add(t)
        for t in sorted(pieces):
            point = tuple(start[i] + t * direction[i] for i in range(3))
            best = max(best, distance_to_K1(point))
            if cap is not None and best >= cap:
                return best
    return best


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def build_all() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    moves = unit_moves()
    bits = side_bits()
    if len(bits) != len(moves):
        raise AssertionError("side bit count does not match the move count")
    composite = composite_linear(moves)
    if composite != MATRIX_A:
        raise AssertionError("E_92 ... E_0 does not equal A")
    documents = [build_generator(k, move, int(bits[k])) for k, move in enumerate(moves)]
    index = {
        "schema": "t73_johnson_pl_generator_index/v1",
        "generator": "scripts/build_t73_johnson_pl_generators.py",
        "generator_count": len(documents),
        "matrix_A": MATRIX_A,
        "composition_order": "psi_A = psi_92 o ... o psi_0; matrices A = E_92 * ... * E_0",
        "composite_linear_matrix": composite,
        "composite_matrix_equals_A": composite == MATRIX_A,
        "side_bits": bits,
        "delta": fs(DELTA),
        "protected_ball_radius": fs(PROTECTED_RADIUS),
        "section_straightening": {
            "file": "gen_093_section_straightening.json",
            "generator": "scripts/build_t73_section_straightening.py",
            "role": (
                "theta, a PL homeomorphism of T^3 supported in the protected cube, equal "
                "to A^{-1} near 0, so that Psi = psi_A o theta is the identity on a ball; "
                "it is not one of the 93 Johnson generators and has its own schema"
            ),
        },
        "files": [
            {
                "file": f"gen_{doc['index']:03d}.json",
                "index": doc["index"],
                "alpha_target": doc["alpha_target"],
                "alpha_prefix": doc["alpha_prefix"],
                "power": doc["power"],
                "side": doc["side"],
                "cells": len(doc["push_Pi"]["tetrahedra"]),
                "generator_sha256": doc["generator_sha256"],
            }
            for doc in documents
        ],
    }
    index["index_sha256"] = canonical_sha(index)
    return documents, index


def dump(document: dict[str, Any]) -> str:
    return json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    documents, index = build_all()
    ok = index["composite_matrix_equals_A"]

    if args.write:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        total = 0
        for document in documents:
            path = OUTPUT_DIR / f"gen_{document['index']:03d}.json"
            path.write_text(dump(document), encoding="utf-8")
            total += path.stat().st_size
        path = OUTPUT_DIR / "index.json"
        path.write_text(dump(index), encoding="utf-8")
        total += path.stat().st_size
        print(f"WROTE={OUTPUT_DIR.relative_to(ROOT)} FILES={len(documents) + 1} BYTES={total}")

    if args.check:
        mismatched = []
        for document in documents:
            path = OUTPUT_DIR / f"gen_{document['index']:03d}.json"
            if not path.exists():
                mismatched.append(f"{path.name}:MISSING")
            elif json.loads(path.read_text(encoding="utf-8")) != document:
                mismatched.append(f"{path.name}:DIFFERS")
        path = OUTPUT_DIR / "index.json"
        if not path.exists():
            mismatched.append("index.json:MISSING")
        elif json.loads(path.read_text(encoding="utf-8")) != index:
            mismatched.append("index.json:DIFFERS")
        print(f"COMMITTED_GENERATORS={'PASS' if not mismatched else 'FAIL ' + ','.join(mismatched)}")
        ok = ok and not mismatched

    print(f"JOHNSON_GENERATORS={len(documents)}")
    print(f"INDUCED_MATRIX_EQUALS_A={'PASS' if index['composite_matrix_equals_A'] else 'FAIL'}")
    print(f"PROTECTED_BALL_RADIUS={index['protected_ball_radius']}")
    print(f"INDEX_SHA256={index['index_sha256']}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
