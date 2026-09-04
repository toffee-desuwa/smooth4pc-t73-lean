#!/usr/bin/env python3
"""Independent checker for the explicit PL Johnson generators of COMMIT 4a.

Everything below is recomputed from the *raw* geometric data of a generator
file: the vertex coordinates, the tetrahedron index tuples, the affine matrices
and translations, the linear matrix of ``E_k`` and the recorded spine polygons.
No boolean status field of the document is read (in particular none of
``*_status``, ``all_*``, ``*_is_*``, ``fixed_pointwise_by_push``,
``embeds_in_torus``); recorded *numbers* such as cell determinants are
recomputed and compared, never trusted.

Checks performed for every generator
------------------------------------
1.  ``E_k`` is integral, unimodular with determinant ``+1``, equals
    ``I + s E_{p,t}`` for the recorded ``(target, prefix, power)``, and the
    recorded inverse really inverts it.
2.  The 24 push cells are non-degenerate and positively oriented; every internal
    triangle is shared by exactly two cells, the remaining 24 triangles form a
    closed connected surface with Euler characteristic 2, and the moved vertex
    lies in all 24 cells and in none of the boundary triangles.  Hence the union
    is the cone from the moved vertex over that surface, i.e. the support.
3.  The affine map recorded on each cell takes the cell's four vertices to the
    recorded image vertices; on every internal face the two adjacent cells agree
    vertex by vertex (continuity); every Jacobian determinant is strictly
    positive.
4.  Every vertex except the moved one is fixed, so the map is the identity on
    the boundary surface and the extension by the identity outside the support
    is continuous.
5.  The image cells are positively oriented, have the *same* boundary triangles,
    and the same total volume: the images tile the same support.
6.  The inverse cells are the image cells, and inverse-matrix times matrix is the
    identity (and the composition fixes every cell vertex).
7.  ``psi`` cells: the domain vertices are ``E_k^{-1}`` of the push vertices, the
    recorded map is ``push_map o E_k``, its determinant is positive, and the
    recorded ``psi^{-1}`` map inverts it.
8.  Protected ball: some ambient coordinate of the whole support (and of every
    single cell) is confined to an open unit interval ``(n, n+1)``; therefore the
    sup-norm distance from the support to ``Z^3`` is at least the corresponding
    clearance, which is compared with ``1/196104``.  Since the support of the
    push is the only place where ``Pi_k`` is not the identity, ``Pi_k`` fixes the
    protected ball pointwise.
9.  The three spine loops of ``K1`` are pushed through ``E_k`` and then through
    the raw cell data by an implementation written independently of the builder
    (brute-force clipping against each translated tetrahedron), and the result is
    compared with the recorded polygons; the endpoint displacements are compared
    with the recorded ``induced_H1_matrix`` and with the matrix of ``E_k``.
10. The planar square triangulation is re-verified from its raw vertices and
    triangles (positive orientation, total area 1, star of ``m`` = the convex
    quadrilateral ``P R1 Q R2``), and the moved vertex's planar image is checked
    against the side bit.

Honest note on ``E_k`` and the protected ball
---------------------------------------------
``Pi_k`` fixes the protected ball ``B(0, 1/196104)`` pointwise.  ``E_k`` does
*not*: it is a linear transvection, so the only points it fixes are the plane
``x_t = 0``.  ``psi_k = Pi_k o E_k`` therefore fixes ``0`` but no ball around it,
and the composite ``psi_A`` is affine with linear part ``A != I`` near ``0``.
This checker reports both facts separately and never conflates them.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
GENERATOR_DIR = ROOT / "geometry" / "t73_johnson_generators"
PROTECTED_RADIUS = Fraction(1, 196104)
EXPECTED_M_PRIME = {
    "prefix-first": (Fraction(1, 4), Fraction(3, 4)),
    "target-first": (Fraction(3, 4), Fraction(1, 4)),
}


# ---------------------------------------------------------------------------
# exact arithmetic helpers (written independently of the builder)
# ---------------------------------------------------------------------------


def F(value) -> Fraction:
    return Fraction(value)


def point(raw) -> tuple[Fraction, Fraction, Fraction]:
    return (F(raw[0]), F(raw[1]), F(raw[2]))


def matrix(raw):
    return [[F(entry) for entry in row] for row in raw]


def diff(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def det(rows) -> Fraction:
    (a0, a1, a2), (b0, b1, b2), (c0, c1, c2) = rows
    return a0 * (b1 * c2 - b2 * c1) - a1 * (b0 * c2 - b2 * c0) + a2 * (b0 * c1 - b1 * c0)


def volume6(cell) -> Fraction:
    return det((diff(cell[1], cell[0]), diff(cell[2], cell[0]), diff(cell[3], cell[0])))


def apply_affine(mat, translation, vector):
    return tuple(
        sum(mat[i][j] * vector[j] for j in range(3)) + translation[i] for i in range(3)
    )


def mat_mul(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def identity3():
    return [[F(1 if i == j else 0) for j in range(3)] for i in range(3)]


def floor_div(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_div(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


class Failure(Exception):
    pass


def require(condition, message):
    if not condition:
        raise Failure(message)


# ---------------------------------------------------------------------------
# an independent push evaluator: brute force over translated tetrahedra
# ---------------------------------------------------------------------------


class RawPush:
    def __init__(self, vertices, tetrahedra, maps):
        self.vertices = vertices
        self.tetrahedra = tetrahedra
        self.maps = maps
        self.low = [min(v[i] for v in vertices) for i in range(3)]
        self.high = [max(v[i] for v in vertices) for i in range(3)]
        self.faces = []
        for cell in tetrahedra:
            planes = []
            corners = [vertices[i] for i in cell]
            for skip in range(4):
                others = corners[:skip] + corners[skip + 1:]
                normal = cross(diff(others[1], others[0]), diff(others[2], others[0]))
                offset = sum(normal[i] * others[0][i] for i in range(3))
                inside = sum(normal[i] * corners[skip][i] for i in range(3))
                if inside < offset:
                    normal = tuple(-c for c in normal)
                    offset = -offset
                planes.append((normal, offset))
            self.faces.append(planes)

    def shifts_for(self, low, high):
        ranges = []
        for i in range(3):
            first = ceil_div(low[i] - self.high[i])
            last = floor_div(high[i] - self.low[i])
            ranges.append(range(first, last + 1))
        return [(x, y, z) for x in ranges[0] for y in ranges[1] for z in ranges[2]]

    def image_of_point(self, p):
        shifts = self.shifts_for(list(p), list(p))
        for shift in shifts:
            local = tuple(p[i] - shift[i] for i in range(3))
            for index, planes in enumerate(self.faces):
                if all(
                    sum(n[i] * local[i] for i in range(3)) >= c for n, c in planes
                ):
                    mat, translation = self.maps[index]
                    moved = apply_affine(mat, translation, local)
                    return tuple(moved[i] + shift[i] for i in range(3))
        return tuple(p)

    def image_of_segment(self, start, end):
        """Breakpoints strictly after ``start``, ending with the image of ``end``."""
        direction = diff(end, start)
        if direction == (0, 0, 0):
            return [self.image_of_point(end)]
        low = [min(start[i], end[i]) for i in range(3)]
        high = [max(start[i], end[i]) for i in range(3)]
        samples: dict[Fraction, tuple] = {}
        for shift in self.shifts_for(low, high):
            for index, planes in enumerate(self.faces):
                lo, hi = Fraction(0), Fraction(1)
                empty = False
                for normal, offset in planes:
                    value = sum(normal[i] * (start[i] - shift[i]) for i in range(3)) - offset
                    slope = sum(normal[i] * direction[i] for i in range(3))
                    if slope == 0:
                        if value < 0:
                            empty = True
                            break
                        continue
                    limit = -value / slope
                    if slope > 0:
                        lo = max(lo, limit)
                    else:
                        hi = min(hi, limit)
                    if lo > hi:
                        empty = True
                        break
                if empty or lo >= hi:
                    continue
                mat, translation = self.maps[index]
                for t in (lo, hi):
                    raw = tuple(start[i] + t * direction[i] for i in range(3))
                    local = tuple(raw[i] - shift[i] for i in range(3))
                    moved = apply_affine(mat, translation, local)
                    moved = tuple(moved[i] + shift[i] for i in range(3))
                    if t in samples and samples[t] != moved:
                        raise Failure("push is discontinuous across a cell wall")
                    samples[t] = moved
        out = []
        for t in sorted(set(samples) | {Fraction(1)}):
            if t == 0:
                continue
            out.append(samples.get(t, tuple(start[i] + t * direction[i] for i in range(3))))
        return out

    def image_of_polyline(self, points):
        out = [self.image_of_point(points[0])]
        for start, end in zip(points, points[1:]):
            out.extend(self.image_of_segment(start, end))
        return simplify(out)


def cross(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def simplify(points):
    out = [points[0]]
    for p in points[1:]:
        if p != out[-1]:
            out.append(p)
    if len(out) < 3:
        return out
    kept = [out[0]]
    for index in range(1, len(out) - 1):
        previous, current, following = kept[-1], out[index], out[index + 1]
        line = diff(following, previous)
        offset = diff(current, previous)
        if line != (0, 0, 0) and cross(offset, line) == (0, 0, 0):
            pivot = next(i for i in range(3) if line[i] != 0)
            if 0 < offset[pivot] / line[pivot] < 1:
                continue
        kept.append(current)
    kept.append(out[-1])
    return kept


# ---------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------


def check_planar(document) -> dict[str, Any]:
    square = document["square_triangulation"]
    vertices = [(F(v[0]), F(v[1])) for v in square["vertices"]]
    triangles = [tuple(tri) for tri in square["triangles"]]

    def area2(tri):
        a, b, c = (vertices[i] for i in tri)
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    require(all(area2(tri) > 0 for tri in triangles), "planar triangle not positively oriented")
    require(sum(area2(tri) for tri in triangles) == 2, "planar triangles do not fill the square")
    counts: dict[tuple[int, int], int] = {}
    for tri in triangles:
        for i in range(3):
            edge = tuple(sorted((tri[i], tri[(i + 1) % 3])))
            counts[edge] = counts.get(edge, 0) + 1
    require(set(counts.values()) <= {1, 2}, "planar edge in more than two triangles")
    require(
        sum(1 for count in counts.values() if count == 1) == 4,
        "planar boundary is not four edges",
    )
    star = [tri for tri in triangles if 0 in tri]
    require(len(star) == 4, "the star of m is not four triangles")
    ring = sorted({v for tri in star for v in tri} - {0})
    require(ring == [1, 2, 3, 4], "the star of m is not the quadrilateral P R1 Q R2")
    quad = [vertices[i] for i in (1, 2, 3, 4)]
    for i in range(4):
        a, b, c = quad[i], quad[(i + 1) % 4], quad[(i + 2) % 4]
        turn = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        require(turn > 0, "P R1 Q R2 is not convex")
    return {"planar_triangles": len(triangles), "star_triangles": len(star)}


def check_generator(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text(encoding="utf-8"))
    report: dict[str, Any] = {"file": path.name, "index": document["index"]}

    target = document["alpha_target"]
    prefix = document["alpha_prefix"]
    third = document["third_axis"]
    power = document["power"]
    side = document["side"]
    require({target, prefix, third} == {0, 1, 2}, "axes are not a permutation")
    require(power in (1, -1), "power is not a unit")
    require(side in EXPECTED_M_PRIME, "unknown side")
    require(document["side_bit"] in (0, 1), "side bit is not a bit")
    require(
        (document["side_bit"] == 0) == (side == "prefix-first"),
        "side bit and side name disagree",
    )

    # 1. the linear part -----------------------------------------------------
    matrix_e = matrix(document["linear_part_E"]["matrix"])
    expected = identity3()
    expected[prefix][target] += power
    require(matrix_e == expected, "E_k is not I + s E_{p,t}")
    require(all(entry.denominator == 1 for row in matrix_e for entry in row), "E_k not integral")
    require(det(matrix_e) == 1, "det E_k is not 1")
    inverse_e = matrix(document["linear_part_E"]["inverse_matrix"])
    require(mat_mul(matrix_e, inverse_e) == identity3(), "recorded E_k inverse is wrong")
    require(
        [F(c) for c in document["linear_part_E"]["translation"]] == [F(0)] * 3,
        "E_k has a translation",
    )

    # 2. the push cells ------------------------------------------------------
    push = document["push_Pi"]
    vertices = [point(v) for v in push["vertices"]]
    image_vertices = [point(v) for v in push["image_vertices"]]
    tetrahedra = [tuple(cell) for cell in push["tetrahedra"]]
    moved = push["moved_vertex"]
    require(len(vertices) == len(image_vertices), "vertex lists differ in length")
    require(len(tetrahedra) == 24, "expected 24 push cells")
    require(all(len(set(cell)) == 4 for cell in tetrahedra), "a cell repeats a vertex")

    determinants = [volume6([vertices[i] for i in cell]) for cell in tetrahedra]
    require(all(value > 0 for value in determinants), "a push cell is not positively oriented")
    image_determinants = [volume6([image_vertices[i] for i in cell]) for cell in tetrahedra]
    require(
        all(value > 0 for value in image_determinants),
        "m' is outside the kernel of the star: an image cell flips orientation",
    )
    require(
        [F(v) for v in push["cell_determinants"]] == determinants,
        "recorded cell determinants disagree",
    )
    require(
        [F(v) for v in push["image_cell_determinants"]] == image_determinants,
        "recorded image cell determinants disagree",
    )

    face_counts: dict[tuple[int, ...], int] = {}
    face_cells: dict[tuple[int, ...], list[int]] = {}
    for index, cell in enumerate(tetrahedra):
        for skip in range(4):
            face = tuple(sorted(cell[:skip] + cell[skip + 1:]))
            face_counts[face] = face_counts.get(face, 0) + 1
            face_cells.setdefault(face, []).append(index)
    require(max(face_counts.values()) <= 2, "a triangle is shared by three cells")
    boundary = sorted(face for face, count in face_counts.items() if count == 1)
    require(len(boundary) == 24, "expected 24 boundary triangles")
    require(
        boundary == sorted(tuple(face) for face in push["boundary_triangles"]),
        "recorded boundary triangles disagree",
    )
    require(all(moved in cell for cell in tetrahedra), "a cell misses the moved vertex")
    require(all(moved not in face for face in boundary), "the moved vertex is on the boundary")

    edge_counts: dict[tuple[int, int], int] = {}
    for face in boundary:
        for i in range(3):
            edge = tuple(sorted((face[i], face[(i + 1) % 3])))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    require(set(edge_counts.values()) == {2}, "the boundary surface is not closed")
    boundary_vertices = {v for face in boundary for v in face}
    chi = len(boundary_vertices) - len(edge_counts) + len(boundary)
    require(chi == 2, f"the boundary surface has Euler characteristic {chi}, not 2")
    adjacency: dict[int, set[int]] = {v: set() for v in boundary_vertices}
    for a, b in edge_counts:
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen = {next(iter(boundary_vertices))}
    stack = list(seen)
    while stack:
        current = stack.pop()
        for other in adjacency[current]:
            if other not in seen:
                seen.add(other)
                stack.append(other)
    require(seen == boundary_vertices, "the boundary surface is disconnected")

    volume_before = sum(determinants, Fraction(0)) / 6
    volume_after = sum(image_determinants, Fraction(0)) / 6
    require(volume_before == volume_after, "the push changes the support volume")
    require(
        F(push["support_volume"]) == volume_before,
        "recorded support volume disagrees",
    )

    # the recorded half-spaces really cut out the support: every cell vertex
    # satisfies them (so the star is inside the convex polytope) and every
    # boundary triangle is tight on one of them (so the star's boundary sphere
    # lies in the polytope's boundary sphere), hence the two coincide.
    halfspaces = [
        (point(item["normal"]), F(item["bound"])) for item in push["support_halfspaces"]
    ]
    require(len(halfspaces) == 6, "expected six support half-spaces")
    for vertex in vertices + image_vertices:
        for normal, bound in halfspaces:
            require(
                sum(normal[i] * vertex[i] for i in range(3)) <= bound,
                "a cell vertex violates a recorded support half-space",
            )
    for face in boundary:
        tight = [
            index
            for index, (normal, bound) in enumerate(halfspaces)
            if all(sum(normal[i] * vertices[v][i] for i in range(3)) == bound for v in face)
        ]
        require(tight, "a boundary triangle is not on the boundary of the support polytope")

    # 3. cell maps, continuity, positivity -----------------------------------
    maps = []
    for index, cell in enumerate(tetrahedra):
        record = push["cell_maps"][index]
        mat = matrix(record["matrix"])
        translation = point(record["translation"])
        for vertex in cell:
            require(
                apply_affine(mat, translation, vertices[vertex]) == image_vertices[vertex],
                f"cell {index} map does not send its vertices to the image vertices",
            )
        jac = det(mat)
        require(jac > 0, f"cell {index} has a non-positive Jacobian")
        require(jac == F(record["jacobian"]), "recorded Jacobian disagrees")
        maps.append((mat, translation))
    for face, cells in face_cells.items():
        if len(cells) != 2:
            continue
        left, right = cells
        for vertex in face:
            require(
                apply_affine(*maps[left], vertices[vertex])
                == apply_affine(*maps[right], vertices[vertex]),
                "the affine maps disagree on a shared face",
            )

    # 4. identity outside the support ----------------------------------------
    for index, (source, image) in enumerate(zip(vertices, image_vertices)):
        if index == moved:
            require(source != image, "the moved vertex does not move")
        else:
            require(source == image, f"vertex {index} is not fixed")
    require(point(push["moved_from"]) == vertices[moved], "moved_from disagrees")
    require(point(push["moved_to"]) == image_vertices[moved], "moved_to disagrees")

    planar_image = (
        image_vertices[moved][target],
        power * image_vertices[moved][prefix],
    )
    require(
        planar_image == EXPECTED_M_PRIME[side],
        f"m' is {planar_image}, not the value required by side {side}",
    )
    require(
        image_vertices[moved][third] == 0,
        "m' leaves the level c = 0",
    )

    # 5/6. inverse cells ------------------------------------------------------
    for index, cell in enumerate(tetrahedra):
        record = push["inverse_cell_maps"][index]
        mat = matrix(record["matrix"])
        translation = point(record["translation"])
        forward = maps[index]
        require(mat_mul(mat, forward[0]) == identity3(), "inverse matrix does not invert")
        for vertex in cell:
            back = apply_affine(mat, translation, image_vertices[vertex])
            require(back == vertices[vertex], "the inverse cell map is not the inverse")
        require(det(mat) > 0, "an inverse cell has a non-positive Jacobian")

    # 7. psi cells ------------------------------------------------------------
    psi = document["psi_cells"]
    domain = [point(v) for v in psi["domain_vertices"]]
    psi_image = [point(v) for v in psi["image_vertices"]]
    require([tuple(c) for c in psi["tetrahedra"]] == tetrahedra, "psi cells use other tetrahedra")
    require(psi_image == image_vertices, "psi image vertices differ from the push image")
    for index, vertex in enumerate(vertices):
        require(
            apply_affine(matrix_e, (F(0), F(0), F(0)), domain[index]) == vertex,
            "psi domain vertices are not E_k^{-1} of the push vertices",
        )
    for index, cell in enumerate(tetrahedra):
        record = psi["cell_maps"][index]
        mat = matrix(record["matrix"])
        translation = point(record["translation"])
        require(mat == mat_mul(maps[index][0], matrix_e), "psi cell map is not push o E")
        require(translation == maps[index][1], "psi cell translation is wrong")
        require(det(mat) > 0, "a psi cell has a non-positive Jacobian")
        for vertex in cell:
            require(
                apply_affine(mat, translation, domain[vertex]) == image_vertices[vertex],
                "psi cell map does not realise psi on its vertices",
            )
        back = psi["inverse_cell_maps"][index]
        back_mat = matrix(back["matrix"])
        back_translation = point(back["translation"])
        require(mat_mul(back_mat, mat) == identity3(), "psi inverse matrix does not invert")
        for vertex in cell:
            require(
                apply_affine(back_mat, back_translation, image_vertices[vertex])
                == domain[vertex],
                "psi inverse cell map is not the inverse",
            )

    # 8. the protected ball ---------------------------------------------------
    lows = [min(v[i] for v in vertices) for i in range(3)]
    highs = [max(v[i] for v in vertices) for i in range(3)]
    require(
        all(highs[i] - lows[i] < 1 for i in range(3)),
        "the support does not embed in the torus",
    )
    clearances = []
    for i in range(3):
        base = floor_div(lows[i])
        clearances.append(
            min(lows[i] - base, base + 1 - highs[i]) if highs[i] < base + 1 else Fraction(0)
        )
    ball_clearance = max(clearances)
    require(
        ball_clearance > PROTECTED_RADIUS,
        "the support of the push meets the protected ball",
    )
    for index, cell in enumerate(tetrahedra):
        corners = [vertices[i] for i in cell]
        best = Fraction(0)
        for i in range(3):
            lo = min(c[i] for c in corners)
            hi = max(c[i] for c in corners)
            base = floor_div(lo)
            if hi < base + 1:
                best = max(best, min(lo - base, base + 1 - hi))
        require(best > PROTECTED_RADIUS, f"cell {index} touches the protected ball")
    axis_clearance = {}
    for axis in range(3):
        value = max(clearances[i] for i in range(3) if i != axis)
        require(value > 0, f"the support meets the axis circle e_{axis}")
        axis_clearance[f"e_{axis}"] = str(value)

    # 9. spine images ---------------------------------------------------------
    raw_push = RawPush(vertices, tetrahedra, maps)
    computed = []
    for axis in range(3):
        end = [F(0), F(0), F(0)]
        end[axis] = F(1)
        loop = [(F(0), F(0), F(0)), tuple(end)]
        linear = [apply_affine(matrix_e, (F(0), F(0), F(0)), p) for p in loop]
        computed.append(raw_push.image_of_polyline(linear))
    for axis in range(3):
        recorded = [point(p) for p in document["spine_images"][f"C_{axis}"]]
        require(recorded == computed[axis], f"recorded spine image C_{axis} disagrees")
    induced = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for column, polyline in enumerate(computed):
        shift = diff(polyline[-1], polyline[0])
        for row in range(3):
            require(shift[row].denominator == 1, "a spine image is not closed in the torus")
            induced[row][column] = int(shift[row])
    require(induced == document["induced_H1_matrix"], "recorded H_1 matrix disagrees")
    require(
        induced == [[int(entry) for entry in row] for row in matrix_e],
        "the spine images do not induce E_k on H_1",
    )
    # the bent path: the diagonal P--Q must become P--m'--Q
    diagonal = computed[target]
    require(len(diagonal) == 5, "the diagonal image is not a four-segment path")
    require(diagonal[2] == image_vertices[moved], "the diagonal does not pass through m'")

    report.update(check_planar(document))
    report.update(
        {
            "cells": len(tetrahedra),
            "boundary_triangles": len(boundary),
            "support_volume": str(volume_before),
            "ball_clearance": str(ball_clearance),
            "ball_margin": str(ball_clearance - PROTECTED_RADIUS),
            "axis_clearance": axis_clearance,
            "spine_image_vertex_counts": [len(p) for p in computed],
            "induced_H1_matrix": induced,
            "status": "PASS",
        }
    )
    return report


THETA_PATH = GENERATOR_DIR / "gen_093_section_straightening.json"


def check_theta(path: Path) -> dict[str, Any]:
    """Independent check of the section-straightening homeomorphism ``theta``.

    Everything is regenerated from the raw node table ``(R_i, M_i)`` and the raw
    combinatorics (14 cube directions, 24 surface triangles, the 72 shell cell
    index tuples and the 24 cone cell index tuples).  No status field of the
    document is read; the recorded determinants and volumes are recomputed and
    compared, never trusted.
    """
    document = json.loads(path.read_text(encoding="utf-8"))
    directions = [point(d) for d in document["surface"]["directions"]]
    triangles = [tuple(t) for t in document["surface"]["triangles"]]
    shell_cells = [tuple(c) for c in document["shell_cells"]]
    cone_cells = [tuple(c) for c in document["cone_cells"]]
    require(len(directions) == 14, "expected 14 cube directions")
    require(len(triangles) == 24, "expected 24 surface triangles")
    require(len(shell_cells) == 72, "expected 72 cells per shell")
    require(len(cone_cells) == 24, "expected 24 cone cells")

    matrix_a = matrix(document["matrix_A"])
    inverse_a = matrix(document["matrix_A_inverse"])
    require(mat_mul(matrix_a, inverse_a) == identity3(), "A^{-1} does not invert A")
    norm = max(sum(abs(entry) for entry in row) for row in inverse_a)
    require(
        norm == document["matrix_A_inverse_infinity_norm"], "recorded |A^{-1}|_inf disagrees"
    )

    nodes = document["nodes"]
    radii = [F(node["radius"]) for node in nodes]
    maps = [matrix(node["matrix"]) for node in nodes]
    require(len(nodes) == document["shell_count"] + 1, "node count does not match shell count")
    require(maps[0] == identity3(), "theta is not the identity on the outer boundary")
    require(maps[-1] == inverse_a, "theta is not A^{-1} on the inner cube")
    require(F(document["r_outer"]) == radii[0], "recorded r_outer disagrees")
    require(F(document["r_inner"]) == radii[-1], "recorded r_inner disagrees")
    for i in range(1, len(radii)):
        require(0 < radii[i] < radii[i - 1], "radii are not strictly decreasing and positive")
    require(2 * radii[0] < 1, "C_out does not embed in the torus")
    for mat in maps:
        require(det(mat) == 1, "a node matrix is not unimodular")

    def surface(radius, mat=None):
        pts = [tuple(radius * c for c in d) for d in directions]
        if mat is None:
            return pts
        return [apply_affine(mat, (F(0), F(0), F(0)), p) for p in pts]

    def oriented(cells, pts):
        out = []
        for cell in cells:
            if volume6([pts[i] for i in cell]) < 0:
                cell = (cell[0], cell[1], cell[3], cell[2])
            out.append(cell)
        return out

    source_volume = Fraction(0)
    image_volume = Fraction(0)
    worst_source = None
    worst_image = None
    checked = 0
    recorded = document["determinants"]["per_shell"]
    require(len(recorded) == document["shell_count"], "per-shell record has the wrong length")
    for i in range(1, len(nodes)):
        src = surface(radii[i]) + surface(radii[i - 1])
        img = surface(radii[i], maps[i]) + surface(radii[i - 1], maps[i - 1])
        cells = oriented(shell_cells, src)
        src_dets = [volume6([src[j] for j in cell]) for cell in cells]
        img_dets = [volume6([img[j] for j in cell]) for cell in cells]
        require(min(src_dets) > 0, f"shell {i} has a non-positive source cell")
        require(min(img_dets) > 0, f"shell {i} has a non-positive image cell")
        require(
            F(recorded[i - 1]["min_source_determinant"]) == min(src_dets)
            and F(recorded[i - 1]["min_image_determinant"]) == min(img_dets),
            f"recorded extremal determinants of shell {i} disagree",
        )
        relative = mat_mul(inverse_of(maps[i - 1]), maps[i])
        spread = max(sum(abs(entry) for entry in row) for row in relative) * radii[i]
        require(spread <= radii[i - 1], f"image shell {i} is not nested in its outer shell")
        source_volume += sum(src_dets, Fraction(0)) / 6
        image_volume += sum(img_dets, Fraction(0)) / 6
        worst_source = min(src_dets) if worst_source is None else min(worst_source, min(src_dets))
        worst_image = min(img_dets) if worst_image is None else min(worst_image, min(img_dets))
        checked += len(cells)

    inner_src = surface(radii[-1]) + [(F(0), F(0), F(0))]
    inner_img = surface(radii[-1], inverse_a) + [(F(0), F(0), F(0))]
    cone = oriented(cone_cells, inner_src)
    cone_src = [volume6([inner_src[j] for j in cell]) for cell in cone]
    cone_img = [volume6([inner_img[j] for j in cell]) for cell in cone]
    require(min(cone_src) > 0 and min(cone_img) > 0, "the inner cone has a non-positive cell")
    source_volume += sum(cone_src, Fraction(0)) / 6
    image_volume += sum(cone_img, Fraction(0)) / 6
    worst_source = min(worst_source, min(cone_src))
    worst_image = min(worst_image, min(cone_img))
    checked += len(cone)

    cube = (2 * radii[0]) ** 3
    require(source_volume == cube, "the source cells do not tile C_out")
    require(image_volume == cube, "the image cells do not have the volume of C_out")
    require(checked == document["cells"]["total"], "recorded cell count disagrees")
    require(
        norm * radii[-1] < radii[0],
        "A^{-1} C_in is not inside the cube on which psi_A is linear",
    )
    require(
        F(document["determinants"]["min_source"]) == worst_source
        and F(document["determinants"]["min_image"]) == worst_image,
        "recorded global extremal determinants disagree",
    )
    return {
        "file": path.name,
        "shells": document["shell_count"],
        "cells_checked": checked,
        "r_outer": str(radii[0]),
        "r_inner": str(radii[-1]),
        "source_volume": str(source_volume),
        "image_volume": str(image_volume),
        "identity_on_outer_boundary": maps[0] == identity3(),
        "linear_A_inverse_on_inner_cube": maps[-1] == inverse_a,
        "inner_cube_image_inside_linear_regime": norm * radii[-1] < radii[0],
        "status": "PASS",
    }


def inverse_of(mat):
    determinant = det(mat)
    if determinant == 0:
        raise Failure("singular node matrix")
    cof = [
        [
            mat[(i + 1) % 3][(j + 1) % 3] * mat[(i + 2) % 3][(j + 2) % 3]
            - mat[(i + 1) % 3][(j + 2) % 3] * mat[(i + 2) % 3][(j + 1) % 3]
            for j in range(3)
        ]
        for i in range(3)
    ]
    return [[Fraction(cof[j][i]) / determinant for j in range(3)] for i in range(3)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", nargs="?", type=Path)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--theta", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    if args.theta:
        target = args.path or THETA_PATH
        try:
            report = check_theta(target)
        except Failure as error:
            print(f"FAIL {target.name}: {error}")
            print("THETA_PL_HOMEOMORPHISM=FAIL")
            return 1
        print(json.dumps(report, indent=2, sort_keys=True))
        print(f"THETA_CELLS_CHECKED={report['cells_checked']}")
        print("THETA_PL_HOMEOMORPHISM=PASS")
        return 0

    if args.all:
        paths = sorted(GENERATOR_DIR.glob("gen_[0-9][0-9][0-9].json"))
    elif args.path:
        paths = [args.path]
    else:
        parser.error("give a generator file or --all")

    failures = []
    reports = []
    for path in paths:
        try:
            reports.append(check_generator(path))
        except Failure as error:
            failures.append(f"{path.name}: {error}")
            reports.append({"file": path.name, "status": "FAIL", "reason": str(error)})

    if args.json:
        print(json.dumps(reports, indent=2, sort_keys=True))
    else:
        for report in reports:
            if report["status"] == "PASS":
                if len(paths) == 1:
                    print(json.dumps(report, indent=2, sort_keys=True))
            else:
                print(f"FAIL {report['file']}: {report['reason']}")
    print(f"GENERATORS_CHECKED={len(paths)}")
    print(f"PL_HOMEOMORPHISM_ALL={'PASS' if not failures else 'FAIL'}")
    print(f"PROTECTED_BALL_FIXED_BY_PUSH_FACTORS={'PASS' if not failures else 'FAIL'}")
    print("PROTECTED_BALL_FIXED_BY_LINEAR_FACTORS=FAIL")
    print(
        "PROTECTED_BALL_FIXED="
        f"{'FAIL' if not failures else 'FAIL'} "
        "reason=E_k is a transvection and moves every point with x_t != 0"
    )
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
