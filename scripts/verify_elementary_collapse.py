#!/usr/bin/env python3
"""Generic elementary-collapse checker, plus the search that produces the
collapse sequences for the four dual-block handlebodies of Lemma P0a.

Definitions
-----------
A simplicial complex is given by a list of *maximal* simplices (sorted tuples of
vertex indices); it is understood to contain every non-empty face of every
listed simplex.  An *elementary collapse* removes a pair ``(sigma, tau)`` where

  * ``tau`` is a simplex currently present,
  * ``sigma`` is a proper face of ``tau`` with ``dim sigma = dim tau - 1``,
  * ``sigma`` is *free*: ``tau`` is the unique proper coface of ``sigma`` in the
    current complex,
  * neither ``sigma`` nor ``tau`` belongs to the protected subcomplex ``K``.

``X`` *collapses to* ``K`` when a sequence of elementary collapses turns ``X``
into exactly ``K``.  By Rourke--Sanderson, a compact PL manifold that collapses
onto a subpolyhedron is a regular neighbourhood of it; a regular neighbourhood
of a rank-3 graph in an orientable 3-manifold is a genus-3 handlebody.

This module recomputes the complex and replays the sequence from scratch.  It
never reads a status field from the certificate; the only inputs it trusts are
the triangulation rules in ``scripts/build_t73_common_heegaard_complex.py`` and
the raw step list.

Freeness is decided by counting *codimension-one* cofaces.  In a complex closed
under faces this is equivalent to counting all proper cofaces: if a present
``rho`` strictly contains ``sigma`` with ``dim rho >= dim sigma + 2`` then every
``sigma + {x}``, ``x`` in ``rho - sigma``, is a present codimension-one coface,
and there are at least two of them.  ``--paranoid`` re-derives the coface set by
brute force on a random sample of steps as an independent cross-check.

Usage
-----
    python3 scripts/verify_elementary_collapse.py --search   # write sequences
    python3 scripts/verify_elementary_collapse.py --check    # replay + verify
    python3 scripts/verify_elementary_collapse.py --self-test
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import random
from collections import defaultdict
from fractions import Fraction
from heapq import heappop, heappush
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "geometry"
SEQUENCE_OUTPUT = GEOMETRY / "t73_handlebody_collapse_sequences.json"

PAIRS = (("H_J0", "K_J0"), ("H_J1", "K_J1"), ("H_AR0", "K_AR0"), ("H_AR1", "K_AR1"))

Simplex = tuple[int, ...]


def canonical_sha(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest().upper()


def load_builder():
    path = ROOT / "scripts" / "build_t73_common_heegaard_complex.py"
    spec = importlib.util.spec_from_file_location("build_t73_common_heegaard_complex", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------
# generic complex machinery
# --------------------------------------------------------------------------


def closure(maximal: Iterable[Simplex]) -> set[Simplex]:
    out: set[Simplex] = set()
    for cell in maximal:
        cell = tuple(cell)
        for size in range(1, len(cell) + 1):
            out.update(itertools.combinations(cell, size))
    return out


def facets(simplex: Simplex) -> tuple[Simplex, ...]:
    return tuple(itertools.combinations(simplex, len(simplex) - 1))


def cofacet_map(present: Iterable[Simplex]) -> dict[Simplex, set[Simplex]]:
    cofacets: dict[Simplex, set[Simplex]] = defaultdict(set)
    for cell in present:
        if len(cell) > 1:
            for facet in facets(cell):
                cofacets[facet].add(cell)
    return cofacets


def brute_force_proper_cofaces(present: set[Simplex], simplex: Simplex) -> list[Simplex]:
    target = set(simplex)
    return [
        cell
        for cell in present
        if len(cell) > len(simplex) and target <= set(cell)
    ]


# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------


def verify_collapse(
    maximal: Iterable[Simplex],
    target_maximal: Iterable[Simplex],
    steps: Sequence[Sequence[Simplex]],
    paranoid_samples: int = 0,
    seed: int = 20260904,
) -> dict[str, Any]:
    """Replay ``steps`` on ``closure(maximal)`` and check the result equals ``closure(target)``."""
    present = closure(maximal)
    protected = closure(target_maximal)
    if not protected <= present:
        return {
            "ok": False,
            "failed_step": None,
            "reason": "the protected subcomplex is not contained in the complex",
        }
    cofacets = cofacet_map(present)

    sampled: set[int] = set()
    if paranoid_samples > 0 and steps:
        rng = random.Random(seed)
        sampled = set(rng.sample(range(len(steps)), min(paranoid_samples, len(steps))))

    for index, step in enumerate(steps):
        if len(step) != 2:
            return {"ok": False, "failed_step": index, "reason": "step is not a pair"}
        sigma = tuple(step[0])
        tau = tuple(step[1])
        if tau not in present:
            return {"ok": False, "failed_step": index, "reason": f"tau {tau} is not present"}
        if sigma not in present:
            return {"ok": False, "failed_step": index, "reason": f"sigma {sigma} is not present"}
        if len(tau) != len(sigma) + 1:
            return {"ok": False, "failed_step": index, "reason": "tau is not a codimension-one coface"}
        if not set(sigma) < set(tau):
            return {"ok": False, "failed_step": index, "reason": "sigma is not a proper face of tau"}
        if sigma in protected:
            return {"ok": False, "failed_step": index, "reason": f"sigma {sigma} lies in K"}
        if tau in protected:
            return {"ok": False, "failed_step": index, "reason": f"tau {tau} lies in K"}
        cofaces = cofacets.get(sigma, set())
        if cofaces != {tau}:
            return {
                "ok": False,
                "failed_step": index,
                "reason": f"sigma {sigma} is not free: cofacets {sorted(cofaces)}",
            }
        if cofacets.get(tau):
            return {
                "ok": False,
                "failed_step": index,
                "reason": f"tau {tau} is not maximal: cofacets {sorted(cofacets[tau])}",
            }
        if index in sampled:
            brute = brute_force_proper_cofaces(present, sigma)
            if brute != [tau]:
                return {
                    "ok": False,
                    "failed_step": index,
                    "reason": f"brute-force coface check disagrees for {sigma}: {sorted(brute)}",
                }
        for facet in facets(tau):
            cofacets[facet].discard(tau)
        present.discard(tau)
        for facet in facets(sigma):
            cofacets[facet].discard(sigma)
        present.discard(sigma)

    if present != protected:
        return {
            "ok": False,
            "failed_step": None,
            "reason": (
                f"remainder has {len(present)} simplices, target has {len(protected)}; "
                f"extra={len(present - protected)} missing={len(protected - present)}"
            ),
        }
    return {"ok": True, "failed_step": None, "reason": "", "steps": len(steps)}


# --------------------------------------------------------------------------
# search
# --------------------------------------------------------------------------


def find_collapse(
    maximal: Iterable[Simplex],
    target_maximal: Iterable[Simplex],
    weight: Callable[[Simplex], float] | None = None,
    seed: int = 0,
    jitter: float = 0.0,
) -> dict[str, Any]:
    """Greedy elementary collapse.

    Priority: collapse the highest-dimensional free pair available, and among
    those the one whose free face is farthest from the protected subcomplex
    (``weight``).  ``jitter`` adds random noise to the priority so that repeated
    restarts explore different orders.
    """
    present = closure(maximal)
    protected = closure(target_maximal)
    if not protected <= present:
        raise ValueError("the protected subcomplex is not contained in the complex")
    cofacets = cofacet_map(present)
    rng = random.Random(seed)

    def score(simplex: Simplex) -> float:
        base = 0.0 if weight is None else weight(simplex)
        if jitter:
            base += rng.uniform(-jitter, jitter)
        return base

    heap: list[tuple[int, float, Simplex]] = []
    pushed: set[Simplex] = set()

    def offer(simplex: Simplex) -> None:
        if simplex in protected or simplex not in present:
            return
        coface_set = cofacets.get(simplex)
        if coface_set is None or len(coface_set) != 1:
            return
        (tau,) = coface_set
        if tau in protected:
            return
        heappush(heap, (-len(simplex), -score(simplex), simplex))
        pushed.add(simplex)

    for simplex in present:
        offer(simplex)

    steps: list[tuple[Simplex, Simplex]] = []
    while heap:
        _, _, sigma = heappop(heap)
        pushed.discard(sigma)
        if sigma in protected or sigma not in present:
            continue
        coface_set = cofacets.get(sigma)
        if coface_set is None or len(coface_set) != 1:
            continue
        (tau,) = coface_set
        if tau in protected or tau not in present or cofacets.get(tau):
            continue
        steps.append((sigma, tau))
        touched: set[Simplex] = set()
        for facet in facets(tau):
            cofacets[facet].discard(tau)
            touched.add(facet)
        present.discard(tau)
        for facet in facets(sigma):
            cofacets[facet].discard(sigma)
            touched.add(facet)
        present.discard(sigma)
        for facet in touched:
            if facet in present:
                offer(facet)

    remainder = present - protected
    return {
        "complete": not remainder,
        "steps": steps,
        "remainder": sorted(remainder),
        "remainder_size": len(remainder),
    }


def search_with_restarts(
    maximal: Iterable[Simplex],
    target_maximal: Iterable[Simplex],
    weight: Callable[[Simplex], float] | None = None,
    restarts: int = 12,
) -> dict[str, Any]:
    maximal = list(maximal)
    target_maximal = list(target_maximal)
    attempts: list[dict[str, Any]] = []
    plans: list[tuple[Callable[[Simplex], float] | None, float, int]] = [(weight, 0.0, 0), (None, 0.0, 0)]
    for attempt in range(restarts):
        plans.append((weight, 1.0 + attempt, 1000 + attempt))
        plans.append((None, 1.0 + attempt, 5000 + attempt))
    for chosen_weight, jitter, seed in plans:
        outcome = find_collapse(
            maximal, target_maximal, weight=chosen_weight, seed=seed, jitter=jitter
        )
        attempts.append(
            {
                "seed": seed,
                "jitter": jitter,
                "weighted": chosen_weight is not None,
                "remainder_size": outcome["remainder_size"],
            }
        )
        if outcome["complete"]:
            outcome["attempts"] = attempts
            return outcome
    outcome["attempts"] = attempts  # type: ignore[name-defined]
    return outcome  # type: ignore[name-defined]


# --------------------------------------------------------------------------
# the four handlebodies
# --------------------------------------------------------------------------


def wrap(value: Fraction, period: int) -> Fraction:
    residue = value % period
    if residue * 2 > period:
        residue -= period
    return residue


def spine_distance_weight(builder, model: dict[str, Any], spine_name: str) -> Callable[[Simplex], float]:
    """Squared torus distance from a T'-simplex barycentre to the axis spine."""
    period = model["period"]
    base = [Fraction(c) for c in model["spines"][spine_name]["base_point"]]
    coordinates = model["Tprime_vertex_barycentre"]

    cache: dict[Simplex, float] = {}

    def weight(simplex: Simplex) -> float:
        cached = cache.get(simplex)
        if cached is not None:
            return cached
        count = len(simplex)
        point = [
            sum((coordinates[v][j] for v in simplex), Fraction(0)) / count for j in range(3)
        ]
        best = None
        for axis in range(3):
            total = Fraction(0)
            for j in range(3):
                if j == axis:
                    continue
                delta = wrap(point[j] - base[j], period)
                total += delta * delta
            if best is None or total < best:
                best = total
        value = float(best)  # type: ignore[arg-type]
        cache[simplex] = value
        return value

    return weight


def handlebody_inputs(builder, model: dict[str, Any], handle_name: str, spine_name: str):
    tprime = model["Tprime_tetrahedra"]
    maximal = [tprime[position] for position in model["handlebodies"][handle_name]["tetrahedron_indices"]]
    spine = model["spines"][spine_name]
    target = [tuple(edge) for edge in spine["subdivided_edges"]]
    isolated = set(spine["subdivided_vertices"]) - {v for edge in target for v in edge}
    target += [(v,) for v in sorted(isolated)]
    return maximal, target


def simplex_index_table(model: dict[str, Any]) -> tuple[dict[Simplex, int], list[Simplex]]:
    table = model["Tprime_simplex_index"]
    inverse = model["Tprime_simplices"]
    return table, inverse


def encode_steps(steps: Sequence[tuple[Simplex, Simplex]], table: dict[Simplex, int]) -> list[list[int]]:
    return [[table[sigma], table[tau]] for sigma, tau in steps]


def decode_steps(encoded: Sequence[Sequence[int]], inverse: Sequence[Simplex]) -> list[tuple[Simplex, Simplex]]:
    return [(tuple(inverse[a]), tuple(inverse[b])) for a, b in encoded]


def run_search(restarts: int = 12) -> dict[str, Any]:
    builder = load_builder()
    model = builder.build_model()
    table, inverse = simplex_index_table(model)
    document: dict[str, Any] = {
        "schema": "t73_handlebody_collapse_sequences/v1",
        "generator": "scripts/verify_elementary_collapse.py --search",
        "triangulation": "geometry/t73_common_torus_triangulation.json",
        "handlebodies": "geometry/t73_handlebody_pair_dual_block.json",
        "encoding": (
            "Each step is [sigma, tau] where the integers index the global list of "
            "simplices of T', ordered by (dimension, sorted tuple of T'-vertex "
            "indices).  T'-vertices are themselves indexed by the order (dimension, "
            "sorted tuple) on the simplices of T, so T'-vertex i < 64 is T-vertex i."
        ),
        "simplex_index_sha256": canonical_sha([list(s) for s in inverse]),
        "simplex_count": len(inverse),
        "collapses": {},
    }
    for handle_name, spine_name in PAIRS:
        maximal, target = handlebody_inputs(builder, model, handle_name, spine_name)
        weight = spine_distance_weight(builder, model, spine_name)
        outcome = search_with_restarts(maximal, target, weight=weight, restarts=restarts)
        entry: dict[str, Any] = {
            "handlebody": handle_name,
            "spine": spine_name,
            "target_maximal": [list(cell) for cell in target],
            "complete": outcome["complete"],
            "step_count": len(outcome["steps"]),
            "attempts": outcome["attempts"],
        }
        if outcome["complete"]:
            entry["steps"] = encode_steps(outcome["steps"], table)
        else:
            entry["steps"] = encode_steps(outcome["steps"], table)
            entry["stuck_remainder"] = [table[cell] for cell in outcome["remainder"]]
            entry["stuck_remainder_size"] = outcome["remainder_size"]
        document["collapses"][handle_name] = entry
    return document


def run_verification(document: dict[str, Any], paranoid_samples: int = 0) -> dict[str, Any]:
    builder = load_builder()
    model = builder.build_model()
    table, inverse = simplex_index_table(model)
    expected_sha = canonical_sha([list(s) for s in inverse])
    results: dict[str, Any] = {
        "simplex_index_sha256_matches": document.get("simplex_index_sha256") == expected_sha,
        "collapses": {},
    }
    for handle_name, spine_name in PAIRS:
        entry = document.get("collapses", {}).get(handle_name)
        if entry is None:
            results["collapses"][handle_name] = {
                "ok": False,
                "reason": "no collapse sequence recorded",
                "step_count": 0,
            }
            continue
        maximal, target = handlebody_inputs(builder, model, handle_name, spine_name)
        recorded_target = [tuple(cell) for cell in entry.get("target_maximal", [])]
        if set(recorded_target) != set(target):
            results["collapses"][handle_name] = {
                "ok": False,
                "reason": "recorded target subcomplex differs from the recomputed spine",
                "step_count": len(entry.get("steps", [])),
            }
            continue
        steps = decode_steps(entry.get("steps", []), inverse)
        outcome = verify_collapse(
            maximal, target, steps, paranoid_samples=paranoid_samples
        )
        outcome["step_count"] = len(steps)
        outcome["complex_simplices"] = len(closure(maximal))
        outcome["target_simplices"] = len(closure(target))
        results["collapses"][handle_name] = outcome
    results["all_ok"] = bool(
        results["simplex_index_sha256_matches"]
        and all(entry["ok"] for entry in results["collapses"].values())
    )
    return results


# --------------------------------------------------------------------------
# self test on tiny examples
# --------------------------------------------------------------------------


def self_test() -> bool:
    ok = True

    # A single triangle collapses to one of its edges.
    triangle = [(0, 1, 2)]
    target = [(0, 1)]
    outcome = find_collapse(triangle, target)
    ok &= outcome["complete"]
    ok &= verify_collapse(triangle, target, outcome["steps"], paranoid_samples=10)["ok"]

    # The 2-sphere (boundary of a tetrahedron) is not collapsible to a point.
    sphere = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    stuck = find_collapse(sphere, [(0,)])
    ok &= not stuck["complete"]

    # A corrupted step must be rejected.
    good = find_collapse(triangle, target)["steps"]
    swapped = [(tau, sigma) for sigma, tau in good]
    ok &= not verify_collapse(triangle, target, swapped)["ok"]

    # A non-free face must be rejected: in the full triangle, (0,1) has two
    # cofaces once nothing has been removed, so (0,) is not free for ((0,),(0,1)).
    ok &= not verify_collapse(triangle, [(2,)], [((0,), (0, 1))])["ok"]

    # A truncated sequence must be rejected by the final equality test.
    ok &= not verify_collapse(triangle, target, good[:-1])["ok"]
    return bool(ok)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def dump(document: dict[str, Any]) -> str:
    return json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"


def print_results(results: dict[str, Any]) -> bool:
    labels = {"H_J0": "COLLAPSE_H_J0", "H_J1": "COLLAPSE_H_J1", "H_AR0": "COLLAPSE_H_AR0", "H_AR1": "COLLAPSE_H_AR1"}
    ok = bool(results["simplex_index_sha256_matches"])
    for handle_name, label in labels.items():
        entry = results["collapses"].get(handle_name, {"ok": False, "reason": "missing"})
        status = "PASS" if entry.get("ok") else "FAIL"
        ok &= bool(entry.get("ok"))
        detail = ""
        if not entry.get("ok"):
            detail = f" step={entry.get('failed_step')} reason={entry.get('reason')}"
        print(f"{label}={status} steps={entry.get('step_count', 0)}{detail}")
    print(f"COLLAPSE_INDEX_SHA={'PASS' if results['simplex_index_sha256_matches'] else 'FAIL'}")
    print(f"COLLAPSE_ALL={'PASS' if ok else 'FAIL'}")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--search", action="store_true", help="find and write collapse sequences")
    parser.add_argument("--check", action="store_true", help="verify the committed sequences")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--paranoid", type=int, default=0, help="brute-force freeness on N sampled steps")
    parser.add_argument("--restarts", type=int, default=12)
    args = parser.parse_args()

    if args.self_test:
        passed = self_test()
        print(f"COLLAPSE_SELF_TEST={'PASS' if passed else 'FAIL'}")
        if not (args.search or args.check):
            return 0 if passed else 1

    if args.search:
        document = run_search(restarts=args.restarts)
        GEOMETRY.mkdir(parents=True, exist_ok=True)
        SEQUENCE_OUTPUT.write_text(dump(document), encoding="utf-8")
        print(f"WROTE={SEQUENCE_OUTPUT.relative_to(ROOT)} BYTES={SEQUENCE_OUTPUT.stat().st_size}")
        for handle_name, entry in sorted(document["collapses"].items()):
            print(
                f"SEARCH_{handle_name}={'COMPLETE' if entry['complete'] else 'STUCK'} "
                f"steps={entry['step_count']} "
                f"remainder={entry.get('stuck_remainder_size', 0)}"
            )

    if args.check or args.search:
        if not SEQUENCE_OUTPUT.exists():
            print("COLLAPSE_ALL=FAIL reason=missing_sequences")
            return 1
        document = json.loads(SEQUENCE_OUTPUT.read_text(encoding="utf-8"))
        results = run_verification(document, paranoid_samples=args.paranoid)
        return 0 if print_results(results) else 1

    if not args.self_test:
        parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
