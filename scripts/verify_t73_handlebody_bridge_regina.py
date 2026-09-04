#!/usr/bin/env python3
"""Independent Regina cross-check of the dual-block P0a handlebody bridge.

The PL proof of Lemma P0a is combinatorial and does not use Regina:

  * ``scripts/build_t73_common_heegaard_complex.py`` builds the common period-4
    Freudenthal torus ``T``, its barycentric subdivision ``T'``, the four
    coordinate spines and the four dual-block regular neighbourhoods
    ``H = N(K;T')`` / ``closure(T' - H)``;
  * ``scripts/verify_elementary_collapse.py`` replays explicit elementary
    collapses ``H^i -> K^i``, which by Rourke--Sanderson identify each ``H^i``
    as a regular neighbourhood of a rank-3 graph in an orientable 3-manifold,
    hence as a genus-3 handlebody.

This script re-runs both of those from scratch and then, independently, hands
each of the four complexes to Regina as a ``Triangulation3`` (one Regina
tetrahedron per ``T'``-tetrahedron, glued along every shared triangle with the
``Perm4`` induced by the vertex labels).  Regina is asked for validity,
orientability, the boundary surface, ``H_1``, and ``recogniseHandlebody()``.

Nothing in this script trusts a status field written by another script: the
handlebodies, the collapses, the union/boundary relations and the translation
``T(v)=v-(1,1,1)`` are all recomputed here.  ``--check`` regenerates the whole
certificate in memory and compares it against the committed JSON.

    ~/ws/venv/bin/python -B scripts/verify_t73_handlebody_bridge_regina.py --write
    ~/ws/venv/bin/python -B scripts/verify_t73_handlebody_bridge_regina.py --check
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import time
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "audit" / "t73_handlebody_bridge_regina.json"
COLLAPSE_SEQUENCES = ROOT / "geometry" / "t73_handlebody_collapse_sequences.json"

ORDER = ("H_J0", "H_J1", "H_AR0", "H_AR1")
SPINE_OF = {"H_J0": "K_J0", "H_J1": "K_J1", "H_AR0": "K_AR0", "H_AR1": "K_AR1"}
RAW_HOMOLOGY_LIMIT = 2000  # raw H_1 on 7776 tetrahedra takes many minutes

# Fields that legitimately vary with the Regina build or its randomised
# simplification heuristics; --check re-asserts their invariants instead of
# comparing them literally.
VOLATILE_FIELDS = ("regina_version", "simplified_tetrahedra", "elapsed_seconds")


def canonical_sha(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest().upper()


def load_script(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------
# Regina
# --------------------------------------------------------------------------


def import_regina():
    try:
        import regina  # type: ignore
    except ImportError:
        return None
    return regina


def build_regina_from_cells(regina, cells: Sequence[Sequence[int]]):
    """One Regina tetrahedron per abstract 4-tuple, joined along shared triangles.

    ``cells`` must be sorted 4-tuples of distinct vertex labels.  Two cells that
    share a triangle are glued with the ``Perm4`` that matches vertex labels, so
    the Regina triangulation carries exactly the simplicial structure of the
    input complex.  Returns ``(triangulation, gluings, boundary_triangles)``.
    """
    cells = [tuple(cell) for cell in cells]
    triangulation = regina.Triangulation3()
    objects = [triangulation.newTetrahedron() for _ in cells]
    owner: dict[tuple[int, ...], tuple[int, int]] = {}
    gluings = 0
    for index, cell in enumerate(cells):
        for i in range(4):
            face = tuple(cell[k] for k in range(4) if k != i)
            if face in owner:
                other_index, j = owner.pop(face)
                other = cells[other_index]
                image = [0, 0, 0, 0]
                image[i] = j
                for k in range(4):
                    if k == i:
                        continue
                    image[k] = other.index(cell[k])
                if sorted(image) != [0, 1, 2, 3]:
                    raise AssertionError("gluing permutation is not a bijection")
                objects[index].join(i, objects[other_index], regina.Perm4(*image))
                gluings += 1
            else:
                owner[face] = (index, i)
    return triangulation, gluings, len(owner)


def build_regina_triangulation(regina, model: dict[str, Any], handle_name: str):
    """Regina triangulation of one dual-block handlebody of the model."""
    tprime = model["Tprime_tetrahedra"]
    cells = [
        tprime[position]
        for position in model["handlebodies"][handle_name]["tetrahedron_indices"]
    ]
    return build_regina_from_cells(regina, cells)


def regina_report(regina, model: dict[str, Any], handle_name: str) -> dict[str, Any]:
    started = time.time()
    triangulation, gluings, boundary_faces = build_regina_triangulation(regina, model, handle_name)
    raw_size = triangulation.size()

    boundary_components = triangulation.countBoundaryComponents()
    boundary_chi = (
        triangulation.boundaryComponent(0).eulerChar() if boundary_components == 1 else None
    )
    report: dict[str, Any] = {
        "raw_tetrahedra": raw_size,
        "internal_gluings": gluings,
        "boundary_triangles": boundary_faces,
        "is_valid": bool(triangulation.isValid()),
        "is_orientable": bool(triangulation.isOrientable()),
        "is_connected": bool(triangulation.isConnected()),
        "boundary_components": boundary_components,
        "boundary_euler_characteristic": int(boundary_chi) if boundary_chi is not None else None,
        "boundary_genus": (2 - int(boundary_chi)) // 2 if boundary_chi is not None else None,
    }

    if raw_size <= RAW_HOMOLOGY_LIMIT:
        raw_h1 = triangulation.homology()
        report["raw_homology_h1"] = raw_h1.str()
        report["raw_homology_rank"] = int(raw_h1.rank())
        report["raw_homology_invariant_factors"] = int(raw_h1.countInvariantFactors())
    else:
        report["raw_homology_h1"] = None
        report["raw_homology_rank"] = None
        report["raw_homology_invariant_factors"] = None
        report["raw_homology_skipped_because"] = (
            f"raw triangulation has {raw_size} > {RAW_HOMOLOGY_LIMIT} tetrahedra; "
            "H_1 is computed on the simplified copy instead"
        )

    simplified = regina.Triangulation3(triangulation)
    if hasattr(simplified, "simplify"):
        simplified.simplify()
        report["simplify_method"] = "simplify"
    elif hasattr(simplified, "intelligentSimplify"):
        simplified.intelligentSimplify()
        report["simplify_method"] = "intelligentSimplify"
    else:
        report["simplify_method"] = None
    report["simplified_tetrahedra"] = simplified.size()

    h1 = simplified.homology()
    report["homology_h1"] = h1.str()
    report["homology_rank"] = int(h1.rank())
    report["homology_invariant_factors"] = int(h1.countInvariantFactors())
    report["homology_is_free_rank_three"] = bool(h1.rank() == 3 and h1.countInvariantFactors() == 0)
    report["simplified_boundary_components"] = simplified.countBoundaryComponents()
    report["simplified_boundary_euler_characteristic"] = (
        int(simplified.boundaryComponent(0).eulerChar())
        if simplified.countBoundaryComponents() == 1
        else None
    )

    if hasattr(simplified, "recogniseHandlebody"):
        report["recognise_handlebody_available"] = True
        report["recognise_handlebody_genus"] = int(simplified.recogniseHandlebody())
    elif hasattr(simplified, "isHandlebody"):
        report["recognise_handlebody_available"] = False
        report["recognise_handlebody_genus"] = 3 if simplified.isHandlebody(3) else -1
    else:
        report["recognise_handlebody_available"] = False
        report["recognise_handlebody_genus"] = None

    report["elapsed_seconds"] = round(time.time() - started, 3)
    report["verdict"] = bool(
        report["is_valid"]
        and report["is_orientable"]
        and report["is_connected"]
        and report["boundary_components"] == 1
        and report["boundary_euler_characteristic"] == -4
        and report["boundary_genus"] == 3
        and report["homology_is_free_rank_three"]
        and report["recognise_handlebody_genus"] == 3
        and (report["raw_homology_rank"] in (None, 3))
        and (report["raw_homology_invariant_factors"] in (None, 0))
    )
    return report


# --------------------------------------------------------------------------
# the certificate
# --------------------------------------------------------------------------


def combinatorial_section(builder, collapser, model: dict[str, Any]) -> dict[str, Any]:
    spines = {
        name: {
            "base_point": spine["base_point"],
            "vertices": len(spine["vertices"]),
            "edges": len(spine["edges"]),
            "components": spine["components"],
            "rank": spine["rank"],
            "full_subcomplex_of_T": spine["full_subcomplex_of_T"],
            "subdivided_full_subcomplex_of_Tprime": spine["subdivided_full_subcomplex_of_Tprime"],
        }
        for name, spine in model["spines"].items()
    }
    handlebodies = {
        name: {
            "spine": report["spine"],
            "tetrahedron_count": report["tetrahedron_count"],
            "simplex_counts": report["simplex_counts"],
            "euler_characteristic": report["euler_characteristic"],
            "boundary_triangle_count": report["boundary_triangle_count"],
            "boundary_euler_characteristic": report["boundary_euler_characteristic"],
            "boundary_genus": report["boundary_genus"],
            "boundary_is_closed_connected_surface": report["boundary_is_closed_connected_surface"],
            "manifold_face_multiplicities": report["manifold_face_multiplicities"],
            "contains_subdivided_spine": report["contains_subdivided_spine"],
        }
        for name, report in model["handlebodies"].items()
    }

    document = json.loads(COLLAPSE_SEQUENCES.read_text(encoding="utf-8"))
    collapse = collapser.run_verification(document)
    collapse_section = {
        name: {
            "ok": entry["ok"],
            "reason": entry.get("reason", ""),
            "failed_step": entry.get("failed_step"),
            "step_count": entry.get("step_count", 0),
            "complex_simplices": entry.get("complex_simplices"),
            "target_simplices": entry.get("target_simplices"),
        }
        for name, entry in collapse["collapses"].items()
    }

    return {
        "triangulation_checks": model["checks"],
        "spines": spines,
        "handlebodies": handlebodies,
        "pair_checks": model["pair_checks"],
        "translation": model["translation"],
        "translation_mutation_100_fails": not builder.translation_maps_pair(model, (-1, 0, 0)),
        "collapse_sequence_index_sha_matches": collapse["simplex_index_sha256_matches"],
        "collapses": collapse_section,
        "collapse_sequences_sha256": hashlib.sha256(
            COLLAPSE_SEQUENCES.read_bytes()
        ).hexdigest().upper(),
    }


def generate(with_regina: bool = True) -> dict[str, Any]:
    builder = load_script("build_t73_common_heegaard_complex")
    collapser = load_script("verify_elementary_collapse")
    model = builder.build_model()

    stable = combinatorial_section(builder, collapser, model)
    stable["matches_committed_ar_torus"] = builder.committed_ar_tetrahedra() == set(
        model["T_simplices_by_dimension"][3]
    )

    regina = import_regina() if with_regina else None
    regina_section: dict[str, Any] = {
        "available": regina is not None,
        "regina_version": regina.versionString() if regina is not None else None,
        "handlebodies": {},
    }
    if regina is not None:
        for name in ORDER:
            regina_section["handlebodies"][name] = regina_report(regina, model, name)

    combinatorial_ok = bool(
        stable["triangulation_checks"]["T_euler_characteristic"] == 0
        and stable["triangulation_checks"]["Tprime_euler_characteristic"] == 0
        and stable["triangulation_checks"]["Tprime_every_triangle_in_two_tetrahedra"]
        and stable["triangulation_checks"]["Tprime_all_vertex_links_are_spheres"]
        and all(spine["rank"] == 3 for spine in stable["spines"].values())
        and all(
            spine["subdivided_full_subcomplex_of_Tprime"] for spine in stable["spines"].values()
        )
        and all(entry["boundary_genus"] == 3 for entry in stable["handlebodies"].values())
        and all(entry["euler_characteristic"] == -2 for entry in stable["handlebodies"].values())
        and all(entry["contains_subdivided_spine"] for entry in stable["handlebodies"].values())
        and all(
            pair["union_is_Tprime"]
            and pair["shared_boundary_triangles"]
            and pair["both_boundaries_genus_three"]
            and pair["both_boundaries_closed_connected"]
            for pair in stable["pair_checks"].values()
        )
        and stable["translation"]["maps_handlebody_pair"]
        and stable["translation_mutation_100_fails"]
        and stable["collapse_sequence_index_sha_matches"]
        and all(entry["ok"] for entry in stable["collapses"].values())
        and stable["matches_committed_ar_torus"]
    )
    regina_ok = bool(
        regina is not None
        and all(regina_section["handlebodies"][name]["verdict"] for name in ORDER)
    )

    certificate: dict[str, Any] = {
        "schema": "t73_handlebody_bridge_regina/v1",
        "generator": "scripts/verify_t73_handlebody_bridge_regina.py",
        "lemma": "lem:P0a (Johnson--AR handlebody bridge), paper/spc4-t73-candidate/main.tex",
        "proof_route": (
            "dual-block regular neighbourhoods in the first barycentric subdivision "
            "plus explicit elementary collapses; Regina is an independent cross-check "
            "and is not used by the PL argument"
        ),
        "inputs": {
            "triangulation": "geometry/t73_common_torus_triangulation.json",
            "johnson_spines": "geometry/t73_johnson_spines.json",
            "ar_spines": "geometry/t73_ar_spines.json",
            "handlebodies": "geometry/t73_handlebody_pair_dual_block.json",
            "collapse_sequences": "geometry/t73_handlebody_collapse_sequences.json",
        },
        "volatile_fields": list(VOLATILE_FIELDS),
        "combinatorial": stable,
        "regina": regina_section,
        "combinatorial_verdict": combinatorial_ok,
        "regina_verdict": regina_ok,
        "p0a_bridge": bool(combinatorial_ok and regina_ok),
    }
    certificate = json.loads(json.dumps(certificate, sort_keys=True, separators=(",", ":")))
    certificate["stable_sha256"] = canonical_sha(strip_volatile(certificate))
    return certificate


def strip_volatile(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: strip_volatile(item)
            for key, item in value.items()
            if key not in VOLATILE_FIELDS and key != "stable_sha256"
        }
    if isinstance(value, list):
        return [strip_volatile(item) for item in value]
    return value


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def print_lines(certificate: dict[str, Any]) -> None:
    combinatorial = certificate["combinatorial"]
    regina_section = certificate["regina"]

    print(f"REGINA_AVAILABLE={regina_section['available']}")
    print(f"REGINA_VERSION={regina_section['regina_version']}")
    for name in ORDER:
        entry = combinatorial["collapses"].get(name, {"ok": False})
        status = "PASS" if entry.get("ok") else "FAIL"
        suffix = "" if entry.get("ok") else f" step={entry.get('failed_step')} reason={entry.get('reason')}"
        print(f"COLLAPSE_{name}={status} steps={entry.get('step_count', 0)}{suffix}")
    ranks = sorted({spine["rank"] for spine in combinatorial["spines"].values()})
    print(f"SPINE_RANK={ranks[0] if len(ranks) == 1 else ranks}")
    print(f"COMMON_BOUNDARY_GENUS={combinatorial['pair_checks']['johnson']['common_boundary_genus']}")
    union_ok = all(
        pair["union_is_Tprime"] and pair["shared_boundary_triangles"]
        for pair in combinatorial["pair_checks"].values()
    )
    print(f"UNION_IS_T3={'PASS' if union_ok else 'FAIL'}")
    print(
        "S_MAPS_HANDLEBODY_PAIR="
        f"{'PASS' if combinatorial['translation']['maps_handlebody_pair'] else 'FAIL'}"
    )
    print(
        "S_MUTATION_100_REJECTED="
        f"{'PASS' if combinatorial['translation_mutation_100_fails'] else 'FAIL'}"
    )
    genus_list = [
        regina_section["handlebodies"].get(name, {}).get("recognise_handlebody_genus")
        for name in ORDER
    ]
    print(f"REGINA_HANDLEBODY_GENUS={json.dumps(genus_list, separators=(',', ':'))}")
    raw_sizes = [regina_section["handlebodies"].get(name, {}).get("raw_tetrahedra") for name in ORDER]
    simplified_sizes = [
        regina_section["handlebodies"].get(name, {}).get("simplified_tetrahedra") for name in ORDER
    ]
    print(f"REGINA_RAW_TETRAHEDRA={json.dumps(raw_sizes, separators=(',', ':'))}")
    print(f"REGINA_SIMPLIFIED_TETRAHEDRA={json.dumps(simplified_sizes, separators=(',', ':'))}")
    print(f"COMBINATORIAL_VERDICT={'PASS' if certificate['combinatorial_verdict'] else 'FAIL'}")
    print(f"REGINA_VERDICT={'PASS' if certificate['regina_verdict'] else 'FAIL'}")
    print(f"P0A_BRIDGE={'PASS' if certificate['p0a_bridge'] else 'FAIL'}")
    print(f"STABLE_SHA256={certificate['stable_sha256']}")


def compare(regenerated: dict[str, Any]) -> list[str]:
    if not CERTIFICATE.exists():
        return [f"{CERTIFICATE.relative_to(ROOT)}:MISSING"]
    committed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    problems: list[str] = []
    if strip_volatile(committed) != strip_volatile(regenerated):
        problems.append("stable_payload_differs")
    if committed.get("stable_sha256") != regenerated["stable_sha256"]:
        problems.append("stable_sha256_differs")
    committed_regina = committed.get("regina", {}).get("handlebodies", {})
    for name in ORDER:
        entry = committed_regina.get(name, {})
        if entry.get("recognise_handlebody_genus") != 3:
            problems.append(f"{name}:committed_genus_not_three")
        fresh = regenerated["regina"]["handlebodies"].get(name, {})
        if fresh.get("recognise_handlebody_genus") != 3:
            problems.append(f"{name}:recomputed_genus_not_three")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--no-regina", action="store_true", help="combinatorial checks only")
    args = parser.parse_args()

    if not COLLAPSE_SEQUENCES.exists():
        print("P0A_BRIDGE=FAIL reason=missing_collapse_sequences")
        return 1

    certificate = generate(with_regina=not args.no_regina)
    if args.write:
        CERTIFICATE.parent.mkdir(parents=True, exist_ok=True)
        CERTIFICATE.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(f"WROTE={CERTIFICATE.relative_to(ROOT)} BYTES={CERTIFICATE.stat().st_size}")

    print_lines(certificate)

    ok = bool(certificate["p0a_bridge"])
    if args.check:
        problems = compare(certificate)
        print(f"COMMITTED_CERTIFICATE={'PASS' if not problems else 'FAIL ' + ','.join(problems)}")
        ok = ok and not problems
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
