"""Tests for the dual-block P0a handlebody bridge.

Covers the builder (``scripts/build_t73_common_heegaard_complex.py``), the
elementary-collapse verifier (``scripts/verify_elementary_collapse.py``), a fast
subset of the Regina cross-check
(``scripts/verify_t73_handlebody_bridge_regina.py``), and mutation tests which
require the checkers to reject corrupted input.

The Regina tests are skipped -- loudly, never silently -- when the interpreter
running the tests cannot import Regina.  Use the Regina virtualenv to exercise
them::

    ~/ws/venv/bin/python -m unittest tests.test_t73_handlebody_bridge -v
"""

from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "geometry"
COLLAPSE_SEQUENCES = GEOMETRY / "t73_handlebody_collapse_sequences.json"
REGINA_CERTIFICATE = ROOT / "audit" / "t73_handlebody_bridge_regina.json"

HANDLE_NAMES = ("H_J0", "H_J1", "H_AR0", "H_AR1")
SPINE_NAMES = ("K_J0", "K_J1", "K_AR0", "K_AR1")
EXPECTED_TETRAHEDRA = {"H_J0": 1440, "H_J1": 7776, "H_AR0": 1440, "H_AR1": 7776}
EXPECTED_STEPS = {"H_J0": 3944, "H_J1": 17688, "H_AR0": 3944, "H_AR1": 17688}


def load_script(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILDER = load_script("build_t73_common_heegaard_complex")
COLLAPSER = load_script("verify_elementary_collapse")
REGINA_BRIDGE = load_script("verify_t73_handlebody_bridge_regina")

try:  # pragma: no cover - environment dependent
    import regina  # type: ignore

    REGINA_SKIP_REASON = ""
except ImportError as error:  # pragma: no cover - environment dependent
    regina = None  # type: ignore[assignment]
    REGINA_SKIP_REASON = (
        "Regina is not importable in this interpreter "
        f"({error}); rerun with the Regina virtualenv, e.g. "
        "~/ws/venv/bin/python -m unittest tests.test_t73_handlebody_bridge"
    )


class DualBlockComplexTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = BUILDER.build_model()
        cls.sequences = json.loads(COLLAPSE_SEQUENCES.read_text(encoding="utf-8"))

    # ---------------------------------------------------------------- model

    def test_common_triangulation_is_a_closed_3_manifold(self) -> None:
        checks = self.model["checks"]
        self.assertEqual(checks["T_vertices"], 64)
        self.assertEqual(checks["T_edges"], 448)
        self.assertEqual(checks["T_triangles"], 768)
        self.assertEqual(checks["T_tetrahedra"], 384)
        self.assertEqual(checks["T_euler_characteristic"], 0)
        self.assertTrue(checks["T_every_triangle_in_two_tetrahedra"])
        self.assertTrue(checks["T_all_vertex_links_are_spheres"])
        self.assertEqual(checks["Tprime_vertices"], 1664)
        self.assertEqual(checks["Tprime_tetrahedra"], 9216)
        self.assertEqual(checks["Tprime_euler_characteristic"], 0)
        self.assertTrue(checks["Tprime_every_triangle_in_two_tetrahedra"])
        self.assertTrue(checks["Tprime_all_vertex_links_are_spheres"])

    def test_matches_the_committed_ar_torus(self) -> None:
        self.assertEqual(
            BUILDER.committed_ar_tetrahedra(),
            set(self.model["T_simplices_by_dimension"][3]),
        )

    def test_spines_are_rank_three_roses(self) -> None:
        for name in SPINE_NAMES:
            spine = self.model["spines"][name]
            with self.subTest(spine=name):
                self.assertEqual(len(spine["vertices"]), 10)
                self.assertEqual(len(spine["edges"]), 12)
                self.assertEqual(spine["components"], 1)
                self.assertEqual(spine["rank"], 3)
                self.assertEqual(len(spine["subdivided_vertices"]), 22)
                self.assertEqual(len(spine["subdivided_edges"]), 24)
                # The regular-neighbourhood theorem is applied to K', which is
                # always full in T'.  K itself is *not* full in T here, and the
                # builder records that honestly rather than claiming otherwise.
                self.assertTrue(spine["subdivided_full_subcomplex_of_Tprime"])
                self.assertFalse(spine["full_subcomplex_of_T"])
                self.assertTrue(spine["non_full_witnesses_in_T"])

    def test_dual_block_handlebodies(self) -> None:
        for name in HANDLE_NAMES:
            report = self.model["handlebodies"][name]
            with self.subTest(handlebody=name):
                self.assertEqual(report["tetrahedron_count"], EXPECTED_TETRAHEDRA[name])
                self.assertEqual(report["euler_characteristic"], -2)
                self.assertTrue(report["manifold_face_multiplicities"])
                self.assertTrue(report["boundary_is_closed_connected_surface"])
                self.assertEqual(report["boundary_euler_characteristic"], -4)
                self.assertEqual(report["boundary_genus"], 3)
                self.assertEqual(report["boundary_triangle_count"], 1056)
                self.assertTrue(report["contains_subdivided_spine"])

    def test_each_pair_fills_the_torus_with_one_common_boundary(self) -> None:
        for pair in ("johnson", "ar"):
            checks = self.model["pair_checks"][pair]
            with self.subTest(pair=pair):
                self.assertTrue(checks["tetrahedra_disjoint"])
                self.assertTrue(checks["union_is_Tprime"])
                self.assertTrue(checks["shared_boundary_triangles"])
                self.assertTrue(checks["both_boundaries_genus_three"])
                self.assertTrue(checks["both_boundaries_closed_connected"])
                self.assertEqual(checks["common_boundary_genus"], 3)

    def test_translation_carries_the_johnson_pair_onto_the_ar_pair(self) -> None:
        translation = self.model["translation"]
        self.assertTrue(translation["is_simplicial_automorphism_of_T"])
        self.assertTrue(all(translation["spine_images"].values()))
        self.assertTrue(all(translation["handlebody_images"].values()))
        self.assertTrue(translation["maps_handlebody_pair"])

    def test_committed_geometry_documents_are_reproducible(self) -> None:
        for path, document in BUILDER.documents(self.model).items():
            with self.subTest(path=path.name):
                self.assertTrue(path.exists(), f"{path} is missing")
                committed = json.loads(path.read_text(encoding="utf-8"))
                self.assertEqual(committed, document)

    # ------------------------------------------------------------- collapse

    def test_collapse_verifier_self_test(self) -> None:
        self.assertTrue(COLLAPSER.self_test())

    def test_committed_collapse_sequences_verify(self) -> None:
        results = COLLAPSER.run_verification(self.sequences, paranoid_samples=25)
        self.assertTrue(results["simplex_index_sha256_matches"])
        for name in HANDLE_NAMES:
            entry = results["collapses"][name]
            with self.subTest(handlebody=name):
                self.assertTrue(entry["ok"], entry.get("reason"))
                self.assertEqual(entry["step_count"], EXPECTED_STEPS[name])
        self.assertTrue(results["all_ok"])

    # -------------------------------------------------------------- helpers

    def _collapse_inputs(self, handle_name: str, spine_name: str):
        maximal, target = COLLAPSER.handlebody_inputs(BUILDER, self.model, handle_name, spine_name)
        steps = COLLAPSER.decode_steps(
            self.sequences["collapses"][handle_name]["steps"], self.model["Tprime_simplices"]
        )
        return maximal, target, steps

    # ------------------------------------------------------------ mutations

    def test_mutation_swapped_collapse_step_is_rejected(self) -> None:
        maximal, target, steps = self._collapse_inputs("H_J0", "K_J0")
        mutated = list(steps)
        sigma, tau = mutated[0]
        mutated[0] = (tau, sigma)
        outcome = COLLAPSER.verify_collapse(maximal, target, mutated)
        self.assertFalse(outcome["ok"])
        self.assertEqual(outcome["failed_step"], 0)

    def test_mutation_non_free_collapse_step_is_rejected(self) -> None:
        maximal, target, steps = self._collapse_inputs("H_J0", "K_J0")
        cofacets = COLLAPSER.cofacet_map(COLLAPSER.closure(maximal))
        not_free = next(step for step in steps if len(cofacets[step[0]]) >= 2)
        outcome = COLLAPSER.verify_collapse(maximal, target, [not_free] + list(steps))
        self.assertFalse(outcome["ok"])
        self.assertEqual(outcome["failed_step"], 0)
        self.assertIn("not free", outcome["reason"])

    def test_mutation_truncated_collapse_sequence_is_rejected(self) -> None:
        maximal, target, steps = self._collapse_inputs("H_AR1", "K_AR1")
        outcome = COLLAPSER.verify_collapse(maximal, target, list(steps)[:-1])
        self.assertFalse(outcome["ok"])
        self.assertIsNone(outcome["failed_step"])
        self.assertIn("remainder", outcome["reason"])

    def test_mutation_collapsing_into_the_spine_is_rejected(self) -> None:
        maximal, target, steps = self._collapse_inputs("H_J0", "K_J0")
        spine_edge = tuple(target[0])
        spine_vertex = (spine_edge[0],)
        outcome = COLLAPSER.verify_collapse(
            maximal, target, [(spine_vertex, spine_edge)] + list(steps)
        )
        self.assertFalse(outcome["ok"])
        self.assertEqual(outcome["failed_step"], 0)
        self.assertIn("lies in K", outcome["reason"])

    def test_mutation_moved_tetrahedron_breaks_the_boundary_checks(self) -> None:
        zero = list(self.model["handlebodies"]["H_J0"]["tetrahedron_indices"])
        one = list(self.model["handlebodies"]["H_J1"]["tetrahedron_indices"])
        moved_zero = BUILDER.handlebody_report(self.model, zero[1:], "K_J0")
        moved_one = BUILDER.handlebody_report(self.model, sorted(one + zero[:1]), "K_J1")
        checks = BUILDER.check_pair(self.model, moved_zero, moved_one)
        self.assertNotEqual(moved_zero["boundary_genus"], 3)
        self.assertFalse(checks["both_boundaries_genus_three"])
        self.assertFalse(checks["both_boundaries_closed_connected"])

    def test_mutation_duplicated_tetrahedron_breaks_the_union_check(self) -> None:
        zero = list(self.model["handlebodies"]["H_J0"]["tetrahedron_indices"])
        one = list(self.model["handlebodies"]["H_J1"]["tetrahedron_indices"])
        report_zero = BUILDER.handlebody_report(self.model, zero, "K_J0")
        report_one = BUILDER.handlebody_report(self.model, sorted(one + zero[:1]), "K_J1")
        checks = BUILDER.check_pair(self.model, report_zero, report_one)
        self.assertFalse(checks["tetrahedra_disjoint"])
        self.assertFalse(checks["union_is_Tprime"])
        self.assertFalse(checks["shared_boundary_triangles"])

    def test_mutation_wrong_translation_is_rejected(self) -> None:
        self.assertTrue(BUILDER.translation_maps_pair(self.model, (-1, -1, -1)))
        for shift in ((-1, 0, 0), (0, -1, -1), (-2, -2, -2), (1, 1, 1)):
            with self.subTest(shift=shift):
                self.assertFalse(BUILDER.translation_maps_pair(self.model, shift))


@unittest.skipUnless(regina is not None, REGINA_SKIP_REASON)
class ReginaCrossCheckTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = BUILDER.build_model()

    def test_regina_recognises_four_genus_three_handlebodies(self) -> None:
        recognised = []
        for name in HANDLE_NAMES:
            triangulation, gluings, boundary_faces = REGINA_BRIDGE.build_regina_triangulation(
                regina, self.model, name
            )
            with self.subTest(handlebody=name):
                self.assertEqual(triangulation.size(), EXPECTED_TETRAHEDRA[name])
                self.assertEqual(boundary_faces, 1056)
                self.assertTrue(triangulation.isValid())
                self.assertTrue(triangulation.isOrientable())
                self.assertTrue(triangulation.isConnected())
                self.assertEqual(triangulation.countBoundaryComponents(), 1)
                self.assertEqual(int(triangulation.boundaryComponent(0).eulerChar()), -4)

                simplified = regina.Triangulation3(triangulation)
                simplified.simplify()
                homology = simplified.homology()
                self.assertEqual(int(homology.rank()), 3)
                self.assertEqual(int(homology.countInvariantFactors()), 0)
                self.assertTrue(hasattr(simplified, "recogniseHandlebody"))
                genus = int(simplified.recogniseHandlebody())
                self.assertEqual(genus, 3)
                recognised.append(genus)
        self.assertEqual(recognised, [3, 3, 3, 3])

    def test_committed_regina_certificate_agrees_with_a_fresh_recomputation(self) -> None:
        self.assertTrue(REGINA_CERTIFICATE.exists(), f"{REGINA_CERTIFICATE} is missing")
        committed = json.loads(REGINA_CERTIFICATE.read_text(encoding="utf-8"))
        for name in HANDLE_NAMES:
            triangulation, _, _ = REGINA_BRIDGE.build_regina_triangulation(
                regina, self.model, name
            )
            simplified = regina.Triangulation3(triangulation)
            simplified.simplify()
            genus = int(simplified.recogniseHandlebody())
            with self.subTest(handlebody=name):
                entry = committed["regina"]["handlebodies"][name]
                self.assertEqual(entry["recognise_handlebody_genus"], genus)
                self.assertEqual(entry["raw_tetrahedra"], triangulation.size())
                self.assertEqual(
                    entry["boundary_euler_characteristic"],
                    int(triangulation.boundaryComponent(0).eulerChar()),
                )

    def _recognise(self, cells):
        triangulation, _, boundary_faces = REGINA_BRIDGE.build_regina_from_cells(regina, cells)
        simplified = regina.Triangulation3(triangulation)
        simplified.simplify()
        return int(simplified.recogniseHandlebody()), boundary_faces

    def test_regina_genus_tracks_the_rank_of_the_spine(self) -> None:
        """Negative control: the same construction on a rank-2 spine gives genus 2.

        Note that deleting a single boundary tetrahedron is *not* a usable
        mutation here: a boundary shelling does not change the homeomorphism
        type, and Regina correctly still reports genus 3.  Changing the spine
        does change the answer, which is what the certificate relies on.
        """
        two_circles = set()
        for axis in (0, 1):
            for step in range(4):
                point = [0, 0, 0]
                point[axis] = step
                two_circles.add(BUILDER.vid(point))
        self.assertEqual(len(two_circles), 7)
        cells = [chain for chain in self.model["Tprime_tetrahedra"] if chain[0] in two_circles]
        genus, _ = self._recognise(cells)
        self.assertEqual(genus, 2)

    def test_regina_rejects_the_closed_torus(self) -> None:
        """Negative control: T' itself is closed, so it is not a handlebody."""
        genus, boundary_faces = self._recognise(self.model["Tprime_tetrahedra"])
        self.assertEqual(boundary_faces, 0)
        self.assertNotEqual(genus, 3)


if __name__ == "__main__":
    unittest.main()
