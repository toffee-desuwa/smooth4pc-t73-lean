from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATOR_DIR = ROOT / "geometry" / "t73_johnson_generators"
MATRIX_A = [[0, 269, 1240], [0, 41, 189], [1, 0, 32]]


def load(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read(index: int) -> dict:
    return json.loads((GENERATOR_DIR / f"gen_{index:03d}.json").read_text(encoding="utf-8"))


def write_temporary(document: dict) -> Path:
    handle = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8")
    handle.write(json.dumps(document))
    handle.close()
    return Path(handle.name)


class GeneratorFileTest(unittest.TestCase):
    def test_ninety_three_generator_files_and_an_index(self):
        files = sorted(GENERATOR_DIR.glob("gen_[0-9][0-9][0-9].json"))
        self.assertEqual(len(files), 93)
        index = json.loads((GENERATOR_DIR / "index.json").read_text(encoding="utf-8"))
        self.assertEqual(index["generator_count"], 93)
        self.assertEqual(index["matrix_A"], MATRIX_A)
        self.assertTrue(index["composite_matrix_equals_A"])
        self.assertEqual(index["protected_ball_radius"], "1/196104")

    def test_linear_parts_multiply_to_A_in_the_recorded_order(self):
        current = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for position in range(93):
            local = read(position)["linear_part_E"]["matrix"]
            current = [
                [sum(local[i][k] * current[k][j] for k in range(3)) for j in range(3)]
                for i in range(3)
            ]
        self.assertEqual(current, MATRIX_A)

    def test_generators_rebuild_bit_for_bit(self):
        builder = load("build_t73_johnson_pl_generators")
        moves = builder.unit_moves()
        bits = builder.side_bits()
        for position in (0, 46, 92):
            rebuilt = builder.build_generator(position, moves[position], int(bits[position]))
            self.assertEqual(rebuilt, read(position))


class VerifierTest(unittest.TestCase):
    def setUp(self):
        self.verifier = load("verify_t73_pl_homeomorphism")

    def test_selected_generators_verify(self):
        for position in (0, 1, 46, 91, 92):
            report = self.verifier.check_generator(GENERATOR_DIR / f"gen_{position:03d}.json")
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["cells"], 24)
            self.assertEqual(report["boundary_triangles"], 24)
            self.assertEqual(report["support_volume"], "3/32")
            self.assertEqual(report["ball_clearance"], "1/8")

    def test_every_generator_verifies(self):
        failures = []
        for path in sorted(GENERATOR_DIR.glob("gen_[0-9][0-9][0-9].json")):
            try:
                self.verifier.check_generator(path)
            except self.verifier.Failure as error:  # pragma: no cover - regression guard
                failures.append(f"{path.name}: {error}")
        self.assertEqual(failures, [])


class MutationTest(unittest.TestCase):
    def setUp(self):
        self.verifier = load("verify_t73_pl_homeomorphism")
        self.builder = load("build_t73_johnson_pl_generators")

    def _rebuild_push_with_new_image(self, document, bad_planar):
        """Rewrite every affine record consistently for a different ``m'``."""
        builder = self.builder
        push = document["push_Pi"]
        target = document["alpha_target"]
        prefix = document["alpha_prefix"]
        third = document["third_axis"]
        power = document["power"]
        vertices = [tuple(Fraction(c) for c in v) for v in push["vertices"]]
        tetrahedra = [tuple(cell) for cell in push["tetrahedra"]]
        moved = push["moved_vertex"]
        image = list(vertices)
        image[moved] = builder.to_ambient(
            bad_planar[0], bad_planar[1], Fraction(0), target, prefix, third, power
        )
        matrix_e = [[Fraction(v) for v in row] for row in document["linear_part_E"]["matrix"]]
        inverse_e = [
            [Fraction(v) for v in row] for row in document["linear_part_E"]["inverse_matrix"]
        ]
        cell_maps, inverse_cell_maps, psi_maps, psi_inverse_maps, dets = [], [], [], [], []
        for cell in tetrahedra:
            source = [vertices[i] for i in cell]
            destination = [image[i] for i in cell]
            mat, translation = builder.affine_from_correspondence(source, destination)
            inverse_mat = builder.inverse3(mat)
            inverse_translation = tuple(
                -builder.matvec(inverse_mat, translation)[i] for i in range(3)
            )
            composite = builder.matmul(mat, matrix_e)
            inverse_composite = builder.matmul(inverse_e, inverse_mat)
            inverse_composite_translation = tuple(
                builder.matvec(inverse_e, inverse_translation)[i] for i in range(3)
            )
            cell_maps.append(
                {
                    "matrix": builder.enc_matrix(mat),
                    "translation": builder.enc_point(translation),
                    "jacobian": builder.fs(builder.det3(mat)),
                }
            )
            inverse_cell_maps.append(
                {
                    "matrix": builder.enc_matrix(inverse_mat),
                    "translation": builder.enc_point(inverse_translation),
                    "jacobian": builder.fs(builder.det3(inverse_mat)),
                }
            )
            psi_maps.append(
                {
                    "matrix": builder.enc_matrix(composite),
                    "translation": builder.enc_point(translation),
                    "jacobian": builder.fs(builder.det3(composite)),
                }
            )
            psi_inverse_maps.append(
                {
                    "matrix": builder.enc_matrix(inverse_composite),
                    "translation": builder.enc_point(inverse_composite_translation),
                    "jacobian": builder.fs(builder.det3(inverse_composite)),
                }
            )
            dets.append(builder.fs(builder.tet_determinant([image[i] for i in cell])))
        push["image_vertices"] = [builder.enc_point(v) for v in image]
        push["moved_to"] = builder.enc_point(image[moved])
        push["m_prime_planar"] = [builder.fs(bad_planar[0]), builder.fs(bad_planar[1])]
        push["image_cell_determinants"] = dets
        push["cell_maps"] = cell_maps
        push["inverse_cell_maps"] = inverse_cell_maps
        document["psi_cells"]["image_vertices"] = push["image_vertices"]
        document["psi_cells"]["cell_maps"] = psi_maps
        document["psi_cells"]["inverse_cell_maps"] = psi_inverse_maps
        return document

    def test_moving_m_prime_outside_the_star_kernel_is_rejected(self):
        document = read(0)
        outside = (Fraction(1, 8), Fraction(1, 8))
        document = self._rebuild_push_with_new_image(document, outside)
        path = write_temporary(document)
        try:
            with self.assertRaises(self.verifier.Failure) as caught:
                self.verifier.check_generator(path)
        finally:
            path.unlink()
        self.assertIn("kernel", str(caught.exception))

    def test_a_valid_alternative_m_prime_is_still_accepted_by_the_kernel_check(self):
        """Control for the previous test: an interior point keeps the kernel check happy."""
        document = read(0)
        interior = (Fraction(9, 16), Fraction(9, 16))
        document = self._rebuild_push_with_new_image(document, interior)
        path = write_temporary(document)
        try:
            with self.assertRaises(self.verifier.Failure) as caught:
                self.verifier.check_generator(path)
        finally:
            path.unlink()
        # the kernel is fine; the side bit is what now fails
        self.assertNotIn("kernel", str(caught.exception))
        self.assertIn("side", str(caught.exception))

    def test_flipping_a_cell_orientation_is_rejected(self):
        document = read(0)
        cell = document["push_Pi"]["tetrahedra"][0]
        document["push_Pi"]["tetrahedra"][0] = [cell[1], cell[0], cell[2], cell[3]]
        path = write_temporary(document)
        try:
            with self.assertRaises(self.verifier.Failure) as caught:
                self.verifier.check_generator(path)
        finally:
            path.unlink()
        self.assertIn("positively oriented", str(caught.exception))

    def test_negating_a_recorded_determinant_is_rejected(self):
        document = read(0)
        value = Fraction(document["push_Pi"]["cell_determinants"][3])
        document["push_Pi"]["cell_determinants"][3] = str(-value)
        path = write_temporary(document)
        try:
            with self.assertRaises(self.verifier.Failure) as caught:
                self.verifier.check_generator(path)
        finally:
            path.unlink()
        self.assertIn("determinants disagree", str(caught.exception))

    def test_changing_a_side_bit_changes_the_spine_image_and_the_sha(self):
        builder = self.builder
        moves = builder.unit_moves()
        bits = builder.side_bits()
        committed = read(7)
        flipped = builder.build_generator(7, moves[7], 1 - int(bits[7]))
        self.assertNotEqual(flipped["side"], committed["side"])
        target = committed["alpha_target"]
        self.assertNotEqual(
            flipped["spine_images"][f"C_{target}"],
            committed["spine_images"][f"C_{target}"],
        )
        self.assertNotEqual(flipped["generator_sha256"], committed["generator_sha256"])
        # the induced map on H_1 is unchanged: both sides lift the same transvection
        self.assertEqual(flipped["induced_H1_matrix"], committed["induced_H1_matrix"])

    def test_changing_a_linear_matrix_entry_breaks_the_composite_and_the_verifier(self):
        current = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for position in range(93):
            document = read(position)
            local = document["linear_part_E"]["matrix"]
            if position == 12:
                local = [row[:] for row in local]
                local[0][2] += 1
            current = [
                [sum(local[i][k] * current[k][j] for k in range(3)) for j in range(3)]
                for i in range(3)
            ]
        self.assertNotEqual(current, MATRIX_A)

        document = read(12)
        document["linear_part_E"]["matrix"][0][2] += 1
        path = write_temporary(document)
        try:
            with self.assertRaises(self.verifier.Failure) as caught:
                self.verifier.check_generator(path)
        finally:
            path.unlink()
        self.assertIn("E_k", str(caught.exception))


class CompositeTest(unittest.TestCase):
    def setUp(self):
        self.compose = load("compose_t73_psi_A")

    def test_committed_composite_document(self):
        document = json.loads(
            (ROOT / "geometry" / "t73_psi_A.json").read_text(encoding="utf-8")
        )
        self.assertEqual(document["generator_count"], 93)
        self.assertEqual(document["composite_linear_matrix"], MATRIX_A)
        self.assertTrue(document["composite_matrix_equals_A"])
        self.assertTrue(document["local_normal_form"]["certificate_holds"])
        self.assertFalse(document["protected_ball"]["psi_A_fixes_the_ball_pointwise"])
        self.assertTrue(document["protected_ball"]["push_factors_fix_the_ball_pointwise"])
        diagnostics = document["handlebody_diagnostics"]
        self.assertEqual(diagnostics["per_generator_pass_count"], 0)
        self.assertEqual(diagnostics["first_prefix_leaving_H_J0"], 1)
        self.assertEqual(diagnostics["setwise_preservation"]["status"], "OPEN")
        self.assertTrue(document["attaching_link_status"].startswith("OPEN"))
        self.assertFalse(document["exact_prefix_transport"]["complete"])

    def test_membership_test_agrees_with_commit3_handlebody(self):
        report = self.compose.crosscheck_membership()
        self.assertEqual(report["mismatches"], 0)
        self.assertEqual(report["tprime_tetrahedra_probed"], 9216)
        self.assertTrue(report["spine_vertex_rule_matches_commit3"])

    def test_short_prefix_transport_reproduces(self):
        generators = self.compose.load_generators()[:6]
        transport = self.compose.prefix_transport(generators, budget=200, progress=False)
        self.assertEqual(transport["max_spine_vertex_counts"][:6], [2, 5, 11, 19, 31, 39])
        self.assertEqual(transport["first_prefix_leaving_H_J0"], 1)
        self.assertFalse(transport["reports"][1]["spine_image_in_H_J0"])
        self.assertTrue(transport["reports"][0]["spine_image_in_H_J0"])

    def test_point_transport_fixes_the_origin_and_is_linear_on_the_protected_ball(self):
        generators = self.compose.load_generators()
        origin = (Fraction(0), Fraction(0), Fraction(0))
        current = origin
        for generator in generators:
            current = generator.apply_point(current)
        self.assertEqual(current, origin)

        probe = (Fraction(1, 400000), Fraction(0), Fraction(0))
        current = probe
        for generator in generators:
            current = generator.apply_point(current)
        expected = tuple(
            sum(Fraction(MATRIX_A[i][j]) * probe[j] for j in range(3)) for i in range(3)
        )
        self.assertEqual(current, expected)


if __name__ == "__main__":
    unittest.main()


THETA_PATH = GENERATOR_DIR / "gen_093_section_straightening.json"


class SectionStraighteningTest(unittest.TestCase):
    def setUp(self):
        self.verifier = load("verify_t73_pl_homeomorphism")
        self.compose = load("compose_t73_psi_A")

    def test_theta_file_is_referenced_from_the_index(self):
        index = json.loads((GENERATOR_DIR / "index.json").read_text(encoding="utf-8"))
        self.assertEqual(
            index["section_straightening"]["file"], "gen_093_section_straightening.json"
        )
        self.assertTrue(THETA_PATH.exists())
        self.assertEqual(len(sorted(GENERATOR_DIR.glob("gen_[0-9][0-9][0-9].json"))), 93)

    def test_theta_verifies_independently(self):
        report = self.verifier.check_theta(THETA_PATH)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["shells"], 186)
        self.assertEqual(report["cells_checked"], 13416)
        self.assertEqual(report["source_volume"], report["image_volume"])
        self.assertTrue(report["identity_on_outer_boundary"])
        self.assertTrue(report["linear_A_inverse_on_inner_cube"])
        self.assertTrue(report["inner_cube_image_inside_linear_regime"])

    def test_theta_rebuilds_bit_for_bit(self):
        builder = load("build_t73_section_straightening")
        self.assertEqual(builder.build(), json.loads(THETA_PATH.read_text(encoding="utf-8")))

    def test_single_shell_straightening_is_impossible_for_this_A(self):
        document = json.loads(THETA_PATH.read_text(encoding="utf-8"))
        attempt = document["single_shell_attempt"]
        self.assertEqual(attempt["result"], "IMPOSSIBLE for this A, at every choice of r_inner")
        for entry in attempt["attempts"]:
            self.assertGreater(entry["non_positive_image_cells"], 0)
            self.assertLess(entry["min_image_determinant_sign"], 1)

    def test_a_live_single_shell_is_rejected_by_the_verifier(self):
        document = json.loads(THETA_PATH.read_text(encoding="utf-8"))
        nodes = document["nodes"]
        document["nodes"] = [nodes[0], nodes[-1]]
        document["shell_count"] = 1
        document["determinants"]["per_shell"] = [
            {"min_source_determinant": "0", "min_image_determinant": "0"}
        ]
        document["cells"]["total"] = 96
        path = write_temporary(document)
        try:
            with self.assertRaises(self.verifier.Failure) as caught:
                self.verifier.check_theta(path)
        finally:
            path.unlink()
        self.assertIn("non-positive image cell", str(caught.exception))

    def test_perturbing_a_node_matrix_is_rejected(self):
        document = json.loads(THETA_PATH.read_text(encoding="utf-8"))
        document["nodes"][-1]["matrix"][0][0] = str(
            Fraction(document["nodes"][-1]["matrix"][0][0]) + 1
        )
        path = write_temporary(document)
        try:
            with self.assertRaises(self.verifier.Failure) as caught:
                self.verifier.check_theta(path)
        finally:
            path.unlink()
        self.assertIn("A^{-1}", str(caught.exception))

    def test_Psi_is_the_identity_on_the_section_ball(self):
        generators = self.compose.load_generators()
        report = self.compose.section_straightening_report(generators)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["failures"], [])
        self.assertGreaterEqual(report["probe_count"], 9)
        self.assertTrue(report["probes_returned_to_themselves"])
        self.assertEqual(report["induced_H1_matrix_of_Psi"], MATRIX_A)
        r_inner = Fraction(report["fixed_ball_radius"])
        self.assertGreater(r_inner, 0)
        self.assertLess(r_inner * report["matrix_A_inverse_infinity_norm"], Fraction(1, 196104))
