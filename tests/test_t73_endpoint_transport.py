from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class EndpointTransportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.b = load("build_t73_endpoint_transport")
        cls.convention, cls.audit = cls.b.compute(sigma=-1, verbose=False)

    def test_oriented_model_is_consistent(self) -> None:
        checks = self.audit["model_checks"]
        self.assertTrue(all(checks["intertwining"].values()))
        self.assertTrue(all(checks["inverses"].values()))
        self.assertTrue(all(checks["braid_relation"].values()))
        self.assertTrue(all(checks["phi"].values()))
        self.assertTrue(all(v is True for v in checks["duality"].values()))
        self.assertEqual(checks["loop_value_at_q1"], -2)

    def test_letterwise_transport_covers_all_orientation_patterns(self) -> None:
        letter = self.audit["letterwise_transport"]
        self.assertTrue(letter["all_letters_pass"])
        self.assertTrue(letter["final_permutation_is_identity"])
        self.assertEqual(letter["letters"], 45360)
        self.assertEqual(set(letter["pattern_counts"]), {"V+V", "V-V", "V+V*", "V-V*", "V*+V", "V*-V", "V*+V*", "V*-V*"})

    def test_cup_and_cap_are_derived_not_hardcoded(self) -> None:
        self.assertEqual(self.audit["u_public_constant_terms"], [[2, 1], [87, -1]])
        self.assertEqual(self.audit["ell_public_constant_terms"], [[2, -1], [87, 1]])
        derived = self.b.derive_endpoint_terms(self.b.CONVENTION)
        self.assertEqual(derived["u_terms"], [[2, 1], [87, -1]])
        self.assertEqual(derived["ell_terms"], [[2, -1], [87, 1]])

    def test_cubic_agrees_three_ways_and_controls_reproduce_erratum(self) -> None:
        d3 = self.audit["delta3"]
        self.assertTrue(d3["agree"])
        self.assertEqual(d3["constant_terms_pipeline"], 2624)
        controls = self.audit["coordinate_controls"]
        self.assertEqual(controls["thxy_u_collar_ell_collar_word"], -59072)
        self.assertEqual(controls["thxy_u_thxy_ell_collar_word"], -2496)
        self.assertEqual(controls["collar_u_collar_ell_collar_word"], 2624)

    def test_committed_files_match(self) -> None:
        committed = json.loads(self.b.CONVENTION.read_text(encoding="utf-8"))
        self.assertEqual(committed, self.convention)
        audit = json.loads(self.b.AUDIT.read_text(encoding="utf-8"))
        volatile = {"elapsed_seconds", "audit_sha256"}
        self.assertEqual({k: v for k, v in audit.items() if k not in volatile}, {k: v for k, v in self.audit.items() if k not in volatile})

    def test_mutations_are_detected(self) -> None:
        verify = load("verify_t73_endpoint_transport")
        detected = verify.mutation_tests(self.b, self.convention, self.audit)
        self.assertTrue(all(detected.values()), detected)

    def test_wrong_pivotal_coefficient_fails_letterwise(self) -> None:
        rec = self.b.load_recompute()
        data = json.loads(self.b.PUBLIC_INPUT.read_text(encoding="utf-8"))
        b44, _ = rec.build_oriented_b44(data)
        b88 = rec.cable_word(b44)
        bad = self.b.letterwise_transport(
            b88[:2000], self.convention["endpoints"], lambda p: self.b.LP.q(-p), {self.b.V: self.b.LP.one(), self.b.VD: self.b.LP.one()}
        )
        self.assertFalse(bad["all_letters_pass"])


if __name__ == "__main__":
    unittest.main()
