#!/usr/bin/env python3
"""Audit the proof state of every load-bearing premise in the T73 paper.

Geometric items are Open until the named acceptance tests pass.  This program
does not hardcode proved: true.

Two layers are reported for every item.  ``state`` is the certificate-internal
verdict: the committed finite models (control strands, local model movies,
model spheres) replay and their generators report PASS.  ``paper_status`` is
the claim boundary of the controlling paper, read from the status table via
``check_t73_claim_boundary.EXPECTED_STATUS``: an input is OPEN there until the
actual Cappell--Shaneson geometry is constructed, whatever the finite model
certificates say.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMMITTED = ROOT / "audit" / "t73_premise_audit.json"

# Which status-table rows of the paper decide the paper status of each item.
PAPER_ROWS = {
    "P0": ("P0a (handlebody bridge)", "P0b (two framed cancellations)", "P0c (MWW cabling framing)", "P0d (finite word)"),
    "C": ("C1 (coefficient bimodule)", "C2 (statewise cocone)"),
    "S": ("S (sphere system, hemisphere maps)",),
    "P3_E11": ("P3/E11",),
    "P3_E12": ("P3/E12",),
    "P3_E13": ("P3/E13",),
}


def paper_statuses() -> dict[str, dict[str, Any]]:
    """Paper status per item from the claim-boundary checker (which reads the paper)."""
    path = ROOT / "scripts" / "check_t73_claim_boundary.py"
    spec = importlib.util.spec_from_file_location("check_t73_claim_boundary", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load check_t73_claim_boundary.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.check()  # the paper must actually carry the expected rows
    expected = module.EXPECTED_STATUS
    result: dict[str, dict[str, Any]] = {}
    for item, rows in PAPER_ROWS.items():
        row_status = {row: expected[row].lstrip("\\") for row in rows}
        result[item] = {
            "paper_status": "PROVED" if all(v == "Discharged" for v in row_status.values()) else "OPEN",
            "paper_rows": row_status,
        }
    return result


def require(text: str, needle: str, source: Path) -> None:
    if needle not in text:
        raise AssertionError(f"{source}: missing premise marker {needle!r}")


def generate() -> dict[str, Any]:
    paper = ROOT / "paper" / "spc4-t73-candidate" / "main.tex"
    paper_text = paper.read_text(encoding="utf-8")
    completion = ROOT / "docs" / "research" / "T73_COMPLETION_AUDIT_2026-09-02.md"
    completion_text = completion.read_text(encoding="utf-8")
    p0_certificate = json.loads(
        (ROOT / "audit" / "t73_p0_johnson_certificate.json").read_text(encoding="utf-8")
    )
    c_witness = json.loads(
        (ROOT / "audit" / "t73_c_comparison_witness.json").read_text(encoding="utf-8")
    )
    s_certificate = json.loads(
        (ROOT / "audit" / "t73_s_relative_moves_certificate.json").read_text(encoding="utf-8")
    )
    p3_certificate = json.loads(
        (ROOT / "audit" / "t73_p3_four_handle.json").read_text(encoding="utf-8")
    )
    e12_certificate = json.loads(
        (ROOT / "audit" / "t73_e12_s4_reduction.json").read_text(encoding="utf-8")
    )
    e13_close = json.loads(
        (ROOT / "audit" / "t73_e13_close.json").read_text(encoding="utf-8")
    )
    e13_identification = json.loads(
        (ROOT / "audit" / "t73_e13_identification.json").read_text(encoding="utf-8")
    )

    for marker in (
        r"\label{thm:P0discharge}",
        r"\label{lem:P0d-link}",
        r"\label{lem:C1}",
        r"\label{thm:Cdischarge}",
        r"\label{thm:Sdischarge}",
        r"\label{hyp:P3}",
        r"that assembly is \Open\ for the Johnson candidate",
        r"\label{sec:final-identifications}",
    ):
        require(paper_text, marker, paper)
    for marker in (
        "| P0 | **PASS**",
        "| C | **PASS**",
        "| S | **PASS**",
        "| P3/E11 | **PASS**",
        "| P3/E12 | **PASS**",
        "| P3/E13 | **PASS**",
    ):
        require(completion_text, marker, completion)

    p0_pass = (
        p0_certificate.get("verdict") == "PASS"
        and p0_certificate.get("P0_status") != "OPEN"
    )
    if p0_pass and p0_certificate.get("checks", {}).get("johnson_ar_affine_bridge") is not True:
        raise AssertionError("P0 cannot be PASS without the Johnson--AR handlebody bridge")
    c_pass = (
        c_witness.get("C_status") == "PASS"
        and c_witness.get("C1_status") == "PASS"
        and c_witness.get("C2_status") == "PASS"
    )
    if c_pass and not p0_pass:
        raise AssertionError("C cannot be PASS without P0")
    s_pass = s_certificate.get("verdict") == "PASS"
    if s_pass and not (p0_pass and c_pass):
        raise AssertionError("S cannot be PASS without P0 and C")
    if s_pass and not s_certificate.get("checks", {}).get("detector_fixed"):
        raise AssertionError("S cannot be PASS unless the detector ball is fixed")
    if s_pass and not s_certificate.get("checks", {}).get("actual_attaching_system_identified"):
        raise AssertionError("S cannot be PASS without an identified attaching system")
    p3_e11_pass = (
        p3_certificate.get("verdict") == "PASS"
        and "PASS" in str(p3_certificate.get("E11_status", ""))
        and p3_certificate.get("closed_manifold", {}).get("identified_with_Sigma_A_0") is False
        and p3_certificate.get("four_handle", {}).get("triangulated_W3") is False
    )
    if p3_e11_pass and not s_pass:
        raise AssertionError("P3/E11 cannot be PASS without S")
    p3_e12_pass = (
        p3_certificate.get("E12_status") == "PASS"
        and p3_certificate.get("e12_s4", {}).get("about_standard_S4_not_candidate")
        and e12_certificate.get("verdict") == "PASS"
        and e12_certificate.get("checks", {}).get("s4_degree_494_zero") is True
        and e12_certificate.get("about_standard_S4_not_candidate") is True
        and e12_certificate.get("lean_s4_reduction_data_inhabited") is False
        and e12_certificate.get("identified_with_X_J") is False
        and e12_certificate.get("p3_certificate_sha256") == p3_certificate.get("certificate_sha256")
    )
    p3_e13_pass = (
        e13_close.get("verdict") == "IDENTIFIED_CS_HANDLE_PICTURE"
        and e13_close.get("E13_status") == "PASS"
        and e13_close.get("checks", {}).get("identified_with_Sigma_A_0") is True
        and e13_close.get("checks", {}).get("lean_cs_topology_data_inhabited") is False
        and e13_close.get("checks", {}).get("uniqueness_of_regular_neighborhoods_used") is False
        and e13_close.get("p3_certificate_sha256") == p3_certificate.get("certificate_sha256")
        and e13_identification.get("verdict") == "IDENTIFIED_CS_HANDLE_PICTURE"
        and e13_identification.get("checks", {}).get("identified_with_Sigma_A_0") is True
        and e13_identification.get("e13_close_sha256") == e13_close.get("certificate_sha256")
        and p3_certificate.get("closed_manifold", {}).get("identified_with_Sigma_A_0") is False
        and p3_certificate.get("e13_determinants", {}).get("identifies_X_J_with_Sigma_A_0")
        is False
    )
    items = {
        "P0": {
            "state": "PASS" if p0_pass else "OPEN",
            "proved": p0_pass,
            "falsified": False,
            "evidence": [
                "scripts/certify_t73_johnson_ar_bridge.py",
                "scripts/reconstruct_t73_p0.py",
                "scripts/certify_t73_spine_star_handlebodies.py",
                "audit/t73_p0a_handlebody_pair.json",
                "scripts/check_t73_p0_pipeline.py",
                "scripts/falsify_t73_linking_from_words.py",
                "scripts/search_t73_johnson_alpha_sides.py",
                "docs/proofs/T73_GAP_FREE_BASIS_RECEIPT.md",
                "tests/test_t73_p0_reconstruction.py",
                "audit/t73_p0_johnson_certificate.json",
                "docs/proofs/T73_GEOMETRIC_PREMISE_STATUS_20260903.md",
            ],
            "blocker": (
                "none"
                if p0_pass
                else "P0 reconstruction input is missing or reconstruct_t73_p0.py rejected it."
            ),
            "control": "The synthetic rational 44-strand control recovers the public word but is deliberately not AR-bound.",
            "falsified_route": "The current 19-step Nielsen representative is falsified as the source of the public 44-channel collar; only that retired route remains rejected.",
            "falsified_compact_lift": "GAP proves the three compact straight spine words are injective but not surjective on F3, so they cannot be the images of a handlebody homeomorphism.",
            "replacement_candidate": "Inner conjugation by x^-1 gives a GAP-certified F3 automorphism with 44 channels, but its exact compact word and geometric owner/framing movie remain open.",
            "replacement_adjudication": "The x^-1 correction is simultaneous inner conjugation, hence only a basepoint change in Out(F3); its two extra word passages are not an embedded 44-channel witness.",
            "johnson_candidate": (
                "The 93-bit Johnson alpha-side lift is a GAP free basis, "
                "matches compact m2, has 44 y-channels, and the six-sweep "
                "reconstruction recovers the public 11340-letter word."
            ),
            "certificate_sha256": p0_certificate["certificate_sha256"],
            "p0a_status": "PASS" if p0_pass else "OPEN",
            "p0b_status": "PASS" if p0_pass else "OPEN",
            "p0c_status": "PASS" if p0_pass else "OPEN",
            "p0d_finite_word_match": True,
        },
        "C": {
            "state": "PASS" if c_pass else "OPEN",
            "proved": c_pass,
            "falsified": False,
            "generator_internal_status": c_witness.get("C_status"),
            "evidence": [
                "scripts/generate_t73_c_comparison_witness.py",
                "scripts/certify_t73_c1_cut_link.py",
                "scripts/certify_t73_c2_comparison.py",
                "audit/t73_c_comparison_witness.json",
                "audit/t73_c1_cut_link.json",
                "audit/t73_c2_comparison.json",
                "Smooth4PC/RepresentableCoefficient.lean",
            ],
            "blocker": (
                "none"
                if c_pass
                else "C1 requires collar-bound product rectangles and C2 comparison maps."
            ),
            "adjudication": (
                "Johnson replacement: C1 product isotopy of P0 reconstruction "
                "strands and C2 action cubes plus RepresentableCoefficient.lean."
            ),
            "certificate_sha256": c_witness["witness_sha256"],
        },
        "S": {
            "state": "PASS" if s_pass else "OPEN",
            "proved": s_pass,
            "falsified": False,
            "generator_internal_status": s_certificate.get("verdict"),
            "evidence": [
                "scripts/certify_t73_s_relative_moves.py",
                "scripts/certify_t73_s_standard_spheres.py",
                "audit/t73_s_relative_moves_certificate.json",
                "audit/t73_s_standard_spheres.json",
            ],
            "relative_geometry_proved": s_pass,
            "blocker": (
                "none"
                if s_pass
                else (
                    "No B-fixing reversed 1-handle picture, or "
                    "actual_standard_sphere_endpoint_foam_computed is false."
                )
            ),
            "adjudication": (
                "Johnson replacement reversed 3-handle picture: belt spheres miss "
                "the P0 cube and the C1 leftover link; HJ Theorem 5.3 is used only "
                "for kernel invariance, not to fix B."
                if s_pass
                else "Paper Lemmas Ssystem and Sendpoint remain slogans."
            ),
            "certificate_sha256": s_certificate["certificate_sha256"],
        },
        "P3_E11": {
            "state": "PASS" if p3_e11_pass else "OPEN",
            "proved": p3_e11_pass,
            "falsified": False,
            "blocker": (
                "none"
                if p3_e11_pass
                else "MWW Proposition 3.4 requires the Johnson replacement 4-handle picture"
            ),
            "adjudication": (
                "Johnson replacement 1-3 cancellations and a PL 4-ball; MWW 3.4 on X_J, not Sigma_A^0."
                if p3_e11_pass
                else "MWW Proposition 3.4 applies only after the four-handle picture exists."
            ),
            "evidence": [
                "scripts/certify_t73_p3_four_handle.py",
                "audit/t73_p3_four_handle.json",
                "MWW Proposition 3.4",
            ],
            "certificate_sha256": p3_certificate.get("certificate_sha256"),
        },
        "P3_E12": {
            "state": "PASS" if p3_e12_pass else "CITED_EXTERNAL",
            "proved": p3_e12_pass,
            "falsified": False,
            "blocker": "none" if p3_e12_pass else "S^4 reduction certificate missing or not about standard S^4",
            "evidence": [
                "scripts/certify_t73_e12_s4.py",
                "audit/t73_e12_s4_reduction.json",
                "scripts/certify_t73_p3_four_handle.py",
                "Smooth4PC/T73S4Control.lean",
                "Smooth4PC/T73S4Inhabitant.lean",
                "MWW Corollary 3.5 unpacked as empty Khovanov and two I^4 glued along S^3",
            ],
            "certificate_sha256": e12_certificate.get("certificate_sha256"),
        },
        "P3_E13": {
            "state": "PASS" if p3_e13_pass else "OPEN",
            "proved": p3_e13_pass,
            "falsified": False,
            "blocker": (
                "none"
                if p3_e13_pass
                else "E13 close certificate missing or does not identify X_J with Sigma_A^0"
            ),
            "evidence": [
                "Smooth4PC/T73Finite.lean",
                "Smooth4PC/T73JohnsonTransvections.lean",
                "Smooth4PC/T73GeometryPack.lean",
                "Iwaki Proposition 2.1 for the matrix criterion",
                "audit/t73_p3_four_handle.json",
                "scripts/certify_t73_e13_close.py",
                "audit/t73_e13_close.json",
                "audit/t73_reduced_link_pd.json",
                "scripts/certify_t73_e13_identification.py",
                "audit/t73_e13_identification.json",
            ],
            "certificate_sha256": e13_close.get("certificate_sha256"),
        },
    }
    paper = paper_statuses()
    for name, entry in items.items():
        entry["paper_status"] = paper[name]["paper_status"]
        entry["paper_rows"] = paper[name]["paper_rows"]
        entry["state_meaning"] = (
            "certificate-internal replay of the committed finite model; "
            "not a statement about the actual Cappell--Shaneson geometry"
        )
    return {
        "schema": "t73_premise_audit/v2",
        "overall": "OPEN",
        "counterexample_claim_proved": False,
        "counterexample_claim_falsified": False,
        "items": items,
        "interpretation": (
            "The committed finite-model certificates for P0, C, S, the MWW "
            "four-handle layer and the E13 CS handle picture replay and report "
            "PASS internally. The controlling paper records P0, C, S, P3/E11 "
            "and P3/E13 as OPEN (actual geometry not constructed) and P0a, P0d, "
            "P3/E12 and the finite detector as discharged. Lean ExternalGeometry "
            "remains uninhabited. The empty-link control S4ReductionData is "
            "inhabited in T73S4Inhabitant.lean. The counterexample claim is not "
            "proved."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    generated = generate()
    if args.write:
        COMMITTED.write_text(
            json.dumps(generated, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(f"WROTE={COMMITTED}")
    if args.check:
        committed = json.loads(COMMITTED.read_text(encoding="utf-8"))
        if committed != generated:
            raise AssertionError("committed premise audit differs from regenerated audit")
        print("T73_PREMISE_AUDIT=OPEN")
        print(f"OVERALL={generated['overall']}")
        print(f"P0={generated['items']['P0']['state']}")
        print(f"C={generated['items']['C']['state']}")
        print(f"S={generated['items']['S']['state']}")
        print(f"P3_E11={generated['items']['P3_E11']['state']}")
        print(f"P3_E12={generated['items']['P3_E12']['state']}")
        print(f"P3_E13={generated['items']['P3_E13']['state']}")
        print(
            "PAPER_STATUS="
            + ",".join(f"{name}:{entry['paper_status']}" for name, entry in generated["items"].items())
        )
        print(f"COUNTEREXAMPLE={generated['counterexample_claim_proved']}")
        return
    if not args.write:
        print(json.dumps(generated, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
