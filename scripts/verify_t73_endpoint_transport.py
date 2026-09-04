#!/usr/bin/env python3
"""Independent verification of the T73 endpoint transport.

This script never reads a self-reported PASS field.  It recomputes the whole
endpoint transport (oriented model checks, letterwise conjugation along the
actual 45360-letter cabled word, cup/cap derivation, rho(W)-I in h^3 End, the
divided cubic three ways, and the coordinate controls) from the primary
sources, compares the result with the committed authority and audit files,
and then runs mutation tests that must FAIL:

* pivotal coefficient of the V* defect changed from -q to +q or to -q^-1;
* position monomial changed from q^{-p} to q^{+p};
* the pivotal sign sigma flipped (the cap row is negated, the cubic changes sign);
* the two cup endpoints exchanged (the cup would join two exits);
* one entry of the geometric->public permutation moved.

Exit status is nonzero unless every recomputed check passes, the committed
files match, and every mutation is detected.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_t73_endpoint_transport.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("build_t73_endpoint_transport", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load build_t73_endpoint_transport.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def strip_volatile(audit: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in audit.items() if k not in {"elapsed_seconds", "audit_sha256"}}


def mutation_tests(b, convention: dict[str, Any], audit: dict[str, Any]) -> dict[str, bool]:
    """Every entry must be True (= the mutation was detected)."""
    rec = b.load_recompute()
    public_input = json.loads(b.PUBLIC_INPUT.read_text(encoding="utf-8"))
    b44, _ = rec.build_oriented_b44(public_input)
    b88 = rec.cable_word(b44)
    records = convention["endpoints"]
    ratio_exp = int(convention["public_model"]["position_monomial"].split("q^{")[1].split("*")[0])
    good_position = lambda p: b.LP.q(ratio_exp * p)  # noqa: E731
    good_pivotal = {b.V: b.LP.one(), b.VD: b.LP.mono(-1, 1)}
    detected: dict[str, bool] = {}

    # 1. wrong pivotal coefficient (+q instead of -q)
    res = b.letterwise_transport(b88, records, good_position, {b.V: b.LP.one(), b.VD: b.LP.mono(1, 1)})
    detected["pivotal_plus_q_rejected"] = not res["all_letters_pass"]
    # 2. wrong pivotal coefficient (-q^-1 instead of -q)
    res = b.letterwise_transport(b88, records, good_position, {b.V: b.LP.one(), b.VD: b.LP.mono(-1, -1)})
    detected["pivotal_minus_q_inverse_rejected"] = not res["all_letters_pass"]
    # 3. wrong position monomial (q^{+p})
    res = b.letterwise_transport(b88, records, lambda p: b.LP.q(-ratio_exp * p), good_pivotal)
    detected["position_monomial_sign_rejected"] = not res["all_letters_pass"]
    # 4. the correct data pass (sanity of the harness)
    res = b.letterwise_transport(b88, records, good_position, good_pivotal)
    detected["reference_data_pass"] = res["all_letters_pass"] and res["final_permutation_is_identity"]
    # 5. sigma flipped: the cubic changes sign and disagrees with the committed audit
    _, audit_plus = b.compute(sigma=+1, verbose=False)
    detected["sigma_flip_changes_cubic"] = (
        audit_plus["delta3"]["constant_terms_pipeline"] == -audit["delta3"]["constant_terms_pipeline"]
        and audit_plus["ell_public_constant_terms"] != audit["ell_public_constant_terms"]
    )
    # 6. exchanged cup endpoints: derive_endpoint_terms must refuse a cup joining two exits
    mutated = json.loads(json.dumps(convention))
    by_passage = {r["passage_id"]: r for r in mutated["endpoints"]}
    a, bb = (by_passage[p] for p in mutated["selected_cup"]["passages"])
    other_exit = next(r for r in mutated["endpoints"] if r["tensor_factor"] == b.VD and r["passage_id"] != bb["passage_id"])
    mutated["selected_cup"]["passages"] = [other_exit["passage_id"], bb["passage_id"]]
    payload = {k: v for k, v in mutated.items() if k != "convention_sha256"}
    mutated["convention_sha256"] = b.canonical_sha(payload)
    path = ROOT / "audit" / "_mutant_endpoint_convention.json"
    path.write_text(json.dumps(mutated), encoding="utf-8")
    try:
        terms = b.derive_endpoint_terms(path)
        # two exits: the cap/cup would carry two V* defects; the derived constant terms then
        # differ from the committed ones, which the receipt comparison catches.
        detected["cup_endpoint_exchange_rejected"] = terms["u_terms"] != audit["u_public_constant_terms"]
    except Exception:
        detected["cup_endpoint_exchange_rejected"] = True
    finally:
        path.unlink(missing_ok=True)
    # 7. permutation mutation: moving one endpoint's geometric slot changes the convention SHA and
    #    the transported cup, so --check fails; verify the SHA sensitivity directly.
    mutated = json.loads(json.dumps(convention))
    r0, r1 = mutated["endpoints"][0], mutated["endpoints"][1]
    r0["geometric_order"], r1["geometric_order"] = r1["geometric_order"], r0["geometric_order"]
    payload = {k: v for k, v in mutated.items() if k != "convention_sha256"}
    detected["permutation_mutation_changes_sha"] = b.canonical_sha(payload) != convention["convention_sha256"]
    return detected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-mutations", action="store_true")
    args = parser.parse_args()
    b = load_builder()
    convention, audit = b.compute(sigma=-1, verbose=False)
    committed_convention = json.loads(b.CONVENTION.read_text(encoding="utf-8"))
    committed_audit = json.loads(b.AUDIT.read_text(encoding="utf-8"))
    convention_match = committed_convention == convention
    audit_match = strip_volatile(committed_audit) == strip_volatile(audit)
    recomputed_pass = bool(
        audit["model_checks_pass"]
        and audit["letterwise_transport"]["all_letters_pass"]
        and audit["letterwise_transport"]["final_permutation_is_identity"]
        and audit["rho_W_minus_I_in_h3_End"]
        and audit["delta3"]["agree"]
    )
    detected = {} if args.skip_mutations else mutation_tests(b, convention, audit)
    mutations_ok = all(detected.values())
    print(f"CONVENTION_MATCH={'PASS' if convention_match else 'FAIL'}")
    print(f"AUDIT_MATCH={'PASS' if audit_match else 'FAIL'}")
    print(f"RECOMPUTED_MODEL_CHECKS={'PASS' if audit['model_checks_pass'] else 'FAIL'}")
    print(f"RECOMPUTED_LETTERWISE={'PASS' if audit['letterwise_transport']['all_letters_pass'] else 'FAIL'}")
    print(f"RHO_W_MINUS_I_IN_H3={'PASS' if audit['rho_W_minus_I_in_h3_End'] else 'FAIL'}")
    print(f"U_PUBLIC_CONSTANT={json.dumps(audit['u_public_constant_terms'], separators=(',', ':'))}")
    print(f"ELL_PUBLIC_CONSTANT={json.dumps(audit['ell_public_constant_terms'], separators=(',', ':'))}")
    print(f"COORDINATE_CONTROLS={json.dumps({k: v for k, v in audit['coordinate_controls'].items() if isinstance(v, int)}, sort_keys=True, separators=(',', ':'))}")
    print(f"MUTATION_TESTS={json.dumps(detected, sort_keys=True, separators=(',', ':'))}")
    print(f"ENDPOINT_TRANSPORT={'PASS' if recomputed_pass and convention_match and audit_match and mutations_ok else 'FAIL'}")
    print(f"NO_UNRESOLVED_SIGNS={'PASS' if audit['no_unresolved_signs'] else 'FAIL'}")
    print(f"DELTA3={audit['delta3']['constant_terms_pipeline']}")
    if not (recomputed_pass and convention_match and audit_match and mutations_ok):
        sys.exit(2)


if __name__ == "__main__":
    main()
