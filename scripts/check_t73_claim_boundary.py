#!/usr/bin/env python3
"""Check that the controlling paper keeps an honest claim boundary.

The English manuscript states the geometric inputs P0, C, S and P3 for the
Johnson replacement as hypotheses and records their status in a status
table.  The checker enforces that:

* P0, C and S are stated as hypotheses whose status is OPEN, not as proved
  theorems;
* the status table rows agree with the expected statuses below;
* the conditional top-level theorem and the Lean boundary sentences are
  present;
* retired or false statements (HJ lemmas "do not appear", undetermined
  endpoint signs, split-unknot inference, "1--3 cancellations") are absent.

Every expected status must be updated deliberately when a later commit
actually establishes an input; the checker reads the paper text, never a
self-reported PASS field.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "spc4-t73-candidate" / "main.tex"
PUBLISHED = ROOT / "paper" / "spc4-t73-candidate" / "sec-published-results.tex"

# Status table rows (label as printed in the paper -> expected status macro).
EXPECTED_STATUS = {
    "P0a (handlebody bridge)": r"\Open",
    "P0b (two framed cancellations)": r"\Open",
    "P0c (MWW cabling framing)": r"\Open",
    "P0d (finite word)": r"\Discharged",
    "C1 (coefficient bimodule)": r"\Open",
    "C2 (statewise cocone)": r"\Open",
    "C3": r"\Unused",
    "S (sphere system, hemisphere maps)": r"\Open",
    "P2/E7": r"\Unused",
    "P3/E11": r"\Open",
    "P3/E12": r"\Discharged",
    "P3/E13": r"\Open",
    "Finite detector": r"\Discharged",
    r"Lean \texttt{ExternalGeometry}": r"\Open",
}

HYPOTHESES = (
    r"\begin{hypothesis}[Embedded candidate presentation, P0]\label{hyp:P0}",
    r"\begin{hypothesis}[Coefficient comparison, C]\label{hyp:P1}",
    r"\begin{hypothesis}[Relative three-handle closure, S]\label{hyp:P2}",
    r"\begin{hypothesis}[Four-handle closure of the replacement picture, P3]\label{hyp:P3}",
)

OPEN_PROPOSITIONS = (
    r"\begin{proposition}[Status of P0 for the Johnson replacement]\label{thm:P0discharge}",
    r"\begin{proposition}[Status of C]\label{thm:Cdischarge}",
    r"\begin{proposition}[Status of S]\label{thm:Sdischarge}",
)


def require(text: str, needle: str, source: Path) -> None:
    if needle not in text:
        raise AssertionError(f"{source}: required claim-boundary text missing: {needle!r}")


def reject(text: str, needle: str, source: Path) -> None:
    if needle in text:
        raise AssertionError(f"{source}: forbidden unconditional claim present: {needle!r}")


def status_rows(text: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in text.splitlines():
        if " & " not in line or not line.rstrip().endswith(r"\\"):
            continue
        cells = [cell.strip() for cell in line.split(" & ")]
        if len(cells) >= 2 and cells[0] in EXPECTED_STATUS:
            rows[cells[0]] = cells[1]
    return rows


def check() -> None:
    paper_text = PAPER.read_text(encoding="utf-8")
    published_text = PUBLISHED.read_text(encoding="utf-8")

    # Conditional top-level theorem (not an unconditional Trace-73 theorem).
    require(paper_text, r"\begin{theorem}[Conditional trace-73 theorem]\label{thm:joined}", PAPER)
    require(paper_text, "That interface is not constructed in Lean, so no counterexample is claimed.", PAPER)
    require(paper_text, r"that assembly is \Open\ for the Johnson candidate", PAPER)
    require(paper_text, r"Hypotheses~\ref{hyp:P0}--\ref{hyp:P2} are \Open.", PAPER)

    # P0/C/S/P3 are hypotheses with OPEN status propositions, not theorems.
    for marker in HYPOTHESES + OPEN_PROPOSITIONS:
        require(paper_text, marker, PAPER)
    for marker in (
        r"\begin{theorem}[Embedded candidate presentation, P0]",
        r"\begin{theorem}[Coefficient comparison, C]",
        r"\begin{theorem}[Relative three-handle closure, S]",
        r"\begin{theorem}[Four-handle closure of the replacement picture, P3]",
        r"\begin{theorem}[Johnson replacement discharges P0]",
        r"\begin{theorem}[Coefficient comparison C]",
        r"\begin{theorem}[Relative three-handle closure S]",
    ):
        reject(paper_text, marker, PAPER)

    # Status table rows.
    rows = status_rows(paper_text)
    for label, expected in EXPECTED_STATUS.items():
        actual = rows.get(label)
        if actual != expected:
            raise AssertionError(
                f"{PAPER}: status row {label!r} is {actual!r}, expected {expected!r}"
            )

    # Labels used by the claim map and the Lean boundary appendix.
    for marker in (
        r"\label{thm:P0discharge}",
        r"\label{lem:P0d-link}",
        r"\label{lem:C1}",
        r"\label{lem:filtered-cubic}",
        r"\label{sec:local-stabilization}",
        r"\label{sec:status-table}",
        r"\label{thm:Cdischarge}",
        r"\label{thm:Sdischarge}",
        r"\label{hyp:P3}",
        r"\label{sec:final-identifications}",
        r"Not the E7",
    ):
        require(paper_text, marker, PAPER)

    # Finite / cited layers that must remain visible in the manuscript.
    require(paper_text, "2624", PAPER)
    require(paper_text, "494", PAPER)
    require(paper_text, r"\cite[Proposition~3.4 and Corollary~3.5]{ManolescuWalkerWedrich2023}", PAPER)
    require(paper_text, "cited not formalized", PAPER)
    require(paper_text, r"scripts/recompute\_t73\_delta3.py --check", PAPER)
    require(paper_text, r"certify\_t73\_e12\_s4.py", PAPER)
    require(paper_text, "frozen Burau computation", PAPER)

    # Horvat--Jablonowski: the lemmas exist in the current version; only
    # Theorem 5.3 is invoked, never the relative Lemma 5.7 to fix B.
    require(paper_text, "Lemma~5.5 (spotted-ball slide lemma) and Lemma~5.7 (relative uniqueness of", PAPER)
    require(paper_text, "does not invoke Lemma~5.7 in its relative form", PAPER)
    require(published_text, "Lemma~5.5 (spotted-ball slide lemma) and Lemma~5.7", PUBLISHED)
    for source, text in ((PAPER, paper_text), (PUBLISHED, published_text)):
        reject(text, "do not appear in", source)
        reject(text, "do not\nappear in that source", source)

    # One-directional Gamma_3 implication; no undetermined endpoint signs;
    # no split-unknot inference; no 1--3 cancellation wording.
    require(paper_text, r"W\in\Gamma_3(P_{44})\ \Longrightarrow", PAPER)
    require(paper_text, "the converse is not asserted", PAPER)
    reject(paper_text, "equivalently, and also", PAPER)
    reject(paper_text, r"\pm e_5", PAPER)
    reject(paper_text, r"\pm e_2^*", PAPER)
    reject(paper_text, "split-unknot\nfactor", PAPER)
    reject(paper_text, "split\n$r_{zx}$ component is retained", PAPER)
    reject(paper_text, "$1$--$3$ cancellations", PAPER)
    reject(paper_text, "1--3 cancellation", PAPER)

    # Forbidden: absent-ledger slogans or an unconditional joined theorem /
    # counterexample claim.
    reject(paper_text, r"\begin{theorem}[Trace-73 theorem]\label{thm:joined}", PAPER)
    reject(paper_text, "gives a counterexample", PAPER)
    reject(paper_text, "The actual reduced PD ledger is absent", PAPER)
    reject(paper_text, "that object is absent", PAPER)
    reject(paper_text, "we have proved Hypotheses", PAPER)
    reject(paper_text, "we prove the geometric\ninputs P0, C, S, and P3", PAPER)


def main() -> None:
    check()
    print("T73_CLAIM_BOUNDARY=OPEN_GEOMETRY")
    print("STATUS=" + ",".join(f"{k}:{v.lstrip(chr(92))}" for k, v in EXPECTED_STATUS.items()))
    print(f"PAPER={PAPER}")


if __name__ == "__main__":
    main()
