from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


AXIOM_THEOREMS = (
    "Smooth4PC.le_add_single",
    "Smooth4PC.reynolds_apply",
    "Smooth4PC.reynolds_comp_of_row_perm",
    "Smooth4PC.card_placement_of_uniform_fibres",
    "Smooth4PC.reynolds_comp_dotted",
    "Smooth4PC.reynolds_comp_undotted",
    "Smooth4PC.CableCocone.Psi_step",
    "Smooth4PC.CableCocone.Psi_comp_of_eq",
    "Smooth4PC.CableCocone.Psi_comp_aux",
    "Smooth4PC.CableCocone.Psi_comp",
    "Smooth4PC.CableCocone.Psi_two_steps",
    "Smooth4PC.row_comp_Psi_of_eq",
    "Smooth4PC.row_descent_along_path",
    "Smooth4PC.extendedRow_eq_comp_of_le",
    "Smooth4PC.extendedRow_eq_of_threshold",
    "Smooth4PC.extendedRow_comp_psi",
)
FOUNDATIONAL_AXIOMS = {"propext", "Quot.sound", "Classical.choice"}
FORBIDDEN_TOKENS = (
    "sorry",
    "admit",
    "axiom",
    "constant",
    "opaque",
    "unsafe",
    "extern",
    "implemented_by",
    "run_tac",
)


def resolve_lake() -> Path:
    configured = os.environ.get("T73_LAKE")
    if configured:
        candidate = Path(configured).expanduser()
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(f"T73_LAKE does not name a file: {candidate}")

    on_path = shutil.which("lake")
    if on_path:
        return Path(on_path)

    toolchain_bin = (
        Path.home()
        / ".elan"
        / "toolchains"
        / "leanprover--lean4---v4.32.1"
        / "bin"
    )
    for name in ("lake.exe", "lake"):
        candidate = toolchain_bin / name
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("lake was not found via T73_LAKE, PATH, or ~/.elan")


def lake_printenv(lake: Path, repo: Path, name: str) -> str:
    # `lake env printenv NAME` reports the search path the workspace would hand
    # to Lean.  Reading it is more robust than re-deriving the package layout
    # from lake-manifest.json, because the mathlib checkout may live outside
    # this repository.
    result = subprocess.run(
        [str(lake), "env", "printenv", name],
        cwd=repo,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"lake env printenv {name} failed: {result.stdout + result.stderr}"
        )
    return result.stdout.strip()


def resolve_lean(lake: Path, repo: Path) -> Path:
    # Lean is invoked directly rather than through `lake env lean`, because the
    # throwaway olean root has to come FIRST in LEAN_PATH and `lake env` only
    # ever appends the inherited LEAN_PATH after its own entries.  Lean commits
    # to the first search-path entry that holds a directory named after the
    # module root (`Smooth4PC`) and reports the olean as missing there instead
    # of continuing the search, so once `.lake/build/lib/lean/Smooth4PC/` exists
    # for any other module of this package, a scratch root appended at the end
    # can never be reached.
    configured = os.environ.get("T73_LEAN")
    if configured:
        candidate = Path(configured).expanduser()
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(f"T73_LEAN does not name a file: {candidate}")

    sysroot = Path(lake_printenv(lake, repo, "LEAN_SYSROOT")).expanduser()
    for name in ("lean.exe", "lean"):
        candidate = sysroot / "bin" / name
        if candidate.is_file():
            return candidate

    on_path = shutil.which("lean")
    if on_path:
        return Path(on_path)
    raise FileNotFoundError("lean was not found via T73_LEAN, LEAN_SYSROOT, or PATH")


def resolve_scratch() -> Path:
    # Directory holding the throwaway olean root for one test run.
    # Deliberately NOT the system temporary directory: under WSL, /tmp is
    # cleared when the distribution restarts after an idle period, which can
    # delete the olean root between the build step and the audit step of a
    # single test and make the audit fail for a reason unrelated to the
    # proofs.  T73_TMP still overrides the location.
    configured = os.environ.get("T73_TMP")
    scratch = (
        Path(configured).expanduser()
        if configured
        else Path.home() / "ws" / "tmp-t73" / "lean2"
    )
    scratch.mkdir(parents=True, exist_ok=True)
    return scratch


class ReynoldsCableCoconeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo = Path(__file__).resolve().parents[1]
        cls.module = cls.repo / "Smooth4PC" / "ReynoldsCableCocone.lean"
        cls.audit = cls.repo / "T73ReynoldsAudit.lean"
        cls.lake = resolve_lake()
        cls.lean = resolve_lean(cls.lake, cls.repo)
        cls.package_lean_path = lake_printenv(cls.lake, cls.repo, "LEAN_PATH")

    def lean_environment(self, olean_root: Path) -> dict[str, str]:
        env = os.environ.copy()
        env.pop("LEAN_PATH", None)
        env.pop("LEAN_SRC_PATH", None)
        paths = [str(olean_root)]
        paths.extend(
            entry for entry in self.package_lean_path.split(os.pathsep) if entry
        )
        env["LEAN_PATH"] = os.pathsep.join(paths)
        return env

    def run_lean(
        self, source: Path, olean_root: Path, output: Path | None = None
    ) -> subprocess.CompletedProcess[str]:
        args = [str(self.lean)]
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
            args.extend(["-o", str(output)])
        args.append(str(source))
        return subprocess.run(
            args,
            cwd=self.repo,
            env=self.lean_environment(olean_root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )

    def assert_axiom_reports(self, output: str) -> None:
        reports = re.findall(
            r"(?m)^'([^']+)' depends on axioms:\s*\[([^\]]*)\]$", output
        )
        self.assertEqual(len(reports), len(AXIOM_THEOREMS), "missing axiom report")
        self.assertEqual({name for name, _ in reports}, set(AXIOM_THEOREMS))
        for _, payload in reports:
            names = {name.strip() for name in payload.split(",") if name.strip()}
            self.assertLessEqual(
                names, FOUNDATIONAL_AXIOMS, f"unexpected axioms: {sorted(names)}"
            )

    def test_axiom_report_gate_rejects_a_synthetically_removed_line(self) -> None:
        complete = "\n".join(
            f"'{name}' depends on axioms: [propext]" for name in AXIOM_THEOREMS
        )
        missing = complete.replace(complete.splitlines()[0] + "\n", "", 1)
        with self.assertRaisesRegex(AssertionError, "missing axiom report"):
            self.assert_axiom_reports(missing)

    def test_axiom_report_gate_rejects_a_nonfoundational_axiom(self) -> None:
        polluted = "\n".join(
            f"'{name}' depends on axioms: [propext, sorryAx]"
            for name in AXIOM_THEOREMS
        )
        with self.assertRaisesRegex(AssertionError, "unexpected axioms"):
            self.assert_axiom_reports(polluted)

    def test_no_hidden_declarations(self) -> None:
        self.assertTrue(
            self.module.is_file(),
            "missing Smooth4PC/ReynoldsCableCocone.lean",
        )
        self.assertTrue(self.audit.is_file(), "missing T73ReynoldsAudit.lean")
        source = self.module.read_text(encoding="utf-8") + self.audit.read_text(
            encoding="utf-8"
        )
        for token in FORBIDDEN_TOKENS:
            self.assertNotRegex(
                source,
                rf"\b{re.escape(token)}\b",
                f"forbidden Lean token: {token}",
            )

    def test_audit_covers_every_theorem_of_the_module(self) -> None:
        declared = set(
            re.findall(
                r"(?m)^theorem\s+([A-Za-z_][A-Za-z0-9_'.]*)",
                self.module.read_text(encoding="utf-8"),
            )
        )
        self.assertEqual(
            {name.removeprefix("Smooth4PC.") for name in AXIOM_THEOREMS},
            declared,
            "the audit must report on exactly the theorems of the module",
        )

    def test_module_builds_and_audit_is_foundational(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="reynolds-cable-", dir=resolve_scratch()
        ) as tmp:
            olean_root = Path(tmp) / "olean"
            module_olean = olean_root / "Smooth4PC" / "ReynoldsCableCocone.olean"
            build = self.run_lean(self.module, olean_root, module_olean)
            self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
            audit = self.run_lean(self.audit, olean_root)
            output = audit.stdout + audit.stderr
            self.assertEqual(audit.returncode, 0, output)
        self.assertNotIn("sorryAx", output)
        self.assertNotIn("Lean.ofReduceBool", output)
        self.assert_axiom_reports(output)


if __name__ == "__main__":
    unittest.main()
