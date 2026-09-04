from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


AXIOM_THEOREMS = (
    "Smooth4PC.pairingCoeff_lt_three_eq_zero",
    "Smooth4PC.pairingCoeff_three_of_startsAtThree",
    "Smooth4PC.pairingCoeff_transport",
    "Smooth4PC.startsAtThree_transport",
    "Smooth4PC.cubic_invariant_under_simultaneous_transport",
    "Smooth4PC.transportVectorSeries_zero_of_id",
    "Smooth4PC.pairingCoeff_three_transportVectorSeries",
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


def resolve_mathlib(repo: Path) -> Path:
    configured = os.environ.get("T73_MATHLIB")
    if configured:
        candidate = Path(configured).expanduser()
    else:
        manifest = json.loads(
            (repo / "lake-manifest.json").read_text(encoding="utf-8")
        )
        entry = next(
            package for package in manifest["packages"] if package["name"] == "mathlib"
        )
        if entry["type"] == "path":
            candidate = repo / entry["dir"]
        elif entry["type"] == "git":
            candidate = (
                repo
                / manifest.get("packagesDir", ".lake/packages")
                / entry["name"]
            )
            if entry.get("subDir"):
                candidate /= entry["subDir"]
        else:
            raise RuntimeError(f"unsupported mathlib source: {entry['type']}")

    if not candidate.is_dir():
        raise FileNotFoundError(
            f"mathlib is not materialized at {candidate}; run `lake update`"
        )
    return candidate.resolve()


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
        else Path.home() / "ws" / "tmp-t73" / "lean"
    )
    scratch.mkdir(parents=True, exist_ok=True)
    return scratch


class FilteredCubicNaturalityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo = Path(__file__).resolve().parents[1]
        cls.module = cls.repo / "Smooth4PC" / "FilteredCubicNaturality.lean"
        cls.audit = cls.repo / "T73FilteredCubicAudit.lean"
        cls.lake = resolve_lake()
        cls.mathlib = resolve_mathlib(cls.repo)

    def lean_environment(self, olean_root: Path) -> dict[str, str]:
        env = os.environ.copy()
        env.pop("LEAN_PATH", None)
        root_packages = self.repo / ".lake" / "packages"
        if root_packages.is_dir():
            package_roots = sorted(root_packages.iterdir())
        else:
            package_roots = [self.mathlib]
            nested_packages = self.mathlib / ".lake" / "packages"
            if nested_packages.is_dir():
                package_roots.extend(sorted(nested_packages.iterdir()))
        paths = [olean_root]
        paths.extend(
            package / ".lake" / "build" / "lib" / "lean"
            for package in package_roots
            if (package / ".lake" / "build" / "lib" / "lean").is_dir()
        )
        env["LEAN_PATH"] = os.pathsep.join(str(path) for path in paths)
        return env

    def run_lean(
        self, source: Path, olean_root: Path, output: Path | None = None
    ) -> subprocess.CompletedProcess[str]:
        args = [str(self.lake), "env", "lean"]
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

    def test_no_hidden_declarations(self) -> None:
        self.assertTrue(
            self.module.is_file(),
            "missing Smooth4PC/FilteredCubicNaturality.lean",
        )
        self.assertTrue(self.audit.is_file(), "missing T73FilteredCubicAudit.lean")
        source = self.module.read_text(encoding="utf-8") + self.audit.read_text(
            encoding="utf-8"
        )
        for token in FORBIDDEN_TOKENS:
            self.assertNotRegex(
                source,
                rf"\b{re.escape(token)}\b",
                f"forbidden Lean token: {token}",
            )

    def test_module_builds_and_audit_is_foundational(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="filtered-cubic-", dir=resolve_scratch()
        ) as tmp:
            olean_root = Path(tmp) / "olean"
            module_olean = (
                olean_root / "Smooth4PC" / "FilteredCubicNaturality.olean"
            )
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
