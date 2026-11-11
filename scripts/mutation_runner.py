"""Deliberate-bug (mutation) runner with an INVALID guard.

Plants one bug at a time from scripts/mutation_sets.py in a throwaway copy of the
repository, runs the full test suite against it, and reports one verdict per bug:

* CAUGHT       some tests fail, and they are a small, specific subset.
* MISSED       every test still passes: no test guards this behaviour.
* INVALID      the result says nothing about the behaviour, because the mutant:
               - has a syntax error, or the package no longer imports;
               - breaks test collection;
               - (app.py only) crashes the app on its start screen;
               - fails more than BROAD_FRACTION of all tests.
* NOT APPLIED  the text to replace was not found exactly once (the code changed;
               update the entry in mutation_sets.py).

The unmodified copy must pass first, or nothing is run. The full suite runs for every
mutant (no -x), so the share of failing tests is always known.

Usage (from the repository root):
    python scripts/mutation_runner.py                   # everything
    python scripts/mutation_runner.py --phase phase5    # one or more phases
    python scripts/mutation_runner.py --module vina_analyser/report.py
    python scripts/mutation_runner.py --name "methods:"  # names containing this text
    python scripts/mutation_runner.py --check           # only check entries still apply
"""

import argparse
import os
import py_compile
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mutation_sets import MUTATIONS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BROAD_FRACTION = 0.30
IGNORE = shutil.ignore_patterns(".venv", "venv", "__pycache__", ".pytest_cache", "figures")
IMPORT_CHECK = (
    "import vina_analyser.parsers, vina_analyser.analysis, vina_analyser.plots, "
    "vina_analyser.inputs, vina_analyser.report"
)
APP_SMOKE = (
    "from streamlit.testing.v1 import AppTest\n"
    "at = AppTest.from_file('app.py', default_timeout=120); at.run()\n"
    "assert not at.exception, at.exception"
)
SUMMARY = re.compile(r"(\d+) (passed|failed|error|errors)")


def _remove(path: Path) -> None:
    def writable_then_retry(function, target, *_):
        os.chmod(target, stat.S_IWRITE)  # git objects are read-only on Windows
        function(target)

    if not path.exists():
        return
    if sys.version_info >= (3, 12):
        shutil.rmtree(path, onexc=writable_then_retry)
    else:
        shutil.rmtree(path, onerror=writable_then_retry)


def _copy(dest: Path) -> None:
    # .git is copied too, so the git-metadata tests run against a real checkout.
    _remove(dest)
    shutil.copytree(ROOT, dest, ignore=IGNORE)


def _counts(output: str) -> dict:
    line = output.strip().splitlines()[-1] if output.strip() else ""
    counts = {"passed": 0, "failed": 0, "error": 0}
    for number, kind in SUMMARY.findall(line):
        counts["error" if kind.startswith("error") else kind] += int(number)
    return counts


def _run(work: Path, *args) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, *args], cwd=work, capture_output=True, text=True)


def _pytest(work: Path, *args) -> subprocess.CompletedProcess:
    return _run(work, "-m", "pytest", "-p", "no:cacheprovider", *args)


def _last_line(text: str) -> str:
    lines = text.strip().splitlines()
    return lines[-1][:90] if lines else "?"


def _apply(text: str, edits) -> tuple[str | None, str]:
    for old, new in edits:
        found = text.count(old)
        if found != 1:
            return None, f"target text found {found} times"
        text = text.replace(old, new)
    return text, ""


def run_one(work: Path, rel: str, mutation) -> tuple:
    _copy(work)
    path = work / rel
    text, problem = _apply(path.read_text(encoding="utf-8"), mutation.edits)
    if text is None:
        return "NOT APPLIED", problem
    path.write_text(text, encoding="utf-8")

    if rel.endswith(".py"):
        try:
            py_compile.compile(str(path), doraise=True)
        except py_compile.PyCompileError as error:
            return "INVALID", "syntax error: " + _last_line(str(error))
    check = _run(work, "-c", IMPORT_CHECK)
    if check.returncode:
        return "INVALID", "import error: " + _last_line(check.stderr)
    if rel == "app.py":
        smoke = _run(work, "-c", APP_SMOKE)
        if smoke.returncode:
            return "INVALID", "app fails to start: " + _last_line(smoke.stderr)
    if _pytest(work, "--collect-only", "-q").returncode:
        return "INVALID", "collection error"

    run = _pytest(work, "-q", "-rfE")
    counts = _counts(run.stdout)
    total = sum(counts.values())
    bad = counts["failed"] + counts["error"]
    failing = [line.split(" - ")[0].split("::", 1)[-1] for line in run.stdout.splitlines()
               if line.startswith(("FAILED ", "ERROR tests/"))]
    detail = f"{bad}/{total} failed" + (": " + ", ".join(failing[:3]) if failing else "")
    if bad == 0:
        return "MISSED", detail
    if bad / total > BROAD_FRACTION:
        return "INVALID", "too broad, " + detail
    return "CAUGHT", detail


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--phase", action="append", help="only this phase (repeatable)")
    parser.add_argument("--module", action="append", help="only mutations of this file (repeatable)")
    parser.add_argument("--name", help="only mutations whose name contains this text")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--scratch", type=Path, help="work folder (default: a new temp folder)")
    parser.add_argument("--check", action="store_true", help="only check that every entry still applies")
    args = parser.parse_args()

    chosen = [(rel, m) for rel, items in MUTATIONS.items() for m in items
              if (not args.module or rel in args.module) and (not args.phase or m.phase in args.phase)
              and (not args.name or args.name.lower() in m.name.lower())]
    if not chosen:
        print("no mutations match those filters")
        return 1

    if args.check:
        stale = 0
        for rel, m in chosen:
            text, problem = _apply((ROOT / rel).read_text(encoding="utf-8"), m.edits)
            if text is None:
                stale += 1
                print(f"NOT APPLIED  [{m.phase}] {rel}: {m.name} ({problem})")
        print(f"{len(chosen) - stale}/{len(chosen)} mutations apply to the current code")
        return 1 if stale else 0

    scratch = args.scratch or Path(tempfile.mkdtemp(prefix="vina_mutants_"))
    try:
        _copy(scratch / "baseline")
        baseline = _pytest(scratch / "baseline", "-q")
        print(f"baseline (no bug planted): {_last_line(baseline.stdout)}", flush=True)
        if baseline.returncode:
            print("the unmodified copy fails, so no mutation result would mean anything")
            return 1
        _remove(scratch / "baseline")

        def task(indexed):
            index, (rel, mutation) = indexed
            work = scratch / f"mutant_{index}"
            verdict, detail = run_one(work, rel, mutation)
            _remove(work)
            print(f"[{mutation.phase}] {verdict:<11} {mutation.name}  ({detail})", flush=True)
            return rel, mutation, verdict, detail

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(task, enumerate(chosen)))
    finally:
        if args.scratch is None:
            _remove(scratch)

    print("\n==== summary ====")
    for phase in dict.fromkeys(m.phase for _, m in chosen):
        print(f"\n{phase}")
        for rel, m, verdict, detail in results:
            if m.phase == phase:
                print(f"  {verdict:<11} {rel}: {m.name}  ({detail})")
    totals = {}
    for _, m, verdict, _ in results:
        key = "guard" if m.phase == "guard" else "real"
        totals.setdefault(key, {}).setdefault(verdict, 0)
        totals[key][verdict] += 1
    print(f"\ntotals: {totals}")
    # A guard mutant must come out INVALID; anything else means the guard is broken.
    guard_broken = any(m.phase == "guard" and verdict != "INVALID" for _, m, verdict, _ in results)
    missed = any(m.phase != "guard" and verdict != "CAUGHT" for _, m, verdict, _ in results)
    return 1 if guard_broken or missed else 0


if __name__ == "__main__":
    sys.exit(main())
