"""Exports: CSV, Excel, figure and results ZIPs, a methods paragraph and run metadata.

Rules for everything written here:

* Only parsed values are reported. A value the logs don't contain is written as
  "not recorded", never filled in from Vina's defaults.
* No absolute paths leave the machine: file names only.
* Text a spreadsheet could execute (=, +, -, @ at the start) is stored as text.
* ZIPs are byte-identical for identical input (fixed timestamps, sorted entries).

No Streamlit code here, so all of it is tested with plain pytest.
"""

import hashlib
import importlib.metadata
import io
import json
import platform
import re
import subprocess
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath, PureWindowsPath

import openpyxl
import pandas as pd

import vina_analyser
from vina_analyser.analysis import IMPLAUSIBLE_RMSD, LABELS, NOTES, RMSD_SUCCESS_CUTOFF
from vina_analyser.parsers import AGREEMENT_TOLERANCE, Issue

# UTF-8 with a byte-order mark: Excel on Windows reads a CSV without one in the local
# code page and turns Δ, Å and the Unicode minus into garbage.
CSV_ENCODING = "utf-8-sig"

MINUS = "−"
NOT_RECORDED = "not recorded"
PACKAGES = ("pandas", "numpy", "matplotlib", "streamlit", "openpyxl")
BUNDLE_FOLDER = "vina_results"
ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)  # the earliest date a zip can hold

# Vina logs carry nothing about how the inputs were prepared, so say so every time:
# leaving it out would read as if the paragraph were complete.
PREPARATION_GAP = (
    "Ligand and receptor preparation (protonation, charges, conversion to PDBQT, removal "
    "of waters and cofactors) is not recorded in Vina log files and must be described separately."
)

LIMITATIONS = [
    "A docking score is not a binding affinity. Vina scores are estimates from an "
    "empirical scoring function, not measured affinities, and say nothing about "
    "activity in a cell or an organism.",
    "The receptor is rigid. Protein flexibility, water molecules and changes in "
    "protonation are not modelled.",
    "No comparison across targets. Scores are ranked and compared only within one "
    "target, grid box and scoring function, never across targets.",
    "Vina's scoring error is roughly 2–3 kcal/mol. Smaller differences between ligands "
    "should not be read as differences in binding.",
    "Scores are not comparable across different scoring functions or very different "
    "grid boxes.",
]

CITATIONS = [
    "Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New "
    "Docking Methods, Expanded Force Field, and Python Bindings. J. Chem. Inf. Model. "
    "2021. DOI 10.1021/acs.jcim.1c00203",
    "Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with "
    "a new scoring function, efficient optimization and multithreading. J. Comput. "
    "Chem. 2010. DOI 10.1002/jcc.21334",
]

# Excel number formats by column header, so 0.3190 never shows as 0.319.
SCORE_FORMAT = "0.000"
NUMBER_FORMATS = {
    LABELS["affinity"]: SCORE_FORMAT,
    LABELS["delta_vs_control"]: SCORE_FORMAT,
    LABELS["mean_affinity"]: SCORE_FORMAT,
    LABELS["sd_affinity"]: SCORE_FORMAT,
    LABELS["range_affinity"]: SCORE_FORMAT,
    LABELS["ligand_efficiency"]: "0.0000",
    LABELS["rmsd_vs_reference"]: "0.00",
    "control_affinity": SCORE_FORMAT,
    "min_affinity": SCORE_FORMAT,
    "max_affinity": SCORE_FORMAT,
    "rmsd_lb": SCORE_FORMAT,
    "rmsd_ub": SCORE_FORMAT,
    "max_log_pdbqt_diff": SCORE_FORMAT,
}

# OWASP's list of characters that make a spreadsheet treat a cell as a formula.
_FORMULA_START = ("=", "+", "-", "@", "\t", "\r")


# ------------------------------------------------------------------------ tables


def display_table(frame: pd.DataFrame) -> pd.DataFrame:
    """A copy with column names replaced by their labels (with units) where known."""
    return frame.rename(columns={column: LABELS[column] for column in frame.columns if column in LABELS})


def _formula_like(value) -> bool:
    return isinstance(value, str) and value.startswith(_FORMULA_START)


def csv_bytes(frame: pd.DataFrame) -> bytes:
    """The table as CSV in UTF-8 with a BOM, without the index, safe to open in Excel.

    Text that a spreadsheet would run as a formula (a ligand typed as "=1+1", for
    example) gets a leading apostrophe. Numbers, negative ones included, are untouched.
    """
    safe = frame.copy()
    for column in safe.columns:
        if safe[column].dtype == object or pd.api.types.is_string_dtype(safe[column]):
            safe[column] = safe[column].map(lambda v: "'" + v if _formula_like(v) else v)
    return safe.to_csv(index=False, lineterminator="\n").encode(CSV_ENCODING)


# ----------------------------------------------------------------------- methods


@dataclass(frozen=True)
class Redocking:
    """A redocking RMSD the user computed against a reference they confirmed."""

    target: str
    ligand: str
    reference_name: str
    reference_sha256: str
    confirmed: bool
    table: pd.DataFrame


@dataclass
class MethodsText:
    text: str
    issues: list[Issue] = field(default_factory=list)


def methods_text(
    runs: pd.DataFrame, redocking: Redocking | None = None, git: dict | None = None
) -> MethodsText:
    """A methods paragraph built only from values parsed from the logs.

    `git` is git_info() for the analyser's own checkout, or None outside one.
    """
    issues: list[Issue] = []

    def setting(column, describe, name):
        values = _values(runs, column)
        if not values:
            return None
        distinct = sorted(set(values.values()), key=str)
        if len(distinct) == 1:
            text = describe(distinct[0])
            missing = len(runs) - len(values)
            return text + (f" (not recorded for {missing} runs)" if missing else "")
        counts = Counter(values.values())
        parts = []
        for value, count in sorted(counts.items(), key=lambda item: (-item[1], str(item[0]))):
            where = [label for label, v in values.items() if v == value]
            parts.append(f"{_plain(value)} ({where[0] if count == 1 else f'{count} runs'})")
        issues.append(
            Issue(
                "warning",
                "methods_value_varies",
                f"{name} differs between runs ({', '.join(parts)}); the methods paragraph "
                "lists every value instead of choosing one.",
            )
        )
        return f"{name} varied between runs: " + ", ".join(parts)

    version = setting("vina_version", lambda v: f"AutoDock Vina {v}", "Vina version")
    if version is None:
        version = "AutoDock Vina (version not recorded)"
    elif "varied" in version:
        version = f"AutoDock Vina ({version})"
    scoring = setting("scoring_function", lambda v: f"using the {v} scoring function", "scoring function")
    scoring = scoring or f"(scoring function: {NOT_RECORDED})"
    if "varied" in scoring:
        scoring = f"({scoring})"

    sentences = [f"Docking was performed with {version} (Eberhardt et al., 2021; Trott and Olson, 2010) {scoring}."]
    sentences.append(_targets_sentence(runs))
    sentences.append(PREPARATION_GAP)

    search = [
        setting("exhaustiveness", lambda v: f"an exhaustiveness of {int(v)}", "exhaustiveness")
        or f"exhaustiveness: {NOT_RECORDED}",
        setting("grid_spacing", lambda v: f"a grid spacing of {_number(v)} Å", "grid spacing")
        or f"grid spacing: {NOT_RECORDED}",
        setting("cpu", lambda v: f"{int(v)} CPU threads", "CPU threads") or f"CPU threads: {NOT_RECORDED}",
    ]
    sentences.append("The search used " + _join(search) + ".")
    sentences.append(_seed_sentence(runs))

    targets = sorted(runs["target"].unique())
    boxes = {target: _box_text(runs[runs["target"] == target], target, issues) for target in targets}
    if len(targets) == 1:
        box = boxes[targets[0]]
        if box.startswith("centre"):
            centre, size = box.split("; size ")
            sentences.append(
                f"The grid box was centred at {centre.removeprefix('centre ')} with a size of {size}."
            )
        else:
            sentences.append(f"Grid box: {box}.")
    else:
        sentences.append("The grid box of each target is listed below.")

    sentences.append(_poses_sentence(runs))
    sentences.append(
        "Scores were taken from the output PDBQT files (three decimals) where available, "
        "otherwise from the logs."
    )
    redock = _redocking_sentence(redocking)
    if redock:
        sentences.append(redock)
    sentences.append(_analyser_sentence(git))

    lines = [" ".join(sentences)]
    if len(targets) > 1:
        lines.append("")
        lines += [f"- {target}: {boxes[target]}" for target in targets]
    lines += ["", "References:", *[f"- {citation}" for citation in CITATIONS], ""]
    lines.append(
        "This paragraph was generated automatically from the Vina log files; check it "
        "before use."
    )
    return MethodsText("\n".join(lines), issues)


def _run_label(row) -> str:
    seed = row["seed"]
    return f"{row['ligand']} on {row['target']}, seed {NOT_RECORDED if pd.isna(seed) else int(seed)}"


def _values(runs: pd.DataFrame, column: str) -> dict:
    """{run label: value} for runs where the value is recorded."""
    return {
        _run_label(row): row[column]
        for _, row in runs.iterrows()
        if column in runs.columns and not pd.isna(row[column])
    }


def _plain(value) -> str:
    if isinstance(value, float):
        return _number(value)
    return str(int(value)) if hasattr(value, "__int__") and not isinstance(value, str) else str(value)


def _number(value: float, decimals: int | None = None) -> str:
    text = f"{value:.{decimals}f}" if decimals is not None else f"{value:g}"
    return text.replace("-", MINUS)


def _join(parts: list[str]) -> str:
    return parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + " and " + parts[-1]


def _basename(path) -> str:
    # Logs may hold Windows or POSIX paths; either way, only the file name is kept.
    return PurePosixPath(PureWindowsPath(str(path)).as_posix()).name


def _targets_sentence(runs: pd.DataFrame) -> str:
    ligands = runs["ligand"].nunique()
    targets = sorted(runs["target"].unique())
    described = []
    for target in targets:
        receptors = sorted({_basename(r) for r in runs.loc[runs["target"] == target, "receptor"].dropna()})
        receptor = " or ".join(receptors) if receptors else f"receptor file {NOT_RECORDED}"
        described.append(f"{target} (receptor {receptor})" if receptors else f"{target} ({receptor})")
    # Journals disallow a sentence that starts with a numeral.
    subject = "A single ligand was" if ligands == 1 else f"A total of {ligands} ligands were"
    return f"{subject} docked into {len(targets)} target{'s' if len(targets) != 1 else ''}: {_join(described)}."


def _analyser_sentence(git: dict | None) -> str:
    credit = f"Results were summarised with Vina Docking Results Analyser v{vina_analyser.__version__}"
    if not git:
        return credit + "."
    dirty = ", with uncommitted changes" if git.get("dirty") else ""
    return f"{credit} (commit {git['commit_short']}{dirty})."


def _seed_sentence(runs: pd.DataFrame) -> str:
    logged = runs[runs["seed_source"] == "log"] if "seed_source" in runs.columns else runs.iloc[0:0]
    if logged.empty:
        return f"Seeds were read from the logs (random seed: {NOT_RECORDED})."
    common = Counter(int(seed) for seed in logged["seed"]).most_common(1)[0][0]
    exceptions = []
    for (target, ligand), group in runs.groupby(["target", "ligand"], sort=True):
        seeds = sorted(int(s) for s, src in zip(group["seed"], group["seed_source"]) if src == "log")
        unlogged = len(group) - len(seeds)
        if seeds == [common] and not unlogged:
            continue
        parts = []
        if seeds:
            parts.append(f"seed{'s' if len(seeds) > 1 else ''} {_join([str(s) for s in seeds])}")
        if unlogged:
            parts.append(f"{unlogged} run{'s' if unlogged > 1 else ''} with no seed recorded in a log")
        exceptions.append(f"{ligand} on {target}: {'; '.join(parts)}")
    sentence = f"Runs used random seed {common} (read from the logs)"
    if exceptions:
        sentence += ", except " + "; ".join(exceptions)
    return sentence + "."


def _box_text(runs: pd.DataFrame, target: str, issues: list) -> str:
    columns = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    boxes = {
        tuple(float(row[c]) for c in columns)
        for _, row in runs.iterrows()
        if not any(pd.isna(row[c]) for c in columns)
    }
    if not boxes:
        return f"grid box {NOT_RECORDED}"
    texts = [
        f"centre ({', '.join(_number(v, 3) for v in box[:3])}) Å; size "
        f"{' × '.join(_number(v, 1) for v in box[3:])} Å"
        for box in sorted(boxes)
    ]
    if len(texts) > 1:
        issues.append(
            Issue(
                "warning",
                "methods_value_varies",
                f"Target {target!r} was docked with {len(texts)} different grid boxes; "
                "all are listed in the methods paragraph.",
            )
        )
        return "grid box varied between runs: " + " | ".join(texts)
    return texts[0]


def _range(values) -> str | None:
    values = [int(v) for v in values if not pd.isna(v)]
    if not values:
        return None
    low, high = min(values), max(values)
    return str(low) if low == high else f"between {low} and {high}"


def _poses_sentence(runs: pd.DataFrame) -> str:
    modes = _range(runs["n_modes_log"])
    poses = _range(runs["n_poses_pdbqt"])
    first = (
        f"The logs list {modes} binding modes per run"
        if modes
        else f"Binding modes in the logs: {NOT_RECORDED}"
    )
    if poses is None:
        return first + "."
    return (
        f"{first}; the output PDBQT files keep only poses within energy_range of the best "
        f"mode (energy_range: {NOT_RECORDED} in the logs), {poses} poses per run."
    )


def _redocking_sentence(redocking: Redocking | None) -> str | None:
    """Only for an RMSD computed against a reference the user confirmed as experimental."""
    if redocking is None or not redocking.confirmed:
        return None
    table = redocking.table
    top = table[(table["pose"] == 1) & table["rmsd_vs_reference"].notna()]
    if top.empty:
        return None
    values = top["rmsd_vs_reference"].astype(float)
    seeds = sorted(int(s) for s in top["seed"].dropna())
    tail = (
        f"against the experimental pose supplied as {redocking.reference_name} (atoms matched "
        "in file order; no superposition or symmetry correction)."
    )
    start = f"The docking protocol was assessed by redocking {redocking.ligand} into {redocking.target}"
    if len(top) == 1:
        return f"{start}, giving a heavy-atom RMSD of {values.iloc[0]:.2f} Å for the top-ranked pose {tail}"
    low, high = values.min(), values.max()
    spread = f"{low:.2f} Å" if f"{low:.2f}" == f"{high:.2f}" else f"{low:.2f}–{high:.2f} Å"
    seed_text = f" (seeds {_join([str(s) for s in seeds])})" if seeds else ""
    return (
        f"{start}: the top-ranked pose of each of {len(top)} runs{seed_text} gave a "
        f"heavy-atom RMSD of {spread} {tail}"
    )


# ------------------------------------------------------------------------- excel


def workbook_bytes(sheets: dict[str, pd.DataFrame], readme: list[str], created: datetime) -> bytes:
    """An .xlsx with a "Read me" sheet first, then one sheet per table.

    `created` is the analysis time, so the workbook and run_metadata.json agree.
    (openpyxl sets the "modified" property to the time of saving; that cannot be fixed
    without editing the file afterwards.)
    """
    book = openpyxl.Workbook()
    first = book.active
    first.title = "Read me"
    lines = [
        f"Vina Docking Results Analyser {vina_analyser.__version__}",
        "",
        "Limitations",
        *LIMITATIONS,
        "",
        *readme,
    ]
    for row, line in enumerate(lines, start=1):
        _write(first.cell(row=row, column=1), line)
    first.column_dimensions["A"].width = 120

    for title, frame in sheets.items():
        sheet = book.create_sheet(title[:31])  # Excel's limit on sheet names
        for column, header in enumerate(frame.columns, start=1):
            _write(sheet.cell(row=1, column=column), str(header))
            width = max(10, min(50, len(str(header)) + 2))
            sheet.column_dimensions[openpyxl.utils.get_column_letter(column)].width = width
        for row_number, row in enumerate(frame.itertuples(index=False), start=2):
            for column, (header, value) in enumerate(zip(frame.columns, row), start=1):
                cell = sheet.cell(row=row_number, column=column)
                _write(cell, value)
                if str(header) in NUMBER_FORMATS and isinstance(cell.value, (int, float)):
                    cell.number_format = NUMBER_FORMATS[str(header)]
        sheet.freeze_panes = "A2"

    when = created.astimezone(timezone.utc).replace(tzinfo=None)
    book.properties.created = when
    book.properties.modified = when
    buffer = io.BytesIO()
    book.save(buffer)
    return buffer.getvalue()


def _write(cell, value) -> None:
    """Write one value: blanks for missing, plain Python types, formulas never."""
    if isinstance(value, (list, tuple)):
        value = ", ".join(str(item) for item in value)
    elif value is None or (not isinstance(value, str) and pd.api.types.is_scalar(value) and pd.isna(value)):
        value = None
    elif hasattr(value, "item"):  # numpy and pandas scalars
        value = value.item()
    cell.value = value
    if _formula_like(value):
        cell.data_type = "s"  # openpyxl would otherwise store "=..." as a formula


# --------------------------------------------------------------------- file names


def safe_name(name: str) -> str:
    """A file-name-safe version of a target or ligand name: no folders, no drive."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("._ ")
    return cleaned or "unnamed"


@dataclass(frozen=True)
class FigureSpec:
    stem: str
    kind: str  # "heatmap", "ranked" or "strip"
    target: str | None


def figure_specs(targets) -> list[FigureSpec]:
    """Every figure in the export, with unique, safe file stems."""
    specs = [FigureSpec("heatmap_best_affinity", "heatmap", None)]
    used = {specs[0].stem}
    for target in targets:
        for prefix, kind in (("ranked", "ranked"), ("poses", "strip")):
            stem = f"{prefix}_{safe_name(target)}"
            candidate, number = stem, 2
            while candidate in used:
                candidate, number = f"{stem}_{number}", number + 1
            used.add(candidate)
            specs.append(FigureSpec(candidate, kind, target))
    return specs


# --------------------------------------------------------------------------- zips


def deterministic_zip(entries: dict[str, bytes]) -> bytes:
    """A ZIP that is byte-identical for identical entries, on every operating system."""
    for name in entries:
        if name.startswith("/") or re.match(r"^[A-Za-z]:", name) or "\\" in name or ".." in name.split("/"):
            raise ValueError(f"Unsafe name for a zip entry: {name!r}")
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name in sorted(entries):
            info = zipfile.ZipInfo(name, date_time=ZIP_TIMESTAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3  # otherwise 0 on Windows and 3 elsewhere
            info.external_attr = 0o644 << 16
            archive.writestr(info, entries[name])
    return buffer.getvalue()


def _captions_text(captions: dict[str, str]) -> bytes:
    blocks = [f"{stem}\n{'-' * len(stem)}\n{text}" for stem, text in captions.items()]
    return ("\n\n".join(blocks) + "\n").encode("utf-8")


def _figure_entries(figures: dict[str, dict], captions: dict[str, str], folder: str = "") -> dict:
    entries = {f"{folder}captions.txt": _captions_text(captions)}
    for stem, files in figures.items():
        for fmt in ("png", "svg"):
            entries[f"{folder}{stem}.{fmt}"] = files[fmt]
    return entries


def figures_zip(figures: dict[str, dict], captions: dict[str, str]) -> bytes:
    """Every figure as PNG and SVG, plus captions.txt."""
    return deterministic_zip(_figure_entries(figures, captions))


def bundle_readme() -> bytes:
    lines = [
        f"Results from the Vina Docking Results Analyser {vina_analyser.__version__}",
        "",
        "results.xlsx       every table, one sheet each, with a Read me sheet",
        "tables/            the same tables as CSV (UTF-8 with BOM)",
        "figures/           300 dpi PNG and SVG figures, with captions.txt",
        "methods.txt        methods paragraph built from the Vina logs (check before use)",
        "run_metadata.json  versions, git commit, input file SHA-256 hashes, choices made",
        "",
        "Limitations",
        *[f"- {line}" for line in LIMITATIONS],
        "",
    ]
    return "\n".join(lines).encode("utf-8")


def bundle_zip(
    *,
    tables: dict[str, pd.DataFrame],
    workbook: bytes,
    methods: str,
    metadata: bytes,
    figures: dict[str, dict],
    captions: dict[str, str],
) -> bytes:
    """Everything in one ZIP, under a single vina_results/ folder."""
    root = f"{BUNDLE_FOLDER}/"
    entries = {
        f"{root}README.txt": bundle_readme(),
        f"{root}results.xlsx": workbook,
        f"{root}methods.txt": methods.encode("utf-8"),
        f"{root}run_metadata.json": metadata,
        **{f"{root}tables/{name}": csv_bytes(frame) for name, frame in tables.items()},
        **_figure_entries(figures, captions, folder=f"{root}figures/"),
    }
    return deterministic_zip(entries)


# ----------------------------------------------------------------------- metadata


_SHA = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")
_REF = re.compile(r"refs/[A-Za-z0-9._/-]+")


def git_info(repo_dir, runner=subprocess.run) -> dict | None:
    """The commit the code was run from, or None when there is no checkout.

    Uses the git command when it works. Where it doesn't (Streamlit Community Cloud may
    have no git binary), the commit is read from .git directly; whether the working
    tree has changes is then unknown, so "dirty" is None.
    """

    def git(*args) -> str | None:
        result = runner(
            ["git", *args], cwd=str(repo_dir), capture_output=True, text=True, timeout=10, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else None

    try:
        full = git("rev-parse", "HEAD")
        short = git("rev-parse", "--short", "HEAD") if full else None
        status = git("status", "--porcelain") if short else None
    except (OSError, subprocess.SubprocessError):
        full = None
    if full and short is not None and status is not None:
        # Uncommitted or untracked files mean the commit alone doesn't describe the code.
        return {"commit": full, "commit_short": short, "dirty": status != ""}
    commit = _commit_from_dot_git(Path(repo_dir) / ".git")
    if commit is None:
        return None
    return {"commit": commit, "commit_short": commit[:7], "dirty": None}


def _commit_from_dot_git(dot_git: Path) -> str | None:
    """HEAD's commit, read from the files git keeps: HEAD, then a loose or packed ref."""
    try:
        if not dot_git.is_dir():
            return None  # no checkout, or a worktree/submodule pointer: not followed
        head = (dot_git / "HEAD").read_text(encoding="utf-8")[:300].strip()
        if _SHA.fullmatch(head):
            return head  # detached HEAD
        if not head.startswith("ref: "):
            return None
        ref = head.removeprefix("ref: ").strip()
        if not _REF.fullmatch(ref) or ".." in ref:
            return None
        loose = dot_git / ref
        if loose.is_file():
            sha = loose.read_text(encoding="utf-8")[:100].strip()
            return sha if _SHA.fullmatch(sha) else None
        packed = dot_git / "packed-refs"
        if packed.is_file():
            for line in packed.read_text(encoding="utf-8").splitlines():
                sha, _, name = line.partition(" ")
                if name == ref and _SHA.fullmatch(sha):
                    return sha
    except (OSError, UnicodeDecodeError):
        return None
    return None


def package_versions() -> dict[str, str]:
    versions = {}
    for package in PACKAGES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = NOT_RECORDED
    return versions


_FIND_GIT = object()


def run_metadata(
    *,
    files: dict[str, bytes],
    names: dict[str, dict],
    runs: pd.DataFrame,
    origin: str,
    control: str | None,
    redocking: Redocking | None,
    analysed_at: datetime,
    repo_dir=None,
    git=_FIND_GIT,
) -> dict:
    """Everything needed to tell which code, packages and inputs produced a result.

    `git` can be passed in (the app caches it); by default it is looked up in repo_dir.
    """
    seeds = {}
    for _, run in runs.iterrows():
        for column in ("pdbqt_file", "log_file"):
            if isinstance(run[column], str):
                seed = None if pd.isna(run["seed"]) else int(run["seed"])
                seeds[run[column]] = (seed, run.get("seed_source"))

    inputs = []
    for name in sorted(files):
        data = files[name]
        entry = {
            "name": name,
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "included": name in names,
        }
        if name in names:
            seed, source = seeds.get(name, (None, None))
            entry.update(
                target=names[name]["target"],
                ligand=names[name]["ligand"],
                seed=seed,
                seed_source=None if source is None or pd.isna(source) else source,
            )
        inputs.append(entry)

    meta = {"app": "Vina Docking Results Analyser", "app_version": vina_analyser.__version__}
    if git is _FIND_GIT:
        git = git_info(repo_dir) if repo_dir is not None else None
    if git:
        meta["git"] = git
    meta.update(
        analysed_at=analysed_at.isoformat(),
        python=platform.python_version(),
        platform=platform.platform(),
        packages=package_versions(),
        origin=origin,
        control=control,
        settings={
            "log_pdbqt_agreement_tolerance": AGREEMENT_TOLERANCE,
            "rmsd_success_cutoff_angstrom": RMSD_SUCCESS_CUTOFF,
            "implausible_rmsd_angstrom": IMPLAUSIBLE_RMSD,
        },
        input_files=inputs,
        rmsd_reference=None
        if redocking is None
        else {
            "file": redocking.reference_name,
            "sha256": redocking.reference_sha256,
            "target": redocking.target,
            "ligand": redocking.ligand,
            "confirmed_experimental_pose": redocking.confirmed,
        },
        notes=[NOTES["scoring_error"], NOTES["affinity"]],
    )
    return meta


def metadata_json(meta: dict) -> bytes:
    return (json.dumps(meta, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
