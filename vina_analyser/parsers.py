"""Read AutoDock Vina output files into tidy pandas tables.

Two sources describe the same docking run:

* the output PDBQT has 3-decimal scores and atoms, but only poses within energy_range;
* the log has every mode, but rounded scores and no atoms.

load_results() merges them. Nothing here touches the disk except load_folder(), so
in-memory uploads can use the same code as local files.
"""

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path, PurePosixPath

import pandas as pd

# Logs seem to keep 4 significant figures, so PDBQT and log values can differ by up to
# 0.005 from rounding alone. Anything over 0.01 is a real disagreement.
AGREEMENT_TOLERANCE = 0.01

# An explicit list of known types, rather than "not hydrogen", means an unexpected type
# is reported instead of being counted as a heavy atom by accident.
HYDROGEN_TYPES = {"H", "HD", "HS"}
HEAVY_TYPES = {
    "C", "A", "N", "NA", "NS", "O", "OA", "OS", "S", "SA", "P", "Se", "SE",
    "F", "Cl", "CL", "Br", "BR", "I", "Si", "SI", "B", "At",
    "Mg", "MG", "Ca", "CA", "Mn", "MN", "Fe", "FE", "Zn", "ZN",
}
# Meeko opens a macrocycle ring for docking and marks the two cut ends. CG0, CG1, ... are
# REAL carbons bonded into the ligand, so they count. G0, G1, ... (named *1, *2) are dummy
# atoms sitting on top of the partner CG carbon, so they don't.
GLUE_CARBON = re.compile(r"CG\d+")
MACROCYCLE_PSEUDO_ATOM = re.compile(r"G\d+")

POSE_COLUMNS = {
    "target": "str",
    "ligand": "str",
    "seed": "Int64",  # Int64 (capital I) allows blanks; plain int64 does not
    "pose": "int64",
    "affinity": "float64",
    "rmsd_lb": "float64",
    "rmsd_ub": "float64",
    "source_file": "str",
    "heavy_atoms": "Int64",
    "has_coordinates": "bool",
}
RUN_COLUMNS = {
    "target": "str",
    "ligand": "str",
    "seed": "Int64",
    "pdbqt_file": "object",
    "log_file": "object",
    "n_poses_pdbqt": "Int64",
    "n_modes_log": "Int64",
    "max_log_pdbqt_diff": "float64",
    "vina_version": "object",
    "scoring_function": "object",
    "receptor": "object",
    "log_ligand": "object",
    "center_x": "float64",
    "center_y": "float64",
    "center_z": "float64",
    "size_x": "float64",
    "size_y": "float64",
    "size_z": "float64",
    "grid_spacing": "float64",
    "exhaustiveness": "Int64",
    "cpu": "Int64",
}


@dataclass(frozen=True)
class Issue:
    """A problem found while parsing. `code` is stable; `message` is for people."""

    severity: str  # "error", "warning" or "info"
    code: str
    message: str
    file: str = ""


@dataclass(frozen=True)
class ParsedName:
    target: str
    ligand: str
    seed: int | None
    kind: str  # "pdbqt" or "log"


@dataclass
class PdbqtPose:
    pose: int
    affinity: float | None = None
    rmsd_lb: float | None = None
    rmsd_ub: float | None = None
    atom_types: list[str] = field(default_factory=list)
    has_flexible_residues: bool = False


@dataclass
class PdbqtData:
    poses: list[PdbqtPose] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)


@dataclass(frozen=True)
class LogMode:
    mode: int
    affinity: float
    rmsd_lb: float
    rmsd_ub: float


@dataclass
class LogData:
    """Header values stay None when absent: Vina defaults are never filled in."""

    modes: list[LogMode] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)
    vina_version: str | None = None
    scoring_function: str | None = None
    receptor: str | None = None
    ligand: str | None = None
    center: tuple[float, float, float] | None = None
    size: tuple[float, float, float] | None = None
    grid_spacing: float | None = None
    exhaustiveness: int | None = None
    cpu: int | None = None
    seed: int | None = None


@dataclass
class ParsedResults:
    poses: pd.DataFrame  # one row per pose
    runs: pd.DataFrame  # one row per docking run
    issues: list[Issue]


def decode_text(data: bytes) -> str:
    """Decode file bytes, honouring a byte-order mark.

    Windows PowerShell 5.1 saves `vina ... > out.log` as UTF-16, which would otherwise
    look like garbage to every regular expression below.
    """
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16")
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        return data.decode("latin-1")  # never fails, and these files are mostly ASCII


def parse_filename(name: str) -> ParsedName | None:
    """Read `<target>__<ligand>[_seedN]_out.pdbqt` or `<target>__<ligand>[_seedN].log`.

    Returns None when the name doesn't clearly fit, so the caller can ask the user.
    """
    basename = PurePosixPath(name.replace("\\", "/")).name
    if basename.endswith("_out.pdbqt"):
        stem, kind = basename.removesuffix("_out.pdbqt"), "pdbqt"
    elif basename.endswith(".log"):
        stem, kind = basename.removesuffix(".log"), "log"
    else:
        return None

    parts = stem.split("__")
    if len(parts) != 2:
        return None
    target, ligand = parts

    seed = None
    seed_match = re.fullmatch(r"(.+)_seed(\d+)", ligand)
    if seed_match:
        ligand, seed = seed_match.group(1), int(seed_match.group(2))

    for part in (target, ligand):
        if not part or part.startswith("_") or part.endswith("_"):
            return None
    return ParsedName(target, ligand, seed, kind)


def _is_heavy(atom_type: str) -> bool:
    return atom_type in HEAVY_TYPES or GLUE_CARBON.fullmatch(atom_type) is not None


def _unknown_atom_types(atom_types: Iterable[str]) -> set[str]:
    return {
        atom_type
        for atom_type in atom_types
        if atom_type not in HYDROGEN_TYPES
        and not _is_heavy(atom_type)
        and not MACROCYCLE_PSEUDO_ATOM.fullmatch(atom_type)
    }


def count_heavy_atoms(atom_types: Iterable[str]) -> int | None:
    """Count real non-hydrogen atoms, or return None if any type is unrecognised."""
    atom_types = list(atom_types)
    if _unknown_atom_types(atom_types):
        return None
    return sum(_is_heavy(atom_type) for atom_type in atom_types)


def parse_pdbqt(text: str) -> PdbqtData:
    """Read the MODEL ... ENDMDL blocks of a Vina output PDBQT."""
    data = PdbqtData()
    seen = set()
    current = None  # the pose being read, or None when between MODEL blocks
    keep_current = False

    for line in text.lstrip("﻿").splitlines():
        if line.startswith("MODEL"):
            if current is not None:
                data.issues.append(_truncated(current.pose))
            fields = line.split()
            if len(fields) < 2 or not fields[1].isdigit():
                data.issues.append(
                    Issue("error", "malformed_model", f"Unreadable MODEL line: {line!r}")
                )
                current, keep_current = None, False
                continue
            current = PdbqtPose(pose=int(fields[1]))
            keep_current = current.pose not in seen
            if keep_current:
                seen.add(current.pose)
                data.poses.append(current)
            else:
                data.issues.append(
                    Issue(
                        "error",
                        "duplicate_model",
                        f"MODEL {current.pose} appears twice; only the first was kept.",
                    )
                )
        elif current is None:
            continue
        elif line.startswith("ENDMDL"):
            current = None
        elif line.startswith("REMARK VINA RESULT:"):
            _read_vina_result(line, current, data)
        elif line.startswith(("ATOM", "HETATM")):
            # Read col 78 to the end of the line, not just cols 78-79: "CG0" has 3 letters.
            current.atom_types.append(line[77:].strip())
        elif line.startswith("BEGIN_RES"):
            current.has_flexible_residues = True

    if current is not None:
        data.issues.append(_truncated(current.pose))
    if not data.poses:
        data.issues.append(
            Issue(
                "error",
                "no_poses",
                "No MODEL blocks found. This doesn't look like a Vina output file.",
            )
        )
    return data


def _truncated(pose: int) -> Issue:
    return Issue("warning", "truncated_pose", f"MODEL {pose} has no ENDMDL; the file may be cut off.")


def _read_vina_result(line: str, pose: PdbqtPose, data: PdbqtData) -> None:
    try:
        pose.affinity, pose.rmsd_lb, pose.rmsd_ub = (
            float(value) for value in line.split(":", 1)[1].split()[:3]
        )
    except ValueError:
        pose.affinity = pose.rmsd_lb = pose.rmsd_ub = None
        data.issues.append(
            Issue("warning", "malformed_score", f"Unreadable score line: {line!r}")
        )


_NUMBER = r"(-?\d+(?:\.\d+)?)"
_XYZ = rf"X\s*{_NUMBER}\s*Y\s*{_NUMBER}\s*Z\s*{_NUMBER}"


def _search(pattern: str, text: str) -> re.Match | None:
    return re.search(pattern, text, flags=re.MULTILINE)


def parse_log(text: str) -> LogData:
    """Read the header values and the mode table from Vina's stdout."""
    text = text.lstrip("﻿").replace("\r\n", "\n")
    data = LogData()

    def group(pattern):
        match = _search(pattern, text)
        return match.group(1).strip() if match else None

    # Anchored to the line start: the citation banner has "# AutoDock Vina 1.2.0: ...".
    data.vina_version = group(r"^AutoDock Vina\s+(\S+)\s*$")
    data.scoring_function = group(r"^Scoring function\s*:\s*(\S+)")
    data.receptor = group(r"^Rigid receptor\s*:\s*(.+)$")
    data.ligand = group(r"^Ligand\s*:\s*(.+)$")

    center = _search(rf"^Grid center\s*:\s*{_XYZ}", text)
    size = _search(rf"^Grid size\s*:\s*{_XYZ}", text)
    data.center = tuple(float(value) for value in center.groups()) if center else None
    data.size = tuple(float(value) for value in size.groups()) if size else None

    spacing = group(rf"^Grid space\s*:\s*{_NUMBER}")
    exhaustiveness = group(r"^Exhaustiveness\s*:\s*(\d+)")
    cpu = group(r"^CPU\s*:\s*(\d+)")
    seed = group(r"random seed:\s*(-?\d+)")
    data.grid_spacing = float(spacing) if spacing else None
    data.exhaustiveness = int(exhaustiveness) if exhaustiveness else None
    data.cpu = int(cpu) if cpu else None
    data.seed = int(seed) if seed else None

    _read_mode_table(text, data)
    return data


def _read_mode_table(text: str, data: LogData) -> None:
    in_table = False
    for line in text.splitlines():
        if not in_table:
            in_table = line.startswith("-----+")
            continue
        fields = line.split()
        if not fields or not fields[0].isdigit():
            break  # a blank line or following text ends the table
        try:
            mode, affinity, rmsd_lb, rmsd_ub = fields
            data.modes.append(LogMode(int(mode), float(affinity), float(rmsd_lb), float(rmsd_ub)))
        except ValueError:
            data.issues.append(
                Issue("warning", "malformed_mode_row", f"Unreadable table row: {line!r}")
            )
    if not data.modes and not data.issues:
        data.issues.append(Issue("error", "no_mode_table", "No mode table found in this log."))


def load_results(
    files: Mapping[str, str],
    names: Mapping[str, Mapping] | None = None,
) -> ParsedResults:
    """Parse and merge files given as {filename: text}.

    `names` lets the user supply names for a file, as
    {filename: {"target": ..., "ligand": ..., "seed": optional}}.
    """
    names = names or {}
    issues = []

    # Step 1: work out what each file claims to be, and group the pairs.
    groups: dict[tuple, dict[str, list[str]]] = {}
    for filename in sorted(files):
        parsed = _resolve_name(filename, names.get(filename))
        if parsed is None:
            issues.append(
                Issue(
                    "error",
                    "unrecognised_name",
                    "Filename doesn't match <target>__<ligand>[_seedN]_out.pdbqt or .log. "
                    "Target and ligand names are needed before it can be used.",
                    filename,
                )
            )
            continue
        key = (parsed.target, parsed.ligand, parsed.seed)
        groups.setdefault(key, {}).setdefault(parsed.kind, []).append(filename)

    # Step 2: parse each run.
    built = []
    for (target, ligand, name_seed), by_kind in groups.items():
        clashes = [found for found in by_kind.values() if len(found) > 1]
        if clashes:
            for filename in sum(clashes, []):
                issues.append(_duplicate_issue(target, ligand, name_seed, filename))
            continue
        built.append(
            _build_run(
                target,
                ligand,
                name_seed,
                by_kind.get("pdbqt", [None])[0],
                by_kind.get("log", [None])[0],
                files,
            )
        )

    # Step 3: two runs can still clash once the log's seed is known.
    counts: dict[tuple, int] = {}
    for run in built:
        counts[run["key"]] = counts.get(run["key"], 0) + 1

    pose_rows, run_rows = [], []
    for run in built:
        if counts[run["key"]] > 1:
            target, ligand, seed = run["key"]
            for filename in (run["run"]["pdbqt_file"], run["run"]["log_file"]):
                if filename:
                    issues.append(_duplicate_issue(target, ligand, seed, filename))
            continue
        pose_rows += run["poses"]
        run_rows.append(run["run"])
        issues += run["issues"]

    sort_keys = ["target", "ligand", "seed"]
    return ParsedResults(
        poses=_frame(pose_rows, POSE_COLUMNS).sort_values([*sort_keys, "pose"], ignore_index=True),
        runs=_frame(run_rows, RUN_COLUMNS).sort_values(sort_keys, ignore_index=True),
        issues=issues,
    )


def _duplicate_issue(target: str, ligand: str, seed: int | None, filename: str) -> Issue:
    return Issue(
        "error",
        "duplicate_run",
        f"More than one file claims {target} / {ligand} / seed {seed}. "
        "None of them were loaded; remove the extra file or rename it.",
        filename,
    )


def _resolve_name(filename: str, override: Mapping | None) -> ParsedName | None:
    parsed = parse_filename(filename)
    if not override or not override.get("target") or not override.get("ligand"):
        return parsed

    if filename.endswith(".pdbqt"):
        kind = "pdbqt"
    elif filename.endswith(".log"):
        kind = "log"
    else:
        return None
    return ParsedName(override["target"], override["ligand"], override.get("seed"), kind)


def _build_run(
    target: str,
    ligand: str,
    name_seed: int | None,
    pdbqt_file: str | None,
    log_file: str | None,
    files: Mapping[str, str],
) -> dict:
    """Merge one PDBQT and/or log into pose rows, a run row and any issues."""
    issues = []

    pdbqt = parse_pdbqt(files[pdbqt_file]) if pdbqt_file else None
    log = parse_log(files[log_file]) if log_file else None
    if pdbqt:
        issues += [replace(issue, file=pdbqt_file) for issue in pdbqt.issues]
    if log:
        issues += [replace(issue, file=log_file) for issue in log.issues]

    if log is None:
        issues.append(
            Issue(
                "info",
                "missing_log",
                "No matching log, so the search settings are not recorded.",
                pdbqt_file,
            )
        )
    if pdbqt is None:
        issues.append(
            Issue(
                "info",
                "missing_pdbqt",
                "No matching output PDBQT, so these modes have no coordinates or atom counts.",
                log_file,
            )
        )

    seed = name_seed
    if log and log.seed is not None:
        seed = log.seed
        if name_seed is not None and name_seed != log.seed:
            issues.append(
                Issue(
                    "warning",
                    "seed_conflict",
                    f"Filename says seed {name_seed}, log says {log.seed}. Using the log.",
                    log_file,
                )
            )

    pdbqt_poses = {pose.pose: pose for pose in pdbqt.poses} if pdbqt else {}
    log_modes = {mode.mode: mode for mode in log.modes} if log else {}
    if pdbqt_poses and log_modes:
        missing = sorted(set(pdbqt_poses) - set(log_modes))
        if missing:
            issues.append(
                Issue(
                    "warning",
                    "pose_not_in_log",
                    f"Poses {missing} are in the PDBQT but not in the log.",
                    pdbqt_file,
                )
            )

    largest_diff = None
    rows = []
    for number in sorted(set(pdbqt_poses) | set(log_modes)):
        pose = pdbqt_poses.get(number)
        mode = log_modes.get(number)
        if pose is not None and mode is not None and pose.affinity is not None:
            diff = max(
                abs(pose.affinity - mode.affinity),
                abs(pose.rmsd_lb - mode.rmsd_lb),
                abs(pose.rmsd_ub - mode.rmsd_ub),
            )
            largest_diff = diff if largest_diff is None else max(largest_diff, diff)
            if diff > AGREEMENT_TOLERANCE:
                issues.append(
                    Issue(
                        "warning",
                        "value_disagreement",
                        f"Pose {number}: PDBQT and log differ by {diff:.3f}. Using the PDBQT.",
                        pdbqt_file,
                    )
                )
        rows.append(_pose_row(target, ligand, seed, number, pose, mode, pdbqt_file, log_file, issues))

    issues += _atom_issues(pdbqt, pdbqt_file) if pdbqt else []

    run_row = {
        "target": target,
        "ligand": ligand,
        "seed": seed,
        "pdbqt_file": pdbqt_file,
        "log_file": log_file,
        "n_poses_pdbqt": len(pdbqt.poses) if pdbqt else None,
        "n_modes_log": len(log.modes) if log else None,
        "max_log_pdbqt_diff": largest_diff,
        **_header_values(log),
    }
    return {"key": (target, ligand, seed), "poses": rows, "run": run_row, "issues": issues}


def _pose_row(target, ligand, seed, number, pose, mode, pdbqt_file, log_file, issues) -> dict:
    row = {"target": target, "ligand": ligand, "seed": seed, "pose": number}

    # The explicit flag means later code never has to infer this from blank cells.
    row["has_coordinates"] = pose is not None
    row["source_file"] = pdbqt_file if pose is not None else log_file

    if pose is not None and pose.affinity is not None:
        source = pose  # PDBQT keeps 3 decimals, so it wins whenever it has a score
    else:
        source = mode
        if pose is not None:
            fallback = "Using the log value." if mode is not None else "No score is available."
            issues.append(
                Issue(
                    "warning",
                    "missing_score",
                    f"Pose {number} has no REMARK VINA RESULT line. {fallback}",
                    pdbqt_file,
                )
            )
    row["affinity"] = source.affinity if source else None
    row["rmsd_lb"] = source.rmsd_lb if source else None
    row["rmsd_ub"] = source.rmsd_ub if source else None

    row["heavy_atoms"] = None
    if pose is not None and pose.atom_types and not pose.has_flexible_residues:
        row["heavy_atoms"] = count_heavy_atoms(pose.atom_types)
    return row


def _atom_issues(pdbqt: PdbqtData, pdbqt_file: str) -> list[Issue]:
    """Report atom-count problems once per file, rather than once per pose."""
    issues = []
    unknown = set()
    for pose in pdbqt.poses:
        unknown |= _unknown_atom_types(pose.atom_types)
    if unknown:
        issues.append(
            Issue(
                "warning",
                "unknown_atom_type",
                f"Unrecognised atom types {sorted(unknown)}: heavy_atoms left blank.",
                pdbqt_file,
            )
        )
    if any(pose.has_flexible_residues for pose in pdbqt.poses):
        issues.append(
            Issue(
                "warning",
                "flexible_residues",
                "Output contains flexible residues, whose atoms aren't part of the ligand: "
                "heavy_atoms left blank.",
                pdbqt_file,
            )
        )
    if any(not pose.atom_types for pose in pdbqt.poses):
        issues.append(
            Issue("warning", "no_atoms", "A pose has no atoms: heavy_atoms left blank.", pdbqt_file)
        )
    counts = {
        count_heavy_atoms(pose.atom_types)
        for pose in pdbqt.poses
        if pose.atom_types and not pose.has_flexible_residues
    } - {None}
    if len(counts) > 1:
        issues.append(
            Issue(
                "warning",
                "inconsistent_heavy_atoms",
                f"Poses of one ligand have different heavy atom counts: {sorted(counts)}.",
                pdbqt_file,
            )
        )
    return issues


def _header_values(log: LogData | None) -> dict:
    center = log.center if log and log.center else (None, None, None)
    size = log.size if log and log.size else (None, None, None)
    return {
        "vina_version": log.vina_version if log else None,
        "scoring_function": log.scoring_function if log else None,
        "receptor": log.receptor if log else None,
        "log_ligand": log.ligand if log else None,
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": size[0],
        "size_y": size[1],
        "size_z": size[2],
        "grid_spacing": log.grid_spacing if log else None,
        "exhaustiveness": log.exhaustiveness if log else None,
        "cpu": log.cpu if log else None,
    }


def _frame(rows: list[dict], columns: Mapping[str, str]) -> pd.DataFrame:
    """Build a table with fixed columns and types, even when there are no rows."""
    frame = pd.DataFrame(rows, columns=list(columns), dtype="object")
    return frame.astype(dict(columns))


def load_folder(*folders: str | Path) -> ParsedResults:
    """Load every .pdbqt and .log directly inside the given folders (not subfolders)."""
    files: dict[str, str] = {}
    for folder in folders:
        for path in sorted(Path(folder).iterdir()):
            if path.suffix not in {".pdbqt", ".log"}:
                continue
            if path.name in files:
                raise ValueError(f"Two files are both named {path.name}; rename one to continue.")
            files[path.name] = decode_text(path.read_bytes())
    return load_results(files)
