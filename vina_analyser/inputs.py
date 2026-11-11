"""Turn uploaded files into named files the parser can use, without guessing.

Two steps sit between an upload and the analysis:

1. read_uploads(): several loose files, or one .zip, become {filename: bytes}.
   Everything stays in memory, nothing uploaded is executed, and a zip is refused
   outright if it is unsafe (path traversal, too large, too many files).
2. recognition_table() / names_from_table(): the user sees what each file was read as
   and can correct it. Every included row is passed to the parser with explicit names,
   so what the table shows is exactly what gets analysed.

No Streamlit code here, so all of it is tested with plain pytest.
"""

import io
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

import pandas as pd

from vina_analyser.parsers import Issue, decode_text, parse_filename, read_log_seed

ACCEPTED_SUFFIXES = (".pdbqt", ".log")

# Kept equal to server.maxUploadSize in .streamlit/config.toml (a test checks this).
MAX_UPLOAD_MB = 50
# Real Vina outputs are tens of kilobytes. These limits are far above any genuine
# dataset and exist to stop a zip that unpacks to gigabytes (a "zip bomb").
MAX_FILE_BYTES = 20 * 1024 * 1024
MAX_TOTAL_BYTES = 200 * 1024 * 1024
MAX_FILES = 2000

RECOGNITION_COLUMNS = [
    "include", "file", "kind", "target", "ligand", "seed", "seed_source", "recognised", "filename_seed",
]
NOT_RECOGNISED = "no: type target and ligand"
NO_SEED = "—"

_DRIVE = re.compile(r"^[A-Za-z]:")


class InputRefused(ValueError):
    """The upload cannot be used, with the reason in the message."""


@dataclass
class Upload:
    files: dict[str, bytes] = field(default_factory=dict)
    issues: list[Issue] = field(default_factory=list)


# --------------------------------------------------------------------- reading


def read_uploads(items: list[tuple[str, bytes]]) -> Upload:
    """Read (name, bytes) pairs: either loose files or exactly one .zip."""
    zips = [name for name, _ in items if name.lower().endswith(".zip")]
    if zips:
        if len(items) > 1:
            raise InputRefused(
                "Upload a .zip on its own, or loose files without a .zip, "
                "so it is clear which files belong to the analysis."
            )
        return _check(_read_zip(items[0][1]))

    upload = Upload()
    for name, data in items:
        name = PurePosixPath(name.replace("\\", "/")).name
        if not name.lower().endswith(ACCEPTED_SUFFIXES):
            upload.issues.append(_skipped(name, "not a .pdbqt or .log file"))
            continue
        if name in upload.files:
            raise InputRefused(f"Two uploaded files are both named {name}; rename one.")
        if len(data) > MAX_FILE_BYTES:
            raise InputRefused(_too_large(name))
        upload.files[name] = data
    return _check(upload)


def _read_zip(data: bytes) -> Upload:
    try:
        archive = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile:
        raise InputRefused("The uploaded file is not a valid zip archive.") from None

    with archive:
        entries = [info for info in archive.infolist() if not info.is_dir()]
        # Check every path before reading anything, so an unsafe archive is refused whole.
        for info in entries:
            if _unsafe(info.filename):
                raise InputRefused(
                    f"The zip contains an unsafe path ({info.filename!r}), "
                    "so it was not opened."
                )
        if len(entries) > MAX_FILES:
            raise InputRefused(f"The zip holds {len(entries)} files; the limit is {MAX_FILES}.")

        upload, origin, total = Upload(), {}, 0
        for info in entries:
            path = info.filename.replace("\\", "/")
            name = PurePosixPath(path).name
            if path.startswith("__MACOSX/") or name.startswith("._"):
                continue  # macOS metadata, not user data
            if not name.lower().endswith(ACCEPTED_SUFFIXES):
                upload.issues.append(_skipped(path, "not a .pdbqt or .log file (nested zips are not opened)"))
                continue
            if info.flag_bits & 0x1:
                raise InputRefused(f"{path} is encrypted; upload an unencrypted zip.")
            if name in origin:
                raise InputRefused(
                    f"Two files in the zip become {name} once folders are removed: "
                    f"{origin[name]} and {path}. Rename one."
                )
            # Read at most one byte past the limit: the size in the zip header can lie.
            with archive.open(info) as handle:
                content = handle.read(MAX_FILE_BYTES + 1)
            if len(content) > MAX_FILE_BYTES:
                raise InputRefused(_too_large(path))
            total += len(content)
            if total > MAX_TOTAL_BYTES:
                raise InputRefused(
                    f"The zip unpacks to more than {MAX_TOTAL_BYTES // 2**20} MB, the limit."
                )
            origin[name] = path
            upload.files[name] = content
    return upload


def _unsafe(path: str) -> bool:
    path = path.replace("\\", "/")
    return path.startswith("/") or bool(_DRIVE.match(path)) or ".." in path.split("/")


def _check(upload: Upload) -> Upload:
    if not upload.files:
        raise InputRefused("No .pdbqt or .log files were found in the upload.")
    return upload


def _skipped(name: str, reason: str) -> Issue:
    return Issue("info", "skipped_file", f"Skipped: {reason}.", name)


def _too_large(name: str) -> str:
    return f"{name} is too large (limit {MAX_FILE_BYTES // 2**20} MB per file)."


def example_files(data_dir: str | Path) -> dict[str, bytes]:
    """The bundled example: every output PDBQT and log in vina_test_data/. Read only."""
    data_dir = Path(data_dir)
    files = {}
    for folder in ("results", "logs"):
        for path in sorted((data_dir / folder).iterdir()):
            if path.suffix in ACCEPTED_SUFFIXES:
                files[path.name] = path.read_bytes()
    return files


# ------------------------------------------------------------ recognition table


def recognition_table(files) -> pd.DataFrame:
    """One row per file: what it will be analysed as, with blanks where nothing is known.

    `files` is {filename: bytes or text}. The seed shown is the one the analysis will
    use: the log header's seed wins over the filename, as in the parser. Only the seed
    line of each log is read here, not the whole file. `filename_seed` stays in the
    table because the parser pairs a PDBQT with its log by filename seed.
    """
    parsed = {name: parse_filename(name) for name in files}
    log_seeds = {}
    for name, content in files.items():
        if name.lower().endswith(".log"):
            text = decode_text(content) if isinstance(content, bytes) else content
            log_seeds[name] = read_log_seed(text)
    # A PDBQT shows the seed of the log it will be paired with, if there is one.
    paired = {
        (p.target, p.ligand, p.seed): log_seeds[name]
        for name, p in parsed.items()
        if p is not None and p.kind == "log" and log_seeds.get(name) is not None
    }

    rows = []
    for name in sorted(files):
        name_parse = parsed[name]
        is_log = name.lower().endswith(".log")
        filename_seed = name_parse.seed if name_parse else None
        if is_log and log_seeds.get(name) is not None:
            seed, source = log_seeds[name], "log header"
        elif not is_log and name_parse and (name_parse.target, name_parse.ligand, filename_seed) in paired:
            seed, source = paired[(name_parse.target, name_parse.ligand, filename_seed)], "matching log"
        elif filename_seed is not None:
            seed, source = filename_seed, "filename"
        else:
            seed, source = None, "not recorded"
        rows.append(
            {
                "include": True,
                "file": name,
                "kind": "log" if is_log else "pdbqt",
                "target": name_parse.target if name_parse else "",
                "ligand": name_parse.ligand if name_parse else "",
                # Text, so a missing seed shows as a dash and never as "None".
                "seed": NO_SEED if seed is None else str(seed),
                "seed_source": source,
                "recognised": "yes" if name_parse else NOT_RECOGNISED,
                "filename_seed": filename_seed,
            }
        )
    table = pd.DataFrame(rows, columns=RECOGNITION_COLUMNS)
    return table.astype(
        {"include": "bool", "target": "str", "ligand": "str", "seed": "str", "filename_seed": "Int64"}
    )


def names_from_table(table: pd.DataFrame) -> tuple[dict[str, dict], list[str]]:
    """Explicit names for every included row, plus the problems that block analysis.

    The filename seed is passed, because the parser uses it to pair a PDBQT with its
    log; the parser then uses the log's seed when there is one, which is the seed the
    table shows.
    """
    names, problems = {}, []
    for row in table.itertuples(index=False):
        if not bool(row.include):
            continue
        target = _clean(row.target)
        ligand = _clean(row.ligand)
        if not target or not ligand:
            problems.append(f"{row.file}: type a target and a ligand, or untick include.")
            continue
        seed = None if pd.isna(row.filename_seed) else int(row.filename_seed)
        names[row.file] = {"target": target, "ligand": ligand, "seed": seed}
    return names, problems


def _clean(value) -> str:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return ""
    return str(value).strip()
