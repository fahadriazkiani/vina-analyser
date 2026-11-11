"""Tests for vina_analyser.report: CSV, Excel, methods text, ZIPs and run metadata.

Real-data expectations come from vina_test_data/ (logs record Vina v1.2.5,
exhaustiveness 8, seeds 42 and 1/2/3, two grid boxes). Made-up variations are built
by copying the parsed tables and changing one value, never by touching the raw data.
"""

import codecs
import hashlib
import importlib.metadata
import io
import json
import re
import shutil
import subprocess
import tomllib
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import openpyxl
import pandas as pd
import pytest

import vina_analyser
from vina_analyser.analysis import LABELS, analyse, parse_reference, redocking_rmsd
from vina_analyser.parsers import load_folder, load_results
from vina_analyser.report import (
    LIMITATIONS,
    PACKAGES,
    Redocking,
    csv_bytes,
    deterministic_zip,
    display_table,
    figure_specs,
    figures_zip,
    git_info,
    methods_text,
    run_metadata,
    safe_name,
    workbook_bytes,
)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "vina_test_data"
REFERENCE = DATA / "ligands" / "imatinib_1iep.pdbqt"
ANALYSED_AT = datetime(2026, 11, 11, 9, 30, 0, tzinfo=timezone.utc)


@pytest.fixture(scope="module")
def real():
    return load_folder(DATA / "results", DATA / "logs")


@pytest.fixture(scope="module")
def raw_files():
    return {
        path.name: path.read_bytes()
        for folder in ("results", "logs")
        for path in sorted((DATA / folder).iterdir())
    }


@pytest.fixture(scope="module")
def redocking(real):
    reference = parse_reference(
        REFERENCE.read_text(encoding="utf-8"), is_experimental_pose=True, filename=REFERENCE.name
    )
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep", runs=real.runs, target="1iep")
    return Redocking(
        target="1iep",
        ligand="imatinib_1iep",
        reference_name=REFERENCE.name,
        reference_sha256=hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
        confirmed=True,
        table=result.table,
    )


# ------------------------------------------------------------------------ CSV


def test_csv_starts_with_the_utf8_byte_order_mark():
    # Excel on Windows only reads a CSV as UTF-8 when it starts with the BOM.
    data = csv_bytes(pd.DataFrame({"a": [1]}))
    assert data.startswith(codecs.BOM_UTF8)
    assert data.count(codecs.BOM_UTF8) == 1


def test_csv_round_trips_delta_angstrom_and_unicode_minus():
    frame = pd.DataFrame({"Δ score (kcal/mol)": ["−1.5"], "RMSD (Å)": [0.91]})
    data = csv_bytes(frame)
    assert "Δ".encode() in data and "Å".encode() in data and "−".encode() in data
    back = pd.read_csv(io.BytesIO(data), encoding="utf-8-sig")
    assert list(back.columns) == list(frame.columns)
    assert back.iloc[0, 0] == "−1.5"


def test_csv_has_no_index_column():
    data = csv_bytes(pd.DataFrame({"a": [1, 2]}))
    assert data.decode("utf-8-sig").splitlines()[0] == "a"


def test_display_table_uses_labels_with_units():
    frame = pd.DataFrame({"target": ["t"], "affinity": [-7.0], "delta_vs_control": [0.5]})
    shown = display_table(frame)
    assert list(shown.columns) == ["target", LABELS["affinity"], LABELS["delta_vs_control"]]
    assert list(frame.columns) == ["target", "affinity", "delta_vs_control"]  # input untouched


@pytest.mark.parametrize("dangerous", ["=1+1", "+1", "-1+cmd", "@SUM(A1)", "\t=1", "\r=1"])
def test_csv_text_that_a_spreadsheet_would_run_is_neutralised(dangerous):
    data = csv_bytes(pd.DataFrame({"ligand": [dangerous], "affinity": [-7.5]}))
    back = pd.read_csv(io.BytesIO(data), encoding="utf-8-sig", keep_default_na=False)
    assert back.loc[0, "ligand"] == "'" + dangerous
    assert back.loc[0, "affinity"] == -7.5  # numbers, including negative ones, untouched


# ------------------------------------------------------------------- methods


def test_methods_states_the_logged_settings(real):
    text = methods_text(real.runs).text
    assert "AutoDock Vina v1.2.5" in text
    assert "exhaustiveness of 8" in text
    assert "grid spacing of 0.375 Å" in text
    assert "vina scoring function" in text


def test_methods_lists_one_box_per_target_as_bullets(real):
    text = methods_text(real.runs).text
    assert "- 1iep: centre (15.190, 53.903, 16.917) Å; size 20.0 × 20.0 × 20.0 Å" in text
    assert "- 5x72: centre (−15.000, 15.000, 129.000) Å; size 30.0 × 24.0 × 24.0 Å" in text


def test_methods_with_one_target_states_the_box_in_the_paragraph(real):
    runs = real.runs[real.runs["target"] == "1iep"]
    text = methods_text(runs).text
    assert "centred at (15.190, 53.903, 16.917) Å with a size of 20.0 × 20.0 × 20.0 Å" in text
    assert "\n- 1iep:" not in text


def test_methods_lists_the_seeds_from_the_logs(real):
    text = methods_text(real.runs).text
    assert "random seed 42" in text
    assert "imatinib_1iep on 1iep: seeds 1, 2, 3 and 42" in text


def test_methods_reports_pose_counts_and_unrecorded_energy_range(real):
    text = methods_text(real.runs).text
    assert "The logs list 9 binding modes per run" in text
    assert "between 3 and 9 poses per run" in text
    assert "energy_range: not recorded in the logs" in text


def test_methods_never_uses_full_paths(real):
    text = methods_text(real.runs).text
    assert "receptors/" not in text and "ligands/" not in text
    assert "1iep_receptor.pdbqt" in text


def test_methods_cites_vina_and_asks_for_a_check(real):
    text = methods_text(real.runs).text
    assert "10.1021/acs.jcim.1c00203" in text and "10.1002/jcc.21334" in text
    assert "check it before use" in text


def test_seed_missing_from_the_log_is_not_recorded_not_guessed():
    log = (DATA / "logs" / "1iep__imatinib_1iep_seed2.log").read_text(encoding="utf-8")
    log = "\n".join(line for line in log.splitlines() if "random seed" not in line)
    pdbqt = (DATA / "results" / "1iep__imatinib_1iep_seed2_out.pdbqt").read_text(encoding="utf-8")
    runs = load_results(
        {"1iep__imatinib_1iep_seed2.log": log, "1iep__imatinib_1iep_seed2_out.pdbqt": pdbqt}
    ).runs
    text = methods_text(runs).text
    assert "random seed: not recorded" in text
    assert "seed 2" not in text and "seed 42" not in text


def test_value_that_differs_within_a_target_is_listed_with_a_warning(real):
    runs = real.runs.copy()
    index = runs.index[(runs["ligand"] == "p59_5x72") & (runs["target"] == "1iep")][0]
    runs.loc[index, "exhaustiveness"] = 32
    result = methods_text(runs)
    assert "exhaustiveness varied between runs: 8 (12 runs), 32 (p59_5x72 on 1iep, seed 42)" in result.text
    assert [issue.code for issue in result.issues] == ["methods_value_varies"]


def test_missing_header_values_say_not_recorded(real):
    runs = real.runs.copy()
    for column in ("vina_version", "exhaustiveness", "grid_spacing", "scoring_function"):
        runs[column] = None
    text = methods_text(runs).text
    assert "AutoDock Vina (version not recorded)" in text
    assert "exhaustiveness: not recorded" in text
    assert "grid spacing: not recorded" in text


def test_redocking_sentence_appears_for_a_confirmed_rmsd(real, redocking):
    text = methods_text(real.runs, redocking=redocking).text
    match = re.search(r"heavy-atom RMSD of (\d\.\d\d)–(\d\.\d\d) Å", text)
    assert match, text
    assert match.groups() == ("0.34", "0.91")
    assert "redocking imatinib_1iep into 1iep" in text
    assert "top-ranked pose of each of 4 runs (seeds 1, 2, 3 and 42)" in text
    assert "atoms matched in file order; no superposition or symmetry correction" in text
    assert "imatinib_1iep.pdbqt" in text


def test_redocking_sentence_for_a_single_run(real, redocking):
    single = Redocking(**{**redocking.__dict__, "table": redocking.table[redocking.table["seed"] == 42]})
    text = methods_text(real.runs, redocking=single).text
    assert "heavy-atom RMSD of 0.91 Å for the top-ranked pose" in text


def test_no_redocking_sentence_without_a_confirmed_rmsd(real, redocking):
    unconfirmed = Redocking(**{**redocking.__dict__, "confirmed": False})
    for value in (None, unconfirmed):
        text = methods_text(real.runs, redocking=value).text
        assert "redocking" not in text.lower()
        assert "RMSD" not in text


def test_no_redocking_sentence_when_no_rmsd_could_be_computed(real, redocking):
    blank = redocking.table.assign(rmsd_vs_reference=float("nan"))
    text = methods_text(real.runs, redocking=Redocking(**{**redocking.__dict__, "table": blank})).text
    assert "redocking" not in text.lower()


# ----------------------------------------------------------------- file names


@pytest.mark.parametrize("name", ["../evil", "a/b", "C:\\x", "..", "", "  "])
def test_safe_name_cannot_escape_a_folder(name):
    safe = safe_name(name)
    assert safe and "/" not in safe and "\\" not in safe and ":" not in safe
    assert safe not in {".", ".."}


def test_figure_stems_never_contain_a_folder_or_parent_reference():
    # Target names come from filenames or the recognition table, so they are untrusted.
    for spec in figure_specs(["../evil", "a/b", "C:\\x"]):
        assert "/" not in spec.stem and "\\" not in spec.stem and ":" not in spec.stem
        assert ".." not in spec.stem


def test_figure_stems_stay_unique_after_cleaning():
    stems = [spec.stem for spec in figure_specs(["a/b", "a_b"])]
    assert len(stems) == len(set(stems)) == 5
    assert stems[0] == "heatmap_best_affinity"


# ----------------------------------------------------------------------- zips


def test_deterministic_zip_is_byte_identical_and_sorted():
    entries = {"b.txt": b"2", "a/c.txt": b"1"}
    first, second = deterministic_zip(entries), deterministic_zip(dict(reversed(entries.items())))
    assert first == second
    with zipfile.ZipFile(io.BytesIO(first)) as archive:
        assert archive.namelist() == ["a/c.txt", "b.txt"]
        assert {info.date_time for info in archive.infolist()} == {(1980, 1, 1, 0, 0, 0)}


def test_deterministic_zip_refuses_unsafe_names():
    for name in ("../x", "/abs", "C:/x"):
        with pytest.raises(ValueError):
            deterministic_zip({name: b""})


def test_figures_zip_holds_png_svg_and_captions():
    figures = {"heatmap_best_affinity": {"png": b"P", "svg": b"S"}}
    captions = {"heatmap_best_affinity": "A caption."}
    with zipfile.ZipFile(io.BytesIO(figures_zip(figures, captions))) as archive:
        assert archive.namelist() == [
            "captions.txt", "heatmap_best_affinity.png", "heatmap_best_affinity.svg",
        ]
        assert "A caption." in archive.read("captions.txt").decode("utf-8")


# ---------------------------------------------------------------------- excel


def sheets_for(real, rmsd=None):
    result = analyse(real)
    sheets = {
        "Best poses": display_table(result.best),
        "Seed summary": display_table(result.seeds),
        "All poses": display_table(real.poses),
        "Runs": real.runs,
        "Issues": pd.DataFrame(columns=["severity", "code", "file", "message"]),
    }
    if rmsd is not None:
        sheets["RMSD"] = display_table(rmsd)
    return sheets


def load(data: bytes):
    return openpyxl.load_workbook(io.BytesIO(data))


def column_cells(sheet, header):
    headers = [cell.value for cell in sheet[1]]
    index = headers.index(header) + 1
    return [row[0] for row in sheet.iter_rows(min_row=2, min_col=index, max_col=index)]


def test_workbook_has_a_read_me_and_one_sheet_per_table(real):
    book = load(workbook_bytes(sheets_for(real), ["Read me line"], ANALYSED_AT))
    assert book.sheetnames == ["Read me", "Best poses", "Seed summary", "All poses", "Runs", "Issues"]
    readme = [row[0] for row in book["Read me"].iter_rows(values_only=True)]
    assert "Read me line" in readme
    assert all(text in readme for text in LIMITATIONS)


def test_workbook_rmsd_sheet_only_when_computed(real, redocking):
    book = load(workbook_bytes(sheets_for(real, redocking.table), [], ANALYSED_AT))
    assert book.sheetnames[-1] == "RMSD"
    cells = column_cells(book["RMSD"], LABELS["rmsd_vs_reference"])
    assert {cell.number_format for cell in cells} == {"0.00"}


def test_workbook_keeps_values_and_fixed_decimals(real):
    book = load(workbook_bytes(sheets_for(real), [], ANALYSED_AT))
    sheet = book["Best poses"]
    scores = column_cells(sheet, LABELS["affinity"])
    assert -13.273 in [cell.value for cell in scores]
    assert {cell.number_format for cell in scores} == {"0.000"}
    efficiency = column_cells(sheet, LABELS["ligand_efficiency"])
    assert {cell.number_format for cell in efficiency} == {"0.0000"}


def test_workbook_stores_formula_like_text_as_text():
    frame = pd.DataFrame({"ligand": ["=1+1", "@SUM(A1)", "plain"], "affinity": [-7.0, -6.0, -5.0]})
    book = load(workbook_bytes({"Best poses": frame}, [], ANALYSED_AT))
    cells = column_cells(book["Best poses"], "ligand")
    assert [cell.value for cell in cells] == ["=1+1", "@SUM(A1)", "plain"]
    assert {cell.data_type for cell in cells} == {"s"}  # "f" would mean a formula


def test_workbook_created_time_is_the_analysis_time(real):
    book = load(workbook_bytes(sheets_for(real), [], ANALYSED_AT))
    assert book.properties.created == ANALYSED_AT.replace(tzinfo=None)


def test_workbook_writes_blanks_for_missing_values():
    # Every real row has a target; a row that is blank everywhere is not written at all.
    frame = pd.DataFrame(
        {"target": ["t", "t"], "seed": pd.array([42, None], dtype="Int64"), "x": [1.0, float("nan")]}
    )
    book = load(workbook_bytes({"Runs": frame}, [], ANALYSED_AT))
    assert [cell.value for cell in column_cells(book["Runs"], "seed")] == [42, None]
    assert [cell.value for cell in column_cells(book["Runs"], "x")] == [1.0, None]


# ------------------------------------------------------------------ metadata


def metadata(real, raw_files, **overrides):
    names = {name: None for name in raw_files}
    from vina_analyser.inputs import names_from_table, recognition_table

    table = recognition_table(raw_files)
    names, _ = names_from_table(table)
    arguments = dict(
        files=raw_files, names=names, runs=real.runs, origin="example", control=None,
        redocking=None, analysed_at=ANALYSED_AT, repo_dir=ROOT,
    )
    arguments.update(overrides)
    return run_metadata(**arguments)


def test_metadata_hashes_every_input_file(real, raw_files):
    meta = metadata(real, raw_files)
    by_name = {entry["name"]: entry for entry in meta["input_files"]}
    assert set(by_name) == set(raw_files)
    for name, data in raw_files.items():
        assert by_name[name]["sha256"] == hashlib.sha256(data).hexdigest()
        assert by_name[name]["bytes"] == len(data)


def test_metadata_records_the_names_and_seed_each_file_was_analysed_as(real, raw_files):
    entry = {e["name"]: e for e in metadata(real, raw_files)["input_files"]}["1iep__imatinib_1iep_out.pdbqt"]
    assert (entry["included"], entry["target"], entry["ligand"], entry["seed"]) == (True, "1iep", "imatinib_1iep", 42)


def test_metadata_marks_excluded_files(real, raw_files):
    names = {name: value for name, value in metadata(real, raw_files)["input_files"][0].items()}
    excluded = next(iter(raw_files))
    from vina_analyser.inputs import names_from_table, recognition_table

    chosen, _ = names_from_table(recognition_table(raw_files))
    chosen.pop(excluded)
    entry = {e["name"]: e for e in metadata(real, raw_files, names=chosen)["input_files"]}[excluded]
    assert entry["included"] is False
    assert "target" not in entry
    assert names  # the first entry exists


def test_metadata_package_versions_match_the_installed_ones(real, raw_files):
    meta = metadata(real, raw_files)
    assert set(meta["packages"]) == set(PACKAGES)
    for package, version in meta["packages"].items():
        assert version == importlib.metadata.version(package)
    assert meta["python"].startswith("3.")


def test_metadata_app_version_matches_pyproject(real, raw_files):
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert vina_analyser.__version__ == project["version"]
    assert metadata(real, raw_files)["app_version"] == project["version"]


def test_metadata_timestamp_is_the_stored_analysis_time(real, raw_files):
    assert metadata(real, raw_files)["analysed_at"] == "2026-11-11T09:30:00+00:00"


def test_metadata_without_rmsd_or_control(real, raw_files):
    meta = metadata(real, raw_files)
    assert meta["control"] is None
    assert meta["rmsd_reference"] is None


def test_metadata_with_rmsd_records_the_confirmed_reference(real, raw_files, redocking):
    meta = metadata(real, raw_files, redocking=redocking, control="imatinib_1iep")
    assert meta["control"] == "imatinib_1iep"
    assert meta["rmsd_reference"] == {
        "file": "imatinib_1iep.pdbqt",
        "sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
        "target": "1iep",
        "ligand": "imatinib_1iep",
        "confirmed_experimental_pose": True,
    }


def test_metadata_is_json_without_absolute_paths(real, raw_files, redocking):
    text = json.dumps(metadata(real, raw_files, redocking=redocking))
    assert str(ROOT) not in text and str(ROOT).replace("\\", "/") not in text
    assert not re.search(r"[A-Za-z]:\\\\", text)


# ------------------------------------------------------------------------ git


def test_git_info_from_this_checkout():
    if shutil.which("git") is None:
        pytest.skip("git is not installed here")
    info = git_info(ROOT)
    full = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    assert info["commit"] == full
    assert re.fullmatch(r"[0-9a-f]{40}", info["commit"])
    assert full.startswith(info["commit_short"]) and 7 <= len(info["commit_short"]) < 40
    assert isinstance(info["dirty"], bool)


def test_git_info_records_a_dirty_working_tree(tmp_path):
    if shutil.which("git") is None:
        pytest.skip("git is not installed here")

    def git(*args):
        subprocess.run(
            ["git", "-c", "user.name=t", "-c", "user.email=t@example.invalid", *args],
            cwd=tmp_path, check=True, capture_output=True,
        )

    git("init", "-q")
    (tmp_path / "f.txt").write_text("1")
    git("add", "f.txt")
    git("commit", "-q", "-m", "c")
    assert git_info(tmp_path)["dirty"] is False
    (tmp_path / "f.txt").write_text("2")
    assert git_info(tmp_path)["dirty"] is True


def test_git_info_outside_a_checkout_is_none(tmp_path):
    assert git_info(tmp_path) is None


def test_git_info_without_git_installed_is_none(tmp_path):
    def missing(*args, **kwargs):
        raise FileNotFoundError("git")

    assert git_info(ROOT, runner=missing) is None


def test_metadata_omits_git_cleanly_when_unavailable(real, raw_files, tmp_path):
    assert "git" not in metadata(real, raw_files, repo_dir=tmp_path)
    if shutil.which("git"):
        assert "commit" in metadata(real, raw_files)["git"]
