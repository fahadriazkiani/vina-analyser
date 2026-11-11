"""Tests for vina_analyser.inputs: turning uploads into named files for the parser.

Uploads are built in memory (zip archives with io.BytesIO), so nothing is written to
disk and vina_test_data/ is only ever read.
"""

import io
import tomllib
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from vina_analyser.inputs import (
    MAX_FILE_BYTES,
    MAX_FILES,
    MAX_UPLOAD_MB,
    RECOGNITION_COLUMNS,
    InputRefused,
    example_files,
    names_from_table,
    read_uploads,
    recognition_table,
)
from vina_analyser.parsers import load_folder, load_results

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "vina_test_data"


def zip_of(entries: dict[str, bytes], **kwargs) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", **kwargs) as archive:
        for name, data in entries.items():
            archive.writestr(name, data)
    return buffer.getvalue()


def real_files() -> dict[str, bytes]:
    return {
        path.name: path.read_bytes()
        for folder in ("results", "logs")
        for path in sorted((DATA / folder).iterdir())
    }


def codes(issues):
    return [issue.code for issue in issues]


# ------------------------------------------------------------------ loose files


def test_loose_files_are_kept_by_name():
    upload = read_uploads([("a__b_out.pdbqt", b"x"), ("a__b.log", b"y")])
    assert upload.files == {"a__b_out.pdbqt": b"x", "a__b.log": b"y"}
    assert upload.issues == []


def test_unsupported_loose_file_is_skipped_and_reported():
    upload = read_uploads([("a__b.log", b"y"), ("notes.docx", b"z")])
    assert list(upload.files) == ["a__b.log"]
    assert codes(upload.issues) == ["skipped_file"]
    assert "notes.docx" == upload.issues[0].file


def test_two_loose_files_with_one_name_are_refused():
    with pytest.raises(InputRefused, match="a__b.log"):
        read_uploads([("a__b.log", b"1"), ("a__b.log", b"2")])


def test_nothing_usable_is_refused():
    with pytest.raises(InputRefused, match="No .pdbqt or .log"):
        read_uploads([("notes.docx", b"z")])


def test_oversized_file_is_refused():
    with pytest.raises(InputRefused, match="too large"):
        read_uploads([("a__b.log", b"x" * (MAX_FILE_BYTES + 1))])


# ------------------------------------------------------------------------- zip


def test_zip_gives_the_same_files_as_loose_uploads():
    files = real_files()
    upload = read_uploads([("all.zip", zip_of(files))])
    assert upload.files == files


def test_zip_folders_are_flattened_to_file_names():
    upload = read_uploads([("run.zip", zip_of({"results/a__b_out.pdbqt": b"x", "logs/a__b.log": b"y"}))])
    assert set(upload.files) == {"a__b_out.pdbqt", "a__b.log"}


def test_zip_mixed_with_other_files_is_refused():
    with pytest.raises(InputRefused, match="on its own"):
        read_uploads([("run.zip", zip_of({"a__b.log": b"y"})), ("c__d.log", b"z")])


def test_two_zip_entries_flattening_to_one_name_are_refused():
    archive = zip_of({"run1/a__b.log": b"1", "run2/a__b.log": b"2"})
    with pytest.raises(InputRefused, match="run1/a__b.log.*run2/a__b.log"):
        read_uploads([("run.zip", archive)])


@pytest.mark.parametrize("name", ["../evil.log", "a/../../evil.log", "/abs/evil.log", "C:/evil.log", "..\\evil.log"])
def test_zip_path_traversal_is_refused(name):
    with pytest.raises(InputRefused, match="unsafe path"):
        read_uploads([("run.zip", zip_of({name: b"x"}))])


def test_nested_zip_is_skipped_not_opened():
    archive = zip_of({"a__b.log": b"y", "inner.zip": zip_of({"c__d.log": b"z"})})
    upload = read_uploads([("run.zip", archive)])
    assert list(upload.files) == ["a__b.log"]
    assert codes(upload.issues) == ["skipped_file"]
    assert "inner.zip" in upload.issues[0].file


def test_macos_metadata_entries_are_skipped_quietly():
    archive = zip_of({"__MACOSX/._a__b.log": b"junk", "._a__b.log": b"junk", "a__b.log": b"y"})
    upload = read_uploads([("run.zip", archive)])
    assert list(upload.files) == ["a__b.log"]


def test_zip_with_too_many_files_is_refused():
    archive = zip_of({f"t__l{index}.log": b"" for index in range(MAX_FILES + 1)})
    with pytest.raises(InputRefused, match="files"):
        read_uploads([("run.zip", archive)])


def test_zip_bomb_is_refused_by_actual_size():
    # Highly compressible: tiny on disk, larger than the per-file limit once opened.
    archive = zip_of({"a__b.log": b"0" * (MAX_FILE_BYTES + 1)}, compression=zipfile.ZIP_DEFLATED)
    assert len(archive) < MAX_FILE_BYTES // 100
    with pytest.raises(InputRefused, match="too large"):
        read_uploads([("run.zip", archive)])


def test_corrupt_zip_is_refused():
    with pytest.raises(InputRefused, match="not a valid zip"):
        read_uploads([("run.zip", b"PK not really")])


# ----------------------------------------------------------- example dataset


def test_example_dataset_has_all_26_files():
    files = example_files(DATA)
    assert len(files) == 26
    assert files == real_files()


# --------------------------------------------------------- recognition table


def test_recognition_table_columns_and_rows():
    table = recognition_table(real_files())
    assert list(table.columns) == RECOGNITION_COLUMNS
    assert len(table) == 26
    assert table["include"].all()
    assert (table["recognised"] == "yes").all()


def test_recognition_reads_target_ligand_and_seed_from_the_name():
    table = recognition_table({"1iep__imatinib_1iep_seed2_out.pdbqt": b"", "5x72__p59_5x72.log": b""})
    rows = table.set_index("file")
    assert rows.loc["1iep__imatinib_1iep_seed2_out.pdbqt", ["kind", "target", "ligand", "seed", "seed_source"]].tolist() == [
        "pdbqt", "1iep", "imatinib_1iep", "2", "filename",
    ]
    assert rows.loc["5x72__p59_5x72.log", "kind"] == "log"
    # An empty log records no seed: shown as a dash, never as "None" or a guess.
    assert rows.loc["5x72__p59_5x72.log", ["seed", "seed_source"]].tolist() == ["—", "not recorded"]


def test_recognition_shows_the_seed_the_analysis_will_use():
    # The log's seed wins in the parser, so the table must show it too, otherwise it
    # contradicts the Summary table (blank here, 42 there).
    rows = recognition_table(real_files()).set_index("file")
    assert rows.loc["1iep__imatinib_1iep.log", ["seed", "seed_source"]].tolist() == ["42", "log header"]
    assert rows.loc["1iep__imatinib_1iep_out.pdbqt", ["seed", "seed_source"]].tolist() == ["42", "matching log"]
    assert rows.loc["1iep__imatinib_1iep_seed2.log", ["seed", "seed_source"]].tolist() == ["2", "log header"]
    assert set(rows["seed_source"]) == {"log header", "matching log"}


def test_seed_column_never_prints_none():
    table = recognition_table({**real_files(), "run7.pdbqt": b"", "1iep__x_out.pdbqt": b""})
    assert not table["seed"].isin(["None", "nan", "<NA>", ""]).any()
    assert table.set_index("file").loc["run7.pdbqt", "seed"] == "—"


def test_shown_seed_matches_the_seed_after_parsing():
    files = {name: data.decode() for name, data in real_files().items()}
    table = recognition_table(files)
    names, _ = names_from_table(table)
    runs = load_results(files, names=names).runs
    parsed = {
        file: str(int(seed))
        for _, run in runs.iterrows()
        for file, seed in ((run["pdbqt_file"], run["seed"]), (run["log_file"], run["seed"]))
        if file
    }
    shown = dict(zip(table["file"], table["seed"]))
    assert shown == parsed


def test_unrecognised_name_is_left_blank_not_guessed():
    table = recognition_table({"run7.pdbqt": b""})
    row = table.iloc[0]
    assert row["target"] == "" and row["ligand"] == ""
    assert row["recognised"].startswith("no")
    assert row["kind"] == "pdbqt"


def test_every_included_row_gets_explicit_names():
    table = recognition_table(real_files())
    names, problems = names_from_table(table)
    assert problems == []
    assert set(names) == set(real_files())
    assert names["1iep__imatinib_1iep_seed2_out.pdbqt"] == {
        "target": "1iep", "ligand": "imatinib_1iep", "seed": 2,
    }
    assert names["1iep__p59_5x72.log"]["seed"] is None


def test_excluded_rows_are_left_out():
    table = recognition_table(real_files())
    table.loc[table["file"] == "1iep__p59_5x72.log", "include"] = False
    names, problems = names_from_table(table)
    assert "1iep__p59_5x72.log" not in names
    assert problems == []


def test_included_row_without_names_is_a_problem():
    table = recognition_table({"run7.pdbqt": b"", "1iep__p59_5x72.log": b""})
    names, problems = names_from_table(table)
    assert len(problems) == 1 and "run7.pdbqt" in problems[0]
    assert "run7.pdbqt" not in names


def test_typed_names_are_stripped_and_used():
    table = recognition_table({"run7.pdbqt": b""})
    table.loc[0, ["target", "ligand"]] = ["  1iep ", " mine"]
    names, problems = names_from_table(table)
    assert problems == []
    assert names["run7.pdbqt"] == {"target": "1iep", "ligand": "mine", "seed": None}


def test_edited_names_of_a_recognised_file_are_what_gets_analysed():
    files = {name: data.decode() for name, data in real_files().items() if "p59" in name and "1iep" in name}
    table = recognition_table(files)
    table["ligand"] = "renamed_ligand"
    names, _ = names_from_table(table)
    results = load_results(files, names=names)
    assert set(results.poses["ligand"]) == {"renamed_ligand"}


def test_renamed_file_with_typed_names_matches_the_original():
    original = load_folder(DATA / "results", DATA / "logs")
    files = {name: data.decode() for name, data in real_files().items()}
    files["mystery.pdbqt"] = files.pop("1iep__p59_5x72_out.pdbqt")
    table = recognition_table(files)
    table.loc[table["file"] == "mystery.pdbqt", ["target", "ligand"]] = ["1iep", "p59_5x72"]
    names, problems = names_from_table(table)
    assert problems == []
    results = load_results(files, names=names)
    expected = original.poses.drop(columns="source_file")
    actual = results.poses.drop(columns="source_file")
    pd.testing.assert_frame_equal(actual, expected)


# ---------------------------------------------------------- server settings


def test_streamlit_upload_limit_matches_the_package_limit():
    config = tomllib.loads((ROOT / ".streamlit" / "config.toml").read_text(encoding="utf-8"))
    assert config["server"]["maxUploadSize"] == MAX_UPLOAD_MB


def test_theme_is_light_and_uses_the_figure_blue():
    # Figures are white-background PNGs, so a dark theme would frame them badly; the
    # primary colour is the Okabe-Ito blue of the bars, leaving red to mean "problem".
    from vina_analyser.plots import BAR_COLOUR

    theme = tomllib.loads((ROOT / ".streamlit" / "config.toml").read_text(encoding="utf-8"))["theme"]
    assert theme["base"] == "light"
    assert theme["primaryColor"].upper() == BAR_COLOUR.upper()
    assert theme["backgroundColor"].upper() == "#FFFFFF"
