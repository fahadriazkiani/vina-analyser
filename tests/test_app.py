"""Tests for app.py, driven by streamlit.testing.v1.AppTest.

AppTest runs the script without a browser. It cannot drive st.file_uploader, so
uploads are injected through st.session_state["input"], the same structure the
uploader handler builds (the upload handling itself is tested in test_inputs.py).
Cell edits are made the way the browser makes them, through the data_editor's
widget state (see edit_table).

Dragging files into the uploader and typing into cells stay on the manual checklist.
"""

from pathlib import Path

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import json

from vina_analyser.analysis import LABELS, NOTES

ROOT = Path(__file__).resolve().parent.parent
APP = str(ROOT / "app.py")
DATA = ROOT / "vina_test_data"
TABS = ["Summary", "Figures", "Reproducibility", "RMSD", "Export"]
CONFIRM = "This reference is an experimental crystal/cryo-EM pose"


def fresh() -> AppTest:
    at = AppTest.from_file(APP, default_timeout=120)
    at.run()
    assert not at.exception, at.exception
    return at


def button(at, key):
    return next(b for b in at.button if b.key == key)


def all_text(at) -> str:
    parts = []
    for kind in ("title", "header", "subheader", "markdown", "caption", "info", "warning", "error", "success"):
        parts += [element.value for element in getattr(at, kind)]
    return "\n".join(str(part) for part in parts)


def frames(container) -> list[pd.DataFrame]:
    return [element.value for element in container.dataframe]


def frame_with(container, column) -> pd.DataFrame:
    """The table that has `column`; fails the test, never skips, if there is none."""
    for frame in frames(container):
        if column in frame.columns:
            return frame
    pytest.fail(f"No table with a {column!r} column is shown.")


@pytest.fixture(scope="module")
def analysed() -> AppTest:
    at = fresh()
    button(at, "load_example").click().run()
    button(at, "run_analysis").click().run()
    assert not at.exception, at.exception
    return at


def tab(at, label):
    return next(t for t in at.tabs if t.label == label)


# ------------------------------------------------------------------- start up


def test_app_starts_without_errors_and_shows_limitations():
    at = fresh()
    text = all_text(at)
    assert "not a binding affinity" in text
    assert "rigid" in text.lower()
    assert "across targets" in text
    assert "2–3 kcal/mol" in text


def test_run_is_disabled_before_any_input():
    at = fresh()
    assert button(at, "run_analysis").disabled


# ---------------------------------------------------------- example dataset


def test_example_dataset_fills_the_recognition_table_before_analysis():
    at = fresh()
    button(at, "load_example").click().run()
    assert not at.exception
    table = next(frame for frame in frames(at) if "recognised" in frame.columns)
    assert len(table) == 26
    assert not button(at, "run_analysis").disabled
    assert len(at.tabs) == 0  # nothing analysed until Run is pressed


def test_analysis_shows_the_five_tabs(analysed):
    assert [t.label for t in analysed.tabs] == TABS


def test_summary_shows_the_known_best_score(analysed):
    summary = next(
        frame for frame in frames(tab(analysed, "Summary")) if LABELS["affinity"] in frame.columns
    )
    row = summary[(summary["target"] == "1iep") & (summary["ligand"] == "imatinib_1iep")]
    assert row[LABELS["affinity"]].iloc[0] == pytest.approx(-13.273)


def test_issues_panel_is_always_shown(analysed):
    # Without a control the example data produces no issues, and the panel says so.
    assert "Issues" in all_text(analysed)
    assert "reported no issues" in all_text(analysed)


def test_no_control_is_chosen_by_default(analysed):
    assert analysed.selectbox(key="control").value == "No control"
    summary = next(
        frame for frame in frames(tab(analysed, "Summary")) if LABELS["affinity"] in frame.columns
    )
    assert LABELS["delta_vs_control"] not in summary.columns


def test_choosing_a_control_adds_the_delta_column():
    at = fresh()
    button(at, "load_example").click().run()
    button(at, "run_analysis").click().run()
    at.selectbox(key="control").set_value("imatinib_1iep").run()
    assert not at.exception
    summary = next(
        frame for frame in frames(tab(at, "Summary")) if LABELS["delta_vs_control"] in frame.columns
    )
    control = summary[summary["ligand"] == "imatinib_1iep"]
    assert (control[LABELS["delta_vs_control"]] == 0).all()
    assert "negative = more favourable" in all_text(at)
    # Imatinib has 4 runs on 1iep, the others 1: a warning, shown and listed.
    assert any("unequal_run_counts" in set(f.get("code", [])) for f in frames(at))
    assert any("the control has 4 run(s)" in w.value for w in at.warning)


def test_reproducibility_tab_shows_the_seed_summary(analysed):
    seeds = frames(tab(analysed, "Reproducibility"))
    assert any(LABELS["mean_affinity"] in frame.columns for frame in seeds)
    assert "2-3 kcal/mol" in all_text(analysed)


def test_figures_tab_shows_images_and_captions(analysed):
    figures = tab(analysed, "Figures")
    assert any("Compare ligands within a single row" in c.value for c in figures.caption)


# ----------------------------------------------------------------------- RMSD


def rmsd_ready() -> AppTest:
    at = fresh()
    button(at, "load_example").click().run()
    button(at, "run_analysis").click().run()
    button(at, "use_example_reference").click().run()
    at.selectbox(key="rmsd_target").set_value("1iep").run()
    at.selectbox(key="rmsd_ligand").set_value("imatinib_1iep").run()
    assert not at.exception, at.exception
    return at


def rmsd_tables(at):
    return [f for f in frames(tab(at, "RMSD")) if LABELS["rmsd_vs_reference"] in f.columns]


def test_confirmation_box_starts_unticked_and_blocks_rmsd():
    at = rmsd_ready()
    box = next(c for c in at.checkbox if c.label.startswith(CONFIRM))
    assert box.value is False
    assert rmsd_tables(at) == []
    assert "Tick the confirmation" in all_text(at)


def test_ticking_the_box_computes_the_known_rmsd():
    at = rmsd_ready()
    next(c for c in at.checkbox if c.label.startswith(CONFIRM)).check().run()
    assert not at.exception, at.exception
    (table,) = rmsd_tables(at)
    assert set(table["target"]) == {"1iep"}
    best = table[table["pose"] == 1].set_index("seed")[LABELS["rmsd_vs_reference"]]
    assert best[42] == pytest.approx(0.91, abs=0.02)
    assert best[2] == pytest.approx(0.34, abs=0.02)
    assert "no symmetry correction" in all_text(at)


def test_rmsd_needs_a_target_as_well_as_a_ligand():
    at = fresh()
    button(at, "load_example").click().run()
    button(at, "run_analysis").click().run()
    button(at, "use_example_reference").click().run()
    at.selectbox(key="rmsd_ligand").set_value("imatinib_1iep").run()
    next(c for c in at.checkbox if c.label.startswith(CONFIRM)).check().run()
    assert rmsd_tables(at) == []
    assert "Choose the target" in all_text(at)


# ------------------------------------------------------ unrecognised input


def test_unrecognised_file_blocks_analysis_until_named():
    at = AppTest.from_file(APP, default_timeout=120)
    text = (DATA / "results" / "1iep__p59_5x72_out.pdbqt").read_bytes()
    at.session_state["input"] = {"origin": "upload", "files": {"run7.pdbqt": text}, "issues": []}
    at.run()
    assert not at.exception
    assert button(at, "run_analysis").disabled
    assert "run7.pdbqt" in all_text(at)


def test_results_are_hidden_when_the_table_no_longer_matches_the_analysis():
    # Stands in for a cell edit after Run (AppTest cannot edit data_editor cells):
    # the stored analysis was made with names that differ from the table now shown.
    at = fresh()
    button(at, "load_example").click().run()
    button(at, "run_analysis").click().run()
    assert len(at.tabs) == 5
    snapshot = dict(at.session_state["analysed"])
    names = {file: dict(entry) for file, entry in snapshot["names"].items()}
    names["1iep__p59_5x72.log"]["ligand"] = "an_older_name"
    at.session_state["analysed"] = {**snapshot, "names": names}
    at.run()
    assert not at.exception
    assert len(at.tabs) == 0
    assert "table has changed" in all_text(at)


# --------------------------------------------------------- manual-check fixes


def captions(container) -> list[str]:
    return [element.value for element in container.caption]


def test_summary_shows_each_note_once_and_no_seed_note(analysed):
    shown = captions(tab(analysed, "Summary"))
    assert shown.count(NOTES["scoring_error"]) == 1
    assert shown.count(NOTES["affinity"]) == 1
    assert NOTES["seed_summary"] not in shown


def test_seed_note_is_on_the_reproducibility_tab_once(analysed):
    shown = captions(tab(analysed, "Reproducibility"))
    assert shown.count(NOTES["seed_summary"]) == 1
    assert shown.count(NOTES["scoring_error"]) == 1


def test_ligand_efficiency_is_shown_with_four_decimals(analysed):
    element = next(
        e for e in tab(analysed, "Summary").dataframe if LABELS["ligand_efficiency"] in e.value.columns
    )
    config = json.loads(element.proto.columns)
    assert config[LABELS["ligand_efficiency"]]["type_config"]["format"] == "%.4f"


def test_runs_preview_shows_the_log_seed_before_analysis():
    at = fresh()
    button(at, "load_example").click().run()
    editor = next(frame for frame in frames(at) if "recognised" in frame.columns).set_index("file")
    assert "seed" not in editor.columns  # the seed used depends on edits: see the preview
    assert editor.loc["1iep__imatinib_1iep_seed2.log", "seed_in_name"] == "2"
    preview = next(frame for frame in frames(at) if "seed used" in frame.columns)
    row = preview[preview["PDBQT file"] == "1iep__imatinib_1iep_out.pdbqt"].iloc[0]
    assert [row["seed used"], row["seed from"]] == ["42", "log header"]
    assert "Runs that will be analysed" in all_text(at)


def test_runs_preview_matches_the_runs_analysed(analysed):
    preview = frame_with(analysed, "seed used")
    runs = frame_with(tab(analysed, "Summary"), "Modes in log")
    assert len(preview) == len(runs) == 13
    key = ["target", "ligand", "PDBQT file"]
    shown = preview.sort_values(key)["seed used"].tolist()
    analysed_seeds = runs.sort_values(key)["seed"]
    assert shown == [str(int(seed)) for seed in analysed_seeds]


# --------------------------------------------------------------------- export


def downloads(at) -> set:
    return {element.key for element in at.get("download_button")}


def test_export_tab_offers_excel_methods_and_metadata(analysed):
    keys = downloads(analysed)
    assert {"download_xlsx", "download_methods", "download_metadata"} <= keys
    assert any("AutoDock Vina v1.2.5" in block.value for block in analysed.code)


def test_methods_in_the_app_credits_this_checkout(analysed):
    from vina_analyser import __version__
    from vina_analyser.report import git_info

    methods = next(block.value for block in analysed.code if "AutoDock Vina" in block.value)
    assert f"Vina Docking Results Analyser v{__version__}" in methods
    assert "is not recorded in Vina log files" in methods
    git = git_info(ROOT)
    assert git is not None, "the tests run from a git checkout, so the commit must be found"
    assert f"(commit {git['commit_short']}" in methods


def test_bundle_is_built_on_request_with_matching_hashes():
    import hashlib
    import io
    import json
    import zipfile

    at = fresh()
    button(at, "load_example").click().run()
    button(at, "run_analysis").click().run()
    assert "download_bundle" not in downloads(at)
    button(at, "prepare_bundle").click().run()
    assert not at.exception, at.exception
    assert {"download_bundle", "download_figures_zip"} <= downloads(at)

    bundle = at.session_state["prepared"]["bundle"]
    with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
        names = set(archive.namelist())
        meta = json.loads(archive.read("vina_results/run_metadata.json"))
    assert {"vina_results/results.xlsx", "vina_results/methods.txt", "vina_results/README.txt"} <= names
    assert "vina_results/figures/heatmap_best_affinity.png" in names
    for entry in meta["input_files"]:
        path = next(DATA.glob(f"*/{entry['name']}"))
        assert entry["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        # The names shown in the recognition table are what the metadata records.
        assert entry["included"] is True
        assert entry["target"] and entry["ligand"]
    seeds = {e["name"]: e["seed"] for e in meta["input_files"]}
    assert seeds["1iep__imatinib_1iep_seed2_out.pdbqt"] == 2


def test_prepared_bundle_is_dropped_when_the_control_changes():
    at = fresh()
    button(at, "load_example").click().run()
    button(at, "run_analysis").click().run()
    button(at, "prepare_bundle").click().run()
    at.selectbox(key="control").set_value("imatinib_1iep").run()
    assert "download_bundle" not in downloads(at)


# -------------------------------------------------------------------- privacy


def test_privacy_statement_is_shown_with_the_uploader():
    at = fresh()  # the start screen, before anything is uploaded
    statement = next(m.value for m in at.markdown if "Your files" in m.value)
    assert "processed in memory" in statement
    assert "never written to disk" in statement
    assert "expire" in statement and "hour" in statement
    assert "anyone with the link" in statement
    assert "own computer" in statement


def test_privacy_statement_comes_right_after_the_upload_header():
    # Where someone decides whether to upload: in section 1, not further down.
    at = fresh()
    order = [e.value for e in at.main if getattr(e, "type", "") in ("markdown", "header")]
    header = next(i for i, v in enumerate(order) if str(v).startswith("1."))
    statement = next(i for i, v in enumerate(order) if "Your files" in str(v))
    assert statement == header + 1


def test_every_cache_of_user_data_expires():
    # The privacy statement promises cached results expire; check every cache says so.
    import ast

    tree = ast.parse(Path(APP).read_text(encoding="utf-8"))
    caches = [
        (node.name, decorator)
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef)
        for decorator in node.decorator_list
        if isinstance(decorator, ast.Call) and ast.unparse(decorator.func) == "st.cache_data"
    ]
    assert len(caches) >= 3
    for name, decorator in caches:
        assert "ttl" in {keyword.arg for keyword in decorator.keywords}, name


# ------------------------------------------------------------ number formats


def formats_of(element) -> dict:
    config = json.loads(element.proto.columns or "{}")
    return {name: c.get("type_config", {}).get("format") for name, c in config.items()}


def test_scores_on_screen_use_the_same_decimals_as_excel(analysed):
    from vina_analyser.report import NUMBER_FORMATS

    printf = {"0.000": "%.3f", "0.0000": "%.4f", "0.00": "%.2f"}
    checked = 0
    for element in analysed.dataframe:
        formats = formats_of(element)
        for column in element.value.columns:
            if column in NUMBER_FORMATS:
                assert formats.get(column) == printf[NUMBER_FORMATS[column]], column
                checked += 1
    assert checked >= 5  # score, LE, mean, SD, range at least
    summary = next(e for e in tab(analysed, "Summary").dataframe if LABELS["affinity"] in e.value.columns)
    assert formats_of(summary)[LABELS["affinity"]] == "%.3f"


# ------------------------------------------------------------- table edits


def edit_table(at, edits: dict[str, dict]) -> None:
    """Edit cells of the recognition table the way the browser does: through the
    data_editor's widget state, keyed by row number."""
    import hashlib

    files = at.session_state["input"]["files"]
    digest = hashlib.sha256()
    for name in sorted(files):
        digest.update(name.encode() + b"\0" + files[name] + b"\0")
    table = frame_with(at, "recognised")
    rows = {list(table["file"]).index(file): change for file, change in edits.items()}
    at.session_state[f"editor_{digest.hexdigest()[:16]}"] = {
        "edited_rows": rows, "added_rows": [], "deleted_rows": [],
    }
    at.run()
    assert not at.exception, at.exception


def preview_row(at, pdbqt):
    preview = frame_with(at, "seed used")
    rows = preview[preview["PDBQT file"] == pdbqt]
    assert len(rows) == 1, pdbqt
    return rows.iloc[0]


def test_preview_follows_an_excluded_log():
    at = fresh()
    button(at, "load_example").click().run()
    edit_table(at, {"1iep__imatinib_1iep_seed2.log": {"include": False}})
    row = preview_row(at, "1iep__imatinib_1iep_seed2_out.pdbqt")
    assert [row["seed used"], row["seed from"], row["log file"]] == ["2", "filename", "—"]


def test_preview_follows_a_retyped_ligand_and_matches_the_analysis():
    at = fresh()
    button(at, "load_example").click().run()
    edit_table(at, {"1iep__imatinib_1iep_out.pdbqt": {"ligand": "renamed"}})
    row = preview_row(at, "1iep__imatinib_1iep_out.pdbqt")
    assert [row["ligand"], row["seed used"], row["log file"]] == ["renamed", "—", "—"]
    preview = frame_with(at, "seed used")
    button(at, "run_analysis").click().run()
    runs = frame_with(tab(at, "Summary"), "Modes in log")
    assert len(preview) == len(runs) == 14  # the log is now a run of its own
    assert sorted(preview["ligand"]) == sorted(runs["ligand"])


def test_heatmap_caption_in_the_app_says_when_scores_are_left_out():
    from vina_analyser.plots import MAX_ANNOTATED_LIGANDS

    pdbqt = (DATA / "results" / "1iep__p59_5x72_out.pdbqt").read_bytes()
    files = {f"1iep__lig{i:02d}_out.pdbqt": pdbqt for i in range(MAX_ANNOTATED_LIGANDS + 1)}
    at = AppTest.from_file(APP, default_timeout=120)
    at.session_state["input"] = {"origin": "upload", "files": files, "issues": []}
    at.run()
    button(at, "run_analysis").click().run()
    assert not at.exception, at.exception
    assert any("best-pose table" in c.value for c in tab(at, "Figures").caption)


# ------------------------------------------------------- reference uploader


def test_reference_uploader_hint_agrees_with_the_instruction(analysed):
    uploader = next(e for e in analysed.get("file_uploader") if e.proto.id.endswith("reference_upload"))
    # The hint under the uploader is generated from these types.
    assert list(uploader.proto.type) == [".pdbqt"]
    instruction = next(m.value for m in tab(analysed, "RMSD").markdown if "Heavy-atom RMSD" in m.value)
    assert "PDBQT references only" in instruction
    # The reason other formats are refused is stated beside the uploader.
    assert "SDF" in instruction and "atom mapping" in instruction
    assert "PDBQT" in uploader.proto.label
