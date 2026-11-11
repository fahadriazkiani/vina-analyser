"""Tests for app.py, driven by streamlit.testing.v1.AppTest.

AppTest runs the script without a browser. It cannot drive st.file_uploader or edit
st.data_editor cells, so:

* uploads are injected through st.session_state["input"], the same structure the
  uploader handler builds (the upload handling itself is tested in test_inputs.py);
* cell edits are covered by the names_from_table tests in test_inputs.py.

The real upload and editing paths are on the manual checklist in PROGRESS.md.
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


def test_recognition_table_shows_the_log_seed():
    at = fresh()
    button(at, "load_example").click().run()
    table = next(frame for frame in frames(at) if "recognised" in frame.columns).set_index("file")
    assert table.loc["1iep__imatinib_1iep.log", "seed"] == "42"
    assert table.loc["1iep__imatinib_1iep.log", "seed_source"] == "log header"
    assert "None" not in set(table["seed"])


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
    if git is not None:
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
