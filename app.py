"""Streamlit interface for the Vina Docking Results Analyser.

UI only: every calculation lives in the vina_analyser package, which is tested on
its own. Run with:  streamlit run app.py
"""

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from vina_analyser import analysis, inputs, plots, report
from vina_analyser.parsers import ParsedResults, decode_text, load_results

ROOT = Path(__file__).resolve().parent
EXAMPLE_DIR = ROOT / "vina_test_data"
EXAMPLE_REFERENCE = EXAMPLE_DIR / "ligands" / "imatinib_1iep.pdbqt"

NO_CONTROL = "No control"
CONFIRM_LABEL = (
    "This reference is an experimental crystal/cryo-EM pose "
    "(not a generated or docked 3D conformer)"
)
# Each cached entry holds parsed files for one dataset. A small limit stops memory
# from growing every time someone uploads a new set of files.
CACHE_ENTRIES = 8
FIGURE_CACHE_ENTRIES = 32
FIGURE_WIDTH_PX = 640

# Cached parses and figures hold users' data, so they expire; the privacy statement
# next to the uploader promises this (a test checks every cache has a ttl).
DATA_CACHE_SECONDS = 3600
GIT_CACHE_SECONDS = 60

PRIVACY_NOTE = (
    "**Your files.** Uploaded files are processed in memory and never written to disk. "
    "Nothing is kept after your session ends, except a temporary in-memory cache of "
    "parsed results and figures that expires within an hour. If this app is running on "
    "a public server, anyone with the link can use it, but other users cannot see your "
    "files or results. For unpublished or confidential results, you can run it on your "
    "own computer instead (see the README)."
)


def limitations_markdown() -> str:
    """The limitations box, from the same list the Excel and ZIP exports use."""
    lines = ["**Read before interpreting any number on this page**", ""]
    for text in report.LIMITATIONS:
        head, _, rest = text.partition(". ")
        lines.append(f"- **{head}.** {rest}" if rest else f"- {text}")
    return "\n".join(lines)


RUN_COLUMNS = {
    "target": "target",
    "ligand": "ligand",
    "seed": "seed",
    "n_poses_pdbqt": "Poses in output PDBQT",
    "n_modes_log": "Modes in log",
    "max_log_pdbqt_diff": "Largest log vs PDBQT difference",
    "pdbqt_file": "PDBQT file",
    "log_file": "Log file",
}


# ---------------------------------------------------------------- cached work


@st.cache_data(max_entries=CACHE_ENTRIES, ttl=DATA_CACHE_SECONDS, show_spinner="Reading the docking files…")
def parse(files: dict[str, str], names: dict[str, dict]) -> ParsedResults:
    return load_results(files, names=names)


@st.cache_data(max_entries=FIGURE_CACHE_ENTRIES, ttl=DATA_CACHE_SECONDS, show_spinner="Drawing the figure…")
def figure_files(kind: str, table: pd.DataFrame, target: str | None, control: str | None) -> dict:
    """PNG and SVG bytes of one figure. The screen shows exactly the PNG offered for download."""
    if kind == "heatmap":
        figure = plots.affinity_heatmap(table)
    elif kind == "ranked":
        figure = plots.ranked_bars(table, target, control=control)
    else:
        figure = plots.pose_strip(table, target)
    try:
        return {"png": plots.figure_bytes(figure, "png"), "svg": plots.figure_bytes(figure, "svg")}
    finally:
        plt.close(figure)  # Streamlit reruns the script often; never leak figures


def content_key(files: dict[str, bytes]) -> str:
    digest = hashlib.sha256()
    for name in sorted(files):
        digest.update(name.encode() + b"\0" + files[name] + b"\0")
    return digest.hexdigest()[:16]


# ---------------------------------------------------------------- state


def set_input(origin: str, files: dict[str, bytes], issues: list) -> None:
    st.session_state["input"] = {"origin": origin, "files": files, "issues": issues}
    st.session_state["input_error"] = None
    st.session_state.pop("analysed", None)
    st.session_state.pop("reference", None)


def handle_uploads() -> None:
    uploaded = st.file_uploader(
        "Vina output files (.pdbqt and .log), or one .zip containing them",
        type=["pdbqt", "log", "zip"],
        accept_multiple_files=True,
        key="uploader",
    )
    signature = tuple((item.name, item.size, item.file_id) for item in uploaded or [])
    if signature == st.session_state.get("upload_signature", ()):
        return
    st.session_state["upload_signature"] = signature
    if not signature:
        current = st.session_state.get("input")
        if current and current["origin"] == "upload":
            st.session_state["input"] = None  # the user removed their files
        return
    try:
        upload = inputs.read_uploads([(item.name, item.getvalue()) for item in uploaded])
    except inputs.InputRefused as error:
        st.session_state["input"] = None
        st.session_state["input_error"] = str(error)
        return
    set_input("upload", upload.files, upload.issues)


# ---------------------------------------------------------------- display helpers


def show_issues(issues: list, heading: bool = True) -> None:
    """Every issue, never collapsed: errors and warnings as banners, all in a table."""
    if heading:
        st.subheader("Issues")
    if not issues:
        if heading:
            st.success("The parser reported no issues for these files.")
        return
    counts = {level: sum(issue.severity == level for issue in issues) for level in ("error", "warning", "info")}
    st.markdown(f"**{counts['error']} errors · {counts['warning']} warnings · {counts['info']} notes**")
    for issue in issues:
        where = f"`{issue.file}`: " if issue.file else ""
        if issue.severity == "error":
            st.error(where + issue.message)
        elif issue.severity == "warning":
            st.warning(where + issue.message)
    st.dataframe(issue_frame(issues), hide_index=True, width="stretch")


def issue_frame(issues: list) -> pd.DataFrame:
    return pd.DataFrame(
        [(i.severity, i.code, i.file, i.message) for i in issues],
        columns=["severity", "code", "file", "message"],
    )


def show_notes(notes) -> None:
    for note in notes:
        st.caption(note)


def show_figure(kind: str, table: pd.DataFrame, caption: str, target=None, control=None) -> None:
    try:
        files = figure_files(kind, table, target, control)
    except ValueError as error:
        st.info(str(error))
        return
    st.image(files["png"], width=FIGURE_WIDTH_PX)
    st.caption(caption)


def best_view(best: pd.DataFrame) -> pd.DataFrame:
    order = [
        "target", "rank", "ligand", "affinity", "ligand_efficiency", "delta_vs_control",
        "n_runs", "control_n_runs", "control_affinity", "seed", "pose", "heavy_atoms", "source_file",
    ]
    return report.display_table(best[[c for c in order if c in best.columns]])


def seeds_view(seeds: pd.DataFrame) -> pd.DataFrame:
    table = seeds.copy()
    table["seeds"] = table["seeds"].map(lambda values: ", ".join(str(v) for v in values))
    return report.display_table(table)


# ---------------------------------------------------------------- page


def main() -> None:
    st.set_page_config(page_title="Vina Docking Results Analyser", layout="wide")
    st.title("Vina Docking Results Analyser")
    with st.container(border=True):
        st.markdown(limitations_markdown())

    # 1. Input --------------------------------------------------------------
    st.header("1. Files")
    st.markdown(PRIVACY_NOTE)
    handle_uploads()
    if st.button("Load example dataset", key="load_example", help="26 real AutoDock Vina 1.2.5 files from vina_test_data/"):
        set_input("example", inputs.example_files(EXAMPLE_DIR), [])
    if st.session_state.get("input_error"):
        st.error(st.session_state["input_error"])

    data = st.session_state.get("input")
    names, problems, key = {}, [], None
    if data:
        if data["origin"] == "example":
            st.info(
                "Example dataset loaded. It is a software test set, not a study: several "
                "runs are deliberate cross-dockings, so its scores mean nothing biologically."
            )
        if data["issues"]:
            show_issues(data["issues"], heading=False)

        # 2. Recognition table ---------------------------------------------
        st.header("2. Check what each file was read as")
        st.markdown(
            "Target, ligand and seed are read from names like "
            "`<target>__<ligand>[_seedN]_out.pdbqt` and `.log`. Correct any row, type names "
            "where the filename was not recognised, or untick **include**. "
            "The runs this produces, with the seed each one will use, are listed below it."
        )
        key = content_key(data["files"])
        edited = st.data_editor(
            inputs.recognition_table(data["files"]),
            key=f"editor_{key}",
            hide_index=True,
            width="stretch",
            disabled=["file", "kind", "seed_in_name", "recognised"],
            column_order=[c for c in inputs.RECOGNITION_COLUMNS if c != "filename_seed"],
            column_config={
                "include": st.column_config.CheckboxColumn("include"),
                "seed_in_name": st.column_config.TextColumn("seed in filename"),
            },
        )
        names, problems = inputs.names_from_table(edited)
        if not names and not problems:
            problems = ["No file is included."]
        for problem in problems:
            st.warning(problem)
        if not problems:
            show_runs_preview(data, names)

    blocked = not data or bool(problems)
    if st.button("Run analysis", key="run_analysis", type="primary", disabled=blocked):
        # The analysis time goes into run_metadata.json and the workbook, so it is
        # stored once here and not recomputed on every rerun.
        st.session_state["analysed"] = {"key": key, "names": names, "at": datetime.now(timezone.utc)}
    if not data:
        st.caption("Upload files or load the example dataset to begin.")
        return

    analysed = st.session_state.get("analysed")
    if not analysed:
        return
    if analysed["key"] != key or analysed["names"] != names:
        st.info("The table has changed since the last analysis. Press **Run analysis** to update the results.")
        return

    # 3. Results ---------------------------------------------------------------
    texts = {name: decode_text(data["files"][name]) for name in names}
    results = parse(texts, names)
    show_results(results, data, analysed)


def show_runs_preview(data: dict, names: dict) -> None:
    """The runs as the analysis will see them: the same cached parse the results use."""
    st.subheader("Runs that will be analysed")
    st.caption(
        "Exactly what gets analysed, recomputed from the table above after every edit. "
        "The seed comes from the log header when there is one, otherwise from the filename."
    )
    texts = {name: decode_text(data["files"][name]) for name in names}
    preview = inputs.runs_preview(parse(texts, names).runs)
    st.dataframe(preview, hide_index=True, width="stretch")


def show_results(results: ParsedResults, data: dict, analysed: dict) -> None:
    st.header("3. Results")
    ligands = sorted(results.poses["ligand"].unique())
    choice = st.selectbox(
        "Control compound (for example a known drug)",
        [NO_CONTROL, *ligands],
        key="control",
        help="Adds the score difference (Δ) to the control within each target. Nothing is chosen for you.",
    )
    control = None if choice == NO_CONTROL else choice
    result = analysis.analyse(results, control=control)
    show_issues(data["issues"] + result.issues)

    summary, figures, reproducibility, rmsd, export = st.tabs(
        ["Summary", "Figures", "Reproducibility", "RMSD", "Export"]
    )
    best = result.best

    with summary:
        if best.empty:
            st.info("No pose has coordinates, so there is nothing to rank. Upload the output PDBQT files.")
        else:
            st.dataframe(
                best_view(best), hide_index=True, width="stretch",
                column_config={
                    # Fixed decimals: otherwise 0.3190 shows as 0.319 beside 0.3587.
                    analysis.LABELS["ligand_efficiency"]: st.column_config.NumberColumn(format="%.4f"),
                },
            )
            # The seed note belongs with the seed table on the Reproducibility tab.
            show_notes(note for note in result.notes if note != analysis.NOTES["seed_summary"])
        st.subheader("Runs")
        st.caption(
            "Vina logs list every mode, but the output PDBQT only keeps poses within "
            "energy_range (3 kcal/mol by default) of the best one, so the two counts can "
            "differ. That is expected, not an error. Where both exist, the 3-decimal "
            "PDBQT score is used."
        )
        st.dataframe(results.runs[list(RUN_COLUMNS)].rename(columns=RUN_COLUMNS), hide_index=True, width="stretch")

    targets = sorted(best["target"].unique())
    with figures:
        if best.empty:
            st.info("No figures: no pose has coordinates.")
        else:
            show_figure("heatmap", best, plots.heatmap_caption())
            target = st.selectbox("Target", targets, key="figure_target")
            show_figure("ranked", best, plots.ranked_bars_caption(target), target=target, control=control)
            show_figure("strip", results.poses, plots.pose_strip_caption(target), target=target)

    with reproducibility:
        seeds = result.seeds
        repeated = seeds[seeds["n_runs"] > 1]
        single = len(seeds) - len(repeated)
        if repeated.empty:
            st.info("No target and ligand pair was run more than once, so there is no spread to report.")
        else:
            st.dataframe(seeds_view(repeated), hide_index=True, width="stretch")
        if single:
            st.caption(f"{single} target and ligand pairs have a single run: no spread available, so they are not listed.")
        st.caption(analysis.NOTES["seed_summary"])
        st.caption(analysis.NOTES["scoring_error"])

    with rmsd:
        redocking = show_rmsd(results, data)

    with export:
        show_export(result, results, redocking, targets, control, data, analysed)


def show_rmsd(results: ParsedResults, data: dict):
    st.markdown(
        "Heavy-atom RMSD of each docked pose against a reference ligand pose. Only "
        "meaningful when the reference is an **experimental** pose of the same ligand in "
        "the **same target and grid box**. PDBQT references only: SDF needs atom "
        "mapping, which is not supported."
    )
    uploaded = st.file_uploader(
        "Reference ligand (PDBQT, heavy atoms in the same order as the docked ligand)",
        type=["pdbqt", "sdf", "mol", "mol2"],  # SDF etc. accepted only to explain the refusal
        key="reference_upload",
    )
    signature = (uploaded.name, uploaded.size, uploaded.file_id) if uploaded else None
    if signature and signature != st.session_state.get("reference_signature"):
        st.session_state["reference_signature"] = signature
        content = uploaded.getvalue()
        st.session_state["reference"] = {"name": uploaded.name, "bytes": content, "text": decode_text(content)}
    if data["origin"] == "example" and st.button(
        "Use the example reference (ligands/imatinib_1iep.pdbqt, crystal coordinates of imatinib in 1IEP)",
        key="use_example_reference",
    ):
        content = EXAMPLE_REFERENCE.read_bytes()
        st.session_state["reference"] = {"name": EXAMPLE_REFERENCE.name, "bytes": content, "text": decode_text(content)}

    reference = st.session_state.get("reference")
    if not reference:
        st.info("Upload a reference ligand to compute RMSD.")
        return None
    st.markdown(f"Reference: `{reference['name']}`")

    targets = sorted(results.poses["target"].unique())
    target = st.selectbox(
        "Target the reference pose belongs to", targets, index=None,
        placeholder="Choose the target", key="rmsd_target",
    )
    ligand_pool = results.poses if target is None else results.poses[results.poses["target"] == target]
    ligand = st.selectbox(
        "Docked ligand the reference describes", sorted(ligand_pool["ligand"].unique()), index=None,
        placeholder="Choose the ligand", key="rmsd_ligand",
    )
    # Keyed to the file content: a new reference always starts unticked.
    reference_key = hashlib.sha256(reference["text"].encode()).hexdigest()[:16]
    confirmed = st.checkbox(CONFIRM_LABEL, value=False, key=f"confirm_{reference_key}")

    missing = []
    if target is None:
        missing.append(
            "Choose the target the reference belongs to. The same reference gives a "
            "meaningful number in its own target and a meaningless one in any other."
        )
    if ligand is None:
        missing.append("Choose the docked ligand the reference describes.")
    if not confirmed:
        missing.append(
            "Tick the confirmation box. No RMSD is computed until you confirm that the "
            "reference is an experimental pose."
        )
    if missing:
        for message in missing:
            st.info(message)
        return None

    try:
        parsed = analysis.parse_reference(reference["text"], is_experimental_pose=confirmed, filename=reference["name"])
        outcome = analysis.redocking_rmsd(results, parsed, ligand=ligand, runs=results.runs, target=target)
    except analysis.ReferenceRefused as error:
        st.error(str(error))
        return None

    show_issues(outcome.issues, heading=False)
    st.dataframe(report.display_table(outcome.table), hide_index=True, width="stretch")
    show_notes(outcome.notes)
    st.caption("Atoms are matched in file order; no superposition and no symmetry correction are applied.")
    return report.Redocking(
        target=target,
        ligand=ligand,
        reference_name=reference["name"],
        reference_sha256=hashlib.sha256(reference["bytes"]).hexdigest(),
        confirmed=confirmed,
        table=outcome.table,
    )


@st.cache_data(ttl=GIT_CACHE_SECONDS, show_spinner=False)
def cached_git_info() -> dict | None:
    # Three quick git calls; cached briefly so they don't run on every widget click.
    return report.git_info(ROOT)


SHEET_NAMES = {
    "best_poses.csv": "Best poses",
    "seed_summary.csv": "Seed summary",
    "all_poses.csv": "All poses",
    "runs.csv": "Runs",
    "issues.csv": "Issues",
    "rmsd.csv": "RMSD",
}


def export_tables(result, results: ParsedResults, redocking) -> dict:
    tables = {
        "best_poses.csv": best_view(result.best) if not result.best.empty else None,
        "seed_summary.csv": seeds_view(result.seeds),
        "all_poses.csv": report.display_table(results.poses),
        "runs.csv": results.runs,
        "issues.csv": issue_frame(result.issues),
        "rmsd.csv": report.display_table(redocking.table) if redocking is not None else None,
    }
    return {name: table for name, table in tables.items() if table is not None}


def show_export(result, results: ParsedResults, redocking, targets, control, data, analysed) -> None:
    tables = export_tables(result, results, redocking)
    methods = report.methods_text(results.runs, redocking=redocking, git=cached_git_info())

    st.subheader("Methods paragraph")
    st.caption(
        'Built only from values in the Vina logs; anything missing says "not recorded". '
        "The redocking sentence appears only after a confirmed RMSD on the RMSD tab."
    )
    st.code(methods.text, language=None, wrap_lines=True)
    show_issues(methods.issues, heading=False)

    readme = ["Notes", *result.notes, "", "Methods", *methods.text.splitlines()]
    workbook = report.workbook_bytes(
        {SHEET_NAMES[name]: table for name, table in tables.items()}, readme, analysed["at"]
    )
    metadata = report.metadata_json(
        report.run_metadata(
            files=data["files"],
            names=analysed["names"],
            runs=results.runs,
            origin=data["origin"],
            control=control,
            redocking=redocking,
            analysed_at=analysed["at"],
            git=cached_git_info(),
        )
    )

    st.subheader("Downloads")
    st.markdown("Tables are CSV in UTF-8 with a byte-order mark, so Excel shows Δ, Å and − correctly.")
    st.download_button(
        "Download results.xlsx (one sheet per table)", workbook, file_name="results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download_xlsx", on_click="ignore",
    )
    st.download_button(
        "Download methods.txt", methods.text.encode("utf-8"), file_name="methods.txt",
        mime="text/plain", key="download_methods", on_click="ignore",
    )
    st.download_button(
        "Download run_metadata.json", metadata, file_name="run_metadata.json",
        mime="application/json", key="download_metadata", on_click="ignore",
    )
    for filename, table in tables.items():
        st.download_button(
            f"Download {filename}", report.csv_bytes(table), file_name=filename,
            mime="text/csv", key=f"download_{filename}", on_click="ignore",
        )

    if result.best.empty:
        return
    specs = report.figure_specs(targets)
    st.subheader("Figures (300 dpi PNG or SVG)")
    labels = {
        "heatmap": lambda spec: "Heatmap of best scores",
        "ranked": lambda spec: f"Ranked scores, {spec.target}",
        "strip": lambda spec: f"All poses, {spec.target}",
    }
    options = {labels[spec.kind](spec): spec for spec in specs}
    spec = options[st.selectbox("Figure", list(options), key="export_figure")]
    files = spec_files(spec, result.best, results.poses, control)
    st.download_button(
        "Download PNG", files["png"], file_name=f"{spec.stem}.png", mime="image/png",
        key="download_png", on_click="ignore",
    )
    st.download_button(
        "Download SVG", files["svg"], file_name=f"{spec.stem}.svg", mime="image/svg+xml",
        key="download_svg", on_click="ignore",
    )

    # Drawing every figure takes a few seconds, so the ZIPs are made on request, and
    # dropped as soon as anything they depend on changes.
    st.subheader("Everything at once")
    signature = (
        content_key(data["files"]),
        tuple(sorted((k, tuple(sorted(v.items()))) for k, v in analysed["names"].items())),
        analysed["at"].isoformat(),
        control,
        None if redocking is None else (redocking.target, redocking.ligand, redocking.reference_sha256),
    )
    prepared = st.session_state.get("prepared")
    if prepared is not None and prepared["signature"] != signature:
        st.session_state.pop("prepared")
        prepared = None
    if prepared is None and st.button("Prepare figures ZIP and full results bundle", key="prepare_bundle"):
        figures, captions = {}, {}
        with st.spinner(f"Drawing {len(specs)} figures for the ZIP…"):
            for item in specs:
                figures[item.stem] = spec_files(item, result.best, results.poses, control)
                captions[item.stem] = spec_caption(item)
            figures_zip = report.figures_zip(figures, captions)
        with st.spinner("Building the results bundle (Excel, tables, figures, methods, metadata)…"):
            bundle = report.bundle_zip(
                tables=tables, workbook=workbook, methods=methods.text,
                metadata=metadata, figures=figures, captions=captions,
            )
        prepared = {"signature": signature, "figures_zip": figures_zip, "bundle": bundle}
        st.session_state["prepared"] = prepared
    if prepared is not None:
        st.download_button(
            "Download figures.zip", prepared["figures_zip"], file_name="figures.zip",
            mime="application/zip", key="download_figures_zip", on_click="ignore",
        )
        st.download_button(
            "Download vina_results.zip (everything)", prepared["bundle"], file_name="vina_results.zip",
            mime="application/zip", key="download_bundle", on_click="ignore", type="primary",
        )


def spec_files(spec, best: pd.DataFrame, poses: pd.DataFrame, control) -> dict:
    table = poses if spec.kind == "strip" else best
    return figure_files(spec.kind, table, spec.target, control)


def spec_caption(spec) -> str:
    if spec.kind == "heatmap":
        return plots.heatmap_caption()
    if spec.kind == "ranked":
        return plots.ranked_bars_caption(spec.target)
    return plots.pose_strip_caption(spec.target)

main()
