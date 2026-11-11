# Project progress: Vina Docking Results Analyser

Handoff notes so a new Claude Code session (or a person) can resume work.
Read this file, then `MASTER_PROMPT_vina_analyser.md`, then `vina_test_data/README.md`.

_Last updated: 2026-11-11, after Phase 5._

## Status

| Phase | Status | Notes |
|---|---|---|
| 1. Parsing (`parsers.py`) | ✅ Done | 117 pose rows (93 with coordinates, 24 log-only), 13 runs |
| 2. Analysis (`analysis.py`) | ✅ Done | Ranking, Δ vs control, LE, seed stats, redocking RMSD |
| 3. Figures (`plots.py`) | ✅ Done | Heatmap, ranked bars, per-target strip plots; reviewed visually |
| 4. Streamlit app (`app.py`) | ✅ Done | Upload/zip, recognition table, 5 tabs, RMSD gate; AppTest-tested |
| 5. Exports & reproducibility (`report.py`) | ✅ Done | Excel, ZIPs, methods text, run_metadata.json with git commit |
| 6. Polish (README, LICENSE, deployment, security) | ⏳ Next | Not started |

Test suite: 384 passed, 0 warnings, after Phase 5 (about 80 s; figure tests dominate).

## Working method (keep using it)
- Plan first, wait for approval, then tests first (show them failing), then code.
- After each phase: full pytest output, a deliberate-bug report (bugs planted in a
  scratch copy), a phase note (built / untested / limitations), and commits.
- Look at generated figures, not just the tests. Tests check structure, not appearance.
- Never modify `vina_test_data/` data files (README corrections are documentation only).

## Decisions made (don't undo without a reason)
1. **Heavy atoms:** exclude H, HD, HS and macrocycle pseudo-atoms `G0, G1…`.
   **`CG0, CG1…` are real carbons and count.** bace1_macrocycle = 38 heavy atoms.
2. **Poses table:** 117 rows, explicit boolean `has_coordinates`
   (False = mode only in the log, outside energy_range).
3. **PDBQT values preferred** over rounded log values; disagreement > 0.01 → warning.
4. **Seed** comes from the log; never assumed (e.g. never assume 42).
5. **No cross-target comparison** anywhere: ranking, Δ vs control and plots are per target.
6. **Control comparison:** same rule for every ligand (best across its runs);
   show `n_runs` for ligand and control; warn when they differ.
7. **Δ direction:** negative = more favourable score than control.
8. **Redocking RMSD:** `is_experimental_pose` is a required argument with no default;
   PDBQT references only (SDF refused: needs atom mapping); in-order matching with
   AutoDock type → element mapping; no superposition, no symmetry correction;
   "within 2.0 Å" column labelled as the common criterion, not proof.
9. **Scoring-error note on outputs:** Vina's error is roughly 2–3 kcal/mol; smaller
   differences must not be read as differences in binding.
10. **Figures:** one shared colour scale in kcal/mol; "compare within a row" caption;
    box only when ≥ 5 poses with coordinates; Okabe–Ito palette; 8.5 cm width;
    300 dpi PNG + SVG; deterministic output; Unicode minus; heatmap cells marked
    "(n=k)" when best of k > 1 runs.
11. **Python:** code must stay compatible with Python 3.11+ (developed on 3.14).
12. **Git:** commits authored by Fahad Riaz Kiani; AI co-author trailer kept for transparency.
13. **Recognition table:** explicit names are passed for EVERY included row; what the
    table shows is exactly what is analysed. Results are hidden if the table changes
    after Run until Run is pressed again. Seed column is read-only (log seed wins).
14. **RMSD in the app:** user must choose target AND ligand; confirmation checkbox is
    unticked by default and keyed to the reference content (new file → unticked again).
    `redocking_rmsd(..., target=...)` filters poses and grid-box checks to that target.
15. **CSV exports:** `utf-8-sig` (BOM) so Excel on Windows shows Δ, Å and −.
16. **Figures:** layout is settled (constrained layout run to a fixed point at 300 dpi)
    and then frozen at build time, so every save of a figure is byte-identical; files
    are written as bytes (no CRLF in SVGs on Windows). Empty tables are refused.
17. **Upload limits:** 50 MB per upload (`.streamlit/config.toml`, kept equal to
    `inputs.MAX_UPLOAD_MB` by a test); zips: 2000 files, 20 MB per file, 200 MB
    unpacked, no `..`/absolute paths, nested zips skipped. Parse cache max 8 entries,
    figure cache max 32.
18. **Methods paragraph:** only parsed log values; "not recorded" otherwise, never
    Vina defaults. Logs do not record `num_modes` or `energy_range`: the observed mode
    count is reported, energy_range is "not recorded". Seeds are quoted only when the
    run's `seed_source` is "log" (a filename seed is never presented as logged). A value
    that differs between runs is listed in full with a warning, never chosen. One
    paragraph; grid boxes as bullets when there is more than one target. File
    basenames only (logs hold full paths). Redocking sentence only for an RMSD computed
    against a confirmed experimental reference; it reports the top-ranked pose of
    each run (range across runs).
19. **Formula injection:** text starting with = + - @ tab CR is stored as text in
    Excel and prefixed with ' in CSV. Numbers are untouched.
20. **Run metadata:** SHA-256 of the original uploaded bytes (excluded files listed
    too), the names/seed each file was analysed as, package versions, analysis time
    (stored when Run is pressed), and git commit (full + short + dirty, untracked files
    count as dirty) when available; the `git` key is omitted otherwise.
21. **ZIPs:** fixed 1980 timestamps, sorted entries, fixed `create_system`, so identical
    input gives identical bytes on every OS. Entry names are validated; target names
    go through `safe_name` before becoming file names.

## Known issues to address later
- Heatmap "(n=4)" label is tight in its cell → put the mark on a second line (Phase 6).
- Large grids (many ligands) untested; long names at 45° use a lot of space.
- Python 3.11 not actually run yet (only grammar + dependency checks) → add CI in Phase 6.
- Licence of the Vina example files to verify in Phase 6 (believed Apache 2.0).
- App figures/tables: wide grids and narrow windows not checked in a browser yet.
- Real file_uploader and data_editor cell editing are NOT covered by AppTest; they are
  on the manual checklist below.
- Figure builds are slower (~0.5-1.5 s each) because the layout is settled at 300 dpi;
  cached in the app.
- results.xlsx "modified" property is set by openpyxl to the save time; only "created"
  equals the analysis time. Not fixable without rewriting the file after saving.
- Spinners on the ZIP build are not visible to AppTest: manual check only.
- The "seed used" shown for a PDBQT in the recognition table comes from the filename
  pairing; retyping a recognised PDBQT's target/ligand does not update it (Summary is
  still correct).
- Score column on screen has no fixed format (Excel does: 0.000).
- README: add a "How this was built" section (AI-assisted, decisions and review by me).
- Repo is private; make public after Phase 6.

## Phase 4 manual checks (not covered by automated tests)
Run with `streamlit run app.py`, then:
1. "Load example dataset" → look at every tab, including the figures in the browser.
2. Upload the 26 files as a multi-select, then as one .zip → identical results.
3. Rename a file (e.g. `run7.pdbqt`) → "needs names", Run disabled; type names → parses.
4. Edit a cell after Run → results hide until Run is pressed again.
5. RMSD: an SDF is refused with the reason; with the PDBQT, nothing is computed until
   target, ligand and the checkbox are all set; 1iep gives 0.91 / 0.34 Å.
6. Pick a control, then go back to "No control".
7. Download a CSV, open it in Excel: Δ, Å and − display correctly.
8. Narrow the window: tables and figures stay usable.

## Phase 5 manual checks
1. Export tab → open results.xlsx in Excel: every sheet, Δ/Å headers, decimals, Read me.
2. "Prepare figures ZIP and full results bundle": the spinner labels appear; unzip both.
3. Read methods.txt as a reviewer: is anything claimed that the logs don't support?
   Compute an RMSD on the RMSD tab and check the redocking sentence appears.
4. `Get-FileHash vina_test_data\logs\1iep__p59_5x72.log` equals the SHA-256 in
   run_metadata.json (lower-case it first).

## Next step: Phase 6 (polish)
Ask for the plan first. README (purpose, screenshot, install/run, input format,
limitations, how to cite Vina), LICENSE (MIT), deployment notes for Streamlit Community
Cloud, input-safety review, CI on Python 3.11, and the known issues above.

## How to resume in a new session
Paste this into Claude Code:

> Read PROGRESS.md, MASTER_PROMPT_vina_analyser.md and vina_test_data/README.md.
> Run `git status`, `git log --oneline -10` and `pytest -q`, and report the results.
> If there are uncommitted changes, show me `git diff --stat` and don't discard
> anything. Then give me the Phase 6 plan and wait for approval.
