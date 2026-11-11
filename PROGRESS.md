# Project progress: Vina Docking Results Analyser

Handoff notes so a new Claude Code session (or a person) can resume work.
Read this file, then `MASTER_PROMPT_vina_analyser.md`, then `vina_test_data/README.md`.

_Last updated: 2026-11-11, after Phase 4._

## Status

| Phase | Status | Notes |
|---|---|---|
| 1. Parsing (`parsers.py`) | ✅ Done | 117 pose rows (93 with coordinates, 24 log-only), 13 runs |
| 2. Analysis (`analysis.py`) | ✅ Done | Ranking, Δ vs control, LE, seed stats, redocking RMSD |
| 3. Figures (`plots.py`) | ✅ Done | Heatmap, ranked bars, per-target strip plots; reviewed visually |
| 4. Streamlit app (`app.py`) | ✅ Done | Upload/zip, recognition table, 5 tabs, RMSD gate; AppTest-tested |
| 5. Exports & reproducibility (`report.py`) | ⏳ Next | CSV (utf-8-sig) already done in Phase 4 |
| 6. Polish (README, LICENSE, deployment, security) | ⏳ | Not started |

Test suite: 318 passed, 0 warnings, after Phase 4 (about 90 s; figure tests dominate).

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

## Next step: Phase 5 (exports & reproducibility)
Ask for the plan first. Excel (one sheet per table), ZIP of all figures, methods
paragraph from parsed log values only ("not recorded" when missing), run_metadata.json.

## How to resume in a new session
Paste this into Claude Code:

> Read PROGRESS.md, MASTER_PROMPT_vina_analyser.md and vina_test_data/README.md.
> Run `git status`, `git log --oneline -10` and `pytest -q`, and report the results.
> If there are uncommitted changes, show me `git diff --stat` and don't discard
> anything. Then give me the Phase 5 plan and wait for approval.
