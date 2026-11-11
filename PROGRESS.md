# Project progress: Vina Docking Results Analyser

Handoff notes so a new Claude Code session (or a person) can resume work.
Read this file, then `MASTER_PROMPT_vina_analyser.md`, then `vina_test_data/README.md`.

_Last updated: 2026-11-12, after Phase 6._

## Status

| Phase | Status | Notes |
|---|---|---|
| 1. Parsing (`parsers.py`) | ✅ Done | 117 pose rows (93 with coordinates, 24 log-only), 13 runs |
| 2. Analysis (`analysis.py`) | ✅ Done | Ranking, Δ vs control, LE, seed stats, redocking RMSD |
| 3. Figures (`plots.py`) | ✅ Done | Heatmap, ranked bars, per-target strip plots; reviewed visually |
| 4. Streamlit app (`app.py`) | ✅ Done | Upload/zip, recognition table, 5 tabs, RMSD gate; AppTest-tested |
| 5. Exports & reproducibility (`report.py`) | ✅ Done | Excel, ZIPs, methods text, run_metadata.json with git commit |
| 6. Polish (README, LICENSE, deployment, security) | ✅ Done | README, MIT + Apache notice, CITATION.cff, CI, runs preview, privacy note |

Test suite: 437 passed, 0 warnings, after Phase 6 (about 90-140 s; figure tests dominate).
Deliberate bugs: 100 in scripts/mutation_sets.py (Phases 2-6); every real one caught,
every guard one INVALID. Phase 6 first missed 3 (app preview after an edit, app
heatmap caption above 16 ligands, `..` in a git ref); each now has a test.
Version 0.1.0. Stays below 1.0.0 until the author has used it on their own data.

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
22. **What is analysed is shown before Run:** the editable table holds facts about
    each file only (its seed in the filename never changes). "Runs that will be
    analysed" below it comes from the same cached `parse(texts, names)` as the
    results, so the seed and PDBQT/log pairing on screen are the ones analysed after
    any retyping or exclusion. (Fixes the stale-seed issue; one parse, one truth.)
23. **Privacy statement next to the uploader**, not only in the README: in memory,
    never on disk, nothing kept after the session except caches, publicly reachable.
    Every `st.cache_data` holding user data has a 1-hour ttl, so that sentence is
    true; a test checks every cache has a ttl.
24. **Git commit without the git command** (Community Cloud): read from `.git/HEAD`
    and a loose or packed ref; `dirty` is then `null` (unknown), and the methods
    credit says only "(commit abc1234)". Refs must be `refs/...` without `..`.
25. **Heatmap legibility:** "(n=k)" on its own line under the score. Figures keep the
    8.5 cm column until the ligands don't fit, then widen (0.85 cm per printed cell,
    0.45 cm per 45° label). Above 16 ligands the scores are not printed and the
    caption points to the best-pose table.
26. **On-screen decimals come from `report.NUMBER_FORMATS`**, the Excel formats, so
    screen and workbook can't drift (score 0.000, LE 0.0000, RMSD 0.00).
27. **Loose uploads get the zip's limits** (2000 files, 200 MB total); Streamlit only
    caps each file.
28. **Tests never assert conditionally and never skip**: tables are found with
    `frame_with()`, which fails with a message; git tests require git.
29. **Version policy:** 0.1.0 is the first public release (tag `v0.1.0`; Zenodo mints a
    DOI per tag). 1.0.0 only after use on real data. CI keeps the Windows job.

## Known issues to address later
- The real file_uploader is NOT covered by AppTest (uploads are injected). Table
  edits ARE: tests set the data_editor's widget state, as the browser does
  (`edit_table` in tests/test_app.py). Dragging files and typing into cells stay
  manual checks.
- Figure builds are slower (~0.5-1.5 s each) because the layout is settled at 300 dpi;
  cached in the app.
- results.xlsx "modified" property is set by openpyxl to the save time; only "created"
  equals the analysis time. Not fixable without rewriting the file after saving.
- Spinners on the ZIP build are not visible to AppTest: manual check only.
- Very large grids: at 60 ligands a heatmap or strip plot is ~30 cm wide, wider than a
  journal page. Readable on screen; for a paper, split the ligand set.
- Text typed into the table (file, target, ligand names) appears in warnings rendered
  as Markdown. Streamlit does not render HTML, and only the person who uploaded the
  files sees it, so this is accepted.
- CI has not run yet: it runs on the first push. Python 3.11 is verified only then.

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

## Deliberate-bug (mutation) runner
`scripts/mutation_runner.py` plants one bug at a time in a throwaway copy of the repo
(in a temp folder, `.git` included), runs the full suite against each, and reports a
verdict per bug. The bugs live beside it in `scripts/mutation_sets.py`, grouped by the
file they change.

- **CAUGHT**: some tests fail, and they are a specific subset.
- **MISSED**: everything still passes, so no test guards that behaviour. Add a test.
- **INVALID**: the result says nothing about the behaviour, because the mutant has a
  syntax error, stops the package importing, breaks test collection, crashes the app
  on its start screen (app.py only), or fails more than 30% of all tests.
- **NOT APPLIED**: the text to replace isn't in the code exactly once any more.
  Re-target the entry to where that behaviour lives now; keep its meaning.

The unmodified copy must pass first, or nothing runs. The suite runs in full for each
mutant (no `-x`), so the share of failing tests is always known. The "guard" phase holds
mutants that must come out INVALID; if one is CAUGHT, the guard itself is broken.

Run from the repository root:

```
python scripts/mutation_runner.py --check                 # seconds: does every entry still apply?
python scripts/mutation_runner.py                         # everything, ~35-40 min on 3 workers
python scripts/mutation_runner.py --phase phase5          # one phase (repeatable)
python scripts/mutation_runner.py --module vina_analyser/report.py
```

It exits non-zero if any real mutant is not CAUGHT or any guard mutant is not INVALID.

**The mutations are a chosen subset per module, by design.** Each one breaks a
behaviour we decided matters (see "Decisions made"), such as comparing within one target,
refusing an SDF reference, or a deterministic ZIP. It is not an exhaustive mutation of
every line, so a clean run means "these decisions are guarded", not "the code is fully
tested". When a phase adds a decision, add a mutation for it. Run `--check` after
refactors.

## Phase 6 manual checks
1. Load the example, untick `1iep__imatinib_1iep_seed2.log`: in "Runs that will be
   analysed" the seed-2 PDBQT now shows seed 2 "from filename" and no log file.
2. Retype the ligand of `1iep__imatinib_1iep_out.pdbqt`: its run loses the log seed
   (—, not recorded) and the log becomes a run of its own.
3. Read the privacy statement under "1. Files" as someone with unpublished data.
4. Summary tab: scores show three decimals (−10.996, not −11).
5. Read README.md top to bottom.

## Next steps (the author's, not code)
1. Push, and check the CI run is green on all four jobs (Actions tab).
2. Make the repository public; connect Zenodo; push the tag `v0.1.0` and make a
   GitHub release from it (Zenodo mints the DOI).
3. Deploy on Streamlit Community Cloud with Python 3.12 (README, "Deploying").
4. Use it on your own data; then decide on 1.0.0.

## How to resume in a new session
Paste this into Claude Code:

> Read PROGRESS.md, MASTER_PROMPT_vina_analyser.md and vina_test_data/README.md.
> Run `git status`, `git log --oneline -10` and `pytest -q`, and report the results.
> If there are uncommitted changes, show me `git diff --stat` and don't discard
> anything. Then tell me which of the "Next steps" are still open.
