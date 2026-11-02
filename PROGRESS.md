# Project progress: Vina Docking Results Analyser

Handoff notes so a new Claude Code session (or a person) can resume work.
Read this file, then `MASTER_PROMPT_vina_analyser.md`, then `vina_test_data/README.md`.

_Last updated: 2026-11-03, after Phase 3._

## Status

| Phase | Status | Notes |
|---|---|---|
| 1. Parsing (`parsers.py`) | ✅ Done | 117 pose rows (93 with coordinates, 24 log-only), 13 runs |
| 2. Analysis (`analysis.py`) | ✅ Done | Ranking, Δ vs control, LE, seed stats, redocking RMSD |
| 3. Figures (`plots.py`) | ✅ Done | Heatmap, ranked bars, per-target strip plots; reviewed visually |
| 4. Streamlit app (`app.py`) | ⏳ Next | Not started |
| 5. Exports & reproducibility (`report.py`) | ⏳ | Not started |
| 6. Polish (README, LICENSE, deployment, security) | ⏳ | Not started |

Test suite: 251 passed at commit `817a8d9`.

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

## Known issues to address later
- Heatmap "(n=4)" label is tight in its cell → put the mark on a second line (Phase 6).
- Large grids (many ligands) untested; long names at 45° use a lot of space.
- Python 3.11 not actually run yet (only grammar + dependency checks) → add CI in Phase 6.
- Licence of the Vina example files to verify in Phase 6 (believed Apache 2.0).
- CSV exports must be UTF-8 (Δ, Å characters) → Phase 5.
- README: add a "How this was built" section (AI-assisted, decisions and review by me).
- Repo is private; make public after Phase 6.

## Next step: Phase 4 (Streamlit app)
Ask for the plan first. Requirements from the master prompt, plus:
- Upload multiple files or one .zip; editable file-recognition table before analysis.
- Tabs: Summary · Figures · Reproducibility · RMSD · Export.
- "Load example dataset" button using `vina_test_data/`.
- Visible limitations box.
- RMSD tab: an unticked confirmation checkbox ("this reference is an experimental
  pose") must be ticked before RMSD is computed.
- Control-compound selector per analysis; warnings shown, not hidden.

## How to resume in a new session
Paste this into Claude Code:

> Read PROGRESS.md, MASTER_PROMPT_vina_analyser.md and vina_test_data/README.md.
> Run `git status`, `git log --oneline -10` and `pytest -q`, and report the results.
> If there are uncommitted changes, show me `git diff --stat` and don't discard
> anything. Then give me the Phase 4 plan and wait for approval.
