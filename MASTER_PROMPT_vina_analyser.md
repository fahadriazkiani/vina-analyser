# Master prompt: Vina Docking Results Analyser

> How to use: put this file and the `vina_test_data/` folder in an empty project folder,
> open Claude Code in that folder, and paste everything below the line.
> Work one phase at a time. Read and approve the plan before any code is written.

---

You are helping me build an open-source tool called **Vina Docking Results Analyser**.
I am a Medical Laboratory Sciences graduate with computational biology experience
(molecular docking, ADMET, network pharmacology) and an intermediate-beginner Python
developer. I know Python, pandas basics, Streamlit and Git. Explain non-obvious decisions
briefly as you go, so I learn from the code.

## 1. Goal
A Streamlit web app plus a reusable Python package that takes AutoDock Vina output files
and produces reproducible, publication-ready tables and figures. It will be a portfolio
project for Master's applications, so **correctness, tests and honest limitations matter
more than features.**

## 2. Working rules (follow strictly)
1. **Plan first.** Before writing code, give me a short plan: file structure, the data
   model, and what each phase delivers. Wait for my approval.
2. **Test-driven, against real data.** `vina_test_data/` contains real Vina 1.2.5 outputs
   and a README with known expected values. Write pytest tests using those values
   *before* implementing each parser or calculation.
3. **Never invent results.** Only say tests pass after running `pytest` and showing me
   the actual output. If something fails, explain the real cause.
4. **Never overwrite or delete my files** without asking. Never modify anything inside
   `vina_test_data/`, because it is raw data.
5. **Keep it simple.** No database, no login, no Docker, no frontend framework beyond
   Streamlit, no extra dependencies without a reason. Allowed: pandas, numpy, matplotlib,
   streamlit, openpyxl, pytest. Ask before adding anything else.
6. **Scientific honesty in the UI.** Docking scores are estimates, not binding affinities
   or proof of activity. Never use wording like "validated hit" or "will inhibit".
7. Small Git commits at the end of each phase, with clear messages.

## 3. Project structure (adjust only with a reason)
```
vina-analyser/
├── vina_analyser/          # importable package, no Streamlit code here
│   ├── __init__.py
│   ├── parsers.py          # read Vina logs and output PDBQT files
│   ├── analysis.py         # ranking, summaries, ligand efficiency, RMSD
│   ├── plots.py            # matplotlib figures
│   └── report.py           # exports + auto-generated methods text
├── app.py                  # Streamlit UI only, calls the package
├── tests/                  # pytest, uses ../vina_test_data
├── requirements.txt        # pinned versions
├── README.md
└── LICENSE                 # MIT
```

## 4. Input formats (verified facts, see vina_test_data/README.md)
- **Output PDBQT**: `MODEL n` … `ENDMDL` blocks. Each has
  `REMARK VINA RESULT:  <affinity>  <rmsd_lb>  <rmsd_ub>`. Atom lines are ATOM/HETATM;
  the AutoDock atom type is the last column (cols 78–79). Hydrogens are `H`/`HD`.
  Macrocycle pseudo-atoms `G0`, `G1`, `CG0`, `CG1`… are NOT real atoms.
- **Vina 1.2 stdout log**: a table after the line starting `-----+`, with columns:
  mode, affinity (kcal/mol), rmsd l.b., rmsd u.b. The log header also records the
  receptor, ligand, grid center, grid size, exhaustiveness and random seed. Parse these
  for the methods text.
- **Known quirks (must have tests):**
  - The log lists all 9 modes, but the PDBQT only contains poses within `energy_range`
    (3 kcal/mol default), e.g. 3 poses for `1iep__imatinib_1iep_out.pdbqt`.
    Report both counts; never treat this as an error.
  - Log values are rounded (−13.27); PDBQT keeps 3 decimals (−13.273). When both are
    present, prefer the PDBQT value and check agreement within 0.01.
- **Naming:** default pattern `<target>__<ligand>[_seedN]_out.pdbqt`. If a filename
  doesn't match, let the user type target/ligand names in an editable table.
  Never guess silently.

## 5. Features by phase

### Phase 1: Parsing (package + tests only, no UI)
- Parse PDBQT and logs into one tidy pandas DataFrame, one row per pose:
  `target, ligand, seed, pose, affinity, rmsd_lb, rmsd_ub, source_file, heavy_atoms`.
- Heavy atom count excludes H, HD and macrocycle pseudo-atoms.
- Tests: the 13 best affinities in the README table, pose counts per file,
  log-vs-PDBQT agreement, and pseudo-atom exclusion.

### Phase 2: Analysis
- Best pose per target × ligand, ranked per target.
- **Reference compound comparison**: user picks a control ligand (for example a known
  drug) and the app shows Δ affinity versus the control per target. Label it as a
  score difference only.
- **Ligand efficiency** = −affinity / heavy atoms (kcal/mol per heavy atom), with a
  note that it is a size-normalised heuristic.
- **Seed reproducibility**: when the same pair has several seeds, report mean, SD and
  range of best affinity (test: imatinib on 1iep, seeds 1/2/3/42).
- **Redocking RMSD (optional)**: user uploads the reference (crystal) ligand PDBQT and
  the app computes heavy-atom RMSD of each pose. Match atoms by order and name, and
  refuse with a clear message if they don't match. Do NOT superimpose; docking poses
  are already in the receptor frame. Show a warning that no symmetry correction is
  applied. Tests: 0.91 Å (seed 42) and 0.34 Å (seed 2) within ±0.02.

### Phase 3: Figures (matplotlib, 300 dpi PNG + SVG)
- Heatmap: targets × ligands, best affinity annotated in cells.
- Ranked horizontal bar chart per target, control compound highlighted.
- Strip/box plot of all pose affinities per ligand.
- Colour-blind-safe palette, clear axis labels with units, no chartjunk.

### Phase 4: Streamlit app
- Upload multiple files or one .zip. Show a file-recognition table (editable names)
  before analysis.
- Tabs: Summary table · Figures · Reproducibility · RMSD · Export.
- A "Load example dataset" button that uses `vina_test_data/`.
- A visible limitations box (docking ≠ binding affinity; rigid receptor; scores are not
  comparable across different scoring functions or very different boxes).

### Phase 5: Exports & reproducibility
- CSV and Excel (one sheet per table) and a ZIP of all figures.
- An auto-generated **Methods paragraph** built only from parsed log values:
  Vina version, exhaustiveness, seed(s), box center/size, number of poses. If a value
  wasn't found, write "not recorded". Never fill in defaults.
- A `run_metadata.json` with app version, Python and package versions, input file names
  with SHA-256 hashes, and the timestamp.

### Phase 6: Polish
- README: purpose, screenshot, install/run steps, input format, limitations,
  how to cite Vina (Eberhardt et al. 2021, J. Chem. Inf. Model., DOI
  10.1021/acs.jcim.1c00203; Trott & Olson 2010, J. Comput. Chem., DOI 10.1002/jcc.21334).
- Deployment notes for Streamlit Community Cloud.
- Input safety: limit upload size, process files in memory, never execute uploaded
  content, and guard against zip path traversal.

## 6. Definition of done for each phase
`pytest` passes and you show the output · I can run it with the commands you give me ·
a short note on what was built, what is untested, and known limitations.

Start with the plan for the whole project and the detailed plan for Phase 1.
