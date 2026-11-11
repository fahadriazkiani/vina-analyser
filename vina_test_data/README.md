# Vina test dataset (for building and testing a docking-results analyser)

**Purpose:** realistic AutoDock Vina output files for developing and testing software.
**This is NOT a scientific study.** Several runs are deliberate cross-dockings (a ligand
docked into a target it was never crystallised with). Do not interpret these scores biologically.

## Provenance
- Docking engine: AutoDock Vina **v1.2.5** (official static Linux binary, GitHub release)
- Input structures: prepared PDBQT files from the official AutoDock Vina repository
  examples: https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example
- Generated: 2026-11-02

## Parameters (identical for every run unless noted)
| Parameter | Value |
|---|---|
| exhaustiveness | 8 |
| num_modes | 9 |
| seed | 42 (repeat runs: seeds 1, 2, 3) |
| cpu | 4 |
| scoring function | vina (default) |
| grid spacing | 0.375 Å (default) |

Grid boxes are in `configs/` (taken unchanged from the official examples):
- `1iep_box.txt`: center (15.190, 53.903, 16.917), size 20×20×20 Å
- `5x72_box.txt`: center (-15.000, 15.000, 129.000), size 30×24×24 Å

## Contents
```
receptors/  1iep_receptor.pdbqt (ABL kinase), 5x72_receptor.pdbqt,
            1iep_receptorH.pdb, 1iep_ligand.sdf (crystal imatinib, for RMSD reference)
ligands/    5 input ligands (PDBQT). imatinib_1iep.pdbqt keeps crystal coordinates.
configs/    Vina box configs
results/    13 Vina output PDBQT files (3-9 poses each, see quirk 1)
logs/       13 Vina stdout logs (contain the affinity / RMSD table)
```
File naming: `<target>__<ligand>[_seedN]_out.pdbqt` and `<target>__<ligand>[_seedN].log`

Ligands: `imatinib_1iep`, `p59_5x72`, `p69_5x72`, `lig_1s63` (from a zinc
metalloprotein example; docked here WITHOUT zinc parameters), `bace1_macrocycle`
(contains macrocycle atom types: real carbons CG0 and pseudo-atoms G0, a useful
parser edge case).

## Known values (use these as test expectations)
Best (mode 1) affinity, kcal/mol:

| target \ ligand | imatinib_1iep | p59_5x72 | p69_5x72 | lig_1s63 | bace1_macrocycle |
|---|---|---|---|---|---|
| 1iep | -13.27 | -8.463 | -7.657 | -8.781 | -6.367 |
| 5x72 | -10.83 | -11 | -10.4 | -9.856 | -7.234 |

Imatinib repeat runs on 1iep: seed1 -13.24, seed2 -13.24, seed3 -13.26.

Redocking RMSD of imatinib (mode 1 vs crystal pose, 37 heavy atoms, same atom order,
no symmetry correction): seed42 0.91 Å, seed1 0.90 Å, seed2 0.34 Å, seed3 0.90 Å.

## Real quirks the parser must handle (verified in these files)
1. **Log and PDBQT pose counts differ.** Logs always list 9 modes, but output PDBQTs only
   contain poses within `energy_range` (default 3 kcal/mol) of the best:
   imatinib on 1iep has 3–4 poses (seeds 42,1 → 3; seeds 2,3 → 4); imatinib on 5x72 has 7;
   all other files have 9.
2. **Precision differs.** Logs round (e.g. -13.27); PDBQT `REMARK VINA RESULT` keeps 3
   decimals (-13.273). Compare with a tolerance of 0.01.
3. **Macrocycle atom types** appear in `*bace1_macrocycle*` files. CG0 atoms are real carbons
   (Meeko's glue carbons, where the ring was opened) and count as heavy atoms. G0 atoms
   (named *1, *2) are pseudo-atoms on top of their partner CG0 and must not be counted.
   bace1_macrocycle has 38 heavy atoms (43 atom lines = 38 heavy + 3 HD + 2 G0).

Note: Vina 1.2.x writes results to stdout, not a `--log` file; logs here are captured stdout.
Each output PDBQT also stores the score in `REMARK VINA RESULT:` lines. These are an
independent source the parser should agree with.

## Licence
The prepared input structures and grid boxes (`receptors/`, `ligands/`, `configs/`) are
taken unchanged from the examples of the official AutoDock Vina repository
(https://github.com/ccsb-scripps/AutoDock-Vina, `example/`), which is distributed under
the Apache License, Version 2.0. A copy of that licence is in `LICENSES/Apache-2.0.txt`
at the repository root; the Vina repository has no NOTICE file. The underlying
structures come from the Protein Data Bank (entries 1IEP and 5X72).

The output files in `results/` and `logs/` were produced for this dataset by running
AutoDock Vina v1.2.5 on those inputs, and are distributed with the rest of this
repository under its MIT licence.
