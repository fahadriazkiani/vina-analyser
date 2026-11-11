"""Turn parsed poses into ranked tables, comparisons and redocking RMSD.

Three rules shape this module:

* Scores are only ever compared inside one target, with one grid box and one scoring
  function. Nothing is subtracted across targets.
* Nothing is estimated. Where an input is missing, the output cell stays blank and an
  issue says why.
* An RMSD is only computed against a reference the user has explicitly confirmed to be
  an experimental pose.
"""

import math
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import pandas as pd

from vina_analyser.parsers import (
    GLUE_CARBON,
    MACROCYCLE_PSEUDO_ATOM,
    Atom,
    Issue,
    ParsedResults,
    heavy_atoms_of,
)

# A pose within this distance of the experimental pose is the usual definition of a
# successful redocking. It is a convention, not a measure of correctness.
RMSD_SUCCESS_CUTOFF = 2.0

# Affinities are compared at the precision the PDBQT itself stores, so that a tie never
# depends on floating-point noise in the last digits.
SCORE_DECIMALS = 3

# A best pose this far from the reference almost always means the wrong reference file
# or the wrong box. Lower-ranked poses are expected to sit far away, since they are
# alternative binding modes, so the check applies to the best pose only.
IMPLAUSIBLE_RMSD = 10.0

LABELS = {
    "affinity": "Vina score (kcal/mol)",
    "rank": "Rank within target",
    "n_runs": "Runs (seeds)",
    "ligand_efficiency": "LE (kcal/mol per heavy atom)",
    "delta_vs_control": "Δ score vs control (kcal/mol)",
    "rmsd_vs_reference": "RMSD vs user-supplied reference (Å)",
    "within_2A": f"Within {RMSD_SUCCESS_CUTOFF} Å (common redocking success criterion)",
    "mean_affinity": "Mean best score (kcal/mol)",
    "sd_affinity": "SD of best score (kcal/mol)",
    "range_affinity": "Range of best score (kcal/mol)",
}

NOTES = {
    # Fixed text, shown with every analysis output.
    "scoring_error": (
        "Vina's scoring error is roughly 2-3 kcal/mol; score differences smaller than "
        "this should not be interpreted as differences in binding."
    ),
    "delta_direction": (
        "In the Δ column, negative = more favourable score than control, and positive "
        "= less favourable."
    ),
    "affinity": (
        "Vina scores are estimates from an empirical scoring function. They are not "
        "measured binding affinities and say nothing about activity in a cell or an "
        "organism."
    ),
    "ligand_efficiency": (
        "Ligand efficiency is a size-normalised heuristic: the score divided by the "
        "number of heavy atoms. It inherits every limitation of the score and is blank "
        "wherever the heavy atom count could not be determined."
    ),
    "delta_vs_control": (
        "A difference between two docking scores for the same target, grid box and "
        "scoring function. It is not a difference in binding affinity, potency or "
        "selectivity, and it is never computed across different targets or boxes."
    ),
    "seed_summary": (
        "Mean, SD and range across repeat runs describe the reproducibility of the "
        "search for one fixed setup. They are not experimental error, and SD is left "
        "blank below three runs."
    ),
    "rmsd_vs_reference": (
        "Heavy-atom RMSD against a reference the user supplied and confirmed to be an "
        "experimental (crystal or cryo-EM) pose. Against a generated 3D conformer the "
        "number is meaningless, because the conformer carries no experimental "
        "coordinates. Atoms are compared in file order with no superposition, since "
        "docked poses are already in the receptor frame, and no symmetry correction is "
        "applied, so symmetry-equivalent atoms can inflate the value."
    ),
    "within_2A": (
        f"A pose within {RMSD_SUCCESS_CUTOFF} Å of the reference is conventionally "
        "called a successful redocking. It is not proof that the pose is correct, and "
        "it depends entirely on the reference being an experimental pose."
    ),
}

BEST_COLUMNS = [
    "target",
    "ligand",
    "n_runs",
    "seed",
    "pose",
    "affinity",
    "heavy_atoms",
    "source_file",
]
RMSD_COLUMNS = {
    "target": "str",
    "ligand": "str",
    "seed": "Int64",
    "pose": "int64",
    "affinity": "float64",
    "has_coordinates": "bool",
    "n_atoms_compared": "Int64",
    "rmsd_vs_reference": "float64",
    "within_2A": "boolean",
    "source_file": "str",
}
SETUP_COLUMNS = [
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "receptor",
    "scoring_function",
]

# AutoDock atom types are not element symbols: "A" is aromatic carbon, "NA" is a
# hydrogen-bond accepting nitrogen (not sodium), and CG* are Meeko's glue carbons.
_ELEMENT_BY_TYPE = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "SI": "Si", "B": "B", "SE": "Se", "AT": "At",
    "MG": "Mg", "CA": "Ca", "MN": "Mn", "FE": "Fe", "ZN": "Zn",
    "H": "H", "HD": "H", "HS": "H",
}

_SDF_MARKER = re.compile(r"^(?:.*V[23]000\s*$|\$\$\$\$\s*$|\s*M\s+END\s*$)", re.MULTILINE)


class ReferenceRefused(ValueError):
    """The supplied reference cannot be used, with the reason in the message."""


@dataclass(frozen=True)
class Reference:
    """A reference ligand pose the user has confirmed to be experimental."""

    atoms: list[Atom]
    filename: str
    is_experimental_pose: bool

    @property
    def heavy_atoms(self) -> list[Atom]:
        return heavy_atoms_of(self.atoms)

    @property
    def n_heavy_atoms(self) -> int:
        return len(self.heavy_atoms)


@dataclass
class ControlComparison:
    table: pd.DataFrame
    issues: list[Issue] = field(default_factory=list)


@dataclass
class RmsdResult:
    table: pd.DataFrame
    issues: list[Issue] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


@dataclass
class AnalysisResult:
    best: pd.DataFrame  # one row per target x ligand, ranked within the target
    seeds: pd.DataFrame  # one row per target x ligand, across repeat runs
    issues: list[Issue] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def output_notes(*keys: str) -> list[str]:
    """The limitations that travel with an output table, scoring error always first."""
    return [NOTES["scoring_error"], NOTES["affinity"], *(NOTES[key] for key in keys)]


def element_of(atom_type: str) -> str | None:
    """Map an AutoDock atom type to a chemical element, or None if it isn't one."""
    if not atom_type:
        return None
    if GLUE_CARBON.fullmatch(atom_type):
        return "C"  # CG0, CG1, ... are real carbons where Meeko opened the ring
    if MACROCYCLE_PSEUDO_ATOM.fullmatch(atom_type):
        return None  # G0, G1, ... are dummy atoms, not elements
    return _ELEMENT_BY_TYPE.get(atom_type.upper())


# ----------------------------------------------------------------- best poses


def best_poses(poses: pd.DataFrame) -> pd.DataFrame:
    """The best pose of each run: lowest score among poses that have coordinates.

    Log-only modes are excluded, however good their score, because a pose without
    coordinates cannot be inspected, exported or used for RMSD.
    """
    with_coordinates = poses[poses["has_coordinates"]]
    ordered = with_coordinates.sort_values(
        ["target", "ligand", "seed", "affinity", "pose"], na_position="last"
    )
    best = ordered.groupby(["target", "ligand", "seed"], dropna=False, sort=True).head(1)
    return best.sort_values(["target", "ligand", "seed"], ignore_index=True)


def best_per_ligand(poses: pd.DataFrame) -> pd.DataFrame:
    """One row per target and ligand: the best run, plus how many runs it came from.

    The same rule applies to every ligand, control included, so a ligand with four
    seeds is not given an advantage over one with a single run. `n_runs` makes the
    difference visible instead of hiding it.
    """
    best = best_poses(poses)
    counts = (
        best.groupby(["target", "ligand"], dropna=False, sort=True)
        .size()
        .rename("n_runs")
        .reset_index()
    )
    # Ties between runs go to the lowest seed, so the choice is reproducible.
    chosen = (
        best.sort_values(["target", "ligand", "affinity", "seed"], na_position="last")
        .groupby(["target", "ligand"], dropna=False, sort=True)
        .head(1)
    )
    merged = chosen.merge(counts, on=["target", "ligand"], how="left")
    merged["n_runs"] = merged["n_runs"].astype("Int64")
    return merged[BEST_COLUMNS].sort_values(["target", "ligand"], ignore_index=True)


def ligand_efficiency(best: pd.DataFrame) -> pd.DataFrame:
    """Add -score / heavy atoms, left blank wherever either input is missing."""
    table = best.copy()
    affinity = pd.to_numeric(table["affinity"], errors="coerce")
    heavy = pd.to_numeric(table["heavy_atoms"], errors="coerce")
    # Guard against a zero count as well as a missing one: never divide by zero.
    heavy = heavy.where(heavy > 0)
    table["ligand_efficiency"] = (-affinity / heavy).astype("float64")
    return table


def rank_per_target(best: pd.DataFrame) -> pd.DataFrame:
    """Rank ligands inside each target, with tied scores sharing a rank.

    Ties take the lowest rank and the next rank is skipped, so two ligands tied at the
    top are both rank 1 and the third is rank 3. Tied rows are ordered by ligand name,
    which keeps the output stable between runs.
    """
    table = best.copy()
    rounded = pd.to_numeric(table["affinity"], errors="coerce").round(SCORE_DECIMALS)
    table["_rounded"] = rounded
    table["rank"] = (
        rounded.groupby(table["target"]).rank(method="min", na_option="keep").astype("Int64")
    )
    table = table.sort_values(
        ["target", "_rounded", "ligand"], na_position="last", ignore_index=True
    )
    return table.drop(columns="_rounded")


# ------------------------------------------------------------ control comparison


def compare_to_control(
    best: pd.DataFrame,
    control: str,
    runs: pd.DataFrame | None = None,
) -> ControlComparison:
    """Add the score difference against a control ligand, within each target only."""
    table = best.copy()
    issues: list[Issue] = []
    for column in ("delta_vs_control", "control_affinity"):
        table[column] = pd.Series(float("nan"), dtype="float64", index=table.index)
    table["control_n_runs"] = pd.Series(pd.NA, dtype="Int64", index=table.index)

    control_rows = table[table["ligand"] == control]
    if control_rows.empty:
        issues.append(
            Issue(
                "error",
                "control_not_found",
                f"Control ligand {control!r} is not in the loaded results, so no "
                "comparison was made.",
            )
        )
        return ControlComparison(table, issues)

    setups = _setup_keys(runs)
    for target, group in table.groupby("target", sort=True):
        control_row = group[group["ligand"] == control]
        if control_row.empty:
            issues.append(
                Issue(
                    "warning",
                    "control_missing_for_target",
                    f"Control ligand {control!r} was not docked against target "
                    f"{target!r}, so that target has no comparison.",
                )
            )
            continue

        control_row = control_row.iloc[0]
        comparable = group.index
        if setups is not None:
            control_setup = setups.get((target, control))
            mismatched = [
                index
                for index in group.index
                if setups.get((target, table.at[index, "ligand"])) != control_setup
            ]
            if mismatched:
                names = sorted({table.at[index, "ligand"] for index in mismatched})
                issues.append(
                    Issue(
                        "warning",
                        "incomparable_setup",
                        f"On target {target!r}, {names} were run with a different grid "
                        "box, receptor or scoring function than the control, so no "
                        "difference was computed for them.",
                    )
                )
                comparable = group.index.difference(mismatched)

        table.loc[comparable, "control_affinity"] = control_row["affinity"]
        table.loc[comparable, "control_n_runs"] = control_row["n_runs"]
        table.loc[comparable, "delta_vs_control"] = (
            pd.to_numeric(table.loc[comparable, "affinity"], errors="coerce")
            - control_row["affinity"]
        )

        unequal = sorted(
            group.loc[
                (group["ligand"] != control) & (group["n_runs"] != control_row["n_runs"]),
                "ligand",
            ]
        )
        if unequal:
            issues.append(
                Issue(
                    "warning",
                    "unequal_run_counts",
                    f"On target {target!r}, the control has {control_row['n_runs']} run(s) "
                    f"but {unequal} do not, so the comparison rests on different amounts "
                    "of sampling.",
                )
            )
    return ControlComparison(table, issues)


def _setup_keys(runs: pd.DataFrame | None) -> dict[tuple, frozenset] | None:
    """Map each target and ligand to the set of search setups it was run with."""
    if runs is None or runs.empty:
        return None
    if any(column not in runs.columns for column in SETUP_COLUMNS):
        return None
    keys: dict[tuple, set] = {}
    for row in runs.itertuples(index=False):
        setup = tuple(getattr(row, column) for column in SETUP_COLUMNS)
        keys.setdefault((row.target, row.ligand), set()).add(setup)
    return {pair: frozenset(setups) for pair, setups in keys.items()}


# --------------------------------------------------------- seed reproducibility


def seed_summary(poses: pd.DataFrame) -> pd.DataFrame:
    """Summarise the best score across repeat runs of the same target and ligand.

    `n_runs` is always reported next to the statistics, and SD is left blank below
    three runs, where it would be misleading rather than informative.
    """
    best = best_poses(poses)
    rows = []
    for (target, ligand), group in best.groupby(["target", "ligand"], dropna=False, sort=True):
        scores = pd.to_numeric(group["affinity"], errors="coerce").dropna()
        n_runs = len(group)
        rows.append(
            {
                "target": target,
                "ligand": ligand,
                "n_runs": n_runs,
                "seeds": sorted(seed for seed in group["seed"] if pd.notna(seed)),
                "mean_affinity": scores.mean() if len(scores) else pd.NA,
                # Sample SD (n - 1), and only once there are three runs to describe.
                "sd_affinity": scores.std(ddof=1) if len(scores) >= 3 else pd.NA,
                "min_affinity": scores.min() if len(scores) else pd.NA,
                "max_affinity": scores.max() if len(scores) else pd.NA,
                "range_affinity": (scores.max() - scores.min()) if len(scores) else pd.NA,
            }
        )
    columns = [
        "target",
        "ligand",
        "n_runs",
        "seeds",
        "mean_affinity",
        "sd_affinity",
        "min_affinity",
        "max_affinity",
        "range_affinity",
    ]
    frame = pd.DataFrame(rows, columns=columns)
    for column in ("mean_affinity", "sd_affinity", "min_affinity", "max_affinity", "range_affinity"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["n_runs"] = frame["n_runs"].astype("Int64")
    return frame


# ------------------------------------------------------------- redocking RMSD


def parse_reference(text: str, *, is_experimental_pose: bool, filename: str = "") -> Reference:
    """Read a reference ligand PDBQT.

    `is_experimental_pose` has no default on purpose. An RMSD is only interpretable
    against an experimental (crystal or cryo-EM) pose, and no filename can establish
    that, so the caller has to state it.
    """
    if not is_experimental_pose:
        raise ReferenceRefused(
            "A redocking RMSD needs a reference that is an experimental (crystal or "
            "cryo-EM) pose. Confirm that this file holds one before continuing: "
            "against a generated 3D conformer the RMSD is meaningless."
        )
    if filename.lower().endswith((".sdf", ".mol", ".mol2")) or _SDF_MARKER.search(text):
        raise ReferenceRefused(
            "SDF references need atom mapping; not supported yet. Supply the reference "
            "as a PDBQT whose heavy atoms are in the same order as the docked ligand."
        )

    models = sum(1 for line in text.splitlines() if line.startswith("MODEL"))
    if models > 1:
        raise ReferenceRefused(
            f"Reference holds {models} models. Provide a single pose, so it is clear "
            "which one the RMSD is measured against."
        )

    atoms = _read_first_model(text)
    if not heavy_atoms_of(atoms):
        raise ReferenceRefused(
            f"Reference {filename or 'file'} contains no atoms that could be read."
        )
    return Reference(atoms=atoms, filename=filename, is_experimental_pose=True)


def _read_first_model(text: str) -> list[Atom]:
    from vina_analyser.parsers import _read_atom  # same fixed-width rules as the results

    class _Sink:
        issues: list = []

    pose = type("_Pose", (), {"atoms": []})()
    pose.atoms = []
    for line in text.lstrip("﻿").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            _read_atom(line, pose, _Sink())
    return pose.atoms


def heavy_atom_rmsd(pose_atoms: Sequence[Atom], reference_atoms: Sequence[Atom]) -> float:
    """Heavy-atom RMSD in file order, with no superposition and no symmetry handling."""
    pose_heavy = heavy_atoms_of(pose_atoms)
    reference_heavy = heavy_atoms_of(reference_atoms)
    if len(pose_heavy) != len(reference_heavy) or not pose_heavy:
        raise ValueError(
            f"Cannot compare {len(pose_heavy)} heavy atoms with {len(reference_heavy)}."
        )
    total = sum(
        (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
        for a, b in zip(pose_heavy, reference_heavy)
    )
    return math.sqrt(total / len(pose_heavy))


def redocking_rmsd(
    results: ParsedResults,
    reference: Reference,
    ligand: str | None = None,
    runs: pd.DataFrame | None = None,
    target: str | None = None,
) -> RmsdResult:
    """Per-pose RMSD against the reference, for poses that have coordinates.

    A reference describes one ligand in one binding site, so callers should name both
    the target and the ligand. The same imatinib reference gives 0.91 Å on 1iep and over
    100 Å on 5x72: only the first number means anything.

    Runs whose ligand cannot match the reference are reported and left blank. If no
    run matches at all and no particular ligand was asked for, the reference is wrong
    for this dataset, so it is refused outright.
    """
    if not reference.is_experimental_pose:  # defensive: Reference cannot be built otherwise
        raise ReferenceRefused("Reference was not confirmed to be an experimental pose.")

    poses = results.poses[results.poses["has_coordinates"]]
    if target is not None:
        if target not in set(results.poses["target"]):
            raise ReferenceRefused(f"Target {target!r} is not in this dataset.")
        poses = poses[poses["target"] == target]
        if runs is not None:
            runs = runs[runs["target"] == target]
    if ligand is not None:
        poses = poses[poses["ligand"] == ligand]

    issues: list[Issue] = []
    rows = []
    matched_any = False
    reference_heavy = reference.heavy_atoms

    for (target, run_ligand, seed), group in poses.groupby(
        ["target", "ligand", "seed"], dropna=False, sort=True
    ):
        sample = results.coordinates.get((target, run_ligand, seed, int(group["pose"].iloc[0])), [])
        problem = _reference_problem(sample, reference_heavy, target, run_ligand)
        if problem is not None:
            issues.append(problem)
        else:
            matched_any = True
            # Names are not used for matching, so a difference is a warning, not a
            # refusal: the elements already line up in order.
            name_warning = _name_warning(sample, reference_heavy, target, run_ligand)
            if name_warning is not None:
                issues.append(name_warning)

        for row in group.itertuples(index=False):
            atoms = results.coordinates.get((target, run_ligand, seed, row.pose), [])
            value = None
            if problem is None and atoms:
                value = heavy_atom_rmsd(atoms, reference_heavy)
                if value > IMPLAUSIBLE_RMSD and row.pose == 1:
                    issues.append(
                        Issue(
                            "warning",
                            "implausible_rmsd",
                            f"{target} / {run_ligand} / seed {seed}: the best pose is "
                            f"{value:.1f} Å from the reference. A value this large "
                            "usually means the reference belongs to a different receptor "
                            "or grid box, rather than that the docking failed.",
                            row.source_file,
                        )
                    )
            rows.append(
                {
                    "target": target,
                    "ligand": run_ligand,
                    "seed": seed,
                    "pose": row.pose,
                    "affinity": row.affinity,
                    "has_coordinates": True,
                    "n_atoms_compared": len(reference_heavy) if value is not None else None,
                    "rmsd_vs_reference": value,
                    "within_2A": (value <= RMSD_SUCCESS_CUTOFF) if value is not None else None,
                    "source_file": row.source_file,
                }
            )

    if not matched_any and ligand is None:
        raise ReferenceRefused(
            f"The reference has {reference.n_heavy_atoms} heavy atoms and no run in this "
            "dataset matches it, so there is nothing it can be compared with."
        )

    issues += _box_issues(reference_heavy, runs)
    frame = pd.DataFrame(rows, columns=list(RMSD_COLUMNS), dtype="object").astype(
        dict(RMSD_COLUMNS)
    )
    return RmsdResult(
        frame.sort_values(["target", "ligand", "seed", "pose"], ignore_index=True),
        issues,
        output_notes("rmsd_vs_reference", "within_2A"),
    )


def _reference_problem(
    atoms: Sequence[Atom], reference_heavy: Sequence[Atom], target: str, ligand: str
) -> Issue | None:
    """Decide whether this run's ligand can be compared with the reference at all."""
    pose_heavy = heavy_atoms_of(atoms)
    if not pose_heavy:
        return Issue(
            "warning", "no_atoms", f"{target} / {ligand}: no atoms to compare.", ""
        )
    if len(pose_heavy) != len(reference_heavy):
        # A different number of heavy atoms means a different molecule, not a bug.
        return Issue(
            "info",
            "reference_not_applicable",
            f"{target} / {ligand} has {len(pose_heavy)} heavy atoms and the reference has "
            f"{len(reference_heavy)}, so the reference does not describe this ligand.",
            "",
        )

    for position, (pose_atom, reference_atom) in enumerate(zip(pose_heavy, reference_heavy), 1):
        pose_element = element_of(pose_atom.atom_type)
        reference_element = element_of(reference_atom.atom_type)
        if pose_element != reference_element:
            return Issue(
                "error",
                "element_mismatch",
                f"{target} / {ligand}: atom {position} is {pose_element} in the pose but "
                f"{reference_element} in the reference, so the atoms are not the same "
                f"molecule in the same order. Pose composition: {_composition(pose_heavy)}.",
                "",
            )
    return None


def _name_warning(
    atoms: Sequence[Atom], reference_heavy: Sequence[Atom], target: str, ligand: str
) -> Issue | None:
    """Warn when the atom names differ, although the elements match in order."""
    for position, (pose_atom, reference_atom) in enumerate(
        zip(heavy_atoms_of(atoms), reference_heavy), 1
    ):
        if pose_atom.name != reference_atom.name:
            return Issue(
                "warning",
                "atom_name_mismatch",
                f"{target} / {ligand}: atom {position} is named {pose_atom.name!r} in the "
                f"pose but {reference_atom.name!r} in the reference. The elements match in "
                "order, so the RMSD was still computed, but check that the two files "
                "describe the same molecule.",
                "",
            )
    return None


def _composition(atoms: Iterable[Atom]) -> str:
    counts: dict[str, int] = {}
    for atom in atoms:
        element = element_of(atom.atom_type) or atom.atom_type
        counts[element] = counts.get(element, 0) + 1
    return " ".join(f"{element}{count}" for element, count in sorted(counts.items()))


def _box_issues(reference_heavy: Sequence[Atom], runs: pd.DataFrame | None) -> list[Issue]:
    """Warn when the reference sits outside a target's grid box: usually a wrong file."""
    if runs is None or runs.empty or not reference_heavy:
        return []
    if any(column not in runs.columns for column in SETUP_COLUMNS):
        return []
    centre = [
        sum(getattr(atom, axis) for atom in reference_heavy) / len(reference_heavy)
        for axis in ("x", "y", "z")
    ]
    issues = []
    for target, group in runs.groupby("target", sort=True):
        row = group.iloc[0]
        box = [row["center_x"], row["center_y"], row["center_z"]]
        size = [row["size_x"], row["size_y"], row["size_z"]]
        if any(pd.isna(value) for value in box + size):
            continue
        if any(abs(centre[axis] - box[axis]) > size[axis] / 2 for axis in range(3)):
            issues.append(
                Issue(
                    "warning",
                    "reference_outside_box",
                    f"The reference centroid lies outside the grid box of target "
                    f"{target!r}, so an RMSD against poses in that box compares two "
                    "different sites.",
                )
            )
    return issues


# ----------------------------------------------------------------- everything


def analyse(results: ParsedResults, control: str | None = None) -> AnalysisResult:
    """Build the ranked best-pose table and the seed summary in one step."""
    best = rank_per_target(ligand_efficiency(best_per_ligand(results.poses)))
    issues = list(results.issues)
    keys = ["ligand_efficiency", "seed_summary"]
    if control is not None:
        comparison = compare_to_control(best, control, runs=results.runs)
        best, extra = comparison.table, comparison.issues
        issues += extra
        keys += ["delta_vs_control", "delta_direction"]
    return AnalysisResult(
        best=best,
        seeds=seed_summary(results.poses),
        issues=issues,
        notes=output_notes(*keys),
    )
