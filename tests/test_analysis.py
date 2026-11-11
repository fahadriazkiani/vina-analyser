"""Tests for vina_analyser.analysis.

Real-data expectations come from vina_test_data/README.md (best affinities, repeat
seeds, redocking RMSD) or from the PDBQT files themselves. Made-up frames cover ties,
missing controls and reference mismatches, which the real data does not contain.
"""

from pathlib import Path

import pandas as pd
import pytest

from vina_analyser.analysis import (
    LABELS,
    NOTES,
    RMSD_SUCCESS_CUTOFF,
    ReferenceRefused,
    analyse,
    best_per_ligand,
    best_poses,
    compare_to_control,
    element_of,
    ligand_efficiency,
    parse_reference,
    rank_per_target,
    redocking_rmsd,
    seed_summary,
)
from vina_analyser.parsers import POSE_COLUMNS, load_folder

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "vina_test_data"
CRYSTAL_IMATINIB = DATA / "ligands" / "imatinib_1iep.pdbqt"

# From vina_test_data/README.md: mode 1 vs the crystal pose, no symmetry correction.
EXPECTED_RMSD = {42: 0.91, 1: 0.90, 2: 0.34, 3: 0.90}

ATOM_PREFIX = "ATOM      1  C1  LIG L   1    "


@pytest.fixture(scope="module")
def real():
    return load_folder(DATA / "results", DATA / "logs")


@pytest.fixture(scope="module")
def reference():
    # The README states this file keeps the crystal coordinates. The user, not the
    # filename, is what confirms that, so the flag is passed explicitly here.
    return parse_reference(
        CRYSTAL_IMATINIB.read_text(),
        filename=CRYSTAL_IMATINIB.name,
        is_experimental_pose=True,
    )


def codes(issues):
    return [issue.code for issue in issues]


def poses_frame(*rows):
    """Build a poses frame with the same columns and dtypes as the parser's."""
    default = {
        "target": "T",
        "ligand": "L",
        "seed": 42,
        "pose": 1,
        "affinity": -7.0,
        "rmsd_lb": 0.0,
        "rmsd_ub": 0.0,
        "source_file": "T__L_out.pdbqt",
        "heavy_atoms": 20,
        "has_coordinates": True,
    }
    built = [{**default, **row} for row in rows]
    return pd.DataFrame(built, columns=list(POSE_COLUMNS)).astype(dict(POSE_COLUMNS))


def runs_frame(*rows):
    default = {
        "target": "T",
        "ligand": "L",
        "seed": 42,
        "center_x": 0.0,
        "center_y": 0.0,
        "center_z": 0.0,
        "size_x": 20.0,
        "size_y": 20.0,
        "size_z": 20.0,
        "receptor": "r.pdbqt",
        "scoring_function": "vina",
    }
    return pd.DataFrame([{**default, **row} for row in rows])


def reference_text(elements=("C",) * 37, coords=None):
    lines = []
    for index, element in enumerate(elements):
        x, y, z = (coords[index] if coords else (float(index), 0.0, 0.0))
        lines.append(f"{ATOM_PREFIX}{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 {element}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------- best poses


def test_best_pose_is_the_lowest_score_of_each_run(real):
    best = best_poses(real.poses)
    assert len(best) == 13  # one row per run
    assert (best["pose"] == 1).all()  # Vina sorts modes, so mode 1 is always best
    row = best[(best["target"] == "1iep") & (best["seed"] == 42) & (best["ligand"] == "imatinib_1iep")]
    assert row["affinity"].iloc[0] == -13.273


def test_best_pose_ignores_log_only_rows():
    # A log-only mode with a better score must never win: it has no coordinates.
    frame = poses_frame(
        {"pose": 1, "affinity": -7.0, "has_coordinates": True},
        {"pose": 2, "affinity": -9.0, "has_coordinates": False, "source_file": "T__L.log"},
    )
    best = best_poses(frame)
    assert len(best) == 1
    assert best["affinity"].iloc[0] == -7.0


def test_run_without_any_coordinates_is_dropped_from_best_poses():
    frame = poses_frame({"has_coordinates": False, "source_file": "T__L.log"})
    assert best_poses(frame).empty


# -------------------------------------------------------- best per ligand, n_runs


def test_best_per_ligand_uses_the_same_rule_for_every_ligand(real):
    best = best_per_ligand(real.poses)
    assert len(best) == 10  # 2 targets x 5 ligands
    row = best[(best["target"] == "1iep") & (best["ligand"] == "imatinib_1iep")].iloc[0]
    # Best across the four seeds, and n_runs records how many runs that came from.
    assert row["affinity"] == -13.273
    assert row["n_runs"] == 4
    assert row["seed"] == 42
    others = best[(best["target"] == "1iep") & (best["ligand"] != "imatinib_1iep")]
    assert set(others["n_runs"]) == {1}


def test_best_per_ligand_breaks_exact_ties_by_lowest_seed():
    frame = poses_frame(
        {"seed": 7, "affinity": -8.0},
        {"seed": 2, "affinity": -8.0},
    )
    row = best_per_ligand(frame).iloc[0]
    assert row["seed"] == 2
    assert row["n_runs"] == 2


# ------------------------------------------------------------- ligand efficiency


def test_ligand_efficiency_known_values(real):
    best = ligand_efficiency(best_per_ligand(real.poses))
    by_ligand = best.set_index(["target", "ligand"])["ligand_efficiency"]
    assert by_ligand[("1iep", "imatinib_1iep")] == pytest.approx(13.273 / 37, abs=1e-9)
    # 38 heavy atoms, using the corrected rule: CG* counts, G* does not.
    assert by_ligand[("1iep", "bace1_macrocycle")] == pytest.approx(6.367 / 38, abs=1e-9)


def test_ligand_efficiency_is_blank_when_heavy_atoms_are_blank_never_estimated():
    frame = poses_frame({"heavy_atoms": None, "affinity": -8.0})
    row = ligand_efficiency(best_per_ligand(frame)).iloc[0]
    assert pd.isna(row["ligand_efficiency"])


def test_blank_heavy_atoms_are_not_filled_in_from_other_ligands():
    # With other ligands present, a missing count must still stay missing: no mean,
    # median or nearest-neighbour substitution.
    frame = poses_frame(
        {"ligand": "known_a", "heavy_atoms": 20, "affinity": -8.0},
        {"ligand": "known_b", "heavy_atoms": 30, "affinity": -9.0},
        {"ligand": "unknown", "heavy_atoms": None, "affinity": -8.5},
    )
    table = ligand_efficiency(best_per_ligand(frame)).set_index("ligand")
    assert pd.isna(table.loc["unknown", "ligand_efficiency"])
    assert table.loc["known_a", "ligand_efficiency"] == pytest.approx(8.0 / 20)


def test_zero_heavy_atoms_gives_a_blank_not_an_error():
    frame = poses_frame({"heavy_atoms": 0, "affinity": -8.0})
    assert pd.isna(ligand_efficiency(best_per_ligand(frame)).iloc[0]["ligand_efficiency"])


def test_ligand_efficiency_is_blank_when_affinity_is_blank():
    frame = poses_frame({"affinity": None})
    row = ligand_efficiency(best_per_ligand(frame)).iloc[0]
    assert pd.isna(row["ligand_efficiency"])


def test_ligand_efficiency_label_calls_it_a_heuristic():
    assert "heavy atom" in LABELS["ligand_efficiency"]
    assert "heuristic" in NOTES["ligand_efficiency"]


# -------------------------------------------------------------------- ranking


def test_ranks_are_within_a_target_only(real):
    ranked = rank_per_target(best_per_ligand(real.poses))
    ranks = ranked.set_index(["target", "ligand"])["rank"]
    assert ranks[("1iep", "imatinib_1iep")] == 1
    assert ranks[("1iep", "lig_1s63")] == 2
    assert ranks[("1iep", "p59_5x72")] == 3
    assert ranks[("1iep", "p69_5x72")] == 4
    assert ranks[("1iep", "bace1_macrocycle")] == 5
    # 5x72 is ranked on its own: p59 wins there, and imatinib is only second.
    assert ranks[("5x72", "p59_5x72")] == 1
    assert ranks[("5x72", "imatinib_1iep")] == 2
    assert ranks[("5x72", "bace1_macrocycle")] == 5
    assert sorted(ranked.loc[ranked["target"] == "1iep", "rank"]) == [1, 2, 3, 4, 5]


def test_tied_scores_share_a_rank_and_the_next_rank_is_skipped():
    frame = poses_frame(
        {"ligand": "best", "affinity": -9.0},
        {"ligand": "zeta", "affinity": -8.0},
        {"ligand": "alpha", "affinity": -8.0},
        {"ligand": "worst", "affinity": -7.0},
    )
    ranked = rank_per_target(ligand_efficiency(best_per_ligand(frame)))
    ranks = ranked.set_index("ligand")["rank"]
    assert [ranks["best"], ranks["alpha"], ranks["zeta"], ranks["worst"]] == [1, 2, 2, 4]


def test_tied_rows_are_ordered_by_ligand_name():
    frame = poses_frame(
        {"ligand": "zeta", "affinity": -8.0},
        {"ligand": "alpha", "affinity": -8.0},
    )
    ranked = rank_per_target(best_per_ligand(frame))
    assert list(ranked["ligand"]) == ["alpha", "zeta"]


def test_ties_are_decided_at_pdbqt_precision():
    # 3 decimals is what the PDBQT itself stores, so a tie cannot hinge on float noise.
    frame = poses_frame(
        {"ligand": "a", "affinity": -8.0001},
        {"ligand": "b", "affinity": -8.0},
    )
    assert list(rank_per_target(best_per_ligand(frame))["rank"]) == [1, 1]


def test_scores_differing_in_the_third_decimal_are_not_tied():
    frame = poses_frame(
        {"ligand": "a", "affinity": -8.002},
        {"ligand": "b", "affinity": -8.0},
    )
    assert list(rank_per_target(best_per_ligand(frame))["rank"]) == [1, 2]


# ------------------------------------------------------------ control comparison


def test_delta_vs_control_is_computed_within_each_target(real):
    best = best_per_ligand(real.poses)
    result = compare_to_control(best, control="imatinib_1iep")
    table = result.table.set_index(["target", "ligand"])
    assert table.loc[("1iep", "imatinib_1iep"), "delta_vs_control"] == 0.0
    assert table.loc[("1iep", "lig_1s63"), "delta_vs_control"] == pytest.approx(4.492, abs=1e-9)
    # On 5x72 the control is beaten, so the sign is negative there.
    assert table.loc[("5x72", "p59_5x72"), "delta_vs_control"] == pytest.approx(-0.168, abs=1e-9)
    assert table.loc[("5x72", "imatinib_1iep"), "control_affinity"] == -10.828
    assert table.loc[("1iep", "lig_1s63"), "control_affinity"] == -13.273


def test_delta_is_never_computed_across_targets():
    # The control exists only for target A, so B must stay blank rather than borrow it.
    frame = poses_frame(
        {"target": "A", "ligand": "ctrl", "affinity": -9.0},
        {"target": "A", "ligand": "x", "affinity": -8.0},
        {"target": "B", "ligand": "x", "affinity": -7.0},
    )
    result = compare_to_control(best_per_ligand(frame), control="ctrl")
    table = result.table.set_index(["target", "ligand"])
    assert table.loc[("A", "x"), "delta_vs_control"] == pytest.approx(1.0)
    assert pd.isna(table.loc[("B", "x"), "delta_vs_control"])
    assert pd.isna(table.loc[("B", "x"), "control_affinity"])
    assert "control_missing_for_target" in codes(result.issues)
    assert "B" in next(i.message for i in result.issues if i.code == "control_missing_for_target")


def test_control_absent_everywhere_is_refused():
    frame = poses_frame({"ligand": "x"})
    result = compare_to_control(best_per_ligand(frame), control="nope")
    assert result.table["delta_vs_control"].isna().all()
    assert "control_not_found" in codes(result.issues)


def test_control_comparison_reports_n_runs_for_both_sides(real):
    result = compare_to_control(best_per_ligand(real.poses), control="imatinib_1iep")
    row = result.table.set_index(["target", "ligand"]).loc[("1iep", "lig_1s63")]
    assert row["n_runs"] == 1
    assert row["control_n_runs"] == 4


def test_warning_when_the_control_has_a_different_number_of_runs(real):
    result = compare_to_control(best_per_ligand(real.poses), control="imatinib_1iep")
    messages = [i.message for i in result.issues if i.code == "unequal_run_counts"]
    # Only 1iep, where the control has 4 runs and every other ligand has 1. On 5x72
    # every ligand has a single run, so there is nothing to warn about there.
    assert len(messages) == 1
    assert "'1iep'" in messages[0]
    assert "lig_1s63" in messages[0]
    assert "'5x72'" not in messages[0]


def test_equal_run_counts_raise_no_warning():
    frame = poses_frame(
        {"ligand": "ctrl", "seed": 1, "affinity": -9.0},
        {"ligand": "ctrl", "seed": 2, "affinity": -9.0},
        {"ligand": "x", "seed": 1, "affinity": -8.0},
        {"ligand": "x", "seed": 2, "affinity": -8.0},
    )
    result = compare_to_control(best_per_ligand(frame), control="ctrl")
    assert "unequal_run_counts" not in codes(result.issues)


def test_delta_is_blank_when_the_grid_box_differs_from_the_control():
    frame = poses_frame(
        {"ligand": "ctrl", "affinity": -9.0},
        {"ligand": "x", "affinity": -8.0},
    )
    runs = runs_frame({"ligand": "ctrl"}, {"ligand": "x", "size_x": 30.0})
    result = compare_to_control(best_per_ligand(frame), control="ctrl", runs=runs)
    table = result.table.set_index("ligand")
    assert pd.isna(table.loc["x", "delta_vs_control"])
    assert "incomparable_setup" in codes(result.issues)


def test_matching_boxes_raise_no_setup_warning(real):
    result = compare_to_control(
        best_per_ligand(real.poses), control="imatinib_1iep", runs=real.runs
    )
    assert "incomparable_setup" not in codes(result.issues)
    assert result.table["delta_vs_control"].notna().all()


def test_delta_label_says_score_difference_not_affinity():
    assert LABELS["delta_vs_control"] == "Δ score vs control (kcal/mol)"
    assert "not" in NOTES["delta_vs_control"]


# --------------------------------------------------------- seed reproducibility


def test_seed_summary_known_values(real):
    summary = seed_summary(real.poses)
    row = summary[(summary["target"] == "1iep") & (summary["ligand"] == "imatinib_1iep")].iloc[0]
    assert row["n_runs"] == 4
    assert row["mean_affinity"] == pytest.approx(-13.25575, abs=1e-5)
    assert row["sd_affinity"] == pytest.approx(0.0152, abs=1e-4)
    assert row["min_affinity"] == -13.273
    assert row["max_affinity"] == -13.242
    assert row["range_affinity"] == pytest.approx(0.031, abs=1e-9)
    assert sorted(row["seeds"]) == [1, 2, 3, 42]


def test_sd_is_blank_below_three_runs_but_values_are_still_shown():
    for seeds in ([1], [1, 2]):
        frame = poses_frame(*({"seed": seed, "affinity": -8.0 - seed / 10} for seed in seeds))
        row = seed_summary(frame).iloc[0]
        assert row["n_runs"] == len(seeds)
        assert pd.isna(row["sd_affinity"])
        assert row["min_affinity"] == pytest.approx(-8.0 - max(seeds) / 10)
        assert row["max_affinity"] == pytest.approx(-8.0 - min(seeds) / 10)
        assert pd.notna(row["range_affinity"])


def test_sd_is_present_from_three_runs():
    frame = poses_frame(
        {"seed": 1, "affinity": -8.0},
        {"seed": 2, "affinity": -8.1},
        {"seed": 3, "affinity": -8.2},
    )
    row = seed_summary(frame).iloc[0]
    assert row["n_runs"] == 3
    assert row["sd_affinity"] == pytest.approx(0.1, abs=1e-9)


def test_seed_summary_covers_every_pair_including_single_runs(real):
    summary = seed_summary(real.poses)
    assert len(summary) == 10
    assert (summary["n_runs"] == 1).sum() == 9


def test_seed_note_does_not_claim_experimental_error():
    assert "search" in NOTES["seed_summary"]
    assert "experimental error" in NOTES["seed_summary"]


# ------------------------------------------------------------ reference parsing


def test_reference_is_refused_without_explicit_confirmation():
    with pytest.raises(ReferenceRefused, match="experimental"):
        parse_reference(reference_text(), filename="x.pdbqt", is_experimental_pose=False)


def test_is_experimental_pose_has_no_default():
    with pytest.raises(TypeError):
        parse_reference(reference_text(), filename="x.pdbqt")


def test_sdf_reference_is_refused_with_a_clear_message():
    sdf = (DATA / "receptors" / "1iep_ligand.sdf").read_text()
    with pytest.raises(ReferenceRefused, match="SDF references need atom mapping"):
        parse_reference(sdf, filename="1iep_ligand.sdf", is_experimental_pose=True)


def test_sdf_content_is_refused_even_with_a_pdbqt_name():
    sdf = (DATA / "receptors" / "1iep_ligand.sdf").read_text()
    with pytest.raises(ReferenceRefused, match="SDF"):
        parse_reference(sdf, filename="mislabelled.pdbqt", is_experimental_pose=True)


def test_reference_without_atoms_is_refused():
    with pytest.raises(ReferenceRefused, match="no atoms"):
        parse_reference("REMARK nothing here\n", filename="x.pdbqt", is_experimental_pose=True)


def test_multi_model_reference_is_refused():
    text = "MODEL 1\n" + reference_text(("C",)) + "ENDMDL\nMODEL 2\n" + reference_text(("C",)) + "ENDMDL\n"
    with pytest.raises(ReferenceRefused, match="single pose"):
        parse_reference(text, filename="x.pdbqt", is_experimental_pose=True)


def test_real_reference_has_37_heavy_atoms(reference):
    assert reference.n_heavy_atoms == 37
    assert reference.is_experimental_pose is True


# ------------------------------------------------------------- element mapping


@pytest.mark.parametrize(
    "atom_type, element",
    [
        ("C", "C"),
        ("A", "C"),  # aromatic carbon
        ("CG0", "C"),  # Meeko glue carbon
        ("CG12", "C"),
        ("N", "N"),
        ("NA", "N"),  # hydrogen-bond accepting nitrogen, not sodium
        ("NS", "N"),
        ("OA", "O"),
        ("OS", "O"),
        ("SA", "S"),
        ("S", "S"),
        ("P", "P"),
        ("F", "F"),
        ("Cl", "Cl"),
        ("CL", "Cl"),
        ("Br", "Br"),
        ("I", "I"),
        ("Zn", "Zn"),
        ("Fe", "Fe"),
        ("HD", "H"),
        ("HS", "H"),
        ("H", "H"),
    ],
)
def test_element_of_maps_autodock_types(atom_type, element):
    assert element_of(atom_type) == element


@pytest.mark.parametrize("atom_type", ["G0", "G12"])
def test_pseudo_atoms_have_no_element(atom_type):
    assert element_of(atom_type) is None


@pytest.mark.parametrize("atom_type", ["XX", "CG", "G", ""])
def test_unknown_types_have_no_element(atom_type):
    assert element_of(atom_type) is None


# ----------------------------------------------------------------- redocking RMSD


def test_redocking_rmsd_matches_known_values(real, reference):
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep")
    table = result.table
    best = table[(table["target"] == "1iep") & (table["pose"] == 1)].set_index("seed")
    for seed, expected in EXPECTED_RMSD.items():
        assert best.loc[seed, "rmsd_vs_reference"] == pytest.approx(expected, abs=0.02)


def test_rmsd_is_only_computed_for_poses_with_coordinates(real, reference):
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep")
    table = result.table
    assert table["has_coordinates"].all()
    assert table["rmsd_vs_reference"].notna().all()
    # The 1iep seed 42 run has 3 poses in the PDBQT but 9 modes in the log.
    assert len(table[(table["target"] == "1iep") & (table["seed"] == 42)]) == 3


def test_within_cutoff_column_flags_the_common_success_criterion(real, reference):
    assert RMSD_SUCCESS_CUTOFF == 2.0
    table = redocking_rmsd(real, reference, ligand="imatinib_1iep").table
    row = table[(table["target"] == "1iep") & (table["seed"] == 42) & (table["pose"] == 1)].iloc[0]
    assert bool(row["within_2A"]) is True
    # Pose 2 of the same run is a different binding mode, about 12 A away.
    other = table[(table["target"] == "1iep") & (table["seed"] == 42) & (table["pose"] == 2)].iloc[0]
    assert other["rmsd_vs_reference"] > 2.0
    assert bool(other["within_2A"]) is False
    assert "criterion" in LABELS["within_2A"]
    assert "not" in NOTES["within_2A"]


def test_within_cutoff_is_blank_when_rmsd_is_blank():
    # A pose whose atoms do not match the reference has no RMSD, so no verdict either.
    reference = parse_reference(
        reference_text(("C",) * 3), filename="r.pdbqt", is_experimental_pose=True
    )
    results = load_folder(DATA / "results", DATA / "logs")
    result = redocking_rmsd(results, reference, ligand="imatinib_1iep")
    assert result.table["rmsd_vs_reference"].isna().all()
    assert result.table["within_2A"].isna().all()


def test_rmsd_label_names_the_user_supplied_reference_and_warns():
    assert LABELS["rmsd_vs_reference"] == "RMSD vs user-supplied reference (Å)"
    note = NOTES["rmsd_vs_reference"]
    assert "conformer" in note
    assert "symmetry" in note
    assert "superposition" in note or "superimpose" in note


def rename_atoms(text):
    """Rewrite the atom-name column, leaving types and coordinates untouched."""
    lines = []
    for index, line in enumerate(text.splitlines()):
        if line.startswith(("ATOM", "HETATM")):
            line = line[:12] + f"Z{index:03d}" + line[16:]
        lines.append(line + "\n")
    return "".join(lines)


def test_atom_name_mismatch_warns_but_rmsd_is_still_computed(real):
    renamed = parse_reference(
        rename_atoms(CRYSTAL_IMATINIB.read_text()),
        filename="renamed.pdbqt",
        is_experimental_pose=True,
    )
    result = redocking_rmsd(real, renamed, ligand="imatinib_1iep")
    table = result.table
    best = table[(table["target"] == "1iep") & (table["seed"] == 42) & (table["pose"] == 1)]
    # Elements still line up, so the value is unchanged from the named reference.
    assert best["rmsd_vs_reference"].iloc[0] == pytest.approx(0.91, abs=0.02)
    issue = next(i for i in result.issues if i.code == "atom_name_mismatch")
    assert issue.severity == "warning"
    # The message names both sides at the first differing position.
    assert "'N'" in issue.message  # imatinib's first heavy atom
    assert "Z0" in issue.message  # the renamed reference


def test_matching_atom_names_raise_no_warning(real, reference):
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep")
    assert "atom_name_mismatch" not in codes(result.issues)


def test_atom_name_mismatch_is_reported_once_per_run(real):
    renamed = parse_reference(
        rename_atoms(CRYSTAL_IMATINIB.read_text()),
        filename="renamed.pdbqt",
        is_experimental_pose=True,
    )
    result = redocking_rmsd(real, renamed, ligand="imatinib_1iep")
    # Four imatinib runs on 1iep plus one on 5x72, not one warning per pose.
    assert len([i for i in result.issues if i.code == "atom_name_mismatch"]) == 5


def test_different_sized_ligand_is_reported_as_not_applicable(real, reference):
    # p59 has 24 heavy atoms against the reference's 37: a different molecule.
    result = redocking_rmsd(real, reference, ligand="p59_5x72")
    assert result.table["rmsd_vs_reference"].isna().all()
    assert "reference_not_applicable" in codes(result.issues)


def test_same_size_but_different_elements_is_an_error(real):
    # 37 atoms as in imatinib, but all carbon, so the elements cannot line up.
    wrong = parse_reference(
        reference_text(("C",) * 37), filename="wrong.pdbqt", is_experimental_pose=True
    )
    result = redocking_rmsd(real, wrong, ligand="imatinib_1iep")
    assert result.table["rmsd_vs_reference"].isna().all()
    assert "element_mismatch" in codes(result.issues)
    message = next(i.message for i in result.issues if i.code == "element_mismatch")
    assert "N" in message  # the message names the first mismatching position


def test_no_matching_run_at_all_is_refused(real):
    wrong = parse_reference(
        reference_text(("C",) * 5), filename="wrong.pdbqt", is_experimental_pose=True
    )
    with pytest.raises(ReferenceRefused, match="no run"):
        redocking_rmsd(real, wrong)


def test_implausible_rmsd_is_warned_about(real, reference):
    # Imatinib docked into 5x72 lands in a different box, about 124 A away.
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep")
    table = result.table.set_index(["target", "seed", "pose"])
    assert table.loc[("5x72", 42, 1), "rmsd_vs_reference"] > 100
    assert "implausible_rmsd" in codes(result.issues)


def test_alternative_binding_modes_are_not_called_implausible(real, reference):
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep")
    messages = [i.message for i in result.issues if i.code == "implausible_rmsd"]
    # Poses 2+ on 1iep are different binding modes about 12 A from the crystal pose,
    # which is expected. Only the cross-docked 5x72 best pose signals a wrong reference.
    assert len(messages) == 1
    assert "5x72" in messages[0]
    assert "best pose" in messages[0]


def test_reference_outside_the_grid_box_is_warned_about(real, reference):
    result = redocking_rmsd(real, reference, runs=real.runs)
    message = next(i.message for i in result.issues if i.code == "reference_outside_box")
    assert "5x72" in message


def test_rmsd_uses_the_corrected_heavy_atom_rule(real, reference):
    # Imatinib has no macrocycle atoms, so this checks the shared helper is used:
    # 37 heavy atoms on both sides, and 3 hydrogens excluded from the 40 atom lines.
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep")
    assert (result.table["n_atoms_compared"] == 37).all()


def test_rmsd_does_not_superimpose(reference):
    # Shifting every atom by 1 A in x must change the RMSD by exactly 1 A.
    shifted = [(atom.x + 1.0, atom.y, atom.z) for atom in reference.atoms]
    moved = parse_reference(
        reference_text([atom.atom_type for atom in reference.atoms], shifted),
        filename="moved.pdbqt",
        is_experimental_pose=True,
    )
    from vina_analyser.analysis import heavy_atom_rmsd

    assert heavy_atom_rmsd(reference.atoms, moved.atoms) == pytest.approx(1.0, abs=1e-9)


# ----------------------------------------------------------------- whole analysis


def test_analyse_puts_the_tables_together(real):
    result = analyse(real, control="imatinib_1iep")
    for _, group in result.best.groupby("target"):
        assert list(group["rank"]) == sorted(group["rank"])
    assert len(result.best) == 10
    assert len(result.seeds) == 10
    assert "ligand_efficiency" in result.best.columns
    assert "delta_vs_control" in result.best.columns
    assert "unequal_run_counts" in codes(result.issues)


def test_analyse_without_a_control_leaves_the_column_out(real):
    result = analyse(real)
    assert "delta_vs_control" not in result.best.columns
    assert "control_not_found" not in codes(result.issues)


def test_analyse_passes_parser_issues_through():
    results = load_folder(DATA / "logs")  # logs only: no coordinates anywhere
    result = analyse(results)
    assert result.best.empty
    assert "missing_pdbqt" in codes(result.issues)


# --------------------------------------------------------- fixed limitations notes


def test_scoring_error_note_is_the_agreed_fixed_text():
    assert NOTES["scoring_error"] == (
        "Vina's scoring error is roughly 2-3 kcal/mol; score differences smaller than "
        "this should not be interpreted as differences in binding."
    )


def test_delta_note_states_which_direction_is_better():
    assert "negative = more favourable score than control" in NOTES["delta_direction"]


def test_analysis_output_carries_the_scoring_error_and_delta_notes(real):
    result = analyse(real, control="imatinib_1iep")
    assert NOTES["scoring_error"] in result.notes
    assert NOTES["delta_direction"] in result.notes
    assert NOTES["ligand_efficiency"] in result.notes
    assert NOTES["seed_summary"] in result.notes


def test_analysis_without_a_control_omits_the_control_notes(real):
    result = analyse(real)
    assert NOTES["scoring_error"] in result.notes
    assert all("control" not in note for note in result.notes)


def test_rmsd_output_carries_its_own_notes(real, reference):
    result = redocking_rmsd(real, reference, ligand="imatinib_1iep")
    assert NOTES["rmsd_vs_reference"] in result.notes
    assert NOTES["within_2A"] in result.notes
    assert NOTES["scoring_error"] in result.notes


# ------------------------------------------------- RMSD restricted to one target


def test_rmsd_for_one_target_only_contains_that_target(real, reference):
    result = redocking_rmsd(real, reference, target="1iep", ligand="imatinib_1iep", runs=real.runs)
    table = result.table
    assert set(table["target"]) == {"1iep"}
    best = table[table["pose"] == 1].set_index("seed")
    for seed, expected in EXPECTED_RMSD.items():
        assert best.loc[seed, "rmsd_vs_reference"] == pytest.approx(expected, abs=0.02)
    # The 124 Å number belongs to 5x72, so asking about 1iep must not produce it.
    assert "implausible_rmsd" not in codes(result.issues)
    assert "reference_outside_box" not in codes(result.issues)


def test_rmsd_for_the_wrong_target_shows_why_the_number_is_meaningless(real, reference):
    result = redocking_rmsd(real, reference, target="5x72", ligand="imatinib_1iep", runs=real.runs)
    assert set(result.table["target"]) == {"5x72"}
    assert result.table.iloc[0]["rmsd_vs_reference"] > 100
    assert {"implausible_rmsd", "reference_outside_box"} <= set(codes(result.issues))


def test_rmsd_for_an_unknown_target_is_refused(real, reference):
    with pytest.raises(ReferenceRefused, match="no_such_target"):
        redocking_rmsd(real, reference, target="no_such_target", ligand="imatinib_1iep")
