"""Tests for vina_analyser.parsers.

Expected values come from vina_test_data/README.md or were checked directly in the
raw files. Real files cover normal behaviour. Small made-up text blocks cover failure
cases that the real files do not contain.
"""

import ast
import tomllib
from pathlib import Path

import pandas as pd
import pytest

from vina_analyser.parsers import (
    count_heavy_atoms,
    decode_text,
    load_folder,
    load_results,
    parse_filename,
    parse_log,
    parse_pdbqt,
)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "vina_test_data"
RESULTS = DATA / "results"
LOGS = DATA / "logs"

# The README table: 10 target x ligand values plus 3 imatinib repeat seeds = 13.
BEST_AFFINITY = {
    ("1iep", "imatinib_1iep", 42): -13.27,
    ("1iep", "p59_5x72", 42): -8.463,
    ("1iep", "p69_5x72", 42): -7.657,
    ("1iep", "lig_1s63", 42): -8.781,
    ("1iep", "bace1_macrocycle", 42): -6.367,
    ("5x72", "imatinib_1iep", 42): -10.83,
    ("5x72", "p59_5x72", 42): -11.0,
    ("5x72", "p69_5x72", 42): -10.4,
    ("5x72", "lig_1s63", 42): -9.856,
    ("5x72", "bace1_macrocycle", 42): -7.234,
    ("1iep", "imatinib_1iep", 1): -13.24,
    ("1iep", "imatinib_1iep", 2): -13.24,
    ("1iep", "imatinib_1iep", 3): -13.26,
}

PDBQT_POSE_COUNT = {
    ("1iep", "imatinib_1iep", 42): 3,
    ("1iep", "imatinib_1iep", 1): 3,
    ("1iep", "imatinib_1iep", 2): 4,
    ("1iep", "imatinib_1iep", 3): 4,
    ("5x72", "imatinib_1iep", 42): 7,
    ("1iep", "p59_5x72", 42): 9,
    ("1iep", "p69_5x72", 42): 9,
    ("1iep", "lig_1s63", 42): 9,
    ("1iep", "bace1_macrocycle", 42): 9,
    ("5x72", "p59_5x72", 42): 9,
    ("5x72", "p69_5x72", 42): 9,
    ("5x72", "lig_1s63", 42): 9,
    ("5x72", "bace1_macrocycle", 42): 9,
}

HEAVY_ATOMS = {
    "imatinib_1iep": 37,
    "p59_5x72": 24,
    "p69_5x72": 24,
    "lig_1s63": 29,
    "bace1_macrocycle": 38,
}

POSE_COLUMNS = [
    "target",
    "ligand",
    "seed",
    "pose",
    "affinity",
    "rmsd_lb",
    "rmsd_ub",
    "source_file",
    "heavy_atoms",
    "has_coordinates",
]

# Copied from a real output file: atom type starts at col 78 (index 77).
ATOM_PREFIX = "ATOM      1  C1  LIG L   1      15.951  52.070  15.374  1.00  0.00     0.135 "


@pytest.fixture(scope="module")
def real():
    return load_folder(RESULTS, LOGS)


def run_id(key):
    target, ligand, seed = key
    return f"{target}-{ligand}-seed{seed}"


def codes(result):
    return [issue.code for issue in result.issues]


def make_pose(number=1, score=(-7.0, 0.0, 0.0), types=("C", "A", "HD")):
    lines = [f"MODEL {number}"]
    if score is not None:
        lines.append(f"REMARK VINA RESULT:  {score[0]}  {score[1]}  {score[2]}")
    lines += [ATOM_PREFIX + atom_type for atom_type in types]
    lines.append("ENDMDL")
    return "\n".join(lines) + "\n"


def make_log(rows=((1, -7.0, 0.0, 0.0),), seed=42):
    table = "\n".join(f"{n}  {a}  {lb}  {ub}" for n, a, lb, ub in rows)
    return (
        "AutoDock Vina v1.2.5\n"
        f"Performing docking (random seed: {seed}) ... \n"
        "mode |   affinity | dist from best mode\n"
        "-----+------------+----------+----------\n"
        f"{table}\n"
    )


# ---------------------------------------------------------------- filename parsing


@pytest.mark.parametrize(
    "name, target, ligand, seed, kind",
    [
        ("1iep__imatinib_1iep_out.pdbqt", "1iep", "imatinib_1iep", None, "pdbqt"),
        ("1iep__imatinib_1iep_seed2_out.pdbqt", "1iep", "imatinib_1iep", 2, "pdbqt"),
        ("5x72__p59_5x72.log", "5x72", "p59_5x72", None, "log"),
        ("1iep__imatinib_1iep_seed3.log", "1iep", "imatinib_1iep", 3, "log"),
        ("some/folder/1iep__lig_1s63.log", "1iep", "lig_1s63", None, "log"),
        ("some\\folder\\1iep__lig_1s63.log", "1iep", "lig_1s63", None, "log"),
    ],
)
def test_parse_filename_matches_the_naming_pattern(name, target, ligand, seed, kind):
    parsed = parse_filename(name)
    assert (parsed.target, parsed.ligand, parsed.seed, parsed.kind) == (target, ligand, seed, kind)


@pytest.mark.parametrize(
    "name",
    [
        "docking_result.pdbqt",  # no target__ligand separator
        "1iep__imatinib.pdbqt",  # output PDBQT must end in _out.pdbqt
        "a__b__c_out.pdbqt",  # ambiguous: refuse rather than guess
        "__ligand_out.pdbqt",  # empty target
        "target___ligand.log",  # three underscores: ambiguous
        "1iep__imatinib.txt",  # unsupported file type
    ],
)
def test_parse_filename_refuses_to_guess(name):
    assert parse_filename(name) is None


# ------------------------------------------------------------------ heavy atoms


def test_heavy_count_excludes_hydrogens_and_g_pseudo_atoms():
    types = ["C", "A", "N", "OA", "Cl", "H", "HD", "HS", "G0", "G1", "CG0", "CG1"]
    assert count_heavy_atoms(types) == 7  # 5 ordinary heavy atoms + CG0 + CG1


def test_glue_carbons_count_as_carbon_but_g_pseudo_atoms_do_not():
    # CG* are real carbons where Meeko opened the ring; G* are dummies on top of them.
    assert count_heavy_atoms(["CG0", "CG1", "CG12"]) == 3
    assert count_heavy_atoms(["G0", "G1", "G12"]) == 0
    assert count_heavy_atoms(["C", "CG0", "G0"]) == 2


def test_truncated_macrocycle_types_are_not_accepted():
    # "CG" is what reading only cols 78-79 would produce from "CG0".
    assert count_heavy_atoms(["CG"]) is None
    assert count_heavy_atoms(["G"]) is None


def test_heavy_count_is_unknown_when_an_atom_type_is_unrecognised():
    assert count_heavy_atoms(["C", "XX"]) is None


@pytest.mark.parametrize("ligand, expected", HEAVY_ATOMS.items())
def test_heavy_atoms_per_ligand_in_real_files(real, ligand, expected):
    rows = real.poses[real.poses["ligand"] == ligand]
    assert set(rows.loc[rows["has_coordinates"], "heavy_atoms"]) == {expected}


def test_macrocycle_atom_types_are_read_whole_and_counted_correctly():
    pose = parse_pdbqt((RESULTS / "1iep__bace1_macrocycle_out.pdbqt").read_text()).poses[0]
    # 43 atom lines = 38 heavy (including 2 CG0 carbons) + 3 HD + 2 G0.
    # "CG0" is 3 characters: reading only cols 78-79 would cut it down to "CG".
    assert len(pose.atom_types) == 43
    assert pose.atom_types.count("CG0") == 2
    assert pose.atom_types.count("G0") == 2
    assert "CG" not in pose.atom_types
    assert count_heavy_atoms(pose.atom_types) == 38


# ------------------------------------------------------------------ table shape


def test_poses_has_the_agreed_columns_and_types(real):
    poses = real.poses
    assert list(poses.columns) == POSE_COLUMNS
    assert poses["has_coordinates"].dtype == bool
    assert poses["seed"].dtype == "Int64"
    assert poses["heavy_atoms"].dtype == "Int64"
    assert poses["affinity"].dtype == "float64"


def test_real_dataset_gives_117_rows_93_with_coordinates(real):
    assert len(real.poses) == 117
    assert real.poses["has_coordinates"].sum() == 93
    assert (~real.poses["has_coordinates"]).sum() == 24


def test_one_row_per_pose(real):
    assert not real.poses.duplicated(["target", "ligand", "seed", "pose"]).any()


def test_real_dataset_raises_no_issues(real):
    assert real.issues == []


def test_has_coordinates_is_consistent_with_other_columns(real):
    poses = real.poses
    with_coords = poses[poses["has_coordinates"]]
    log_only = poses[~poses["has_coordinates"]]
    assert with_coords["source_file"].str.endswith("_out.pdbqt").all()
    assert with_coords["heavy_atoms"].notna().all()
    assert log_only["source_file"].str.endswith(".log").all()
    assert log_only["heavy_atoms"].isna().all()


# ------------------------------------------------------------- known real values


@pytest.mark.parametrize("key, expected", BEST_AFFINITY.items(), ids=map(run_id, BEST_AFFINITY))
def test_best_affinity_matches_readme(real, key, expected):
    target, ligand, seed = key
    poses = real.poses
    row = poses[
        (poses["target"] == target)
        & (poses["ligand"] == ligand)
        & (poses["seed"] == seed)
        & (poses["pose"] == 1)
    ]
    assert len(row) == 1
    assert row["affinity"].iloc[0] == pytest.approx(expected, abs=0.01)
    assert bool(row["has_coordinates"].iloc[0]) is True


@pytest.mark.parametrize(
    "key, exact",
    [
        (("1iep", "imatinib_1iep", 42), -13.273),
        (("1iep", "imatinib_1iep", 1), -13.244),
        (("5x72", "p59_5x72", 42), -10.996),
        (("5x72", "p69_5x72", 42), -10.397),
    ],
    ids=lambda value: run_id(value) if isinstance(value, tuple) else str(value),
)
def test_pdbqt_precision_is_preferred_over_rounded_log(real, key, exact):
    target, ligand, seed = key
    poses = real.poses
    row = poses[
        (poses["target"] == target)
        & (poses["ligand"] == ligand)
        & (poses["seed"] == seed)
        & (poses["pose"] == 1)
    ]
    assert row["affinity"].iloc[0] == exact


def test_log_only_rows_keep_the_log_value(real):
    poses = real.poses
    row = poses[
        (poses["target"] == "1iep")
        & (poses["ligand"] == "imatinib_1iep")
        & (poses["seed"] == 42)
        & (poses["pose"] == 4)
    ]
    # README: mode 4 is more than 3 kcal/mol from the best, so it is only in the log.
    assert bool(row["has_coordinates"].iloc[0]) is False
    assert row["source_file"].iloc[0] == "1iep__imatinib_1iep.log"
    assert row["affinity"].iloc[0] > -13.273 + 3


# --------------------------------------------------------------- pose count quirk


@pytest.mark.parametrize(
    "key, expected", PDBQT_POSE_COUNT.items(), ids=map(run_id, PDBQT_POSE_COUNT)
)
def test_both_pose_counts_are_reported(real, key, expected):
    target, ligand, seed = key
    runs = real.runs
    row = runs[(runs["target"] == target) & (runs["ligand"] == ligand) & (runs["seed"] == seed)]
    assert len(row) == 1
    assert row["n_poses_pdbqt"].iloc[0] == expected
    assert row["n_modes_log"].iloc[0] == 9


def test_runs_table_has_one_row_per_docking_run(real):
    assert len(real.runs) == 13
    assert real.runs["n_poses_pdbqt"].sum() == 93
    assert real.runs["n_modes_log"].sum() == 117


def test_log_only_modes_are_always_the_lowest_ranked(real):
    with_dropped = 0
    for _, run in real.poses.groupby(["target", "ligand", "seed"]):
        kept = run.loc[run["has_coordinates"], "pose"]
        dropped = run.loc[~run["has_coordinates"], "pose"]
        if dropped.empty:
            continue  # every mode kept: nothing to order
        with_dropped += 1
        assert kept.max() < dropped.min()
    # Imatinib's five runs (four on 1iep, one on 5x72) drop modes; without them this
    # test would check nothing.
    assert with_dropped == 5


# --------------------------------------------------------- log vs PDBQT agreement


@pytest.mark.parametrize("log_path", sorted(LOGS.glob("*.log")), ids=lambda p: p.stem)
def test_log_and_pdbqt_agree_independently(log_path):
    # Compares the two parsers directly, so the merge step cannot hide a difference.
    log = parse_log(log_path.read_text())
    pdbqt = parse_pdbqt((RESULTS / f"{log_path.stem}_out.pdbqt").read_text())
    modes = {mode.mode: mode for mode in log.modes}
    assert len(pdbqt.poses) > 0
    for pose in pdbqt.poses:
        mode = modes[pose.pose]
        assert pose.affinity == pytest.approx(mode.affinity, abs=0.01)
        assert pose.rmsd_lb == pytest.approx(mode.rmsd_lb, abs=0.01)
        assert pose.rmsd_ub == pytest.approx(mode.rmsd_ub, abs=0.01)


def test_runs_table_records_the_largest_difference(real):
    assert real.runs["max_log_pdbqt_diff"].max() <= 0.01
    assert real.runs["max_log_pdbqt_diff"].max() > 0  # rounding exists, so never exactly 0


# ------------------------------------------------------------------ log headers


def test_log_header_values():
    log = parse_log((LOGS / "5x72__imatinib_1iep.log").read_text())
    # The citation banner also says "AutoDock Vina 1.2.0"; that must not be taken.
    assert log.vina_version == "v1.2.5"
    assert log.scoring_function == "vina"
    assert log.receptor == "receptors/5x72_receptor.pdbqt"
    assert log.ligand == "ligands/imatinib_1iep.pdbqt"
    assert log.center == (-15.0, 15.0, 129.0)
    assert log.size == (30.0, 24.0, 24.0)
    assert log.grid_spacing == 0.375
    assert log.exhaustiveness == 8
    assert log.cpu == 4
    assert log.seed == 42


def test_missing_header_values_stay_none_never_defaults():
    log = parse_log(make_log().replace("AutoDock Vina v1.2.5\n", ""))
    assert log.vina_version is None
    assert log.exhaustiveness is None
    assert log.center is None
    assert log.grid_spacing is None
    assert log.cpu is None
    assert len(log.modes) == 1


def test_seed_comes_from_the_log_for_files_without_a_seed_in_the_name(real):
    runs = real.runs
    assert (runs["seed"] == 42).sum() == 10
    assert sorted(runs["seed"].tolist()) == [1, 2, 3] + [42] * 10


def test_seed_is_blank_without_a_log_never_assumed():
    text = (RESULTS / "1iep__imatinib_1iep_out.pdbqt").read_text()
    result = load_results({"1iep__imatinib_1iep_out.pdbqt": text})
    assert result.poses["seed"].isna().all()
    assert len(result.poses) == 3
    assert result.poses["has_coordinates"].all()
    assert "missing_log" in codes(result)


def test_log_without_pdbqt_gives_only_log_only_rows():
    text = (LOGS / "1iep__imatinib_1iep.log").read_text()
    result = load_results({"1iep__imatinib_1iep.log": text})
    assert len(result.poses) == 9
    assert not result.poses["has_coordinates"].any()
    assert result.poses["heavy_atoms"].isna().all()
    assert "missing_pdbqt" in codes(result)


# ------------------------------------------------------------ made-up failure cases


def test_windows_line_endings():
    pdbqt = make_pose().replace("\n", "\r\n")
    pose = parse_pdbqt(pdbqt).poses[0]
    assert pose.atom_types == ["C", "A", "HD"]
    log = parse_log(make_log().replace("\n", "\r\n"))
    assert log.seed == 42
    assert log.vina_version == "v1.2.5"
    assert len(log.modes) == 1


def test_disagreement_keeps_pdbqt_value_and_warns():
    result = load_results(
        {
            "t__l_out.pdbqt": make_pose(score=(-7.0, 0.0, 0.0)),
            "t__l.log": make_log(rows=((1, -7.5, 0.0, 0.0),)),
        }
    )
    assert result.poses["affinity"].iloc[0] == -7.0
    assert "value_disagreement" in codes(result)


def test_rmsd_disagreement_also_warns():
    result = load_results(
        {
            "t__l_out.pdbqt": make_pose(number=2, score=(-7.0, 1.0, 2.0)),
            "t__l.log": make_log(rows=((2, -7.0, 1.5, 2.0),)),
        }
    )
    assert result.poses["rmsd_lb"].iloc[0] == 1.0
    assert "value_disagreement" in codes(result)


def test_small_rounding_difference_does_not_warn():
    result = load_results(
        {
            "t__l_out.pdbqt": make_pose(score=(-13.273, 0.0, 0.0)),
            "t__l.log": make_log(rows=((1, -13.27, 0.0, 0.0),)),
        }
    )
    assert result.issues == []


def test_missing_vina_result_uses_log_value_and_warns():
    result = load_results(
        {
            "t__l_out.pdbqt": make_pose(score=None),
            "t__l.log": make_log(rows=((1, -7.5, 0.0, 0.0),)),
        }
    )
    row = result.poses.iloc[0]
    assert row["affinity"] == -7.5
    assert bool(row["has_coordinates"]) is True
    assert "missing_score" in codes(result)


def test_missing_vina_result_without_log_stays_blank():
    result = load_results({"t__l_out.pdbqt": make_pose(score=None)})
    assert pd.isna(result.poses["affinity"].iloc[0])
    assert "missing_score" in codes(result)


def test_unmatched_filename_loads_nothing_and_is_reported():
    result = load_results({"docking_result.pdbqt": make_pose()})
    assert result.poses.empty
    assert result.runs.empty
    assert list(result.poses.columns) == POSE_COLUMNS
    assert codes(result) == ["unrecognised_name"]
    assert result.issues[0].file == "docking_result.pdbqt"


def test_user_supplied_names_rescue_an_unmatched_filename():
    result = load_results(
        {"docking_result.pdbqt": make_pose()},
        names={"docking_result.pdbqt": {"target": "ABL", "ligand": "drug", "seed": 7}},
    )
    row = result.poses.iloc[0]
    assert (row["target"], row["ligand"], row["seed"]) == ("ABL", "drug", 7)
    assert "unrecognised_name" not in codes(result)


def test_seed_conflict_uses_the_log_and_warns():
    result = load_results(
        {
            "t__l_seed5_out.pdbqt": make_pose(),
            "t__l_seed5.log": make_log(seed=99),
        }
    )
    assert result.poses["seed"].iloc[0] == 99
    assert "seed_conflict" in codes(result)


def test_duplicate_run_is_refused_not_merged():
    # Different filenames, but both logs record the same target, ligand and seed.
    result = load_results(
        {
            "t__l.log": make_log(seed=42),
            "t__l_seed42.log": make_log(seed=42),
        }
    )
    assert result.poses.empty
    assert "duplicate_run" in codes(result)


def test_same_filename_in_two_folders_is_refused():
    result = load_results({"a/t__l.log": make_log(), "b/t__l.log": make_log()})
    assert result.poses.empty
    assert "duplicate_run" in codes(result)


def test_unknown_atom_type_blanks_heavy_atoms_and_warns():
    result = load_results({"t__l_out.pdbqt": make_pose(types=("C", "XX"))})
    assert pd.isna(result.poses["heavy_atoms"].iloc[0])
    assert "unknown_atom_type" in codes(result)
    assert "XX" in next(i.message for i in result.issues if i.code == "unknown_atom_type")


def test_flexible_residues_blank_heavy_atoms_and_warn():
    pose = make_pose().replace(
        "ENDMDL",
        f"BEGIN_RES TYR A 100\n{ATOM_PREFIX}C\nEND_RES TYR A 100\nENDMDL",
    )
    result = load_results({"t__l_out.pdbqt": pose})
    assert pd.isna(result.poses["heavy_atoms"].iloc[0])
    assert "flexible_residues" in codes(result)


def test_differing_heavy_counts_within_a_run_warn():
    text = make_pose(1, types=("C", "A")) + make_pose(2, types=("C",))
    result = load_results({"t__l_out.pdbqt": text})
    assert "inconsistent_heavy_atoms" in codes(result)


def test_truncated_pose_is_reported():
    text = make_pose().replace("ENDMDL\n", "")
    parsed = parse_pdbqt(text)
    assert len(parsed.poses) == 1
    assert "truncated_pose" in [issue.code for issue in parsed.issues]


def test_pdbqt_with_no_models_is_reported():
    # An input ligand has atoms but no MODEL blocks, so it is not a docking result.
    text = (DATA / "ligands" / "imatinib_1iep.pdbqt").read_text()
    parsed = parse_pdbqt(text)
    assert parsed.poses == []
    assert "no_poses" in [issue.code for issue in parsed.issues]


def test_log_without_table_is_reported():
    parsed = parse_log("AutoDock Vina v1.2.5\nno results here\n")
    assert parsed.modes == []
    assert "no_mode_table" in [issue.code for issue in parsed.issues]


def test_malformed_table_row_is_reported_and_others_kept():
    text = make_log(rows=((1, -7.0, 0.0, 0.0),)) + "2  not_a_number  1.0  2.0\n3  -6.0  1.0  2.0\n"
    parsed = parse_log(text)
    assert [mode.mode for mode in parsed.modes] == [1, 3]
    assert "malformed_mode_row" in [issue.code for issue in parsed.issues]


def test_pose_in_pdbqt_but_not_log_is_reported():
    result = load_results(
        {
            "t__l_out.pdbqt": make_pose(1) + make_pose(2, score=(-6.0, 1.0, 2.0)),
            "t__l.log": make_log(rows=((1, -7.0, 0.0, 0.0),)),
        }
    )
    assert len(result.poses) == 2
    assert "pose_not_in_log" in codes(result)


def test_hetatm_lines_count_as_atoms():
    text = make_pose(types=("C",)).replace("ATOM  ", "HETATM")
    assert parse_pdbqt(text).poses[0].atom_types == ["C"]


def test_unreadable_score_line_is_reported():
    text = make_pose().replace("-7.0", "not_a_number")
    parsed = parse_pdbqt(text)
    assert parsed.poses[0].affinity is None
    assert "malformed_score" in [issue.code for issue in parsed.issues]


def test_unreadable_model_line_is_reported():
    parsed = parse_pdbqt(make_pose().replace("MODEL 1", "MODEL x"))
    assert "malformed_model" in [issue.code for issue in parsed.issues]
    assert parsed.poses == []


def test_repeated_model_number_keeps_the_first_and_is_reported():
    text = make_pose(1, score=(-7.0, 0.0, 0.0)) + make_pose(1, score=(-5.0, 0.0, 0.0))
    parsed = parse_pdbqt(text)
    assert [pose.affinity for pose in parsed.poses] == [-7.0]
    assert "duplicate_model" in [issue.code for issue in parsed.issues]


def test_pose_with_no_atoms_blanks_heavy_atoms_and_warns():
    result = load_results({"t__l_out.pdbqt": make_pose(types=())})
    assert pd.isna(result.poses["heavy_atoms"].iloc[0])
    assert "no_atoms" in codes(result)


def test_incomplete_user_names_do_not_rescue_a_file():
    result = load_results(
        {"docking_result.pdbqt": make_pose()},
        names={"docking_result.pdbqt": {"target": "ABL", "ligand": ""}},
    )
    assert result.poses.empty
    assert codes(result) == ["unrecognised_name"]


# ---------------------------------------------------------------- reading files


def test_decode_text_handles_powershell_utf16_redirect():
    # In Windows PowerShell 5.1, `vina ... > out.log` writes UTF-16 with a byte-order mark.
    text = make_log().replace("\n", "\r\n")
    assert parse_log(decode_text(text.encode("utf-16"))).seed == 42


def test_decode_text_handles_utf8_with_and_without_bom():
    assert decode_text("abc".encode("utf-8")) == "abc"
    assert decode_text("abc".encode("utf-8-sig")) == "abc"


def test_decode_text_never_fails_on_other_encodings():
    assert decode_text("Å".encode("latin-1")) == "Å"


def test_load_folder_rejects_two_files_with_one_name(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    (tmp_path / "a" / "t__l.log").write_text(make_log())
    (tmp_path / "b" / "t__l.log").write_text(make_log())
    with pytest.raises(ValueError, match="t__l.log"):
        load_folder(tmp_path / "a", tmp_path / "b")


# --------------------------------------------------------------- Python 3.11 floor


def test_requires_python_is_3_11():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert config["project"]["requires-python"] == ">=3.11"


@pytest.mark.parametrize(
    "path",
    sorted([*(ROOT / "vina_analyser").glob("*.py"), *(ROOT / "tests").glob("*.py")]),
    ids=lambda p: p.name,
)
def test_source_parses_with_python_3_11_grammar(path):
    # A partial check: it rejects newer grammar such as `type X = int`, but it does not
    # catch newer standard-library calls. No real 3.11 interpreter is used here.
    ast.parse(path.read_text(encoding="utf-8"), feature_version=(3, 11))


# ---------------------------------------------------- where the seed came from


def test_runs_record_where_the_seed_came_from():
    runs = load_folder(DATA / "results", DATA / "logs").runs
    assert set(runs["seed_source"]) == {"log"}


def test_seed_from_the_filename_is_marked_as_such():
    # Log with its seed line removed: the parser falls back to the filename seed,
    # and the run row must say so, so nothing downstream reports it as logged.
    log = (DATA / "logs" / "1iep__imatinib_1iep_seed2.log").read_text(encoding="utf-8")
    log = "\n".join(line for line in log.splitlines() if "random seed" not in line)
    pdbqt = (DATA / "results" / "1iep__imatinib_1iep_seed2_out.pdbqt").read_text(encoding="utf-8")
    runs = load_results(
        {"1iep__imatinib_1iep_seed2.log": log, "1iep__imatinib_1iep_seed2_out.pdbqt": pdbqt}
    ).runs
    assert runs["seed"].tolist() == [2]
    assert runs["seed_source"].tolist() == ["filename"]


def test_missing_seed_has_no_source():
    pdbqt = (DATA / "results" / "1iep__p59_5x72_out.pdbqt").read_text(encoding="utf-8")
    runs = load_results({"1iep__p59_5x72_out.pdbqt": pdbqt}).runs
    assert runs["seed"].isna().all()
    assert runs["seed_source"].isna().all()
