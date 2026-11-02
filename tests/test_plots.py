"""Tests for vina_analyser.plots.

Figures are checked for structure, not appearance: how many bars, which labels, which
bar carries the control colour, how many markers are hollow. Two things about the files
themselves are checked, because tests depend on them: the resolution and that two saves
of one figure are byte-identical.
"""

import struct
from pathlib import Path

import matplotlib
import matplotlib.pyplot
import pandas as pd
import pytest

from vina_analyser.analysis import NOTES, best_per_ligand, rank_per_target
from vina_analyser.parsers import POSE_COLUMNS, load_folder
from vina_analyser.plots import (
    BAR_COLOUR,
    COLOURBAR_LABEL,
    CONTROL_COLOUR,
    FIGURE_WIDTH_CM,
    MIN_POSES_FOR_BOX,
    affinity_heatmap,
    heatmap_caption,
    pose_strip,
    pose_strip_caption,
    ranked_bars,
    ranked_bars_caption,
    save_figure,
    short_labels,
)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "vina_test_data"


@pytest.fixture(autouse=True)
def close_figures():
    yield
    matplotlib.pyplot.close("all")


@pytest.fixture(scope="module")
def real():
    return load_folder(DATA / "results", DATA / "logs")


@pytest.fixture(scope="module")
def best(real):
    return rank_per_target(best_per_ligand(real.poses))


def poses_frame(*rows):
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


def many_poses(target, ligand, n, has_coordinates=True, start=-9.0):
    return [
        {
            "target": target,
            "ligand": ligand,
            "pose": index + 1,
            "affinity": start + index * 0.1,
            "has_coordinates": has_coordinates,
        }
        for index in range(n)
    ]


def luminance(rgba):
    red, green, blue = rgba[:3]
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def png_dpi(path):
    """Read the pHYs chunk: pixels per metre, which is how PNG stores resolution."""
    data = path.read_bytes()
    index = data.index(b"pHYs")
    x_per_metre, y_per_metre, unit = struct.unpack(">IIB", data[index + 4 : index + 13])
    assert unit == 1, "pHYs unit is not metres"
    return round(x_per_metre * 0.0254), round(y_per_metre * 0.0254)


# ----------------------------------------------------------------- environment


def test_backend_is_agg_so_no_display_is_needed():
    assert matplotlib.get_backend().lower() == "agg"


def test_figures_are_single_column_width():
    assert FIGURE_WIDTH_CM == pytest.approx(8.5)
    figure = affinity_heatmap(poses_to_best(poses_frame()))
    width_cm = figure.get_size_inches()[0] * 2.54
    assert width_cm == pytest.approx(8.5, abs=0.01)
    assert all(size <= 9 for size in _font_sizes(figure))


def _font_sizes(figure):
    return [text.get_fontsize() for text in figure.findobj(matplotlib.text.Text)]


def poses_to_best(poses):
    return rank_per_target(best_per_ligand(poses))


# -------------------------------------------------------------------- heatmap


def test_heatmap_has_one_cell_per_target_and_ligand(best):
    figure = affinity_heatmap(best)
    image = figure.axes[0].images[0]
    assert image.get_array().shape == (2, 5)  # 2 targets x 5 ligands


def test_heatmap_annotates_every_known_value(best):
    figure = affinity_heatmap(best)
    axes = figure.axes[0]
    annotations = {text.get_text() for text in axes.texts}
    assert "-13.27" in annotations
    assert "-6.37" in annotations  # bace1 on 1iep, rounded from -6.367
    assert "-11.00" in annotations  # p59 on 5x72, the PDBQT's -10.996
    assert len(axes.texts) == 10


def test_heatmap_uses_one_shared_scale_in_real_units(best):
    figure = affinity_heatmap(best)
    image = figure.axes[0].images[0]
    assert image.get_clim() == pytest.approx((-13.273, -6.367))
    # A single scale across the whole grid: no per-row normalisation.
    assert image.get_array().shape == (2, 5)


def test_colourbar_is_labelled_with_units_and_direction(best):
    figure = affinity_heatmap(best)
    assert COLOURBAR_LABEL == "Vina score (kcal/mol), more negative = more favourable"
    assert figure.axes[-1].get_ylabel() == COLOURBAR_LABEL


def test_more_favourable_scores_are_darker(best):
    figure = affinity_heatmap(best)
    image = figure.axes[0].images[0]
    best_colour = image.cmap(image.norm(-13.273))
    worst_colour = image.cmap(image.norm(-6.367))
    assert luminance(best_colour) < luminance(worst_colour)


def test_missing_pair_is_a_blank_cell_not_a_zero():
    poses = poses_frame(
        {"target": "A", "ligand": "x", "affinity": -8.0},
        {"target": "A", "ligand": "y", "affinity": -7.0},
        {"target": "B", "ligand": "x", "affinity": -9.0},
    )
    figure = affinity_heatmap(poses_to_best(poses))
    array = figure.axes[0].images[0].get_array()
    assert array.mask.sum() == 1  # B / y was never docked
    assert len(figure.axes[0].texts) == 3


def test_heatmap_caption_restricts_comparison_to_one_row():
    caption = heatmap_caption()
    assert "within a single row" in caption
    assert "one target" in caption
    assert NOTES["scoring_error"] in caption


# ---------------------------------------------------------------- ranked bars


def test_one_bar_per_ligand_ordered_best_first(best):
    figure = ranked_bars(best, "1iep")
    axes = figure.axes[0]
    assert len(axes.patches) == 5
    labels = [label.get_text() for label in axes.get_yticklabels()]
    assert labels == ["imatinib_1iep", "lig_1s63", "p59_5x72", "p69_5x72", "bace1_macrocy…"]
    # The best score sits at the top of the figure.
    assert axes.get_ylim()[0] > axes.get_ylim()[1]


def test_bars_are_labelled_with_their_scores(best):
    figure = ranked_bars(best, "1iep")
    texts = {text.get_text() for text in figure.axes[0].texts}
    assert "-13.27" in texts
    assert "-6.37" in texts


def test_control_bar_is_highlighted_and_named_in_the_legend(best):
    figure = ranked_bars(best, "1iep", control="imatinib_1iep")
    axes = figure.axes[0]
    colours = [patch.get_facecolor() for patch in axes.patches]
    control_colour = matplotlib.colors.to_rgba(CONTROL_COLOUR)
    assert colours[0] == control_colour  # imatinib ranks first on 1iep
    assert all(colour == matplotlib.colors.to_rgba(BAR_COLOUR) for colour in colours[1:])
    assert "control" in axes.get_legend().get_texts()[0].get_text()


def test_control_is_highlighted_wherever_it_ranks(best):
    figure = ranked_bars(best, "5x72", control="imatinib_1iep")
    colours = [patch.get_facecolor() for patch in figure.axes[0].patches]
    # On 5x72 imatinib is only second, so the second bar carries the highlight.
    assert colours[1] == matplotlib.colors.to_rgba(CONTROL_COLOUR)
    assert colours[0] == matplotlib.colors.to_rgba(BAR_COLOUR)


def test_without_a_control_no_bar_is_highlighted(best):
    figure = ranked_bars(best, "1iep")
    colours = [patch.get_facecolor() for patch in figure.axes[0].patches]
    assert all(colour == matplotlib.colors.to_rgba(BAR_COLOUR) for colour in colours)
    assert figure.axes[0].get_legend() is None


def test_control_missing_from_this_target_is_not_an_error():
    poses = poses_frame(
        {"target": "A", "ligand": "ctrl", "affinity": -9.0},
        {"target": "B", "ligand": "x", "affinity": -8.0},
    )
    figure = ranked_bars(poses_to_best(poses), "B", control="ctrl")
    colours = [patch.get_facecolor() for patch in figure.axes[0].patches]
    assert colours == [matplotlib.colors.to_rgba(BAR_COLOUR)]
    assert figure.axes[0].get_legend() is None


def test_unknown_target_is_refused_clearly(best):
    with pytest.raises(ValueError, match="nosuch"):
        ranked_bars(best, "nosuch")


def test_ranked_bars_caption_names_the_target_and_the_uncertainty():
    caption = ranked_bars_caption("1iep")
    assert "1iep" in caption
    assert NOTES["scoring_error"] in caption


# ----------------------------------------------------------------- pose strip


def test_pose_strip_shows_every_pose_of_one_target(real):
    figure = pose_strip(real.poses, "1iep")
    filled, hollow = _marker_counts(figure)
    # 1iep has 72 of the 117 rows: 50 with coordinates and 22 log-only.
    assert filled + hollow == 72
    assert hollow == 22


def test_pose_strip_does_not_pool_poses_across_targets():
    # The same ligand name in two targets must not share a column: scores from
    # different boxes are not comparable, which is the rule everywhere else.
    poses = poses_frame(
        *many_poses("A", "shared", 6, start=-9.0),
        *many_poses("B", "shared", 6, start=-4.0),
    )
    figure = pose_strip(poses, "A")
    offsets = figure.axes[0].collections[0].get_offsets()
    scores = [point[1] for point in offsets]
    assert len(scores) == 6
    assert max(scores) < -8.0  # only target A's values, never B's -4.0 range


def test_log_only_poses_are_hollow_and_explained(real):
    figure = pose_strip(real.poses, "1iep")
    axes = figure.axes[0]
    hollow = _collection(axes, "log only")
    assert len(hollow.get_facecolors()) == 0 or hollow.get_facecolors()[0][3] == 0.0
    assert len(hollow.get_edgecolors()) > 0
    legend_labels = [text.get_text() for text in axes.get_legend().get_texts()]
    assert any("log only" in label for label in legend_labels)
    assert any("no coordinates" in label for label in legend_labels)


def test_box_is_drawn_only_with_enough_poses():
    assert MIN_POSES_FOR_BOX == 5
    poses = poses_frame(
        *many_poses("T", "plenty", 5),
        *many_poses("T", "too_few", 4),
    )
    figure = pose_strip(poses, "T")
    boxes = [patch for patch in figure.axes[0].patches if patch.get_gid().startswith("box:")]
    assert [patch.get_gid() for patch in boxes] == ["box:plenty"]


def test_log_only_poses_do_not_count_towards_the_box():
    poses = poses_frame(
        *many_poses("T", "mostly_log", 3),
        *many_poses("T", "mostly_log", 6, has_coordinates=False),
    )
    figure = pose_strip(poses, "T")
    assert not [patch for patch in figure.axes[0].patches if patch.get_gid().startswith("box:")]


def test_every_real_ligand_gets_a_box(real):
    for target, expected in (("1iep", 5), ("5x72", 5)):
        figure = pose_strip(real.poses, target)
        boxes = [p for p in figure.axes[0].patches if p.get_gid().startswith("box:")]
        assert len(boxes) == expected


def test_pose_strip_caption_explains_markers_and_the_box_rule():
    caption = pose_strip_caption("1iep")
    assert "1iep" in caption
    assert "hollow" in caption.lower()  # the caption opens a sentence with "Hollow"
    assert "log" in caption.lower()
    assert str(MIN_POSES_FOR_BOX) in caption
    assert NOTES["scoring_error"] in caption


def test_pose_strip_jitter_is_reproducible(real):
    first = pose_strip(real.poses, "1iep").axes[0].collections[0].get_offsets()
    second = pose_strip(real.poses, "1iep").axes[0].collections[0].get_offsets()
    assert (first == second).all()


def _collection(axes, label_fragment):
    for collection in axes.collections:
        if label_fragment in (collection.get_label() or ""):
            return collection
    raise AssertionError(f"no collection labelled {label_fragment!r}")


def _marker_counts(figure):
    axes = figure.axes[0]
    filled = len(_collection(axes, "pose").get_offsets())
    hollow = len(_collection(axes, "log only").get_offsets())
    return filled, hollow


# --------------------------------------------------------------- tick labels


def test_short_names_are_left_alone():
    assert short_labels(["abc", "defg"]) == ["abc", "defg"]


def test_long_names_are_truncated_with_an_ellipsis():
    label = short_labels(["a_very_long_ligand_name_indeed"])[0]
    assert label.endswith("…")
    assert len(label) <= 14


def test_names_that_truncate_alike_stay_distinguishable():
    names = ["compound_series_alpha", "compound_series_beta"]
    labels = short_labels(names)
    assert len(set(labels)) == 2
    assert labels[0] != labels[1]


def test_many_similar_names_all_stay_unique():
    names = [f"identical_prefix_entry_{index}" for index in range(12)]
    labels = short_labels(names)
    assert len(set(labels)) == 12


def test_a_name_that_is_a_prefix_of_another_stays_distinct():
    labels = short_labels(["abcdefghijklmnop", "abcdefghijklmnopqrst"])
    assert len(set(labels)) == 2


def test_identical_names_keep_identical_labels():
    assert short_labels(["same", "same"]) == ["same", "same"]


# -------------------------------------------------------------------- saving


def test_save_writes_a_png_and_an_svg(best, tmp_path):
    png, svg = save_figure(affinity_heatmap(best), tmp_path / "heatmap")
    assert png.name == "heatmap.png"
    assert svg.name == "heatmap.svg"
    assert png.stat().st_size > 0
    assert svg.stat().st_size > 0


def test_png_is_300_dpi(best, tmp_path):
    png, _ = save_figure(affinity_heatmap(best), tmp_path / "heatmap")
    assert png_dpi(png) == (300, 300)


def test_two_saves_of_one_figure_are_byte_identical(best, tmp_path):
    figure = ranked_bars(best, "1iep", control="imatinib_1iep")
    first_png, first_svg = save_figure(figure, tmp_path / "first")
    second_png, second_svg = save_figure(figure, tmp_path / "second")
    assert first_png.read_bytes() == second_png.read_bytes()
    assert first_svg.read_bytes() == second_svg.read_bytes()


def test_two_builds_of_the_same_figure_are_byte_identical(real, best, tmp_path):
    # Rebuilding from the same data must also match, so the jitter and the SVG id
    # salt cannot drift between runs.
    first, _ = save_figure(pose_strip(real.poses, "1iep"), tmp_path / "a")
    second, _ = save_figure(pose_strip(real.poses, "1iep"), tmp_path / "b")
    assert first.read_bytes() == second.read_bytes()


def test_svg_carries_no_timestamp(best, tmp_path):
    _, svg = save_figure(affinity_heatmap(best), tmp_path / "heatmap")
    text = svg.read_text(encoding="utf-8")
    assert "<dc:date>" not in text
    assert "Date" not in text
