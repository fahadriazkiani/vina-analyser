"""Publication-ready matplotlib figures for parsed and analysed docking results.

Design rules:

* Agg backend, so figures build with no display and no Streamlit involvement.
* Colour-blind-safe (Okabe-Ito) colours, single-column width, small readable fonts.
* Deterministic output: a fixed SVG id salt, a fixed jitter seed and no timestamps,
  so saving the same figure twice gives byte-identical files and tests can rely on it.
* Nothing is compared across targets. The heatmap says so in its caption, the bar
  chart covers one target, and the strip plot has one panel per target.
"""

import functools
from pathlib import Path

MINUS = "−"  # typographic minus, matching matplotlib's own axis ticks

import matplotlib

matplotlib.use("Agg")  # before pyplot: no display is ever needed

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vina_analyser.analysis import NOTES  # noqa: E402

# A fixed salt makes the element ids inside an SVG reproducible between runs.
matplotlib.rcParams["svg.hashsalt"] = "vina-analyser"

FIGURE_WIDTH_CM = 8.5  # journal single column
CM_PER_INCH = 2.54
MIN_POSES_FOR_BOX = 5
LABEL_LIMIT = 14
JITTER_SEED = 0
DPI = 300

# Okabe-Ito: distinguishable with every common form of colour vision deficiency.
BAR_COLOUR = "#0072B2"  # blue
CONTROL_COLOUR = "#E69F00"  # orange
POINT_COLOUR = "#0072B2"
BOX_COLOUR = "#999999"
MISSING_COLOUR = "#F2F2F2"
# Blues reversed: the most negative (most favourable) score is the darkest.
HEATMAP_CMAP = "Blues_r"

COLOURBAR_LABEL = "Vina score (kcal/mol), more negative = more favourable"
SCORE_LABEL = "Vina score (kcal/mol)"
FILLED_LABEL = "pose in output PDBQT (has coordinates)"
HOLLOW_LABEL = "log only: no coordinates (outside energy_range)"

RC = {
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",  # keep text as text in the SVG, not outlines
    "svg.hashsalt": "vina-analyser",
}


def _styled(function):
    """Build the whole figure under our rcParams.

    Text objects take their size when they are created, so labels, titles and the
    colour bar all have to be made inside the context, not just the axes.
    """

    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        with plt.rc_context(RC):
            figure = function(*args, **kwargs)
            _settle_layout(figure)
            return figure

    return wrapper


MAX_LAYOUT_PASSES = 10


def _settle_layout(figure) -> None:
    """Run constrained layout until the axes stop moving, then freeze it.

    Each save runs one more layout pass starting from where the last one left off.
    The heatmap with its horizontal colour bar needs about three passes to settle, so
    without this, saving the same figure twice gave slightly different files. Passes
    are run at the export resolution, because text extents depend on it.
    """
    import io

    def positions():
        return [tuple(axes.get_position().bounds) for axes in figure.axes]

    before = positions()
    for _ in range(MAX_LAYOUT_PASSES):
        figure.savefig(io.BytesIO(), format="png", dpi=DPI)
        after = positions()
        if after == before:
            break
        before = after
    figure.set_layout_engine("none")


def _new_figure(height_cm: float):
    """A single-column figure with small fonts and no overlapping labels."""
    size = (FIGURE_WIDTH_CM / CM_PER_INCH, height_cm / CM_PER_INCH)
    return plt.subplots(figsize=size, layout="constrained")


def short_labels(names, limit: int = LABEL_LIMIT) -> list[str]:
    """Shorten tick labels to `limit` characters, keeping distinct names distinct.

    Truncation is widened until no two different names share a label, so a figure can
    never show two ligands under one name.
    """
    names = list(names)
    longest = max((len(name) for name in names), default=0)
    width = limit
    labels = [_truncate(name, width) for name in names]
    while _collides(labels, names) and width < longest:
        width += 1
        labels = [_truncate(name, width) for name in names]
    if _collides(labels, names):  # unreachable for distinct names; kept as a guard
        labels = _with_suffixes(labels, names)
    return labels


def _truncate(name: str, width: int) -> str:
    return name if len(name) <= width else name[: width - 1] + "…"


def _collides(labels: list[str], names: list[str]) -> bool:
    """True when one label stands for two different names (identical names are fine)."""
    owner: dict[str, str] = {}
    for label, name in zip(labels, names):
        if owner.setdefault(label, name) != name:
            return True
    return False


def _with_suffixes(labels: list[str], names: list[str]) -> list[str]:
    counts: dict[str, int] = {}
    result = []
    for label, name in zip(labels, names):
        counts[label] = counts.get(label, 0) + 1
        result.append(label if counts[label] == 1 else f"{label[:-2]}~{counts[label]}")
    return result


# -------------------------------------------------------------------- heatmap


@_styled
def affinity_heatmap(best: pd.DataFrame):
    """Best score per target and ligand, on one shared colour scale in kcal/mol."""
    grid = best.pivot_table(index="target", columns="ligand", values="affinity", aggfunc="min")
    grid = grid.sort_index().sort_index(axis=1)
    if grid.empty:
        # An empty grid has zero-width axes, which matplotlib can only warn about.
        raise ValueError("The best-pose table has no scores, so there is nothing to plot.")
    values = np.ma.masked_invalid(grid.to_numpy(dtype=float))
    # Where a cell is the best of several seeds, say so: the score alone would hide
    # that one ligand was sampled more than the others.
    runs = None
    if "n_runs" in best.columns:
        runs = (
            best.pivot_table(index="target", columns="ligand", values="n_runs", aggfunc="max")
            .reindex(index=grid.index, columns=grid.columns)
        )

    height = 3.4 + 0.7 * len(grid.index)
    figure, axes = _new_figure(height)
    # A pair that was never docked shows as a blank cell, never as a zero.
    colours = plt.get_cmap(HEATMAP_CMAP).with_extremes(bad=MISSING_COLOUR)
    image = axes.imshow(values, cmap=colours, aspect="auto")
    # Explicit limits: each cell spans one unit around its index, whatever the grid
    # shape. One target (one row) is the common case, so it must not rely on autoscaling.
    n_rows, n_columns = values.shape
    axes.set_xlim(-0.5, n_columns - 0.5)
    axes.set_ylim(n_rows - 0.5, -0.5)  # first target at the top
    low, high = image.get_clim()
    if low == high:  # one score, or several identical ones: give the colour scale a width
        image.set_clim(low - 0.5, high + 0.5)

    axes.set_xticks(range(len(grid.columns)), short_labels(grid.columns), rotation=45, ha="right")
    axes.set_yticks(range(len(grid.index)), list(grid.index))
    for spine in axes.spines.values():
        spine.set_visible(False)
    axes.tick_params(length=0)

    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            if values.mask.ndim and values.mask[row, column]:
                continue
            value = float(values[row, column])
            label = format_score(value)
            if runs is not None:
                count = runs.iat[row, column]
                if pd.notna(count) and int(count) > 1:
                    label += f" (n={int(count)})"
            axes.text(
                column,
                row,
                label,
                ha="center",
                va="center",
                fontsize=6,
                color=_readable_on(image.cmap(image.norm(value))),
            )

    # Horizontal: the label is long, and upright it was clipped at the figure edge.
    bar = figure.colorbar(image, ax=axes, orientation="horizontal", fraction=0.12, pad=0.06)
    bar.set_label(COLOURBAR_LABEL, fontsize=6)
    bar.ax.tick_params(labelsize=6)
    figure.set_label(heatmap_caption())
    return figure


def format_score(value: float, decimals: int = 2) -> str:
    """Format a score with the Unicode minus sign the axis ticks already use."""
    return f"{value:.{decimals}f}".replace("-", MINUS)


def _tick_step(span: float) -> float:
    """A whole-number tick step, so every ranked chart is spaced alike."""
    for step in (1, 2, 4, 5, 10, 20, 50):
        if span / step <= 5:
            return step
    return 100


def _readable_on(rgba) -> str:
    luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
    return "white" if luminance < 0.5 else "black"


def heatmap_caption() -> str:
    return (
        "Best Vina score for each target and ligand, in kcal/mol, on one shared colour "
        "scale. Compare ligands within a single row (one target) only: different "
        "targets use different receptors and grid boxes, so their scores are not "
        "comparable. A cell marked (n=…) shows the best score across that many repeat "
        "runs with different seeds, so it had more sampling than an unmarked cell, "
        "which comes from a single run. Blank cells were not docked. " + NOTES["scoring_error"]
    )


# ---------------------------------------------------------------- ranked bars


@_styled
def ranked_bars(best: pd.DataFrame, target: str, control: str | None = None):
    """Ligands ranked by best score for one target, with the control highlighted."""
    rows = best[best["target"] == target]
    if rows.empty:
        raise ValueError(f"No results for target {target!r}; nothing to plot.")
    rows = rows.sort_values(["affinity", "ligand"], na_position="last")

    height = 2.0 + 0.5 * len(rows)
    figure, axes = _new_figure(height)
    positions = range(len(rows))
    scores = rows["affinity"].astype(float).to_list()
    ligands = rows["ligand"].to_list()
    highlighted = control is not None and control in ligands

    for position, (ligand, score) in enumerate(zip(ligands, scores)):
        is_control = highlighted and ligand == control
        colour = CONTROL_COLOUR if is_control else BAR_COLOUR
        axes.barh(
            position,
            score,
            color=colour,
            label="control compound" if is_control else None,
            height=0.7,
        )
        # Just inside the tip of the bar: outside it, the minus sign was cut off.
        # The text colour follows the bar, since the control's orange is light.
        axes.text(
            score + abs(score) * 0.02,
            position,
            format_score(score),
            va="center",
            ha="left",
            fontsize=6,
            color=_readable_on(matplotlib.colors.to_rgba(colour)),
        )

    axes.set_yticks(list(positions), short_labels(ligands))
    axes.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(_tick_step(abs(min(scores)))))
    axes.set_xlabel(SCORE_LABEL)
    axes.set_title(f"Target {target}", loc="left")
    axes.invert_yaxis()  # best score at the top
    axes.margins(x=0.02)
    if highlighted:
        axes.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0))
    figure.set_label(ranked_bars_caption(target))
    return figure


def ranked_bars_caption(target: str) -> str:
    return (
        f"Ligands ranked by their best Vina score against target {target}, most "
        "favourable at the top. Ligands with equal scores share a rank and are ordered "
        "by name. " + NOTES["scoring_error"]
    )


# ----------------------------------------------------------------- pose strip


@_styled
def pose_strip(poses: pd.DataFrame, target: str):
    """Every pose score for one target, one column per ligand.

    One panel per target on purpose: pooling a ligand's poses across targets would mix
    different receptors and grid boxes, which is exactly what the rest of the package
    refuses to do.
    """
    rows = poses[poses["target"] == target]
    if rows.empty:
        raise ValueError(f"No poses for target {target!r}; nothing to plot.")
    ligands = sorted(rows["ligand"].unique())
    index_of = {ligand: position for position, ligand in enumerate(ligands)}

    height = 6.4
    figure, axes = _new_figure(height)
    generator = np.random.default_rng(JITTER_SEED)  # reseeded per call, so reproducible

    filled_x, filled_y, hollow_x, hollow_y = [], [], [], []
    for row in rows.sort_values(["ligand", "pose"]).itertuples(index=False):
        if pd.isna(row.affinity):
            continue
        position = index_of[row.ligand] + float(generator.uniform(-0.18, 0.18))
        if row.has_coordinates:
            filled_x.append(position)
            filled_y.append(row.affinity)
        else:
            hollow_x.append(position)
            hollow_y.append(row.affinity)

    # Boxes first, so the markers sit on top of them.
    for ligand in ligands:
        with_coordinates = rows[(rows["ligand"] == ligand) & rows["has_coordinates"]]
        scores = with_coordinates["affinity"].dropna().astype(float).to_list()
        if len(scores) < MIN_POSES_FOR_BOX:
            continue
        drawn = axes.boxplot(
            [scores],
            positions=[index_of[ligand]],
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 0.8},
            boxprops={"facecolor": "white", "edgecolor": BOX_COLOUR, "linewidth": 0.8},
            whiskerprops={"color": BOX_COLOUR, "linewidth": 0.8},
            capprops={"color": BOX_COLOUR, "linewidth": 0.8},
            zorder=1,
        )
        for box in drawn["boxes"]:
            box.set_gid(f"box:{ligand}")

    axes.scatter(
        filled_x,
        filled_y,
        s=9,
        color=POINT_COLOUR,
        alpha=0.85,
        linewidths=0,
        label=FILLED_LABEL,
        zorder=3,
    )
    axes.scatter(
        hollow_x,
        hollow_y,
        s=11,
        facecolors="none",
        edgecolors=POINT_COLOUR,
        linewidths=0.6,
        label=HOLLOW_LABEL,
        zorder=2,
    )

    axes.set_xticks(range(len(ligands)), short_labels(ligands), rotation=45, ha="right")
    axes.set_ylabel(SCORE_LABEL)
    axes.invert_yaxis()  # more favourable (more negative) at the top, as in the bars
    axes.set_xlim(-0.6, len(ligands) - 0.4)

    handles = [handle for handle in axes.collections if len(handle.get_offsets())]
    if handles:
        axes.legend(handles=handles, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0))
    figure.set_label(pose_strip_caption(target))
    return figure


def pose_strip_caption(target: str) -> str:
    return (
        f"Every pose score for target {target}, one column per ligand, most "
        "favourable at the top. Repeat runs of "
        "the same ligand with different seeds share a column. Filled markers "
        "are poses present in the output PDBQT. Hollow markers are modes that appear "
        "only in the log, outside Vina's energy_range, so they have no coordinates and "
        "their scores are the log's rounded values. A box (median and quartiles) is "
        f"drawn only for ligands with at least {MIN_POSES_FOR_BOX} poses with "
        "coordinates in this panel, and summarises those poses alone; otherwise the "
        "points are shown without a box. " + NOTES["scoring_error"]
    )


# -------------------------------------------------------------------- saving


# Dropping Software and Date keeps repeated saves byte-identical.
_SAVE_OPTIONS = {
    "png": {"format": "png", "dpi": DPI, "metadata": {"Software": None}},
    "svg": {"format": "svg", "metadata": {"Date": None}},
}


def figure_bytes(figure, fmt: str) -> bytes:
    """The figure as PNG or SVG bytes, identical to what save_figure writes.

    The app shows and offers exactly these bytes, so what is on screen is what is
    downloaded.
    """
    import io

    buffer = io.BytesIO()
    figure.savefig(buffer, **_SAVE_OPTIONS[fmt])
    return buffer.getvalue()


def save_figure(figure, stem) -> tuple[Path, Path]:
    """Write `<stem>.png` at 300 dpi and `<stem>.svg`, both free of timestamps."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    svg = stem.with_suffix(".svg")
    # Written as bytes: given a path, matplotlib writes SVG in text mode, which turns
    # line endings into CRLF on Windows and makes the file differ between systems.
    png.write_bytes(figure_bytes(figure, "png"))
    svg.write_bytes(figure_bytes(figure, "svg"))
    return png, svg
