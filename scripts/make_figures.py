"""Write the example figures into figures/ so they can be looked at.

Run from the project root:  python scripts/make_figures.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vina_analyser.analysis import analyse  # noqa: E402
from vina_analyser.parsers import load_folder  # noqa: E402
from vina_analyser.plots import (  # noqa: E402
    affinity_heatmap,
    heatmap_caption,
    pose_strip,
    pose_strip_caption,
    ranked_bars,
    ranked_bars_caption,
    save_figure,
)

CONTROL = "imatinib_1iep"


def main(out_dir=ROOT / "figures"):
    results = load_folder(ROOT / "vina_test_data" / "results", ROOT / "vina_test_data" / "logs")
    best = analyse(results, control=CONTROL).best
    written, captions = [], []

    written += save_figure(affinity_heatmap(best), out_dir / "heatmap_best_affinity")
    captions.append(("heatmap_best_affinity", heatmap_caption()))

    for target in sorted(best["target"].unique()):
        written += save_figure(
            ranked_bars(best, target, control=CONTROL), out_dir / f"ranked_{target}"
        )
        captions.append((f"ranked_{target}", ranked_bars_caption(target)))
        written += save_figure(pose_strip(results.poses, target), out_dir / f"poses_{target}")
        captions.append((f"poses_{target}", pose_strip_caption(target)))

    caption_file = out_dir / "captions.txt"
    caption_file.write_text(
        "\n\n".join(f"{name}\n{'-' * len(name)}\n{text}" for name, text in captions),
        encoding="utf-8",
    )
    for path in [*written, caption_file]:
        print(path)


if __name__ == "__main__":
    main()
