"""Exports. Phase 4 part: CSV of the tables shown in the app.

Phase 5 adds Excel, the figure ZIP, the methods paragraph and run_metadata.json.
"""

import pandas as pd

from vina_analyser.analysis import LABELS

# UTF-8 with a byte-order mark: Excel on Windows reads a CSV without one in the local
# code page and turns Δ, Å and the Unicode minus into garbage.
CSV_ENCODING = "utf-8-sig"


def display_table(frame: pd.DataFrame) -> pd.DataFrame:
    """A copy with column names replaced by their labels (with units) where known."""
    return frame.rename(columns={column: LABELS[column] for column in frame.columns if column in LABELS})


def csv_bytes(frame: pd.DataFrame) -> bytes:
    """The table as CSV bytes in UTF-8 with a BOM, without the pandas index."""
    return frame.to_csv(index=False, lineterminator="\n").encode(CSV_ENCODING)
