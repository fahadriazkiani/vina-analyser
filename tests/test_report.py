"""Tests for vina_analyser.report (Phase 4 part: CSV export of on-screen tables)."""

import codecs
import io

import pandas as pd

from vina_analyser.analysis import LABELS
from vina_analyser.report import csv_bytes, display_table


def test_csv_starts_with_the_utf8_byte_order_mark():
    # Excel on Windows only reads a CSV as UTF-8 when it starts with the BOM.
    data = csv_bytes(pd.DataFrame({"a": [1]}))
    assert data.startswith(codecs.BOM_UTF8)
    assert data.count(codecs.BOM_UTF8) == 1


def test_csv_round_trips_delta_angstrom_and_unicode_minus():
    frame = pd.DataFrame({"Δ score (kcal/mol)": ["−1.5"], "RMSD (Å)": [0.91]})
    data = csv_bytes(frame)
    assert "Δ".encode() in data and "Å".encode() in data and "−".encode() in data
    back = pd.read_csv(io.BytesIO(data), encoding="utf-8-sig")
    assert list(back.columns) == list(frame.columns)
    assert back.iloc[0, 0] == "−1.5"


def test_csv_has_no_index_column():
    data = csv_bytes(pd.DataFrame({"a": [1, 2]}))
    assert data.decode("utf-8-sig").splitlines()[0] == "a"


def test_display_table_uses_labels_with_units():
    frame = pd.DataFrame({"target": ["t"], "affinity": [-7.0], "delta_vs_control": [0.5]})
    shown = display_table(frame)
    assert list(shown.columns) == ["target", LABELS["affinity"], LABELS["delta_vs_control"]]
    assert list(frame.columns) == ["target", "affinity", "delta_vs_control"]  # input untouched
