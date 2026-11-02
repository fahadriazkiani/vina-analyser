"""Shared fixtures for the test suite."""

import hashlib
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).resolve().parent.parent / "vina_test_data"


def _hash_all_data_files():
    return {
        path.relative_to(DATA_DIR).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(DATA_DIR.rglob("*"))
        if path.is_file()
    }


@pytest.fixture(scope="session", autouse=True)
def raw_data_is_unchanged():
    """Fail the run if any test modifies the raw data (working rule 4)."""
    before = _hash_all_data_files()
    yield
    assert _hash_all_data_files() == before, "vina_test_data/ was modified during the tests"
