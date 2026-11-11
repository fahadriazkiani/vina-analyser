"""The README, CITATION.cff and LICENSE must agree with the code they describe."""

import re
import tomllib
from pathlib import Path

import vina_analyser
from vina_analyser.report import CITATIONS, LIMITATIONS

ROOT = Path(__file__).resolve().parent.parent
README = (ROOT / "README.md").read_text(encoding="utf-8")
CITATION = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
AUTHOR = "Fahad Riaz"


def test_readme_carries_every_limitation_the_app_shows():
    # Same wording and layout as the limitations box in the app.
    for text in LIMITATIONS:
        head, _, rest = text.partition(". ")
        line = f"- **{head}.** {rest}" if rest else f"- {text}"
        assert line in README.splitlines(), head


def test_readme_cites_both_vina_papers():
    for citation in CITATIONS:
        doi = citation.rsplit("DOI ", 1)[1]
        assert doi in README


def test_readme_links_and_images_point_to_files_that_exist():
    targets = re.findall(r"\]\(([^)#\s]+)", README)
    local = [t for t in targets if not re.match(r"[a-z]+:", t)]
    assert local, "the README links at least the licence and the screenshot"
    for target in local:
        assert (ROOT / target).exists(), target


def test_readme_gives_the_run_command_and_python_version():
    assert "streamlit run app.py" in README
    assert "pip install -r requirements.txt" in README
    assert "3.11" in README


def test_version_is_the_same_everywhere():
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    version = pyproject["project"]["version"]
    assert vina_analyser.__version__ == version
    assert re.search(rf"^version: \"?{re.escape(version)}\"?$", CITATION, re.MULTILINE)


def test_author_name_is_the_same_everywhere():
    assert f"Copyright (c) 2026 {AUTHOR}\n" in (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert re.search(r"family-names: \"?Riaz\"?$", CITATION, re.MULTILINE)
    assert re.search(r"given-names: \"?Fahad\"?$", CITATION, re.MULTILINE)
    assert AUTHOR in README


def test_citation_file_names_the_licence_and_repository():
    assert re.search(r"^license: MIT$", CITATION, re.MULTILINE)
    assert "repository-code: \"https://github.com/fahadriazkiani/vina-analyser\"" in CITATION
    assert re.search(r"^cff-version: 1\.2\.0$", CITATION, re.MULTILINE)


CONCEPT_DOI = "10.5281/zenodo.24146185"  # all versions; resolves to the latest release
VERSION_DOI = "10.5281/zenodo.24146186"  # v0.1.0


def test_readme_shows_the_doi_and_ci_badges_at_the_top():
    top = README.split("\n## ", 1)[0]
    assert f"[![DOI](https://zenodo.org/badge/1443434015.svg)](https://doi.org/{CONCEPT_DOI})" in top
    assert "actions/workflows/tests.yml/badge.svg" in top


def test_readme_cites_the_tool_with_the_concept_doi_and_the_release_doi():
    section = README.split("## Citing this tool\n", 1)[1].split("\n## ", 1)[0]
    assert CONCEPT_DOI in section and "latest release" in section
    assert VERSION_DOI in section and "v0.1.0" in section


def test_citation_file_carries_the_concept_doi():
    assert re.search(rf"^doi: \"?{re.escape(CONCEPT_DOI)}\"?$", CITATION, re.MULTILINE)
