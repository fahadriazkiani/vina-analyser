"""Deliberate bugs for scripts/mutation_runner.py, grouped by the file they change.

Each entry is Mutation(phase, name, edits). `edits` is a list of (old, new) pairs; each
`old` must occur exactly once in the current file, or the runner reports NOT APPLIED.

These are a chosen subset per module, by design: each mutation plants one bug in a
behaviour the project has decided matters (see "Decisions made" in PROGRESS.md), and
asks whether a specific test notices. They are not an exhaustive mutation of every line.

The "guard" phase holds broken mutants that the runner must report as INVALID. If one
of those ever comes out CAUGHT, the guard itself is broken.
"""

from typing import NamedTuple


class Mutation(NamedTuple):
    phase: str
    name: str
    edits: list[tuple[str, str]]


MUTATIONS: dict[str, list[Mutation]] = {
    'vina_analyser/analysis.py': [
        Mutation('guard', 'Broken indentation (must be INVALID)', [
            ('RMSD_SUCCESS_CUTOFF = 2.0\n',
             'RMSD_SUCCESS_CUTOFF = 2.0\n   oops = 1\n'),
        ]),
        Mutation('phase2', 'Rank across targets, not within', [
            ('rounded.groupby(table["target"]).rank(',
             'rounded.rank('),
        ]),
        Mutation('phase2', "Delta against the control's best score on ANY target", [
            ('        control_row = group[group["ligand"] == control]\n',
             '        control_row = table[table["ligand"] == control].sort_values("affinity")\n'),
        ]),
        Mutation('phase2', 'RMSD computed while is_experimental_pose is False', [
            ('    if not is_experimental_pose:\n',
             '    if False:\n'),
        ]),
        Mutation('phase2', 'redocking_rmsd ignores the target argument', [
            ('        poses = poses[poses["target"] == target]\n',
             ''),
        ]),
        Mutation('phase2', 'SDF reference accepted', [
            ('    if filename.lower().endswith((".sdf", ".mol", ".mol2")) or _SDF_MARKER.search(text):',
             '    if False:'),
        ]),
        Mutation('phase2', 'Blank heavy_atoms filled with the median for LE', [
            ('    heavy = heavy.where(heavy > 0)\n',
             '    heavy = heavy.where(heavy > 0)\n    heavy = heavy.fillna(heavy.median())\n'),
        ]),
        Mutation('phase2', 'LE sign flipped', [
            ('(-affinity / heavy)',
             '(affinity / heavy)'),
        ]),
        Mutation('phase2', 'Best pose may be a log-only mode', [
            ('    with_coordinates = poses[poses["has_coordinates"]]\n',
             '    with_coordinates = poses\n'),
        ]),
        Mutation('phase2', 'Tied scores get different ranks', [
            ('.rank(method="min"',
             '.rank(method="first"'),
        ]),
        Mutation('phase2', 'SD reported with only two runs', [
            ('if len(scores) >= 3 else pd.NA',
             'if len(scores) >= 2 else pd.NA'),
        ]),
        Mutation('phase2', 'Delta direction flipped', [
            ('            - control_row["affinity"]\n        )\n',
             '            - control_row["affinity"]\n        ) * -1\n'),
        ]),
        Mutation('phase2', 'Unequal run counts not warned about', [
            ('        if unequal:\n',
             '        if False:\n'),
        ]),
        Mutation('phase2', 'Best run = lowest seed, not best score', [
            ('best.sort_values(["target", "ligand", "affinity", "seed"], na_position="last")',
             'best.sort_values(["target", "ligand", "seed"], na_position="last")'),
        ]),
        Mutation('phase2', 'Implausible RMSD not warned about', [
            ('if value > IMPLAUSIBLE_RMSD and row.pose == 1:',
             'if False:'),
        ]),
        Mutation('phase2', 'Element mismatch accepted', [
            ('        if pose_element != reference_element:\n',
             '        if False:\n'),
        ]),
        Mutation('phase2', 'RMSD includes hydrogens', [
            ('    pose_heavy = heavy_atoms_of(pose_atoms)\n    reference_heavy = heavy_atoms_of(reference_atoms)\n',
             '    pose_heavy = list(pose_atoms)\n    reference_heavy = list(reference_atoms)\n'),
        ]),
        Mutation('phase2', 'Different grid box still compared with the control', [
            ('            if mismatched:\n',
             '            if False:\n'),
        ]),
        Mutation('phase4', 'RMSD target choice ignored', [
            ('        poses = poses[poses["target"] == target]\n',
             ''),
        ]),
    ],
    'vina_analyser/plots.py': [
        Mutation('guard', 'Import of a missing name (must be INVALID)', [
            ('from vina_analyser.analysis import NOTES',
             'from vina_analyser.analysis import NOTES, DOES_NOT_EXIST'),
        ]),
        Mutation('phase3', 'Missing pair drawn as a zero score', [
            ('values = np.ma.masked_invalid(grid.to_numpy(dtype=float))',
             'values = np.ma.masked_invalid(grid.fillna(0).to_numpy(dtype=float))'),
        ]),
        Mutation('phase3', 'Heatmap colours reversed (favourable = light)', [
            ('HEATMAP_CMAP = "Blues_r"',
             'HEATMAP_CMAP = "Blues"'),
        ]),
        Mutation('phase3', '(n=k) mark on single-run cells', [
            ('if pd.notna(count) and int(count) > 1:',
             'if pd.notna(count) and int(count) >= 1:'),
        ]),
        Mutation('phase3', 'Colour bar label loses units and direction', [
            ('COLOURBAR_LABEL = "Vina score (kcal/mol), more negative = more favourable"',
             'COLOURBAR_LABEL = "Vina score"'),
        ]),
        Mutation('phase3', 'Heatmap caption drops the one-row rule', [
            ('"scale. Compare ligands within a single row (one target) only: different "',
             '"scale. different "'),
        ]),
        Mutation('phase3', 'Control bar not highlighted', [
            ('colour = CONTROL_COLOUR if is_control else BAR_COLOUR',
             'colour = BAR_COLOUR'),
        ]),
        Mutation('phase3', 'Bars ordered worst first', [
            ('    axes.invert_yaxis()  # best score at the top\n',
             '    pass  # best score at the top\n'),
        ]),
        Mutation('phase3', 'Box drawn for fewer than 5 poses', [
            ('if len(scores) < MIN_POSES_FOR_BOX:',
             'if len(scores) < 1:'),
        ]),
        Mutation('phase3', 'Log-only poses counted towards the box', [
            ('with_coordinates = rows[(rows["ligand"] == ligand) & rows["has_coordinates"]]',
             'with_coordinates = rows[rows["ligand"] == ligand]'),
        ]),
        Mutation('phase3', 'Strip plot pools poses across targets', [
            ('    rows = poses[poses["target"] == target]\n    if rows.empty:\n        raise ValueError(f"No poses',
             '    rows = poses if target in set(poses["target"]) else poses.iloc[0:0]\n    if rows.empty:\n        raise ValueError(f"No poses'),
        ]),
        Mutation('phase3', 'Hyphen instead of the Unicode minus', [
            ('return f"{value:.{decimals}f}".replace("-", MINUS)',
             'return f"{value:.{decimals}f}"'),
        ]),
        Mutation('phase3', 'Jitter not seeded', [
            ('np.random.default_rng(JITTER_SEED)',
             'np.random.default_rng()'),
        ]),
        Mutation('phase3', 'SVG keeps its timestamp', [
            ('"svg": {"format": "svg", "metadata": {"Date": None}},',
             '"svg": {"format": "svg"},'),
        ]),
        Mutation('phase3', 'PNG at 150 dpi', [
            ('DPI = 300',
             'DPI = 150'),
        ]),
        Mutation('phase3', 'Truncated labels may collide', [
            ('    while _collides(labels, names) and width < longest:\n',
             '    while False:\n'),
            ('    if _collides(labels, names):  # unreachable',
             '    if False:  # unreachable'),
        ]),
        Mutation('phase4', 'Empty heatmap not refused', [
            ('    if grid.empty:\n',
             '    if False:\n'),
        ]),
        Mutation('phase4', 'Layout not settled before saving', [
            ('            _settle_layout(figure)\n',
             ''),
        ]),
        Mutation('phase4', 'SVG written in text mode again', [
            ('    svg.write_bytes(figure_bytes(figure, "svg"))',
             '    figure.savefig(svg, **_SAVE_OPTIONS["svg"])'),
        ]),
    ],
    'vina_analyser/parsers.py': [
        Mutation('guard', 'Valid code that breaks almost everything (must be INVALID)', [
            ('    if data.startswith((b"\\xff\\xfe", b"\\xfe\\xff")):',
             '    return ""\n    if data.startswith((b"\\xff\\xfe", b"\\xfe\\xff")):'),
        ]),
        Mutation('phase5', 'Filename seed reported as logged', [
            ('seed_source = "filename" if name_seed is not None else None',
             'seed_source = "log" if name_seed is not None else None'),
        ]),
    ],
    'app.py': [
        Mutation('phase5', "Methods: app doesn't pass its git commit", [
            ('methods_text(results.runs, redocking=redocking, git=cached_git_info())',
             'methods_text(results.runs, redocking=redocking)'),
        ]),
        Mutation('guard', 'App crashes on start-up (must be INVALID)', [
            ('import streamlit as st\n',
             'import streamlit as st\nst.undefined_call()\n'),
        ]),
        Mutation('phase4', 'RMSD computed with the box unticked', [
            ('    if not confirmed:\n        missing.append(',
             '    if False:\n        missing.append('),
        ]),
        Mutation('phase4', 'Confirmation checkbox ticked by default', [
            ('value=False, key=f"confirm_',
             'value=True, key=f"confirm_'),
        ]),
        Mutation('phase4', "App doesn't pass target to RMSD", [
            ('runs=results.runs, target=target)',
             'runs=results.runs)'),
        ]),
        Mutation('phase4', 'Control defaults to the first ligand', [
            ('[NO_CONTROL, *ligands]',
             '[*ligands, NO_CONTROL]'),
        ]),
        Mutation('phase4', 'Run enabled despite problems', [
            ('blocked = not data or bool(problems)',
             'blocked = not data'),
        ]),
        Mutation('phase4', 'Stale results shown after edits', [
            ('if analysed["key"] != key or analysed["names"] != names:',
             'if analysed["key"] != key:'),
        ]),
        Mutation('phase4', 'Warnings hidden (table only)', [
            ('            st.warning(where + issue.message)',
             '            pass'),
        ]),
        Mutation('phase5', 'Prepared bundle survives a control change', [
            ('    if prepared is not None and prepared["signature"] != signature:',
             '    if False:'),
        ]),
        Mutation('phase5', 'App metadata ignores the edited names', [
            ('            names=analysed["names"],',
             '            names={},'),
        ]),
    ],
    'vina_analyser/inputs.py': [
        Mutation('phase4', 'Names passed only for unrecognised rows', [
            ('        names[row.file] = {',
             '        if row.recognised != "yes": names[row.file] = {'),
        ]),
        Mutation('phase4', 'Seed from table dropped', [
            ('seed = None if pd.isna(row.filename_seed) else int(row.filename_seed)',
             'seed = None'),
        ]),
        Mutation('phase4', 'Blank names allowed through', [
            ('if not target or not ligand:',
             'if False:'),
        ]),
        Mutation('phase4', 'Zip path traversal allowed', [
            ('or ".." in path.split("/")',
             ''),
        ]),
        Mutation('phase4', 'Zip bomb: header size trusted', [
            ('content = handle.read(MAX_FILE_BYTES + 1)',
             'content = handle.read(1)'),
        ]),
        Mutation('phase4', 'Zip mixed with loose files accepted', [
            ('        if len(items) > 1:\n',
             '        if False:\n'),
        ]),
    ],
    'vina_analyser/report.py': [
        Mutation('phase5', 'Methods: preparation gap left out', [
            ('    sentences.append(PREPARATION_GAP)\n',
             ''),
        ]),
        Mutation('phase5', 'Methods: analyser commit left out', [
            ('return f"{credit} (commit {git[\'commit_short\']}{dirty})."',
             'return credit + "."'),
        ]),
        Mutation('phase5', 'Methods: uncommitted changes not mentioned', [
            ('dirty = ", with uncommitted changes" if git.get("dirty") else ""',
             'dirty = ""'),
        ]),
        Mutation('phase5', 'Methods: sentence starts with a numeral', [
            ('f"A total of {ligands} ligands were"',
             'f"{ligands} ligands were"'),
        ]),
        Mutation('phase4', 'Limitations box loses the rigid-receptor point', [
            ('"The receptor is rigid. Protein flexibility',
             '"The receptor. Protein flexibility'),
        ]),
        Mutation('phase4', 'CSV without BOM', [
            ('CSV_ENCODING = "utf-8-sig"',
             'CSV_ENCODING = "utf-8"'),
        ]),
        Mutation('phase5', 'Methods quotes the filename seed', [
            ('logged = runs[runs["seed_source"] == "log"] if "seed_source" in runs.columns else runs.iloc[0:0]',
             'logged = runs[runs["seed"].notna()]'),
        ]),
        Mutation('phase5', "Methods falls back to Vina's default exhaustiveness", [
            ('or f"exhaustiveness: {NOT_RECORDED}",',
             'or "an exhaustiveness of 8",'),
        ]),
        Mutation('phase5', 'Varying value: first one silently chosen', [
            ('        if len(distinct) == 1:\n',
             '        if len(distinct) >= 1:\n'),
        ]),
        Mutation('phase5', 'Full receptor path leaks into methods', [
            ('receptors = sorted({_basename(r) for r',
             'receptors = sorted({str(r) for r'),
        ]),
        Mutation('phase5', 'Redocking sentence without confirmation', [
            ('    if redocking is None or not redocking.confirmed:\n',
             '    if redocking is None:\n'),
        ]),
        Mutation('phase5', 'Redocking reports one run as the whole range', [
            ('    low, high = values.min(), values.max()\n',
             '    low = high = values.iloc[0]\n'),
        ]),
        Mutation('phase5', 'Box bullets missing for several targets', [
            ('        lines += [f"- {target}: {boxes[target]}" for target in targets]\n',
             ''),
        ]),
        Mutation('phase5', 'Excel formula guard removed', [
            ('        cell.data_type = "s"  # openpyxl',
             '        pass  # openpyxl'),
        ]),
        Mutation('phase5', 'CSV formula guard removed', [
            ('safe[column] = safe[column].map(lambda v: "\'" + v if _formula_like(v) else v)',
             'pass'),
        ]),
        Mutation('phase5', "'-' missing from the formula characters", [
            ('_FORMULA_START = ("=", "+", "-", "@", "\\t", "\\r")',
             '_FORMULA_START = ("=", "+", "@", "\\t", "\\r")'),
        ]),
        Mutation('phase5', 'LE loses its 4 decimals in Excel', [
            ('LABELS["ligand_efficiency"]: "0.0000",',
             'LABELS["ligand_efficiency"]: "0.000",'),
        ]),
        Mutation('phase5', "Workbook created time is 'now'", [
            ('    book.properties.created = when\n',
             '    book.properties.created = datetime.now()\n'),
        ]),
        Mutation('phase5', 'ZIP entries get the current time', [
            ('info = zipfile.ZipInfo(name, date_time=ZIP_TIMESTAMP)',
             'import time; info = zipfile.ZipInfo(name, date_time=time.localtime()[:6])'),
        ]),
        Mutation('phase5', 'ZIP entry names not checked', [
            ('            raise ValueError(f"Unsafe name for a zip entry: {name!r}")',
             '            pass'),
        ]),
        Mutation('phase5', 'Target names used raw in file names', [
            ('stem = f"{prefix}_{safe_name(target)}"',
             'stem = f"{prefix}_{target}"'),
        ]),
        Mutation('phase5', 'Hash of decoded text, not original bytes', [
            ('"sha256": hashlib.sha256(data).hexdigest(),',
             '"sha256": hashlib.sha256(data.decode("latin-1").strip().encode()).hexdigest(),'),
        ]),
        Mutation('phase5', 'Excluded files left out of the metadata', [
            ('    for name in sorted(files):\n        data = files[name]\n',
             '    for name in sorted(names):\n        data = files[name]\n'),
        ]),
        Mutation('phase5', 'Git dirty flag ignored', [
            ('"dirty": status != ""',
             '"dirty": False'),
        ]),
        Mutation('phase5', 'Git failure crashes metadata', [
            ('    except (OSError, subprocess.SubprocessError):\n        return None\n',
             '    except ZeroDivisionError:\n        return None\n'),
        ]),
        Mutation('phase5', "Metadata timestamp is 'now'", [
            ('analysed_at=analysed_at.isoformat(),',
             'analysed_at=datetime.now(timezone.utc).isoformat(),'),
        ]),
    ],
    '.streamlit/config.toml': [
        Mutation('phase4', 'Upload limit config drifts', [
            ('maxUploadSize = 50',
             'maxUploadSize = 200'),
        ]),
    ],
}
