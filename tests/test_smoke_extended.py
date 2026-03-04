"""Extended smoke tests — CLI tools not covered by test_smoke.py.

Tools tested here:
  - discover_experiment_files.py  (pure Python, imports lcms_analyzer)
  - cas_resolver.py               (pure Python + PubChem HTTP)
  - scheme_aligner.py             (requires RDKit)
  - prepare_reaction_scheme.py    (--no-polish path = pure Python)
  - ole_extractor.py              (requires olefile)
  - procedure_writer.py           (pure Python + LCMS parsing, CSV fallback)
  - eln_enrichment.py             (library module, import + match_csv_to_scheme)
  - reactant_heuristic.py         (CLI with cdxml/smiles subcommands)
"""

import glob
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest

# Paths -------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PYTHON = sys.executable  # same interpreter running pytest
TEST_DATA = os.environ.get(
    "CHEM_TEST_DATA",
    os.path.join(os.path.dirname(PROJECT_ROOT), "chem-test-data"),
)


def _run(args, **kwargs):
    """Run a subprocess, return CompletedProcess."""
    return subprocess.run(
        args,
        capture_output=True,
        text=True,
        cwd=PROJECT_ROOT,
        timeout=120,
        **kwargs,
    )


def _assert_valid_cdxml(path):
    """Assert file exists, is non-empty, and parses as valid XML with a CDXML root."""
    assert os.path.isfile(path), f"output file missing: {path}"
    assert os.path.getsize(path) > 0, f"output file is empty: {path}"
    tree = ET.parse(path)
    root = tree.getroot()
    assert root.tag == "CDXML", f"unexpected root tag: {root.tag}"


# =========================================================================
# discover_experiment_files.py — experiment file discovery
# =========================================================================

class TestDiscoverExperimentFiles:
    """Smoke tests for discover_experiment_files.py."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.input_dir = os.path.join(
            TEST_DATA, "procedurefilltest", "KL-7001-incomplete"
        )
        if not os.path.isdir(self.input_dir):
            pytest.skip("KL-7001-incomplete test data not found")

    def test_text_output_contains_file_categories(self):
        r = _run([
            PYTHON, "discover_experiment_files.py",
            "--input-dir", self.input_dir,
            "--experiment", "KL-7001-004",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        text = r.stdout
        # Text report should mention the main file categories
        assert "CSV" in text, "output should mention CSV files"
        assert "LCMS" in text, "output should mention LCMS files"
        assert "KL-7001-004" in text, "output should echo experiment name"

    def test_json_output_has_expected_structure(self, tmp_path):
        out = os.path.join(str(tmp_path), "discover.json")
        r = _run([
            PYTHON, "discover_experiment_files.py",
            "--input-dir", self.input_dir,
            "--experiment", "KL-7001-004",
            "--json",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert os.path.isfile(out)

        with open(out) as f:
            data = json.load(f)
        assert data["experiment"] == "KL-7001-004"
        files = data["files"]
        assert "csv" in files
        assert "lcms" in files
        assert "cdx" in files
        assert "rxn" in files
        assert "nmr" in files
        # KL-7001-004 has at least a CSV
        assert len(files["csv"]) >= 1, "should find at least one CSV file"

    def test_lcms_entries_have_categories(self, tmp_path):
        out = os.path.join(str(tmp_path), "discover.json")
        r = _run([
            PYTHON, "discover_experiment_files.py",
            "--input-dir", self.input_dir,
            "--experiment", "KL-7001-004",
            "--json",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        with open(out) as f:
            data = json.load(f)
        lcms = data["files"]["lcms"]
        if lcms:  # may be empty if LCMS PDFs aren't in expected location
            for entry in lcms:
                assert "path" in entry
                assert "category" in entry
                assert "sort_key" in entry
                assert entry["category"] in (
                    "tracking", "workup", "purification", "final", "other"
                )


# =========================================================================
# cas_resolver.py — PubChem CAS resolution (network-dependent)
# =========================================================================

class TestCasResolver:
    """Smoke tests for cas_resolver.py.

    Tests marked @pytest.mark.network require internet access to PubChem.
    Run with: pytest -m network   (to include network tests)
    Skip with: pytest -m "not network"
    """

    @pytest.mark.network
    def test_resolves_aspirin_cas(self):
        """Resolve aspirin (CAS 50-78-2)."""
        r = _run([PYTHON, "-m", "cdxml_toolkit.cas_resolver", "50-78-2", "--pretty"])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        data = json.loads(r.stdout)
        assert data["cas"] == "50-78-2"
        assert data["smiles"], "should return a non-empty SMILES"
        assert data["mw"] is not None, "should return molecular weight"
        assert data["formula"], "should return a molecular formula"

    @pytest.mark.network
    def test_batch_resolution(self):
        """Resolve two CAS numbers: aspirin (50-78-2) and caffeine (58-08-2)."""
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.cas_resolver",
            "50-78-2", "58-08-2",
            "--pretty",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        data = json.loads(r.stdout)
        assert isinstance(data, list)
        assert len(data) == 2
        # Both should have SMILES
        for entry in data:
            assert "smiles" in entry

    @pytest.mark.network
    def test_invalid_cas_returns_error(self):
        """Invalid CAS check digit should exit with code 1."""
        # 999-99-9 has wrong check digit (should be 5, not 9)
        r = _run([PYTHON, "-m", "cdxml_toolkit.cas_resolver", "999-99-9"])
        assert r.returncode == 1

    def test_no_args_exits_with_error(self):
        """Running with no CAS numbers should show usage error (no network needed)."""
        r = _run([PYTHON, "-m", "cdxml_toolkit.cas_resolver"])
        assert r.returncode != 0


# =========================================================================
# scheme_aligner.py — MCS-based structure alignment (RDKit-dependent)
# =========================================================================

_rdkit_available = False
try:
    from rdkit import Chem
    _rdkit_available = True
except ImportError:
    pass


class TestSchemeAligner:
    """Smoke tests for scheme_aligner.py (requires RDKit)."""

    @pytest.fixture(autouse=True)
    def setup(self):
        if not _rdkit_available:
            pytest.skip("RDKit not available")
        # Buchwald-output scheme.cdxml has a reaction scheme with <step>
        self.cdxml = os.path.join(TEST_DATA, "Buchwald-output scheme.cdxml")
        if not os.path.isfile(self.cdxml):
            pytest.skip("Buchwald-output scheme.cdxml not found")

    def test_aligns_and_produces_valid_cdxml(self, tmp_path):
        out = os.path.join(str(tmp_path), "aligned.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_aligner",
            self.cdxml,
            "-o", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        _assert_valid_cdxml(out)

        # Verify output has fragment elements (reaction structures preserved)
        tree = ET.parse(out)
        frags = tree.findall(".//fragment")
        assert len(frags) >= 2, "aligned output should contain multiple fragments"

    def test_stdout_reports_alignment_progress(self):
        r = _run([PYTHON, "-m", "cdxml_toolkit.scheme_aligner", self.cdxml])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        text = r.stdout
        # Should report fragments and step info
        assert "Fragment" in text or "Step" in text, (
            "output should report alignment progress"
        )
        assert "Output:" in text, "should print output path"


# =========================================================================
# prepare_reaction_scheme.py — no-polish path (pure Python only)
# =========================================================================

class TestPrepareReactionScheme:
    """Smoke tests for prepare_reaction_scheme.py.

    Only tests the --no-polish + --cdxml path, which calls reaction_cleanup.py
    (pure Python). The full pipeline requires ChemDraw COM.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        # Use a CDXML with a reaction scheme
        self.cdxml = os.path.join(TEST_DATA, "Buchwald-output scheme.cdxml")
        if not os.path.isfile(self.cdxml):
            pytest.skip("Buchwald-output scheme.cdxml not found")

    def test_no_polish_produces_valid_cdxml(self, tmp_path):
        out = os.path.join(str(tmp_path), "prepared.cdxml")
        r = _run([
            PYTHON, "prepare_reaction_scheme.py",
            "--cdxml", self.cdxml,
            "-o", out,
            "--no-polish",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        _assert_valid_cdxml(out)

    def test_json_output_reports_layout_step(self, tmp_path):
        out = os.path.join(str(tmp_path), "prepared.cdxml")
        r = _run([
            PYTHON, "prepare_reaction_scheme.py",
            "--cdxml", self.cdxml,
            "-o", out,
            "--no-polish",
            "--json",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        data = json.loads(r.stdout)
        assert "steps_completed" in data
        assert "layout" in data["steps_completed"]
        assert data.get("error") is None

    def test_missing_file_returns_error(self):
        r = _run([
            PYTHON, "prepare_reaction_scheme.py",
            "--cdxml", "nonexistent_file.cdxml",
        ])
        assert r.returncode != 0


# =========================================================================
# ole_extractor.py — OLE ChemDraw extraction from Office files
# =========================================================================

_olefile_available = False
try:
    import olefile
    _olefile_available = True
except ImportError:
    pass


class TestOleExtractor:
    """Smoke tests for ole_extractor.py (requires olefile package)."""

    @pytest.fixture(autouse=True)
    def setup(self):
        if not _olefile_available:
            pytest.skip("olefile package not available")
        self.pptx = os.path.join(TEST_DATA, "Sample-weeklyreport-slide.pptx")
        if not os.path.isfile(self.pptx):
            pytest.skip("Sample-weeklyreport-slide.pptx not found")

    def test_extracts_from_pptx(self, tmp_path):
        out_dir = os.path.join(str(tmp_path), "extracted")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.ole_extractor",
            self.pptx,
            "-o", out_dir,
            "--format", "cdx",  # CDX only — avoids COM dependency for conversion
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        text = r.stdout
        assert "OLE Extractor" in text

    def test_output_summary_mentions_chemdraw(self):
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.ole_extractor",
            self.pptx,
            "--format", "cdx",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        text = r.stdout
        # Should report either found objects or no objects
        assert "ChemDraw" in text, (
            "output should mention ChemDraw (either found or not found)"
        )

    def test_missing_file_returns_error(self):
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.ole_extractor",
            "nonexistent_file.pptx",
        ])
        assert r.returncode != 0


# =========================================================================
# procedure_writer.py — lab book entry assembler
# =========================================================================

KL7001_DIR = os.path.join(
    TEST_DATA, "procedurefilltest", "KL-7001-incomplete"
)


class TestProcedureWriter:
    """Smoke tests for procedure_writer.py.

    Uses KL-7001-004 experiment data (Mitsunobu). ChemScript is optional —
    mass resolution falls back to CSV MW values.
    """

    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.tmp = tmp_path
        if not os.path.isdir(KL7001_DIR):
            pytest.skip("KL-7001-incomplete test data not found")

    def test_runs_and_exits_zero(self):
        out = os.path.join(str(self.tmp), "procedure.txt")
        r = _run([
            PYTHON, "procedure_writer.py",
            "--input-dir", KL7001_DIR,
            "--experiment", "KL-7001-004",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert os.path.isfile(out)
        assert os.path.getsize(out) > 0

    def test_output_contains_procedure_section(self):
        out = os.path.join(str(self.tmp), "procedure.txt")
        r = _run([
            PYTHON, "procedure_writer.py",
            "--input-dir", KL7001_DIR,
            "--experiment", "KL-7001-004",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        with open(out) as f:
            text = f.read()
        assert "PROCEDURE" in text, "output should contain PROCEDURE section"

    def test_output_contains_characterization_section(self):
        out = os.path.join(str(self.tmp), "procedure.txt")
        r = _run([
            PYTHON, "procedure_writer.py",
            "--input-dir", KL7001_DIR,
            "--experiment", "KL-7001-004",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        with open(out) as f:
            text = f.read()
        assert "CHARACTERIZATION" in text, (
            "output should contain CHARACTERIZATION section"
        )

    def test_output_contains_notes_section(self):
        out = os.path.join(str(self.tmp), "procedure.txt")
        r = _run([
            PYTHON, "procedure_writer.py",
            "--input-dir", KL7001_DIR,
            "--experiment", "KL-7001-004",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        with open(out) as f:
            text = f.read()
        assert "NOTES" in text, "output should contain NOTES section"

    def test_output_contains_mass_spec_data(self):
        out = os.path.join(str(self.tmp), "procedure.txt")
        r = _run([
            PYTHON, "procedure_writer.py",
            "--input-dir", KL7001_DIR,
            "--experiment", "KL-7001-004",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        with open(out) as f:
            text = f.read()
        # LCMS data uses adduct notation: [M+H]+, [M-H]-, [M+Na]+, etc.
        assert "[M" in text, (
            "output should contain mass adduct data (e.g. [M+H]+) from LCMS"
        )

    def test_missing_experiment_produces_placeholders(self):
        """Nonexistent experiment exits 0 but produces placeholder text."""
        r = _run([
            PYTHON, "procedure_writer.py",
            "--input-dir", KL7001_DIR,
            "--experiment", "KL-NONEXISTENT-999",
        ])
        assert r.returncode == 0
        # Should still produce the three sections, but with placeholder text
        assert "PROCEDURE" in r.stdout
        assert "not available" in r.stdout.lower()


# =========================================================================
# eln_enrichment.py — ELN CSV → scheme annotation (library module)
# =========================================================================

class TestElnEnrichment:
    """Smoke tests for eln_enrichment.py.

    This is a library module (no CLI entry point). Tests verify:
      - Module imports successfully (pure Python top-level imports)
      - Public data structures can be instantiated
      - match_csv_to_scheme runs with real test data (requires reagent_db)
    """

    def test_module_imports(self):
        """Top-level imports (cdxml_utils, constants, text_formatting) work."""
        import cdxml_toolkit.eln_enrichment as eln_enrichment
        assert hasattr(eln_enrichment, "match_csv_to_scheme")
        assert hasattr(eln_enrichment, "enrich_phase_a")
        assert hasattr(eln_enrichment, "enrich_phase_b")

    def test_data_structures_instantiate(self):
        """MatchedReagent and EnrichmentData dataclasses work."""
        from cdxml_toolkit.eln_enrichment import MatchedReagent, EnrichmentData

        mr = MatchedReagent(
            csv_name="Cs2CO3",
            csv_equiv="2.0",
            csv_mass="2.15 g",
            csv_is_substrate=False,
            csv_mw=325.82,
            scheme_element_id="100",
            scheme_position="above_arrow",
            scheme_display="Cs2CO3",
        )
        assert mr.csv_name == "Cs2CO3"

        ed = EnrichmentData()
        assert ed.matches == []
        assert ed.substrate is None

    def test_match_csv_to_scheme_with_real_data(self):
        """match_csv_to_scheme runs on KL-CC-001 data without crashing.

        Full MW matching requires ChemScript + RDKit; name matching via
        reagent_db should still work without them.
        """
        csv_path = os.path.join(TEST_DATA, "KL-CC-001", "KL-CC-001.csv")
        cdxml_path = os.path.join(
            TEST_DATA, "KL-CC-001", "KL-CC-001-refactor-test.cdxml"
        )
        if not os.path.isfile(csv_path):
            pytest.skip("KL-CC-001.csv not found")
        if not os.path.isfile(cdxml_path):
            pytest.skip("KL-CC-001-refactor-test.cdxml not found")

        from cdxml_toolkit.eln_enrichment import match_csv_to_scheme

        tree = ET.parse(cdxml_path)
        root = tree.getroot()
        result = match_csv_to_scheme(root, csv_path, verbose=False)

        # Should return an EnrichmentData (may have matches or not,
        # depending on whether ChemScript/RDKit are available)
        assert result is not None
        assert hasattr(result, "matches")
        assert isinstance(result.matches, list)


# =========================================================================
# reactant_heuristic.py — reagent classification (role lookup + RDKit MCS)
# =========================================================================

class TestReactantHeuristic:
    """Smoke tests for reactant_heuristic.py.

    The cdxml subcommand requires ChemScript for SMILES extraction and
    RDKit for MCS. The --help flag and smiles subcommand (with RDKit)
    are tested where possible.
    """

    def test_help_exits_zero(self):
        r = _run([PYTHON, "-m", "cdxml_toolkit.reactant_heuristic", "--help"])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert "cdxml" in r.stdout, "help should mention cdxml subcommand"
        assert "smiles" in r.stdout, "help should mention smiles subcommand"

    def test_cdxml_help_exits_zero(self):
        r = _run([PYTHON, "-m", "cdxml_toolkit.reactant_heuristic", "cdxml", "--help"])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert "--input" in r.stdout or "-i" in r.stdout

    def test_smiles_help_exits_zero(self):
        r = _run([PYTHON, "-m", "cdxml_toolkit.reactant_heuristic", "smiles", "--help"])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert "--reagents" in r.stdout
        assert "--product" in r.stdout

    @pytest.mark.skipif(not _rdkit_available, reason="RDKit not available")
    def test_smiles_mode_classifies_reagents(self):
        """Classify cyclohexanol + benzoic acid → ester."""
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reactant_heuristic", "smiles",
            "--reagents", "OC1CCCCC1", "OC(=O)c1ccccc1",
            "--product", "O=C(c1ccccc1)OC1CCCCC1",
            "--names", "cyclohexanol", "benzoic acid",
            "--pretty",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        data = json.loads(r.stdout)
        # Output is a dict with "reagents" list
        assert isinstance(data, dict)
        reagents = data["reagents"]
        assert isinstance(reagents, list)
        assert len(reagents) == 2
        for entry in reagents:
            assert "name" in entry or "smiles" in entry
            assert "classification" in entry

    def test_cdxml_mode_with_real_file(self):
        """Run cdxml mode on Buchwald scheme (may degrade without ChemScript)."""
        cdxml = os.path.join(TEST_DATA, "Buchwald-output scheme.cdxml")
        if not os.path.isfile(cdxml):
            pytest.skip("Buchwald-output scheme.cdxml not found")

        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reactant_heuristic", "cdxml",
            "-i", cdxml,
            "--pretty",
        ])
        # Should succeed even without ChemScript (graceful degradation)
        assert r.returncode == 0, f"stderr: {r.stderr}"
        data = json.loads(r.stdout)
        # Output is a dict with "reagents" list
        assert isinstance(data, dict), "output should be a JSON dict"
        reagents = data["reagents"]
        assert isinstance(reagents, list)
        assert len(reagents) >= 1, "should classify at least one reagent"


# ===================================================================
# scheme_merger.py
# ===================================================================

class TestSchemeMerger:
    """Smoke tests for scheme_merger.py — auto-detect, parallel, and validation."""

    REF_FILE = os.path.join(
        TEST_DATA, "KL-7001-004-with-009-scheme.cdxml")

    def test_help(self):
        r = _run([PYTHON, "-m", "cdxml_toolkit.scheme_merger", "--help"])
        assert r.returncode == 0
        assert "--mode" in r.stdout

    def test_parse_scheme_library(self):
        """Import and parse the reference file via library API."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        from cdxml_toolkit.scheme_merger import parse_scheme
        ps = parse_scheme(self.REF_FILE)
        assert len(ps.fragments) == 3
        assert len(ps.run_arrows) == 2
        assert ps.run_arrow_data[0].sm_mass_text == "50.0 mg"
        assert ps.run_arrow_data[1].sm_mass_text == "2.00 g"

    def test_parallel_merge_self(self, tmp_path):
        """Parallel merge of the same file with itself produces valid output."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        out = str(tmp_path / "merged.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_merger", "--mode", "parallel",
            self.REF_FILE, self.REF_FILE,
            "-o", out, "-v",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert os.path.isfile(out)
        # Parse output and verify structure
        tree = ET.parse(out)
        root = tree.getroot()
        page = root.find(".//page")
        arrows = [el for el in page if el.tag == "arrow"]
        frags = [el for el in page if el.tag == "fragment"]
        assert len(frags) == 3, "should keep 3 fragments"
        assert len(arrows) == 5, "1 main + 4 run arrows (2 from each input)"

    def test_parallel_merge_no_equiv(self, tmp_path):
        """--no-equiv removes equiv labels."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        out = str(tmp_path / "merged-noeq.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_merger", "--mode", "parallel",
            self.REF_FILE, self.REF_FILE,
            "-o", out, "--no-equiv",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        # Verify standalone equiv labels are gone
        tree = ET.parse(out)
        root = tree.getroot()
        page = root.find(".//page")
        import re
        equiv_re = re.compile(r'^\s*\(\d+\.?\d*\s*eq\.\)\s*$')
        for t_el in page.iter("t"):
            for s in t_el.iter("s"):
                if s.text and equiv_re.match(s.text):
                    pytest.fail(f"Found equiv label after --no-equiv: {s.text}")

    def test_parallel_merge_equiv_range(self, tmp_path):
        """--equiv-range runs without error."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        out = str(tmp_path / "merged-range.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_merger", "--mode", "parallel",
            self.REF_FILE, self.REF_FILE,
            "-o", out, "--equiv-range",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

    def test_deprecated_parallel_flag(self, tmp_path):
        """Legacy --parallel flag still works with deprecation warning."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        out = str(tmp_path / "merged-legacy.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_merger", "--parallel",
            self.REF_FILE, self.REF_FILE,
            "-o", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert "deprecated" in r.stderr.lower()
        assert os.path.isfile(out)

    def test_auto_detect_parallel(self, tmp_path):
        """Auto-detect mode correctly identifies parallel merge (same file)."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        out = str(tmp_path / "auto-merged.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_merger",
            self.REF_FILE, self.REF_FILE,
            "-o", out, "-v",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert os.path.isfile(out)
        # Should produce same result as explicit parallel
        tree = ET.parse(out)
        page = tree.getroot().find(".//page")
        arrows = [el for el in page if el.tag == "arrow"]
        assert len(arrows) == 5, "auto-detect should produce parallel merge"

    def test_classify_pair_library(self):
        """classify_pair correctly identifies parallel schemes."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        from cdxml_toolkit.scheme_merger import parse_scheme, classify_pair
        ps = parse_scheme(self.REF_FILE)
        assert classify_pair(ps, ps) == "parallel"

    def test_auto_detect_library(self):
        """auto_detect produces valid MergePlan for identical inputs."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        from cdxml_toolkit.scheme_merger import parse_scheme, auto_detect
        ps = parse_scheme(self.REF_FILE)
        plan = auto_detect([ps, ps])
        assert len(plan.parallel_groups) == 1
        assert plan.parallel_groups[0] == [0, 1]
        assert len(plan.sequential_chain) == 1

    def test_parallel_strict_rejects_mismatch(self):
        """parallel_merge with strict=True rejects different reactions."""
        if not os.path.isfile(self.REF_FILE):
            pytest.skip("reference CDXML not found")
        from cdxml_toolkit.scheme_merger import parse_scheme, parallel_merge, ParsedScheme
        import copy
        ps1 = parse_scheme(self.REF_FILE)
        # Create a modified copy with different product SMILES
        ps2 = copy.deepcopy(ps1)
        for pid in ps2.product_ids:
            ps2.fragment_smiles[pid] = "CCCCCC"  # fake SMILES
        try:
            parallel_merge([ps1, ps2], strict=True)
            pytest.fail("Should have raised ValueError for mismatched products")
        except ValueError as e:
            assert "different products" in str(e).lower()


# ===================================================================
# query_status.py — pipeline status introspection
# ===================================================================

class TestQueryStatus:
    """Smoke tests for query_status.py."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.input_dir = os.path.join(
            TEST_DATA, "procedurefilltest", "KL-7001-incomplete"
        )
        if not os.path.isdir(self.input_dir):
            pytest.skip("KL-7001-incomplete test data not found")

    def test_list_mode_shows_experiments(self):
        r = _run([
            PYTHON, "query_status.py",
            "--input-dir", self.input_dir,
            "--list",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert "Experiments" in r.stdout
        assert "KL-7001-004" in r.stdout

    def test_list_mode_json(self):
        r = _run([
            PYTHON, "query_status.py",
            "--input-dir", self.input_dir,
            "--list", "--json",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        data = json.loads(r.stdout)
        assert "experiments" in data
        assert isinstance(data["experiments"], list)
        assert len(data["experiments"]) >= 1

    def test_experiment_detail_json(self, tmp_path):
        out = os.path.join(str(tmp_path), "status.json")
        r = _run([
            PYTHON, "query_status.py",
            "--input-dir", self.input_dir,
            "--experiments", "KL-7001-004",
            "--json",
            "--output", out,
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert os.path.isfile(out)

        with open(out) as f:
            data = json.load(f)
        assert isinstance(data, list)
        assert len(data) == 1
        exp = data[0]
        assert exp["experiment"] == "KL-7001-004"
        assert "input" in exp
        assert "csv" in exp["input"]
        assert "lcms_count" in exp["input"]
        assert "nmr_count" in exp["input"]


# =========================================================================
# reaction_parser.py — unified reaction semantic layer
# =========================================================================

class TestReactionParser:
    """Smoke tests for reaction_parser.py.

    Tests CLI interface with real test data (KL-CC-001 Buchwald coupling).
    """

    CDXML_FILE = os.path.join(TEST_DATA, "KL-CC-001", "KL-CC-001.cdxml")
    CSV_FILE = os.path.join(TEST_DATA, "KL-CC-001", "KL-CC-001.csv")
    RXN_FILE = os.path.join(TEST_DATA, "KL-CC-001", "KL-CC-001.rxn")

    @pytest.fixture(autouse=True)
    def setup(self):
        if not os.path.isfile(self.CDXML_FILE):
            pytest.skip("KL-CC-001.cdxml not found")

    def test_help_exits_zero(self):
        r = _run([PYTHON, "-m", "cdxml_toolkit.reaction_parser", "--help"])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert "--csv" in r.stdout
        assert "--rxn" in r.stdout
        assert "--pretty" in r.stdout

    def test_cdxml_only(self, tmp_path):
        """Parse CDXML only — verify JSON output structure."""
        out = os.path.join(str(tmp_path), "reaction.json")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            self.CDXML_FILE,
            "-o", out,
            "--pretty",
            "--no-rxnmapper", "--no-rxn-insight", "--no-network",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert os.path.isfile(out)
        assert os.path.getsize(out) > 0

        with open(out) as f:
            data = json.load(f)
        assert data["version"] == "1.3"
        assert isinstance(data["species"], list)
        assert len(data["species"]) >= 1, "should find at least one species"
        # Every species should have required fields
        for sp in data["species"]:
            assert "id" in sp
            assert "name" in sp
            assert "role" in sp

    def test_cdxml_with_csv(self, tmp_path):
        """Parse CDXML + CSV — verify CSV fields populated."""
        if not os.path.isfile(self.CSV_FILE):
            pytest.skip("KL-CC-001.csv not found")
        out = os.path.join(str(tmp_path), "reaction.json")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            self.CDXML_FILE,
            "--csv", self.CSV_FILE,
            "-o", out,
            "--pretty",
            "--no-rxnmapper", "--no-rxn-insight", "--no-network",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        with open(out) as f:
            data = json.load(f)
        species = data["species"]
        assert len(species) >= 2, "Buchwald should have multiple species"

        # At least one species should have CSV matching data
        csv_matched = [sp for sp in species if sp.get("csv_name")]
        assert len(csv_matched) >= 1, "CSV matching should find at least one match"

    def test_cdxml_with_csv_and_rxn(self, tmp_path):
        """Parse all three inputs — verify SM and DP identification."""
        if not os.path.isfile(self.CSV_FILE):
            pytest.skip("KL-CC-001.csv not found")
        if not os.path.isfile(self.RXN_FILE):
            pytest.skip("KL-CC-001.rxn not found")
        out = os.path.join(str(tmp_path), "reaction.json")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            self.CDXML_FILE,
            "--csv", self.CSV_FILE,
            "--rxn", self.RXN_FILE,
            "-o", out,
            "--pretty",
            "--no-rxnmapper", "--no-rxn-insight", "--no-network",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        with open(out) as f:
            data = json.load(f)
        species = data["species"]

        # Should identify SM and DP
        sms = [sp for sp in species if sp.get("is_sm")]
        dps = [sp for sp in species if sp.get("is_dp")]
        assert len(sms) >= 1, "should identify at least one SM"
        assert len(dps) >= 1, "should identify at least one DP"

    def test_json_output_file(self, tmp_path):
        """Verify JSON file is well-formed and contains metadata."""
        out = os.path.join(str(tmp_path), "reaction.json")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            self.CDXML_FILE,
            "-o", out,
            "--no-rxnmapper", "--no-rxn-insight", "--no-network",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        with open(out) as f:
            data = json.load(f)
        # Metadata should be present
        assert "metadata" in data
        assert "parser_version" in data["metadata"]
        assert "timestamp" in data["metadata"]
        # Input files should be recorded
        assert "input_files" in data

    def test_stdout_output_when_no_output_flag(self):
        """Without -o, JSON should go to stdout."""
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            self.CDXML_FILE,
            "--no-rxnmapper", "--no-rxn-insight", "--no-network",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        data = json.loads(r.stdout)
        assert data["version"] == "1.3"

    def test_missing_file_returns_error(self):
        """Nonexistent input file should fail gracefully."""
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            "nonexistent_file.cdxml",
        ])
        assert r.returncode != 0

    @pytest.mark.skipif(not _rdkit_available, reason="RDKit not available")
    def test_species_have_masses(self, tmp_path):
        """Species with SMILES should have computed masses and adducts."""
        out = os.path.join(str(tmp_path), "reaction.json")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            self.CDXML_FILE,
            "-o", out,
            "--pretty",
            "--no-rxnmapper", "--no-rxn-insight", "--no-network",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        with open(out) as f:
            data = json.load(f)
        # Find species with SMILES — they should have masses
        with_smiles = [sp for sp in data["species"] if sp.get("smiles")]
        assert len(with_smiles) >= 1
        for sp in with_smiles:
            assert sp.get("exact_mass", 0) > 0, (
                f"species {sp['id']} has SMILES but no exact_mass"
            )
            assert sp.get("mw", 0) > 0, (
                f"species {sp['id']} has SMILES but no MW"
            )
            assert sp.get("adducts"), (
                f"species {sp['id']} has SMILES but no adducts"
            )


# =========================================================================
# scheme_maker.py — experimental CDXML scheme builder from reaction JSON
# =========================================================================

class TestSchemeMaker:
    """Smoke tests for scheme_maker.py (experimental).

    Requires RDKit for SMILES → 2D coord generation.
    Uses reaction_parser to create input JSON, then tests scheme_maker.
    """

    CDXML_FILE = os.path.join(TEST_DATA, "KL-CC-001", "KL-CC-001.cdxml")
    CSV_FILE = os.path.join(TEST_DATA, "KL-CC-001", "KL-CC-001.csv")

    @pytest.fixture(autouse=True)
    def setup(self):
        if not os.path.isfile(self.CDXML_FILE):
            pytest.skip("KL-CC-001.cdxml not found")
        if not _rdkit_available:
            pytest.skip("RDKit required for scheme_maker")

    def _make_json(self, tmp_path, use_csv=True):
        """Helper: run reaction_parser to create input JSON."""
        out = os.path.join(str(tmp_path), "reaction.json")
        args = [
            PYTHON, "-m", "cdxml_toolkit.reaction_parser",
            self.CDXML_FILE,
            "-o", out, "--pretty",
            "--no-rxnmapper", "--no-rxn-insight", "--no-network",
        ]
        if use_csv and os.path.isfile(self.CSV_FILE):
            args.extend(["--csv", self.CSV_FILE])
        r = _run(args)
        assert r.returncode == 0, f"reaction_parser failed: {r.stderr}"
        return out

    def test_help_exits_zero(self):
        r = _run([PYTHON, "-m", "cdxml_toolkit.scheme_maker", "--help"])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        assert "--approach" in r.stdout
        assert "--align-mode" in r.stdout
        assert "--no-run-arrow" in r.stdout

    def test_basic_scheme(self, tmp_path):
        """Build scheme from JSON — verify CDXML output is valid."""
        json_path = self._make_json(tmp_path)
        out = os.path.join(str(tmp_path), "scheme.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_maker",
            json_path, "-o", out,
            "--align-mode", "none",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        _assert_valid_cdxml(out)

    def test_scheme_has_fragments(self, tmp_path):
        """Output CDXML should contain fragment elements."""
        json_path = self._make_json(tmp_path)
        out = os.path.join(str(tmp_path), "scheme.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_maker",
            json_path, "-o", out,
            "--align-mode", "none",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        tree = ET.parse(out)
        frags = tree.findall(".//fragment")
        assert len(frags) >= 2, "should have at least reactant + product fragments"

    def test_scheme_has_arrow(self, tmp_path):
        """Output should have a reaction arrow."""
        json_path = self._make_json(tmp_path)
        out = os.path.join(str(tmp_path), "scheme.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_maker",
            json_path, "-o", out,
            "--align-mode", "none",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        tree = ET.parse(out)
        arrows = tree.findall(".//arrow")
        assert len(arrows) >= 1, "should have at least one arrow"

    def test_scheme_has_conditions_text(self, tmp_path):
        """Output should have condition text elements."""
        json_path = self._make_json(tmp_path, use_csv=True)
        out = os.path.join(str(tmp_path), "scheme.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_maker",
            json_path, "-o", out,
            "--align-mode", "none",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"

        tree = ET.parse(out)
        page = tree.find(".//page")
        assert page is not None
        # Find <t> elements directly under <page> (conditions text)
        texts = [el for el in page if el.tag == "t"]
        assert len(texts) >= 1, "should have conditions text"

    def test_no_run_arrow_flag(self, tmp_path):
        """--no-run-arrow should skip run arrow."""
        json_path = self._make_json(tmp_path)
        out = os.path.join(str(tmp_path), "scheme.cdxml")
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_maker",
            json_path, "-o", out,
            "--align-mode", "none", "--no-run-arrow",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        _assert_valid_cdxml(out)

    def test_missing_file_returns_error(self):
        """Nonexistent input file should fail gracefully."""
        r = _run([PYTHON, "-m", "cdxml_toolkit.scheme_maker", "nonexistent.json"])
        assert r.returncode != 0

    def test_default_output_path(self, tmp_path):
        """Without -o, output should be {stem}-scheme.cdxml."""
        json_path = self._make_json(tmp_path)
        r = _run([
            PYTHON, "-m", "cdxml_toolkit.scheme_maker",
            json_path,
            "--align-mode", "none",
        ])
        assert r.returncode == 0, f"stderr: {r.stderr}"
        expected = os.path.join(str(tmp_path), "reaction-scheme.cdxml")
        assert os.path.isfile(expected), (
            f"expected default output at {expected}"
        )
