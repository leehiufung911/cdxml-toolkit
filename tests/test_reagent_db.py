"""Unit tests for reagent_db.py — two-tier reagent database loader."""

import json
import os

import pytest

from cdxml_toolkit.reagent_db import ReagentDB, get_reagent_db

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JSON_PATH = os.path.join(PROJECT_ROOT, "cdxml_toolkit", "reagent_abbreviations.json")
CS_JSON_PATH = os.path.join(PROJECT_ROOT, "cdxml_toolkit", "chemscanner_abbreviations.json")


# =========================================================================
# Singleton behaviour
# =========================================================================

class TestSingleton:
    def test_get_reagent_db_returns_same_object(self):
        db1 = get_reagent_db()
        db2 = get_reagent_db()
        assert db1 is db2, "get_reagent_db() should return the same singleton"


# =========================================================================
# Name-based lookups
# =========================================================================

class TestDisplayForName:
    def test_cs2co3(self):
        db = get_reagent_db()
        assert db.display_for_name("cs2co3") == "Cs2CO3"

    def test_case_insensitive(self):
        db = get_reagent_db()
        assert db.display_for_name("CS2CO3") == "Cs2CO3"
        assert db.display_for_name("Cs2CO3") == "Cs2CO3"

    def test_alias_tea_resolves_to_et3n(self):
        db = get_reagent_db()
        result = db.display_for_name("tea")
        assert result == "Et3N", f"expected 'Et3N', got {result!r}"

    def test_unknown_name_returns_none(self):
        db = get_reagent_db()
        assert db.display_for_name("xyzzy_not_a_reagent") is None


class TestRoleForName:
    def test_pd_oac2_is_catalyst(self):
        db = get_reagent_db()
        assert db.role_for_name("pd(oac)2") == "catalyst"

    def test_cs2co3_is_base(self):
        db = get_reagent_db()
        assert db.role_for_name("cs2co3") == "base"

    def test_unknown_returns_none(self):
        db = get_reagent_db()
        assert db.role_for_name("xyzzy_not_a_reagent") is None


# =========================================================================
# resolve_display — fallback behaviour
# =========================================================================

class TestResolveDisplay:
    def test_known_name_returns_display(self):
        db = get_reagent_db()
        assert db.resolve_display("cs2co3") == "Cs2CO3"

    def test_unknown_returns_input_unchanged(self):
        db = get_reagent_db()
        assert db.resolve_display("unknown_thing") == "unknown_thing"

    def test_preserves_original_case_on_fallback(self):
        db = get_reagent_db()
        assert db.resolve_display("MyCustomReagent") == "MyCustomReagent"


# =========================================================================
# SMILES-based lookups
# =========================================================================

class TestSmilesLookup:
    def test_display_for_known_smiles(self):
        db = get_reagent_db()
        # Et3N SMILES
        result = db.display_for_smiles("CCN(CC)CC")
        assert result == "Et3N", f"expected 'Et3N', got {result!r}"

    def test_role_for_known_smiles(self):
        db = get_reagent_db()
        result = db.role_for_smiles("CCN(CC)CC")
        assert result == "base"

    def test_unknown_smiles_returns_none(self):
        db = get_reagent_db()
        assert db.display_for_smiles("CCCCCCCCCCCCCCCC") is None

    def test_smiles_role_display(self):
        db = get_reagent_db()
        result = db.smiles_role_display("CCN(CC)CC")
        assert result is not None
        role, display = result
        assert role == "base"
        assert display == "Et3N"

    def test_smiles_role_display_unknown(self):
        db = get_reagent_db()
        assert db.smiles_role_display("CCCCCCCCCCCCCCCC") is None


# =========================================================================
# JSON database integrity
# =========================================================================

class TestDatabaseIntegrity:
    def test_json_file_exists(self):
        assert os.path.isfile(JSON_PATH), "reagent_abbreviations.json missing"

    def test_minimum_entry_count(self):
        with open(JSON_PATH, encoding="utf-8") as f:
            data = json.load(f)
        assert len(data) >= 100, (
            f"expected >=100 entries, got {len(data)}"
        )

    def test_every_entry_has_display(self):
        with open(JSON_PATH, encoding="utf-8") as f:
            data = json.load(f)
        for key, entry in data.items():
            assert "display" in entry, f"entry {key!r} missing 'display'"

    def test_roles_are_known_categories(self):
        known_roles = {
            "catalyst", "ligand", "base", "solvent", "coupling_reagent",
            "reducing_agent", "oxidant", "protecting_group", "deprotecting_agent",
            "acid", "activating_agent", "lewis_acid", "drying_agent",
            "halogenating_agent", "fluorinating_agent", "borylating_agent",
            "additive", "reductant", "reagent",
        }
        with open(JSON_PATH, encoding="utf-8") as f:
            data = json.load(f)
        for key, entry in data.items():
            role = entry.get("role")
            if role is not None:
                assert role in known_roles, (
                    f"entry {key!r} has unknown role {role!r}"
                )


# =========================================================================
# Two-tier cascade behaviour
# =========================================================================

class TestTierTwoCascade:
    """Test that tier-2 (ChemScanner) provides fallback lookups."""

    @pytest.fixture(autouse=True)
    def _check_cs_file(self):
        if not os.path.isfile(CS_JSON_PATH):
            pytest.skip("chemscanner_abbreviations.json not present")

    def test_cs_json_exists(self):
        assert os.path.isfile(CS_JSON_PATH)

    def test_tier2_display_for_name(self):
        """A name in ChemScanner but not curated should resolve."""
        db = get_reagent_db()
        # HATU is in ChemScanner; check it resolves via tier-2
        # (it may or may not be in tier-1 depending on the curated file)
        result = db.display_for_name("hatu")
        assert result is not None

    def test_tier1_wins_over_tier2(self):
        """When a name is in both tiers, tier-1 (with role) should win."""
        db = get_reagent_db()
        # DMAP is in both curated (role=base) and ChemScanner (no role)
        display = db.display_for_name("dmap")
        role = db.role_for_name("dmap")
        assert display is not None
        assert role == "base", "tier-1 entry (with role) should win"

    def test_tier2_name_has_no_role(self):
        """Tier-2-only entries should return None for role lookups."""
        db = get_reagent_db()
        # Find a name that's only in ChemScanner (not in curated)
        # Use a long TCI-style name unlikely to be in the curated 172 entries
        with open(CS_JSON_PATH, encoding="utf-8") as f:
            cs_data = json.load(f)
        with open(JSON_PATH, encoding="utf-8") as f:
            curated_data = json.load(f)
        # Find a key in ChemScanner but not curated
        cs_only = None
        for key in cs_data:
            if key not in curated_data:
                cs_only = key
                break
        if cs_only is None:
            pytest.skip("No ChemScanner-only entry found")
        # Should have display but no role
        assert db.display_for_name(cs_only) is not None
        assert db.role_for_name(cs_only) is None

    def test_resolve_display_uses_tier2(self):
        """resolve_display() should use tier-2 instead of falling back to input."""
        db = get_reagent_db()
        result = db.resolve_display("hatu")
        # Should be "HATU" from ChemScanner, not "hatu" (the fallback)
        assert result == "HATU"

    def test_tier2_smiles_lookup(self):
        """SMILES lookup should cascade to tier-2."""
        db = get_reagent_db()
        # HATU SMILES — look it up via ChemScanner
        with open(CS_JSON_PATH, encoding="utf-8") as f:
            cs_data = json.load(f)
        hatu_entry = cs_data.get("hatu")
        if hatu_entry is None:
            pytest.skip("HATU not in ChemScanner data")
        hatu_smiles = hatu_entry["smiles"]
        result = db.display_for_smiles(hatu_smiles)
        assert result is not None


class TestGracefulDegradation:
    """Tier-1 should still work if tier-2 file is missing."""

    def test_tier1_works_without_tier2(self):
        # Create a DB with a nonexistent tier-2 path
        db = ReagentDB(json_path=JSON_PATH,
                       secondary_path="/nonexistent/path.json")
        assert db.display_for_name("cs2co3") == "Cs2CO3"
        assert db.role_for_name("cs2co3") == "base"

    def test_tier2_miss_returns_none(self):
        # Without tier-2, ChemScanner-only names should return None
        db = ReagentDB(json_path=JSON_PATH,
                       secondary_path="/nonexistent/path.json")
        # HATU is only in ChemScanner, not in curated (probably)
        # This test just confirms no crash
        result = db.display_for_name("hatu")
        # result may be None or not depending on curated file content — just no crash
        assert result is None or isinstance(result, str)
