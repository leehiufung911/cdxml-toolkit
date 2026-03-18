"""
LLM-assisted molecule construction via IUPAC name manipulation.

Provides composable tools designed for use with an LLM orchestrator.
The LLM translates natural language descriptions of molecules into tool calls
that manipulate IUPAC names, which are then validated and converted to
structures.

Architecture::

    NL description --> LLM orchestrator --> tool calls --> IUPAC name --> CDXML

The key insight: IUPAC names are a lossless text representation of molecules.
Instead of having an LLM generate SMILES or manipulate CDXML directly, we let
the LLM do "name surgery" — assembling, modifying, and validating IUPAC names
using grounded tools.  This avoids hallucinated SMILES while leveraging LLMs'
strength with natural language.

Layer 2 — Name manipulation tools:
    resolve_to_smiles   — Resolve any chemical identifier to SMILES
    get_prefix_form     — Get IUPAC substituent prefix for a group
    assemble_name       — Build IUPAC name from parent + substituents
    modify_name         — Add/swap/remove substituents in an existing name
    validate_name       — Check if an IUPAC name resolves to a valid molecule
    name_to_structure   — Convert validated name to CDXML
    enumerate_names     — List alternative IUPAC name forms for a molecule

Layer 3 — Graph manipulation tools (for structural transformations):
    list_reactions      — List available named reaction templates
    apply_reaction      — Apply a reaction template (Suzuki, Buchwald, etc.)
    deprotect           — Remove protecting groups (Boc, Fmoc, Cbz, etc.)

Meta:
    get_tool_definitions — Export all tool schemas for LLM function calling

Usage (Python)::

    from cdxml_toolkit.naming.mol_builder import (
        get_prefix_form, assemble_name, validate_name, name_to_structure,
    )

    pf = get_prefix_form("CF3")
    # {'prefix': 'trifluoromethyl', 'source': 'table', 'ok': True}

    result = assemble_name("pyridine", [
        {"locant": "2", "prefix": "chloro"},
        {"locant": "3", "prefix": pf["prefix"]},
    ])
    # {'name': '2-chloro-3-(trifluoromethyl)pyridine', 'valid': True,
    #  'smiles': '...', 'ok': True}

    cdxml = name_to_structure(result["name"])
    # {'cdxml': '<?xml ...', 'ok': True}
"""

import json
import logging
import os
import re
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lazy singletons — avoid import-time cost for heavy dependencies
# ---------------------------------------------------------------------------

_cs_instance = None
_cs_failed = False


def _get_cs():
    """Lazily obtain a ChemScriptBridge instance (or None)."""
    global _cs_instance, _cs_failed
    if _cs_failed:
        return None
    if _cs_instance is not None:
        return _cs_instance
    try:
        from cdxml_toolkit.chemdraw.chemscript_bridge import ChemScriptBridge
        _cs_instance = ChemScriptBridge()
        return _cs_instance
    except Exception as exc:
        logger.debug("ChemScript unavailable: %s", exc)
        _cs_failed = True
        return None


def _rdkit_canonical(smiles: str) -> Optional[str]:
    """Canonical SMILES via RDKit, or None."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


# ---------------------------------------------------------------------------
# Prefix lookup table — covers common med-chem substituents
# ---------------------------------------------------------------------------

# Maps group identifiers (abbreviations, names, formulae) to IUPAC prefix
# form suitable for direct insertion into a substituted name.
# The table is checked case-insensitively.
_PREFIX_TABLE: Dict[str, str] = {
    # --- Halogens ---
    "f": "fluoro", "cl": "chloro", "br": "bromo", "i": "iodo",
    "fluorine": "fluoro", "chlorine": "chloro",
    "bromine": "bromo", "iodine": "iodo",

    # --- Oxygen ---
    "oh": "hydroxy", "ome": "methoxy", "oet": "ethoxy",
    "oac": "acetyloxy", "obn": "benzyloxy", "oph": "phenoxy",
    "methoxy": "methoxy", "ethoxy": "ethoxy", "hydroxy": "hydroxy",
    "ocf3": "trifluoromethoxy", "oipr": "isopropoxy",

    # --- Nitrogen ---
    "nh2": "amino", "nhme": "methylamino", "nme2": "dimethylamino",
    "nhac": "acetamido", "no2": "nitro", "n3": "azido",
    "amino": "amino", "nitro": "nitro", "azido": "azido",

    # --- Simple carbon ---
    "me": "methyl", "et": "ethyl", "pr": "propyl", "npr": "propyl",
    "ipr": "propan-2-yl", "bu": "butyl", "nbu": "butyl",
    "tbu": "tert-butyl", "sbu": "sec-butyl", "ibu": "isobutyl",
    "methyl": "methyl", "ethyl": "ethyl",
    "vinyl": "ethenyl", "allyl": "prop-2-en-1-yl",
    "isopropyl": "propan-2-yl",

    # --- Cycloalkyl ---
    "cyclopropyl": "cyclopropyl", "cyclobutyl": "cyclobutyl",
    "cyclopentyl": "cyclopentyl", "cyclohexyl": "cyclohexyl",
    "cyclopropane": "cyclopropyl", "cyclobutane": "cyclobutyl",
    "cyclopentane": "cyclopentyl", "cyclohexane": "cyclohexyl",

    # --- Aryl ---
    "ph": "phenyl", "bn": "benzyl", "bz": "benzoyl",
    "phenyl": "phenyl", "benzyl": "benzyl",

    # --- Functional groups (abbreviations) ---
    "cn": "cyano", "cho": "formyl", "cooh": "carboxy",
    "co2h": "carboxy", "-cooh": "carboxy", "-co2h": "carboxy",
    "come": "acetyl", "ac": "acetyl",
    "conh2": "carbamoyl", "-conh2": "carbamoyl",
    "coome": "methoxycarbonyl", "co2me": "methoxycarbonyl",
    "meo2c": "methoxycarbonyl", "meoco": "methoxycarbonyl",
    "cooch3": "methoxycarbonyl", "-cooch3": "methoxycarbonyl",
    "-coome": "methoxycarbonyl", "-co2me": "methoxycarbonyl",
    "cooet": "ethoxycarbonyl", "co2et": "ethoxycarbonyl",
    "eto2c": "ethoxycarbonyl", "etoco": "ethoxycarbonyl",
    "cooipr": "isopropoxycarbonyl", "co2ipr": "isopropoxycarbonyl",
    "cootbu": "tert-butoxycarbonyl", "co2tbu": "tert-butoxycarbonyl",
    "-cho": "formyl",

    # --- Functional group descriptors (natural language → prefix) ---
    "methyl ester": "methoxycarbonyl",
    "me ester": "methoxycarbonyl",
    "ome ester": "methoxycarbonyl",
    "ethyl ester": "ethoxycarbonyl",
    "et ester": "ethoxycarbonyl",
    "isopropyl ester": "isopropoxycarbonyl",
    "tert-butyl ester": "tert-butoxycarbonyl",
    "aldehyde": "formyl",
    "ketone": "oxo",
    "carboxylic acid": "carboxy",
    "nitrile": "cyano",
    "amide": "carbamoyl",
    "primary amide": "carbamoyl",
    "alcohol": "hydroxy",
    "hydroxyl": "hydroxy",
    "thiol": "sulfanyl",
    "mercaptan": "sulfanyl",
    "sulfonic acid": "sulfo",
    "sulfonamide": "sulfamoyl",

    # --- Fluorocarbons ---
    "cf3": "trifluoromethyl", "chf2": "difluoromethyl",
    "ccl3": "trichloromethyl",

    # --- Sulphur ---
    "sh": "sulfanyl", "sme": "methylsulfanyl",
    "so2me": "methanesulfonyl", "ms": "methanesulfonyl",
    "so2nh2": "sulfamoyl",

    # --- Heterocycles as substituent prefix ---
    "morpholine": "morpholino", "morpholinyl": "morpholino",
    "morpholino": "morpholino",
    "piperidine": "piperidin-1-yl", "piperidinyl": "piperidin-1-yl",
    "piperazine": "piperazin-1-yl", "piperazinyl": "piperazin-1-yl",
    "pyrrolidine": "pyrrolidin-1-yl", "pyrrolidinyl": "pyrrolidin-1-yl",
    "pyridine": "pyridinyl", "pyridinyl": "pyridinyl",
    "pyrimidine": "pyrimidinyl",
    "thiophene": "thiophen-2-yl", "thienyl": "thiophen-2-yl",
    "furan": "furan-2-yl", "furyl": "furan-2-yl",
    "pyrrole": "pyrrol-1-yl",
    "imidazole": "imidazolyl", "imidazolyl": "imidazolyl",
    "thiazole": "thiazolyl", "thiazolyl": "thiazolyl",
    "oxazole": "oxazolyl", "oxazolyl": "oxazolyl",
    "indole": "indolyl",
}

# IUPAC multiplying prefixes for identical substituents
_MULTIPLIERS = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_query(query: str, use_network: bool = True) -> Optional[Dict]:
    """Multi-tier resolution chain: reagent DB → formula → ChemScript → PubChem.

    Returns {"smiles": ..., "source": ...} or None.
    """
    from rdkit import Chem

    clean = query.strip()
    if not clean:
        return None

    # Tier 1: Reagent DB
    try:
        from cdxml_toolkit.resolve.reagent_db import get_reagent_db
        db = get_reagent_db()
        entry = db.entry_for_name(clean.lower())
        if entry:
            smi = entry.get("smiles")
            if isinstance(smi, list):
                smi = smi[0]
            if smi and Chem.MolFromSmiles(smi):
                return {"smiles": _rdkit_canonical(smi), "source": "reagent_db"}
    except Exception:
        pass

    # Tier 2: Condensed formula
    try:
        from cdxml_toolkit.resolve.condensed_formula import resolve_condensed_formula
        smi = resolve_condensed_formula(clean)
        if smi:
            canon = _rdkit_canonical(smi)
            if canon:
                return {"smiles": canon, "source": "formula"}
    except Exception:
        pass

    # Tier 3: ChemScript (name → SMILES)
    cs = _get_cs()
    if cs is not None:
        try:
            smi = cs.write_data(clean, "smiles", source_format="name")
            if smi and Chem.MolFromSmiles(smi):
                return {"smiles": _rdkit_canonical(smi), "source": "chemscript"}
        except Exception:
            pass

    # Tier 4: PubChem (online)
    if use_network:
        try:
            from cdxml_toolkit.resolve.cas_resolver import resolve_name_to_smiles
            smi = resolve_name_to_smiles(clean)
            if smi:
                canon = _rdkit_canonical(smi)
                if canon:
                    return {"smiles": canon, "source": "pubchem"}
        except Exception:
            pass

    return None


def _name_to_smiles_cs(name: str) -> Optional[str]:
    """Resolve an IUPAC name to SMILES via ChemScript."""
    cs = _get_cs()
    if cs is None:
        return None
    try:
        smi = cs.write_data(name, "smiles", source_format="name")
        if smi:
            return _rdkit_canonical(smi)
    except Exception:
        pass
    return None


def _smiles_to_name_cs(smiles: str) -> Optional[str]:
    """Get IUPAC name for a SMILES string via ChemScript."""
    cs = _get_cs()
    if cs is None:
        return None
    try:
        return cs.get_name(smiles)
    except Exception:
        return None


def _is_complex_prefix(prefix: str) -> bool:
    """Check if a prefix needs parentheses when inserted into a name.

    Complex prefixes contain hyphens, digits, or commas that would be
    ambiguous without enclosing parentheses.
    """
    # Already parenthesised
    if prefix.startswith("(") and prefix.endswith(")"):
        return False
    # Contains internal structure that needs brackets
    return bool(re.search(r"[\d,]", prefix)) and "-" in prefix


def _locant_sort_key(loc: str):
    """Sort locants: numeric before alphabetic, ascending."""
    m = re.match(r"(\d+)(.*)", loc)
    if m:
        return (0, int(m.group(1)), m.group(2))
    return (1, 0, loc)


def _prefix_alpha_key(prefix: str) -> str:
    """IUPAC alphabetical sort key: ignore leading locants/multipliers.

    ``"1,1-difluoroethyl"`` → ``"difluoroethyl"``
    ``"tert-butyl"`` → ``"tert-butyl"``
    """
    stripped = re.sub(r"^[\d,]+-", "", prefix)
    return stripped.lower()


def _try_validate(name: str, use_network: bool = True) -> Optional[str]:
    """Try to resolve a name to canonical SMILES by any available means.

    Returns canonical SMILES or None.
    """
    from rdkit import Chem

    # ChemScript (most reliable for IUPAC names)
    smi = _name_to_smiles_cs(name)
    if smi:
        return smi

    # PubChem fallback (for common names)
    if use_network:
        try:
            from cdxml_toolkit.resolve.cas_resolver import resolve_name_to_smiles
            smi = resolve_name_to_smiles(name)
            if smi:
                return _rdkit_canonical(smi)
        except Exception:
            pass

    return None


# ---------------------------------------------------------------------------
# Tool 1: resolve_to_smiles
# ---------------------------------------------------------------------------

def resolve_to_smiles(query: str, use_network: bool = True) -> Dict[str, Any]:
    """Resolve a chemical identifier to its canonical SMILES string.

    Accepts common names, IUPAC names, abbreviations, condensed formulae,
    and CAS numbers.  Uses a 4-tier resolution chain:
    reagent DB → condensed formula → ChemScript → PubChem.

    Args:
        query: Chemical identifier.  Examples: ``"aspirin"``,
               ``"PhB(OH)2"``, ``"2-chloropyridine"``, ``"534-17-8"``.
        use_network: Allow PubChem lookup (requires internet).

    Returns:
        Dict with keys ``ok``, ``smiles``, ``source``.
        On failure: ``ok=False`` with an ``error`` message.

    Example::

        >>> resolve_to_smiles("Et3N")
        {'ok': True, 'smiles': 'CCN(CC)CC', 'source': 'formula'}
    """
    result = _resolve_query(query, use_network=use_network)
    if result:
        return {"ok": True, "smiles": result["smiles"], "source": result["source"]}
    return {"ok": False, "error": f"Could not resolve '{query}' to a structure."}


# ---------------------------------------------------------------------------
# Tool 2: get_prefix_form
# ---------------------------------------------------------------------------

def get_prefix_form(group: str) -> Dict[str, Any]:
    """Get the IUPAC substituent prefix form for a chemical group.

    Given a group name, abbreviation, or formula, returns the prefix
    string suitable for insertion into an IUPAC name.

    Uses a curated lookup table for common groups (fast, offline), then
    falls back to ChemScript-based naming with the Se-probe for anything
    not in the table.

    Args:
        group: Group identifier.  Examples: ``"CF3"``, ``"morpholine"``,
               ``"NO2"``, ``"cyclopropyl"``, ``"OMe"``.

    Returns:
        Dict with keys ``ok``, ``prefix``, ``source``.
        ``source`` is ``"table"`` for lookup hits, ``"probe"`` for
        ChemScript probe, or ``"passthrough"`` if the input was already
        a valid prefix form.

    Examples::

        >>> get_prefix_form("CF3")
        {'ok': True, 'prefix': 'trifluoromethyl', 'source': 'table'}
        >>> get_prefix_form("morpholine")
        {'ok': True, 'prefix': 'morpholino', 'source': 'table'}
    """
    clean = group.strip()
    if not clean:
        return {"ok": False, "error": "Empty group."}

    # --- Table lookup (case-insensitive) ---
    key = clean.lower()
    if key in _PREFIX_TABLE:
        return {"ok": True, "prefix": _PREFIX_TABLE[key], "source": "table"}

    # --- Check if it's already a valid prefix ---
    # If appending it to "benzene" gives a valid name, it's a prefix.
    test_name = f"1-{clean}benzene" if not clean[0].isdigit() else f"{clean}benzene"
    smi = _try_validate(test_name)
    if smi:
        return {"ok": True, "prefix": clean, "source": "passthrough"}

    # --- Se-probe via name_fragment_as_substituent ---
    # Resolve group to SMILES, add [*] attachment, call the decomposer.
    resolved = _resolve_query(clean, use_network=True)
    if resolved:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(resolved["smiles"])
        if mol:
            # Build [*]-fragment SMILES by attaching dummy at the most
            # likely bonding position (first atom in canonical SMILES).
            # For many simple groups this is correct.
            edit = Chem.RWMol(mol)
            dummy_idx = edit.AddAtom(Chem.Atom(0))  # [*]
            edit.AddBond(0, dummy_idx, Chem.BondType.SINGLE)
            try:
                Chem.SanitizeMol(edit)
                frag_smi = Chem.MolToSmiles(edit.GetMol())
                from .name_decomposer import name_fragment_as_substituent
                prefix = name_fragment_as_substituent(frag_smi, verbose=False)
                if prefix:
                    return {"ok": True, "prefix": prefix, "source": "probe"}
            except Exception:
                pass

    return {
        "ok": False,
        "error": f"Could not determine prefix form for '{group}'.",
    }


# ---------------------------------------------------------------------------
# Tool 3: assemble_name
# ---------------------------------------------------------------------------

def assemble_name(parent: str,
                  substituents: List[Dict[str, str]],
                  validate: bool = True,
                  use_network: bool = True) -> Dict[str, Any]:
    """Assemble an IUPAC name from a parent and substituent list.

    Handles alphabetical ordering, multiplicative prefixes (di-, tri-),
    and parenthesisation of complex substituents.  Optionally validates
    the assembled name by resolving it to SMILES.

    Args:
        parent: Parent ring or chain name (e.g. ``"pyridine"``,
                ``"benzene"``, ``"pentane"``).
        substituents: List of dicts, each with ``"locant"`` (str) and
                      ``"prefix"`` (str).  Example::

                          [{"locant": "2", "prefix": "chloro"},
                           {"locant": "3", "prefix": "methyl"}]
        validate: If True, resolve the assembled name and confirm validity.
        use_network: Allow PubChem for validation.

    Returns:
        Dict with ``ok``, ``name``, and (if validated) ``valid``, ``smiles``.

    Example::

        >>> assemble_name("pyridine", [
        ...     {"locant": "2", "prefix": "chloro"},
        ...     {"locant": "5", "prefix": "nitro"},
        ... ])
        {'ok': True, 'name': '2-chloro-5-nitropyridine', 'valid': True,
         'smiles': '...'}
    """
    if not parent:
        return {"ok": False, "error": "Parent name is required."}
    if not substituents:
        # Bare parent — still valid
        if validate:
            smi = _try_validate(parent, use_network=use_network)
            if smi:
                return {"ok": True, "name": parent, "valid": True, "smiles": smi}
            return {"ok": True, "name": parent, "valid": False, "smiles": None}
        return {"ok": True, "name": parent}

    # --- Group identical prefixes for multipliers ---
    from collections import defaultdict
    groups: Dict[str, List[str]] = defaultdict(list)
    for sub in substituents:
        prefix = sub.get("prefix", "").strip()
        locant = sub.get("locant", "").strip()
        if prefix:
            groups[prefix].append(locant)

    # --- Build prefix fragments, sorted alphabetically by prefix ---
    fragments = []
    for prefix in sorted(groups.keys(), key=_prefix_alpha_key):
        locants = sorted(groups[prefix], key=_locant_sort_key)
        locant_str = ",".join(loc for loc in locants if loc)
        n = len(locants)

        # Format the prefix with optional multiplier
        if n > 1 and prefix in _MULTIPLIERS:
            mult = _MULTIPLIERS.get(n, str(n))
        elif n > 1:
            mult = _MULTIPLIERS.get(n, str(n))
        else:
            mult = ""

        # Parenthesise complex prefixes
        needs_parens = _is_complex_prefix(prefix)
        pfx = f"({prefix})" if needs_parens else prefix

        if mult:
            part = f"{locant_str}-{mult}{pfx}" if locant_str else f"{mult}{pfx}"
        else:
            part = f"{locant_str}-{pfx}" if locant_str else pfx

        fragments.append(part)

    # --- Assemble final name ---
    name = "-".join(fragments) + parent

    result: Dict[str, Any] = {"ok": True, "name": name}

    if validate:
        smi = _try_validate(name, use_network=use_network)
        result["valid"] = smi is not None
        result["smiles"] = smi
        if not smi:
            # Try without parentheses as alternative
            alt_frags = []
            for prefix in sorted(groups.keys(), key=_prefix_alpha_key):
                locants = sorted(groups[prefix], key=_locant_sort_key)
                locant_str = ",".join(loc for loc in locants if loc)
                n = len(locants)
                mult = _MULTIPLIERS.get(n, "") if n > 1 else ""
                part = f"{locant_str}-{mult}{prefix}" if locant_str else f"{mult}{prefix}"
                alt_frags.append(part)
            alt_name = "-".join(alt_frags) + parent
            if alt_name != name:
                alt_smi = _try_validate(alt_name, use_network=use_network)
                if alt_smi:
                    result["name"] = alt_name
                    result["valid"] = True
                    result["smiles"] = alt_smi

    return result


# ---------------------------------------------------------------------------
# Tool 4: modify_name
# ---------------------------------------------------------------------------

def modify_name(name: str,
                operation: str,
                target: Optional[str] = None,
                replacement: Optional[str] = None,
                locant: Optional[str] = None,
                validate: bool = True,
                use_network: bool = True) -> Dict[str, Any]:
    """Modify an IUPAC name by swapping, adding, or removing a substituent.

    Operations:

    - ``"swap"``: Replace *target* prefix with *replacement*.
      E.g. swap "nitro" → "amino" in "4-nitropyridine" → "4-aminopyridine".

    - ``"add"``: Insert *replacement* at *locant*.
      E.g. add "methyl" at "3" to "2-chloropyridine" → "2-chloro-3-methylpyridine".

    - ``"remove"``: Delete the *target* prefix.
      E.g. remove "chloro" from "2-chloro-3-methylpyridine" → "3-methylpyridine".

    For ``"swap"``, the name is re-alphabetised automatically.

    Args:
        name: The IUPAC name to modify.
        operation: ``"swap"``, ``"add"``, or ``"remove"``.
        target: Prefix to replace (swap) or remove (remove).
        replacement: New prefix (swap) or prefix to insert (add).
        locant: Position for insertion (add only).
        validate: Resolve the result to confirm validity.
        use_network: Allow PubChem for validation.

    Returns:
        Dict with ``ok``, ``name``, ``valid``, ``smiles``.

    Examples::

        >>> modify_name("4-nitropyridine", "swap",
        ...             target="nitro", replacement="amino")
        {'ok': True, 'name': '4-aminopyridine', ...}

        >>> modify_name("2-chloropyridine", "add",
        ...             replacement="methyl", locant="3")
        {'ok': True, 'name': '2-chloro-3-methylpyridine', ...}
    """
    if operation == "swap":
        return _modify_swap(name, target, replacement, validate, use_network)
    elif operation == "add":
        return _modify_add(name, replacement, locant, validate, use_network)
    elif operation == "remove":
        return _modify_remove(name, target, validate, use_network)
    else:
        return {"ok": False, "error": f"Unknown operation '{operation}'. "
                "Use 'swap', 'add', or 'remove'."}


def _parse_name_components(name: str) -> Optional[Dict]:
    """Best-effort parse of a substituted IUPAC name into components.

    Splits a name like ``"2-chloro-5-(trifluoromethyl)pyridine"`` into::

        {"parent": "pyridine",
         "substituents": [{"locant": "2", "prefix": "chloro"},
                          {"locant": "5", "prefix": "trifluoromethyl"}]}

    Uses the aligned namer's ring system list for parent detection.
    """
    try:
        from .aligned_namer import _KNOWN_RINGS
        rings = _KNOWN_RINGS
    except ImportError:
        rings = set()

    # Also try common chain parents
    chains = [
        "icosane", "nonadecane", "octadecane", "heptadecane", "hexadecane",
        "pentadecane", "tetradecane", "tridecane", "dodecane", "undecane",
        "decane", "nonane", "octane", "heptane", "hexane", "pentane",
        "butane", "propane", "ethane", "methane",
        "icosanoic acid", "nonadecanoic acid", "octadecanoic acid",
    ]
    all_parents = sorted(
        list(rings) + chains, key=len, reverse=True
    )

    # Find the parent: longest known name that matches the tail
    parent = None
    prefix_part = ""
    for p in all_parents:
        if name.endswith(p):
            prefix_part = name[:-len(p)]
            parent = p
            break

    # Fallback: if no known parent, try splitting at the last segment
    # that doesn't start with a digit
    if parent is None:
        # Try to identify parent as the last non-prefixed segment
        # Pattern: everything after the last "-" that isn't a locant-prefix pair
        parts = name.rsplit("-", 1)
        if len(parts) == 2 and not re.match(r"^\d", parts[1]):
            parent = parts[1]
            prefix_part = parts[0] + "-"
        else:
            parent = name
            prefix_part = ""

    if not prefix_part.strip("-"):
        return {"parent": parent, "substituents": []}

    # Parse prefix_part into (locant, prefix) pairs
    prefix_str = prefix_part.rstrip("-")
    substituents = []

    # Pattern: locant(s)-[multiplier][(]prefix[)] or locant(s)-[multiplier]prefix
    # Walk through segments
    segments = _split_prefix_segments(prefix_str)
    for seg in segments:
        parsed = _parse_single_prefix(seg)
        if parsed:
            substituents.extend(parsed)

    return {"parent": parent, "substituents": substituents}


def _split_prefix_segments(prefix_str: str) -> List[str]:
    """Split a prefix string into individual prefix segments.

    Handles parenthesised prefixes correctly:
    ``"2-chloro-3-(trifluoromethyl)"`` → ``["2-chloro", "3-(trifluoromethyl)"]``
    """
    segments = []
    current = ""
    depth = 0
    for ch in prefix_str:
        if ch == "(":
            depth += 1
            current += ch
        elif ch == ")":
            depth -= 1
            current += ch
        elif ch == "-" and depth == 0:
            if current:
                # Check: is this a separator between segments, or within one?
                # A segment boundary is after a prefix (lowercase letter or ')').
                # Within a segment: after a locant (digit) or multiplier.
                if current and (current[-1].isalpha() and current[-1].islower()
                                or current[-1] == ")"):
                    segments.append(current)
                    current = ""
                else:
                    current += ch
            else:
                current += ch
        else:
            current += ch
    if current:
        segments.append(current)
    return segments


def _parse_single_prefix(segment: str) -> Optional[List[Dict[str, str]]]:
    """Parse a single prefix segment like '2-chloro' or '2,4-dichloro'.

    Returns list of {"locant": ..., "prefix": ...} dicts.
    """
    # Handle multiplied: 2,4-dichloro
    m = re.match(
        r"^([\d,]+)-(?:di|tri|tetra|penta|hexa)"
        r"[\(\[]?([a-zA-Z][\w,\-]*?)[\)\]]?$",
        segment,
    )
    if m:
        locants = m.group(1).split(",")
        prefix = m.group(2)
        return [{"locant": loc, "prefix": prefix} for loc in locants]

    # Handle parenthesised: 3-(trifluoromethyl)
    m = re.match(r"^(\d+\w?)-\((.+)\)$", segment)
    if m:
        return [{"locant": m.group(1), "prefix": m.group(2)}]

    # Handle simple: 2-chloro
    m = re.match(r"^(\d+\w?)-([a-zA-Z][\w\-]*)$", segment)
    if m:
        return [{"locant": m.group(1), "prefix": m.group(2)}]

    # No locant: just a prefix (e.g., "amino" without locant)
    if re.match(r"^[a-zA-Z]", segment):
        return [{"locant": "", "prefix": segment}]

    return None


def _modify_swap(name, target, replacement, validate, use_network):
    """Swap one prefix for another and re-assemble."""
    if not target or not replacement:
        return {"ok": False, "error": "Both 'target' and 'replacement' required for swap."}

    parsed = _parse_name_components(name)
    if parsed is None:
        return {"ok": False, "error": f"Could not parse name '{name}'."}

    subs = parsed["substituents"]
    found = False
    for sub in subs:
        if sub["prefix"] == target:
            sub["prefix"] = replacement
            found = True
    if not found:
        return {
            "ok": False,
            "error": f"Prefix '{target}' not found in '{name}'.",
            "found_prefixes": [s["prefix"] for s in subs],
        }

    return assemble_name(parsed["parent"], subs, validate=validate,
                         use_network=use_network)


def _modify_add(name, prefix, locant, validate, use_network):
    """Add a new substituent to an existing name."""
    if not prefix:
        return {"ok": False, "error": "'replacement' (prefix to add) is required."}
    if not locant:
        return {"ok": False, "error": "'locant' is required for add operation."}

    parsed = _parse_name_components(name)
    if parsed is None:
        return {"ok": False, "error": f"Could not parse name '{name}'."}

    parsed["substituents"].append({"locant": locant, "prefix": prefix})
    return assemble_name(parsed["parent"], parsed["substituents"],
                         validate=validate, use_network=use_network)


def _modify_remove(name, target, validate, use_network):
    """Remove a substituent from a name."""
    if not target:
        return {"ok": False, "error": "'target' prefix is required for remove."}

    parsed = _parse_name_components(name)
    if parsed is None:
        return {"ok": False, "error": f"Could not parse name '{name}'."}

    original_len = len(parsed["substituents"])
    parsed["substituents"] = [
        s for s in parsed["substituents"] if s["prefix"] != target
    ]
    if len(parsed["substituents"]) == original_len:
        return {
            "ok": False,
            "error": f"Prefix '{target}' not found in '{name}'.",
            "found_prefixes": [s["prefix"] for s in parsed["substituents"]],
        }

    return assemble_name(parsed["parent"], parsed["substituents"],
                         validate=validate, use_network=use_network)


# ---------------------------------------------------------------------------
# Tool 5: validate_name
# ---------------------------------------------------------------------------

def validate_name(name: str,
                  use_network: bool = True) -> Dict[str, Any]:
    """Validate an IUPAC name and return its SMILES if valid.

    Attempts to resolve the name to a structure using ChemScript
    (preferred) or PubChem (fallback).  Returns whether the name is
    valid and the canonical SMILES.

    Args:
        name: IUPAC name to validate.
        use_network: Allow PubChem lookup.

    Returns:
        Dict with ``ok``, ``valid``, ``smiles``, ``name``.

    Example::

        >>> validate_name("2-chloropyridine")
        {'ok': True, 'valid': True, 'smiles': 'Clc1ccccn1', 'name': '2-chloropyridine'}
    """
    smi = _try_validate(name, use_network=use_network)
    if smi:
        # Also get the canonical IUPAC name if ChemScript is available
        canonical = _smiles_to_name_cs(smi)
        return {
            "ok": True,
            "valid": True,
            "smiles": smi,
            "name": name,
            "canonical_name": canonical,
        }
    return {"ok": True, "valid": False, "smiles": None, "name": name}


# ---------------------------------------------------------------------------
# Tool 6: name_to_structure
# ---------------------------------------------------------------------------

def name_to_structure(name: str,
                      output_format: str = "cdxml") -> Dict[str, Any]:
    """Convert a chemical name to a structure in the requested format.

    Resolves the name, generates 2D coordinates, and returns the
    structure as a string (CDXML, SMILES, or MOL).

    Args:
        name: IUPAC or common name.
        output_format: ``"cdxml"`` (default), ``"smiles"``, or ``"mol"``.

    Returns:
        Dict with ``ok`` and the structure data (key matches format name).

    Example::

        >>> result = name_to_structure("2-chloropyridine")
        >>> result["ok"]
        True
        >>> result["cdxml"][:20]
        '<?xml version="1.0"'
    """
    fmt = output_format.lower()

    if fmt == "smiles":
        smi = _try_validate(name)
        if smi:
            return {"ok": True, "smiles": smi}
        return {"ok": False, "error": f"Could not resolve '{name}'."}

    # For CDXML and MOL, prefer ChemScript (gives ACS-styled 2D)
    cs = _get_cs()
    if cs is not None:
        try:
            if fmt == "cdxml":
                cdxml = cs.name_to_cdxml(name)
                return {"ok": True, "cdxml": cdxml}
            elif fmt == "mol":
                mol_data = cs.write_data(name, "mol", source_format="name")
                return {"ok": True, "mol": mol_data}
        except Exception:
            pass

    # Fallback: resolve to SMILES, then generate structure via RDKit
    smi = _try_validate(name)
    if not smi:
        return {"ok": False, "error": f"Could not resolve '{name}'."}

    if fmt == "cdxml":
        # Try ChemScript with SMILES input
        if cs is not None:
            try:
                cdxml = cs.smiles_to_cdxml(smi)
                return {"ok": True, "cdxml": cdxml}
            except Exception:
                pass
        return {"ok": False,
                "error": "CDXML output requires ChemScript. "
                         f"Name resolved to SMILES: {smi}",
                "smiles": smi}

    if fmt == "mol":
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smi)
        if mol:
            AllChem.Compute2DCoords(mol)
            return {"ok": True, "mol": Chem.MolToMolBlock(mol)}
        return {"ok": False, "error": "RDKit could not generate MOL block."}

    return {"ok": False, "error": f"Unknown format '{fmt}'. Use cdxml, smiles, or mol."}


# ---------------------------------------------------------------------------
# Tool 7: enumerate_names
# ---------------------------------------------------------------------------

def enumerate_names(identifier: str,
                    use_network: bool = True) -> Dict[str, Any]:
    """Enumerate alternative IUPAC name forms for a molecule.

    Given a chemical name or SMILES, returns the canonical IUPAC name plus
    alternative forms where substituents appear as different prefixes or
    where a different parent ring/chain is chosen.  This is essential for
    name surgery: it lets you see functional groups as swappable prefixes.

    For example, ``"1-(4-bromophenyl)ethan-1-one"`` (a ketone in suffix
    form) generates alternatives including ``"1-acetyl-4-bromobenzene"``
    where the ketone appears as the prefix ``"acetyl"`` — now swappable
    via ``modify_name``.

    Args:
        identifier: Chemical name, SMILES, abbreviation, or any
            identifier accepted by ``resolve_to_smiles``.
        use_network: Allow PubChem for resolution.

    Returns:
        Dict with:

        - ``ok``: bool
        - ``canonical_name``: the ChemDraw canonical IUPAC name
        - ``smiles``: canonical SMILES
        - ``names``: list of dicts, each with ``name`` (str),
          ``valid`` (bool), ``strategy`` (str), and ``prefixes``
          (list of prefix strings visible in that name form).
          The canonical name is always the first entry.

    Example::

        >>> result = enumerate_names("1-(4-bromophenyl)ethan-1-one")
        >>> for n in result["names"]:
        ...     print(n["name"], n["prefixes"])
        1-(4-bromophenyl)ethan-1-one  ['(4-bromophenyl)']
        1-acetyl-4-bromobenzene       ['acetyl', 'bromo']
        ...
    """
    # Resolve to SMILES — try direct SMILES parse first, then name resolution
    from rdkit import Chem as _Chem
    _test_mol = _Chem.MolFromSmiles(identifier)
    if _test_mol is not None:
        smiles = _Chem.MolToSmiles(_test_mol)
    else:
        resolved = _resolve_query(identifier, use_network=use_network)
        if not resolved:
            return {"ok": False,
                    "error": f"Could not resolve '{identifier}' to a structure."}
        smiles = resolved["smiles"]

    # Run decomposition
    try:
        from .name_decomposer import decompose_name
        result = decompose_name(smiles, verbose=False, timeout=30.0)
    except Exception as exc:
        return {"ok": False,
                "error": f"Decomposition failed: {exc}",
                "smiles": smiles}

    if result.errors:
        return {"ok": False,
                "error": "; ".join(result.errors),
                "smiles": smiles}

    canon = result.canonical_name
    if not canon:
        return {"ok": False,
                "error": "Could not determine canonical name.",
                "smiles": smiles}

    # Build the output list, canonical first
    names = []

    # Parse prefixes from each name form
    canon_parsed = _parse_name_components(canon)
    canon_prefixes = ([s["prefix"] for s in canon_parsed["substituents"]]
                      if canon_parsed else [])
    names.append({
        "name": canon,
        "valid": True,
        "strategy": "canonical",
        "prefixes": canon_prefixes,
    })

    # Add valid alternatives
    seen = {canon}
    for alt in result.alternatives:
        if not alt.valid:
            continue
        if alt.name in seen:
            continue
        seen.add(alt.name)

        parsed = _parse_name_components(alt.name)
        prefixes = ([s["prefix"] for s in parsed["substituents"]]
                    if parsed else [])
        names.append({
            "name": alt.name,
            "valid": True,
            "strategy": alt.strategy,
            "prefixes": prefixes,
        })

    return {
        "ok": True,
        "canonical_name": canon,
        "smiles": smiles,
        "names": names,
    }


# ---------------------------------------------------------------------------
# Layer 3: Graph manipulation — reaction templates
# ---------------------------------------------------------------------------

# Hand-curated templates for common med-chem transformations that an LLM
# will recognise by name.  These supplement the larger ring-forming
# collection loaded from reactions_datamol.json.
_CLASSIC_TEMPLATES: Dict[str, Dict[str, Any]] = {
    "suzuki_coupling": {
        "description": "Suzuki coupling: aryl halide + boronic acid to biaryl",
        "smarts": "[c:1][Br,I].[#6:2][B]([OH])[OH]>>[c:1]-[#6:2]",
        "n_reactants": 2,
        "substrate_hint": "aryl bromide or iodide",
        "reagent_hint": "boronic acid",
        "conditions": ["Pd(dppf)Cl2", "K2CO3", "dioxane/H2O", "80 °C"],
        "category": "coupling",
    },
    "buchwald_amination": {
        "description":
            "Buchwald-Hartwig amination: aryl halide + amine to aryl amine",
        "smarts": "[c:1][Cl,Br,I].[NX3;H2,H1:2]>>[c:1]-[N:2]",
        "n_reactants": 2,
        "substrate_hint": "aryl halide",
        "reagent_hint": "primary or secondary amine",
        "conditions": ["Pd2(dba)3", "XPhos", "Cs2CO3", "toluene", "100 °C"],
        "category": "coupling",
    },
    "snar": {
        "description":
            "Nucleophilic aromatic substitution: activated aryl halide + "
            "nucleophile",
        "smarts": "[c:1][F,Cl].[NX3;H2,H1:2]>>[c:1]-[N:2]",
        "n_reactants": 2,
        "substrate_hint": "electron-poor aryl fluoride or chloride",
        "reagent_hint": "amine nucleophile",
        "conditions": ["DIPEA", "DMSO or NMP", "80-120 °C"],
        "category": "coupling",
    },
    "amide_coupling": {
        "description": "Amide bond formation: carboxylic acid + amine",
        "smarts":
            "[C:1](=[O:2])[OH].[NX3;H2,H1:3]>>[C:1](=[O:2])-[N:3]",
        "n_reactants": 2,
        "substrate_hint": "carboxylic acid",
        "reagent_hint": "primary or secondary amine",
        "conditions": ["HATU", "DIPEA", "DMF", "rt"],
        "category": "coupling",
    },
    "reductive_amination": {
        "description":
            "Reductive amination: aldehyde or ketone + amine to amine",
        "smarts": "[C:1](=[O:2]).[NX3;H2,H1:3]>>[C:1]-[N:3]",
        "n_reactants": 2,
        "substrate_hint": "aldehyde or ketone",
        "reagent_hint": "primary or secondary amine",
        "conditions": ["NaBH(OAc)3", "AcOH", "DCE", "rt"],
        "category": "functional_group",
    },
    "nitro_reduction": {
        "description": "Nitro group reduction to amine (ArNO2 to ArNH2)",
        "smarts": "[c:1][N+](=[O])[O-]>>[c:1]N",
        "n_reactants": 1,
        "substrate_hint": "aromatic nitro compound",
        "reagent_hint": None,
        "conditions": ["SnCl2·2H2O", "EtOH", "80 °C"],
        "category": "functional_group",
    },
    "ester_hydrolysis": {
        "description": "Ester hydrolysis to carboxylic acid",
        "smarts":
            "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH]",
        "n_reactants": 1,
        "substrate_hint": "ester",
        "reagent_hint": None,
        "conditions": ["LiOH", "THF/H2O", "rt"],
        "category": "functional_group",
    },
    "n_alkylation": {
        "description": "N-Alkylation: amine + alkyl halide",
        "smarts":
            "[NX3;H2,H1:1].[C:2][Cl,Br,I]>>[N:1]-[C:2]",
        "n_reactants": 2,
        "substrate_hint": "amine",
        "reagent_hint": "alkyl halide",
        "conditions": ["K2CO3", "DMF", "60 °C"],
        "category": "coupling",
    },
    "sonogashira_coupling": {
        "description":
            "Sonogashira coupling: aryl halide + terminal alkyne",
        "smarts":
            "[c:1][Br,I].[CH:2]#[C:3]>>[c:1]-[C:2]#[C:3]",
        "n_reactants": 2,
        "substrate_hint": "aryl bromide or iodide",
        "reagent_hint": "terminal alkyne",
        "conditions": [
            "PdCl2(PPh3)2", "CuI", "Et3N", "THF", "rt",
        ],
        "category": "coupling",
    },
    "heck_reaction": {
        "description":
            "Heck reaction: aryl halide + alkene to substituted alkene",
        "smarts":
            "[c:1][Br,I].[CH:2]=[CH2:3]>>[c:1]/[CH:2]=[CH2:3]",
        "n_reactants": 2,
        "substrate_hint": "aryl halide",
        "reagent_hint": "terminal alkene",
        "conditions": [
            "Pd(OAc)2", "P(o-tol)3", "Et3N", "DMF", "100 °C",
        ],
        "category": "coupling",
    },
    "alcohol_oxidation": {
        "description": "Alcohol oxidation to aldehyde or ketone",
        "smarts":
            "[C:1][OH:2]>>[C:1]=[O:2]",
        "n_reactants": 1,
        "substrate_hint": "primary or secondary alcohol",
        "reagent_hint": None,
        "conditions": ["Dess-Martin periodinane", "DCM", "rt"],
        "category": "functional_group",
    },
    "grignard_addition": {
        "description":
            "Grignard / organometallic addition to aldehyde or ketone",
        "smarts":
            "[C:1](=[O:2])[#6:4].[#6:3][Mg]>>[C:1]([OH:2])([#6:4])-[#6:3]",
        "n_reactants": 2,
        "substrate_hint": "aldehyde or ketone",
        "reagent_hint": "Grignard reagent (RMgBr SMILES)",
        "conditions": ["THF", "-78 °C to rt"],
        "category": "functional_group",
    },
}


# ---------------------------------------------------------------------------
# Dynamic loading of ring-forming heterocyclic templates from datamol
# ---------------------------------------------------------------------------

_datamol_cache: Optional[Dict[str, Dict[str, Any]]] = None


def _load_datamol_templates() -> Dict[str, Dict[str, Any]]:
    """Load ring-forming heterocyclic reaction templates from datamol JSON.

    Reads ``reactions_datamol.json`` (127 curated reaction templates from
    the datamol project, Apache 2.0) and extracts those tagged as
    heterocycle formation / cyclization.  Each entry is converted to our
    standard template format with snake_case keys.

    Returns:
        Dict mapping template name to template dict.
    """
    global _datamol_cache
    if _datamol_cache is not None:
        return _datamol_cache

    json_path = os.path.join(os.path.dirname(__file__), "reactions_datamol.json")
    if not os.path.exists(json_path):
        logger.warning("reactions_datamol.json not found — ring-forming "
                        "templates unavailable")
        _datamol_cache = {}
        return _datamol_cache

    with open(json_path, encoding="utf-8") as fh:
        raw = json.load(fh)

    ring_tags = {"heterocycle formation", "cyclization", "ring formation"}
    templates: Dict[str, Dict[str, Any]] = {}

    for key, entry in raw.items():
        tags = set(entry.get("tags", []))
        if not tags & ring_tags:
            continue
        syn_smarts = entry.get("syn_smarts", "")
        if not syn_smarts:
            continue

        # Derive template name: JSON key is already kebab-case
        # Convert to snake_case for consistency
        tname = key.replace("-", "_")

        # Count reactant fragments (separated by '.')
        reactant_part = syn_smarts.split(">>")[0] if ">>" in syn_smarts else ""
        # Count top-level dots (outside brackets)
        n_reactants = 1
        depth = 0
        for ch in reactant_part:
            if ch == "[":
                depth += 1
            elif ch == "]":
                depth -= 1
            elif ch == "." and depth == 0:
                n_reactants += 1

        templates[tname] = {
            "description": entry.get("description", entry.get("long_name", key)),
            "long_name": entry.get("long_name", ""),
            "smarts": syn_smarts,
            "n_reactants": n_reactants,
            "substrate_hint": ", ".join(entry.get("rhs_classes", [])),
            "reagent_hint": (", ".join(entry.get("rhs_classes", [])[1:])
                             if n_reactants > 1 and len(entry.get("rhs_classes", [])) > 1
                             else None),
            "conditions": [],  # literature conditions vary
            "category": "heterocycle_formation",
            "tags": list(tags & ring_tags),
            "source": "datamol",
        }

    _datamol_cache = templates
    logger.debug("Loaded %d ring-forming templates from datamol", len(templates))
    return _datamol_cache


# Merged registry: classic hand-written + datamol ring-forming
_merged_templates: Optional[Dict[str, Dict[str, Any]]] = None


def _get_reaction_templates() -> Dict[str, Dict[str, Any]]:
    """Return the merged reaction template registry (lazy-loaded).

    Classic hand-written templates (couplings, functional group transforms)
    are merged with ~60 ring-forming heterocyclic templates from datamol.
    Classic templates take priority on name collisions.
    """
    global _merged_templates
    if _merged_templates is not None:
        return _merged_templates

    datamol = _load_datamol_templates()
    merged = dict(datamol)  # datamol first, classic overrides
    merged.update(_CLASSIC_TEMPLATES)
    _merged_templates = merged
    return _merged_templates


# ---------------------------------------------------------------------------
# Tool 7: list_reactions
# ---------------------------------------------------------------------------

def list_reactions(category: Optional[str] = None) -> Dict[str, Any]:
    """List available named reaction templates.

    Returns a summary of each reaction: name, description, number of
    reactants required, and typical conditions.  Use this to find the
    right template before calling ``apply_reaction``.

    Args:
        category: Optional filter.  One of ``"coupling"``,
            ``"functional_group"``, or ``"heterocycle_formation"``.
            If *None*, all templates are returned.

    Returns:
        Dict with ``ok``, ``reactions`` (list of summaries), and
        ``categories`` (list of available category names).

    Example::

        >>> result = list_reactions()
        >>> for r in result["reactions"]:
        ...     print(r["name"], "—", r["description"])
        suzuki_coupling — Suzuki coupling: aryl halide + boronic acid to biaryl
        ...
        >>> result = list_reactions(category="heterocycle_formation")
    """
    templates = _get_reaction_templates()
    rxns = []
    cats = set()
    for name, tmpl in templates.items():
        cat = tmpl.get("category", "other")
        cats.add(cat)
        if category and cat != category:
            continue
        rxns.append({
            "name": name,
            "description": tmpl["description"],
            "n_reactants": tmpl["n_reactants"],
            "substrate_hint": tmpl["substrate_hint"],
            "reagent_hint": tmpl.get("reagent_hint"),
            "conditions": tmpl.get("conditions", []),
            "category": cat,
        })
    return {"ok": True, "reactions": rxns, "categories": sorted(cats)}


# ---------------------------------------------------------------------------
# Tool 8: apply_reaction
# ---------------------------------------------------------------------------

def apply_reaction(reaction_name: str,
                   substrate: str,
                   reagent: Optional[str] = None) -> Dict[str, Any]:
    """Apply a named reaction template to transform a substrate.

    Takes a substrate SMILES (and optionally a reagent SMILES for
    bimolecular reactions) and returns the product(s).

    Args:
        reaction_name: Template name from ``list_reactions()``.
        substrate: SMILES of the main substrate.
        reagent: SMILES of the coupling partner (for 2-reactant rxns).
                 Can also be a chemical name or abbreviation.

    Returns:
        Dict with ``ok``, ``products`` (list of product dicts with
        ``smiles`` and ``name`` keys), and ``conditions``.

    Example::

        >>> apply_reaction("nitro_reduction", "c1ccc([N+](=O)[O-])cc1")
        {'ok': True,
         'products': [{'smiles': 'Nc1ccccc1', 'name': 'aniline'}],
         'conditions': ['SnCl2·2H2O', 'EtOH', '80 °C']}
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Look up template
    templates = _get_reaction_templates()
    tmpl = templates.get(reaction_name)
    if tmpl is None:
        available = ", ".join(sorted(templates.keys()))
        return {
            "ok": False,
            "error": f"Unknown reaction '{reaction_name}'. "
                     f"Available: {available}",
        }

    # Parse reaction SMARTS
    try:
        rxn = AllChem.ReactionFromSmarts(tmpl["smarts"])
    except Exception as exc:
        return {"ok": False, "error": f"Invalid reaction SMARTS: {exc}"}

    # Parse substrate
    sub_mol = Chem.MolFromSmiles(substrate)
    if sub_mol is None:
        # Maybe it's a name — try resolving
        resolved = _resolve_query(substrate)
        if resolved:
            sub_mol = Chem.MolFromSmiles(resolved["smiles"])
        if sub_mol is None:
            return {"ok": False, "error": f"Could not parse substrate '{substrate}'."}

    # Handle reagent for bimolecular reactions
    if tmpl["n_reactants"] == 2:
        if not reagent:
            return {
                "ok": False,
                "error": f"Reaction '{reaction_name}' requires a reagent. "
                         f"Expected: {tmpl.get('reagent_hint', 'coupling partner')}.",
            }
        rea_mol = Chem.MolFromSmiles(reagent)
        if rea_mol is None:
            # Try resolving as name
            resolved = _resolve_query(reagent)
            if resolved:
                rea_mol = Chem.MolFromSmiles(resolved["smiles"])
            if rea_mol is None:
                return {"ok": False,
                        "error": f"Could not parse reagent '{reagent}'."}
        reactants = (sub_mol, rea_mol)
    else:
        reactants = (sub_mol,)

    # Run reaction
    try:
        product_sets = rxn.RunReactants(reactants)
    except Exception as exc:
        return {"ok": False, "error": f"Reaction failed: {exc}"}

    if not product_sets:
        # Try swapped reactant order for bimolecular reactions
        if tmpl["n_reactants"] == 2:
            try:
                product_sets = rxn.RunReactants((reactants[1], reactants[0]))
            except Exception:
                pass
        if not product_sets:
            return {
                "ok": False,
                "error": "No products formed. Check that the substrate "
                         f"matches: {tmpl['substrate_hint']}.",
            }

    # Collect unique products
    seen = set()
    products = []
    for prod_tuple in product_sets:
        for prod in prod_tuple:
            try:
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                if smi not in seen:
                    seen.add(smi)
                    name = _smiles_to_name_cs(smi)
                    products.append({"smiles": smi, "name": name})
            except Exception:
                continue

    if not products:
        return {"ok": False, "error": "Products could not be sanitised."}

    return {
        "ok": True,
        "products": products,
        "conditions": tmpl["conditions"],
        "reaction": reaction_name,
    }


# ---------------------------------------------------------------------------
# Tool 9: deprotect
# ---------------------------------------------------------------------------

def deprotect(smiles: str) -> Dict[str, Any]:
    """Remove common protecting groups from a molecule.

    Uses RDKit's built-in deprotection library (25 templates covering
    Boc, Fmoc, Cbz, TBS, THP, Bn, Ac, PMB, Tr, and more).

    Args:
        smiles: SMILES of the protected molecule.

    Returns:
        Dict with ``ok``, ``product_smiles``, ``product_name``,
        and ``removed`` (list of protecting group abbreviations removed).

    Example::

        >>> deprotect("O=C(OC(C)(C)C)Nc1ccccc1")  # Boc-aniline
        {'ok': True, 'product_smiles': 'Nc1ccccc1', 'product_name': 'aniline',
         'removed': ['Boc']}
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        resolved = _resolve_query(smiles)
        if resolved:
            mol = Chem.MolFromSmiles(resolved["smiles"])
        if mol is None:
            return {"ok": False, "error": f"Could not parse '{smiles}'."}

    original_smi = Chem.MolToSmiles(mol)

    try:
        from rdkit.Chem import rdDeprotect
        result = rdDeprotect.Deprotect(mol)
    except ImportError:
        return {"ok": False, "error": "rdDeprotect not available in this RDKit build."}
    except Exception as exc:
        return {"ok": False, "error": f"Deprotection failed: {exc}"}

    product_smi = Chem.MolToSmiles(result)

    if product_smi == original_smi:
        return {
            "ok": True,
            "product_smiles": product_smi,
            "product_name": _smiles_to_name_cs(product_smi),
            "removed": [],
            "note": "No protecting groups detected.",
        }

    # Identify which PGs were removed by checking each template
    removed = []
    try:
        deprots = rdDeprotect.GetDeprotections()
        for d in deprots:
            rxn_sma = d.reaction_smarts
            from rdkit.Chem import AllChem
            rxn = AllChem.ReactionFromSmarts(rxn_sma)
            try:
                prods = rxn.RunReactants((mol,))
                if prods:
                    for ptuple in prods:
                        for p in ptuple:
                            try:
                                Chem.SanitizeMol(p)
                                if Chem.MolToSmiles(p) != original_smi:
                                    removed.append(d.abbreviation)
                                    break
                            except Exception:
                                continue
                        if removed and removed[-1] == d.abbreviation:
                            break
            except Exception:
                continue
    except Exception:
        pass

    name = _smiles_to_name_cs(product_smi)
    return {
        "ok": True,
        "product_smiles": product_smi,
        "product_name": name,
        "removed": removed,
    }


# ---------------------------------------------------------------------------
# Tool definitions for LLM function calling
# ---------------------------------------------------------------------------

def get_tool_definitions() -> List[Dict[str, Any]]:
    """Return tool schemas suitable for LLM function calling (Claude/OpenAI).

    Each tool definition follows the Anthropic tool-use format::

        {"name": "...", "description": "...", "input_schema": {...}}

    The LLM orchestrator should register these as available tools and
    call the corresponding Python functions based on the LLM's output.

    Returns:
        List of tool definition dicts.
    """
    return [
        {
            "name": "resolve_to_smiles",
            "description": (
                "Resolve a chemical identifier (name, abbreviation, formula, "
                "or CAS number) to a canonical SMILES string.  Use this to "
                "look up any chemical you need to work with.\n\n"
                "Examples of valid queries:\n"
                '  - Common names: "aspirin", "morpholine", "HATU"\n'
                '  - IUPAC names: "2-chloropyridine", "4-methylbenzoic acid"\n'
                '  - Formulae: "PhB(OH)2", "Et3N", "CF3COOH"\n'
                '  - CAS numbers: "534-17-8"\n'
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Chemical identifier to resolve.",
                    },
                },
                "required": ["query"],
            },
        },
        {
            "name": "get_prefix_form",
            "description": (
                "Get the IUPAC substituent prefix form for a chemical group "
                "so it can be used in assemble_name.  Returns the prefix "
                "string (e.g. 'trifluoromethyl' for 'CF3', 'morpholino' for "
                "'morpholine').\n\n"
                "Use this when you know what group to attach but need its "
                "correct IUPAC prefix name.\n\n"
                "Examples:\n"
                '  - "CF3" -> "trifluoromethyl"\n'
                '  - "NO2" -> "nitro"\n'
                '  - "OMe" -> "methoxy"\n'
                '  - "morpholine" -> "morpholino"\n'
                '  - "cyclopropane" -> "cyclopropyl"\n'
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "group": {
                        "type": "string",
                        "description": (
                            "Group to look up: abbreviation ('CF3', 'OMe'), "
                            "name ('morpholine'), or formula ('CHF2')."
                        ),
                    },
                },
                "required": ["group"],
            },
        },
        {
            "name": "assemble_name",
            "description": (
                "Build an IUPAC name from a parent ring/chain and a list of "
                "substituents.  Handles alphabetical ordering and multiplying "
                "prefixes (di-, tri-) automatically.  Validates the assembled "
                "name by resolving it to a structure.\n\n"
                "Example:\n"
                "  parent: 'pyridine'\n"
                "  substituents: [\n"
                '    {"locant": "2", "prefix": "chloro"},\n'
                '    {"locant": "5", "prefix": "trifluoromethyl"}\n'
                "  ]\n"
                "  -> '2-chloro-5-(trifluoromethyl)pyridine'\n"
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "parent": {
                        "type": "string",
                        "description": (
                            "Parent ring or chain name "
                            "(e.g. 'pyridine', 'benzene', 'pentane')."
                        ),
                    },
                    "substituents": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "locant": {
                                    "type": "string",
                                    "description": "Position number (e.g. '2', '3').",
                                },
                                "prefix": {
                                    "type": "string",
                                    "description": (
                                        "IUPAC prefix (e.g. 'chloro', 'methyl'). "
                                        "Use get_prefix_form first if unsure."
                                    ),
                                },
                            },
                            "required": ["locant", "prefix"],
                        },
                        "description": "List of substituents with positions.",
                    },
                },
                "required": ["parent", "substituents"],
            },
        },
        {
            "name": "modify_name",
            "description": (
                "Modify an existing IUPAC name by swapping, adding, or "
                "removing a substituent.  The name is re-alphabetised and "
                "validated automatically.\n\n"
                "Operations:\n"
                "  - 'swap': Replace target prefix with replacement.\n"
                "    Example: swap 'nitro' -> 'amino' in '4-nitropyridine'\n"
                "  - 'add': Insert replacement prefix at locant.\n"
                "    Example: add 'methyl' at '3' to '2-chloropyridine'\n"
                "  - 'remove': Delete the target prefix.\n"
                "    Example: remove 'chloro' from '2-chloro-3-methylpyridine'\n"
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "The IUPAC name to modify.",
                    },
                    "operation": {
                        "type": "string",
                        "enum": ["swap", "add", "remove"],
                        "description": "Type of modification.",
                    },
                    "target": {
                        "type": "string",
                        "description": "Prefix to replace (swap) or remove (remove).",
                    },
                    "replacement": {
                        "type": "string",
                        "description": "New prefix (swap) or prefix to insert (add).",
                    },
                    "locant": {
                        "type": "string",
                        "description": "Position for insertion (add only).",
                    },
                },
                "required": ["name", "operation"],
            },
        },
        {
            "name": "validate_name",
            "description": (
                "Check whether an IUPAC name is valid by attempting to "
                "resolve it to a molecular structure.  Returns the canonical "
                "SMILES if valid.  Use this to verify names before generating "
                "structures.\n\n"
                "Example:\n"
                '  "2-chloro-3-(trifluoromethyl)pyridine" -> valid, SMILES\n'
                '  "2-chloro-99-methylpyridine" -> invalid\n'
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "IUPAC name to validate.",
                    },
                },
                "required": ["name"],
            },
        },
        {
            "name": "name_to_structure",
            "description": (
                "Convert a validated chemical name to a structure file "
                "(CDXML for ChemDraw, or SMILES/MOL).  This is the final "
                "step: call this after assembling and validating the name.\n\n"
                "Output formats:\n"
                '  - "cdxml": ChemDraw XML (requires ChemScript)\n'
                '  - "smiles": canonical SMILES string\n'
                '  - "mol": MDL MOL block with 2D coordinates\n'
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Chemical name to convert.",
                    },
                    "output_format": {
                        "type": "string",
                        "enum": ["cdxml", "smiles", "mol"],
                        "description": "Output format (default: cdxml).",
                    },
                },
                "required": ["name"],
            },
        },
        {
            "name": "enumerate_names",
            "description": (
                "List alternative IUPAC name forms for a molecule.  "
                "Given a name or SMILES, returns the canonical name plus "
                "alternative forms that express the same molecule using "
                "different parent rings/chains and substituent prefixes.\n\n"
                "IMPORTANT: Call this BEFORE modify_name when doing name "
                "surgery on functional groups that appear as suffixes in "
                "the canonical name (ketones '-one', alcohols '-ol', "
                "amines '-amine', acids '-oic acid', etc.).  The "
                "alternatives expose these groups as swappable prefixes.\n\n"
                "Example:\n"
                "  '1-(4-bromophenyl)ethan-1-one'  (ketone as suffix)\n"
                "  -> alternatives include '1-acetyl-4-bromobenzene'\n"
                "     where the ketone is now the prefix 'acetyl'\n"
                "  -> you can then swap 'acetyl' for another prefix\n\n"
                "Each name form includes a 'prefixes' list showing "
                "which substituent prefixes are visible and swappable.\n"
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "identifier": {
                        "type": "string",
                        "description": (
                            "Chemical name, SMILES, abbreviation, or any "
                            "identifier.  Will be resolved to a structure."
                        ),
                    },
                },
                "required": ["identifier"],
            },
        },
        # --- Layer 3: Graph manipulation tools ---
        {
            "name": "list_reactions",
            "description": (
                "List available named reaction templates.  Returns the "
                "name, description, number of reactants, and typical "
                "conditions for each reaction.  Call this to find the "
                "right template before using apply_reaction.\n\n"
                "Categories:\n"
                "  - 'coupling': Suzuki, Buchwald, SNAr, amide, "
                "Sonogashira, Heck, N-alkylation\n"
                "  - 'functional_group': nitro reduction, ester hydrolysis, "
                "alcohol oxidation, reductive amination, Grignard\n"
                "  - 'heterocycle_formation': ~60 ring-forming reactions "
                "including Huisgen triazole, Fischer indole, Paal-Knorr "
                "pyrrole, Hantzsch pyridine/thiazole, benzimidazole, "
                "benzoxazole, Pictet-Spengler, Biginelli, and many more\n\n"
                "Use the optional category filter to narrow results.\n"
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "category": {
                        "type": "string",
                        "enum": [
                            "coupling",
                            "functional_group",
                            "heterocycle_formation",
                        ],
                        "description": (
                            "Optional: filter by category. "
                            "Omit to list all reactions."
                        ),
                    },
                },
                "required": [],
            },
        },
        {
            "name": "apply_reaction",
            "description": (
                "Apply a named reaction template to a substrate molecule.  "
                "For two-component reactions (e.g. Suzuki, Buchwald), "
                "provide both substrate and reagent SMILES.  For single-"
                "component reactions (e.g. nitro reduction), only the "
                "substrate is needed.\n\n"
                "This tool covers ~70 reactions including:\n"
                "  - Classic couplings (Suzuki, Buchwald, Heck, etc.)\n"
                "  - Functional group transforms (reductions, oxidations)\n"
                "  - Ring-forming heterocyclic reactions (Fischer indole, "
                "Huisgen triazole, Paal-Knorr pyrrole, Hantzsch thiazole, "
                "benzimidazole synthesis, Pictet-Spengler, etc.)\n\n"
                "Use list_reactions() to find the right template name.\n\n"
                "The substrate and reagent can be SMILES strings or "
                "chemical names/abbreviations (they will be resolved "
                "automatically).\n\n"
                "Returns the product SMILES, IUPAC name, and suggested "
                "reaction conditions.\n\n"
                "Examples:\n"
                '  - apply_reaction("nitro_reduction", '
                '"c1ccc([N+](=O)[O-])cc1")\n'
                '  - apply_reaction("suzuki_coupling", '
                '"c1ccc(Br)cc1", "c1ccc(B(O)O)cc1")\n'
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "reaction_name": {
                        "type": "string",
                        "description": (
                            "Reaction template name from list_reactions "
                            "(e.g. 'suzuki_coupling', 'nitro_reduction')."
                        ),
                    },
                    "substrate": {
                        "type": "string",
                        "description": (
                            "SMILES or name of the main substrate."
                        ),
                    },
                    "reagent": {
                        "type": "string",
                        "description": (
                            "SMILES or name of the coupling partner "
                            "(for 2-reactant reactions only)."
                        ),
                    },
                },
                "required": ["reaction_name", "substrate"],
            },
        },
        {
            "name": "deprotect",
            "description": (
                "Remove common protecting groups from a molecule.  "
                "Uses 25 built-in deprotection templates covering:\n"
                "  Boc, Fmoc, Cbz (amines)\n"
                "  TBS/TBDMS, THP, Bn, Ac, PMB, TMS (alcohols)\n"
                "  Acetal/Ketal (carbonyls)\n\n"
                "Accepts SMILES or a chemical name.  Returns the "
                "deprotected product and which PGs were removed.\n\n"
                "Example:\n"
                '  deprotect("O=C(OC(C)(C)C)Nc1ccccc1")  # Boc-aniline\n'
                '  -> product: aniline, removed: [Boc]\n'
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "smiles": {
                        "type": "string",
                        "description": (
                            "SMILES or chemical name of the "
                            "protected molecule."
                        ),
                    },
                },
                "required": ["smiles"],
            },
        },
        # --- Reaction JSON summary ---
        {
            "name": "reaction_summary",
            "description": (
                "Load a reaction JSON file and return a slim summary "
                "with only the fields you need.  Use this instead of "
                "reading the full JSON, which contains bulky geometry "
                "and mass data.\n\n"
                "Default fields (per species): id, name, role, "
                "role_detail, smiles, display_text, formula, mw.\n"
                "Default top-level: experiment, conditions.\n"
                "Default eln_data: product_yield, reaction_type.\n\n"
                "Request additional fields by name when needed:\n"
                "  - LCMS: add species fields ['exact_mass', 'adducts']\n"
                "  - Procedure: add species fields ['csv_mass', "
                "'csv_equiv', 'csv_volume'] and eln fields "
                "['procedure_plain', 'product_obtained', 'sm_mass']\n"
                "  - Scheme drawing: defaults are sufficient\n"
                "  - Pass ['*'] to any field list for all fields.\n\n"
                "Available species fields:\n"
                "  id, name, role, role_detail, smiles, smiles_neutral, "
                "classification_method, is_sm, is_dp, is_substrate, "
                "is_solvent, exact_mass, exact_mass_full, mw, formula, "
                "adducts, source, source_id, csv_equiv, csv_mass, "
                "csv_name, csv_volume, csv_supplier, display_text, "
                "original_geometry\n\n"
                "Available top-level fields:\n"
                "  version, experiment, input_files, reaction_smiles, "
                "reaction_class, reaction_name, "
                "classification_confidence, warnings, metadata, "
                "conditions\n\n"
                "Available eln_data fields:\n"
                "  sm_mass, product_obtained, product_yield, "
                "procedure_text, procedure_plain, reaction_type, "
                "start_date, labbook_name, solvents, solvent_details\n"
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "json_path": {
                        "type": "string",
                        "description": "Path to the reaction JSON file.",
                    },
                    "species_fields": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": (
                            "Per-species fields to include. Omit for "
                            "defaults. Pass ['*'] for all fields."
                        ),
                    },
                    "top_fields": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": (
                            "Top-level fields to include. Omit for "
                            "defaults. Pass ['*'] for all fields."
                        ),
                    },
                    "eln_fields": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": (
                            "eln_data sub-fields to include. Omit for "
                            "defaults. Pass ['*'] for all. Pass [] to "
                            "omit eln_data entirely."
                        ),
                    },
                },
                "required": ["json_path"],
            },
        },
    ]
