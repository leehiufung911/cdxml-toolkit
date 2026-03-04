#!/usr/bin/env python3
"""
reactant_heuristic.py — Classify reaction reagents as atom-contributing
or non-contributing using role lookup + RDKit MCS.

Two input modes:
  cdxml   Parse a CDXML reaction file; extract fragments + text from <step>
  smiles  Accept reagent SMILES + product SMILES directly on the CLI

Examples:
  python reactant_heuristic.py cdxml -i reaction.cdxml --pretty
  python reactant_heuristic.py smiles --reagents "C1COCCN1" "c1cc2scnc2Br" \\
         --product "c1cc2scnc2N1CCOCC1" --pretty
"""

import argparse
import json
import os
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple
from xml.etree import ElementTree as ET

from .constants import CDXML_FOOTER, CDXML_MINIMAL_HEADER


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ReagentInfo:
    """Information about a single reagent being classified."""
    source_id: str = ""
    source_type: str = ""          # "fragment", "text", "smiles_input"
    name: Optional[str] = None
    smiles: Optional[str] = None
    position: str = ""             # "reactant", "above_arrow", "below_arrow"
    classification: str = ""       # "atom_contributing", "non_contributing", "unclassified"
    classification_method: str = ""  # "role_lookup", "mcs", "rxnmapper", "fm_type", "metal_check", "inorganic_check"
    mcs_ratio: Optional[float] = None
    rxnmapper_confidence: Optional[float] = None  # RXNMapper mapping confidence (0-1)
    role: Optional[str] = None     # "catalyst", "ligand", "base", "solvent", etc.


# ---------------------------------------------------------------------------
# Role Lookup  (Tier 1) — via shared reagent database
# ---------------------------------------------------------------------------

from .reagent_db import get_reagent_db

# Transition metals commonly used as catalysts (by atomic number)
CATALYST_METALS = {46, 28, 29, 77, 45, 44, 78, 76, 79}
# Pd=46, Ni=28, Cu=29, Ir=77, Rh=45, Ru=44, Pt=78, Os=76, Au=79


# ---------------------------------------------------------------------------
# CDXML Parsing Helpers  (adapted from cdxml_combiner.py)
# ---------------------------------------------------------------------------

def _get_page(root: ET.Element) -> ET.Element:
    """Find the <page> element in a CDXML root."""
    page = root.find("page")
    if page is None:
        raise SystemExit("ERROR: no <page> element in CDXML")
    return page


def _count_heavy_atoms(frag: ET.Element) -> int:
    """Count non-hydrogen atoms in a fragment."""
    count = 0
    for n in frag.iter("n"):
        if n.get("NodeType") in ("ExternalConnectionPoint", "Fragment",
                                  "Unspecified"):
            continue
        count += 1
    return count


def _get_text_content(el: ET.Element) -> str:
    """Extract concatenated text from all <s> children of a <t> element."""
    parts = []
    for s in el.iter("s"):
        if s.text:
            parts.append(s.text.strip())
    return " ".join(parts).strip()


def _get_fm_molecule_type(el: ET.Element) -> Optional[int]:
    """Read the Findmolecule MOLECULE TYPE objecttag.
    Values: 0=molecule, 1=solvent, 2=condition text, 3=product."""
    for ot in el.iter("objecttag"):
        if ot.get("Name") == "FM MOLECULE TYPE":
            try:
                return int(ot.get("Value", ""))
            except ValueError:
                return None
    return None


def _attrs_to_str(el: ET.Element) -> str:
    parts = []
    for k, v in el.attrib.items():
        v = v.replace("&", "&amp;").replace('"', "&quot;").replace("<", "&lt;")
        parts.append(f'{k}="{v}"')
    return " ".join(parts)


def _element_to_string(el: ET.Element) -> str:
    tag = el.tag
    attrs = _attrs_to_str(el)
    children = list(el)
    text = el.text or ""
    if attrs:
        open_tag = f"<{tag} {attrs}"
    else:
        open_tag = f"<{tag}"
    if not children and not text.strip():
        return f"{open_tag}/>"
    result = f"{open_tag}>"
    if text.strip():
        safe = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        result += safe
    for child in children:
        result += _element_to_string(child)
    result += f"</{tag}>"
    return result


def _fragment_to_cdxml(frag: ET.Element) -> str:
    """Wrap a single <fragment> in a minimal CDXML document."""
    return (
        CDXML_MINIMAL_HEADER + "\n<page id=\"1\">\n"
        + _element_to_string(frag)
        + "\n</page>\n" + CDXML_FOOTER
    )


# ---------------------------------------------------------------------------
# SMILES Extraction
# ---------------------------------------------------------------------------

# Lazy ChemScript singleton
_cs_instance = None
_cs_tried = False


def _get_chemscript():
    """Return a ChemScriptBridge instance (lazy singleton), or None."""
    global _cs_instance, _cs_tried
    if _cs_tried:
        return _cs_instance
    _cs_tried = True
    try:
        from .chemscript_bridge import ChemScriptBridge
        _cs_instance = ChemScriptBridge()
    except Exception as e:
        print(f"  [warn] ChemScript not available: {e}", file=sys.stderr)
        _cs_instance = None
    return _cs_instance


def _fragment_to_smiles(frag: ET.Element) -> Optional[str]:
    """Convert a CDXML <fragment> to SMILES via ChemScript."""
    cs = _get_chemscript()
    if cs is None:
        return None
    cdxml_str = _fragment_to_cdxml(frag)
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".cdxml", mode="w",
                                          delete=False, encoding="utf-8") as f:
            f.write(cdxml_str)
            tmp_path = f.name
        smiles = cs.write_data(tmp_path, "smiles")
        return smiles.strip() if smiles else None
    except Exception as e:
        print(f"  [warn] fragment→SMILES failed: {e}", file=sys.stderr)
        return None
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.unlink(tmp_path)


def _text_to_smiles(text_content: str) -> Optional[str]:
    """Resolve a reagent name to SMILES.

    Resolution chain (first success wins):
      1. py2opsin name->SMILES (offline, handles IUPAC/systematic names)
      2. PubChem name->SMILES via cas_resolver (online, fallback)
    """
    # --- Try OPSIN first (offline) ---
    smiles = _opsin_name_to_smiles(text_content)
    if smiles:
        return smiles

    # --- Fall back to PubChem (online) ---
    try:
        from .cas_resolver import resolve_name_to_smiles
        return resolve_name_to_smiles(text_content)
    except Exception as e:
        print(f"  [warn] name->SMILES failed for '{text_content}': {e}",
              file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# OPSIN name resolution (offline)
# ---------------------------------------------------------------------------

_opsin_available: Optional[bool] = None


def _opsin_name_to_smiles(name: str) -> Optional[str]:
    """Try to resolve a chemical name to SMILES via OPSIN (offline).

    OPSIN handles systematic/IUPAC names and many common names well
    (e.g. "cesium carbonate", "triethylamine", "sodium tert-butoxide").
    Fails on abbreviations (BINAP, Pd2dba3) and some organometallics.

    Requires Java (JRE) on PATH and the py2opsin package.
    """
    global _opsin_available
    if _opsin_available is False:
        return None
    try:
        import warnings
        from py2opsin import py2opsin
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            result = py2opsin(name)
        if result:
            _opsin_available = True
            print(f"  OPSIN: '{name}' -> {result}", file=sys.stderr)
            return result
        _opsin_available = True
        return None
    except FileNotFoundError:
        # Java not found
        if _opsin_available is None:
            print("  [info] OPSIN unavailable (Java not on PATH)", file=sys.stderr)
        _opsin_available = False
        return None
    except ImportError:
        if _opsin_available is None:
            print("  [info] py2opsin not installed", file=sys.stderr)
        _opsin_available = False
        return None
    except Exception as e:
        print(f"  [info] OPSIN name->SMILES failed for '{name}': {e}",
              file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Tier 1 — Role Lookup
# ---------------------------------------------------------------------------

def _contains_catalyst_metal(smiles: str) -> bool:
    """Check if a molecule contains a transition-metal catalyst atom."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return any(a.GetAtomicNum() in CATALYST_METALS for a in mol.GetAtoms())
    except Exception:
        return False


def _is_inorganic(smiles: str) -> bool:
    """Heuristic: molecule has no carbons, or only 1 C with ≥4 heavy atoms
    (likely carbonate, cyanide, etc.)."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        total = mol.GetNumHeavyAtoms()
        if carbons == 0:
            return True
        if carbons == 1 and total >= 4:
            return True
        return False
    except Exception:
        return False


def role_lookup(smiles: Optional[str], name: Optional[str]
                ) -> Optional[Tuple[str, str]]:
    """Tier 1 classification.  Returns (role, method) or None."""
    db = get_reagent_db()

    # 1. SMILES-based lookup (exact canonical match)
    if smiles:
        role = db.role_for_smiles(smiles)
        if role:
            return (role, "role_lookup")

    # 1b. Stereo-agnostic SMILES lookup.  RDKit-only SMILES extraction
    #     often omits E/Z on double bonds (e.g. DEAD's N=N) because
    #     frag_to_mol doesn't set bond stereo from 2D coordinates.
    if smiles:
        role = _role_for_smiles_no_stereo(smiles, db)
        if role:
            return (role, "role_lookup_no_stereo")

    # 2. Name-based lookup
    if name:
        role = db.role_for_name(name)
        if role:
            return (role, "role_lookup")

    # 3. Metal-containing → catalyst
    if smiles and _contains_catalyst_metal(smiles):
        return ("catalyst", "metal_check")

    # 4. Inorganic salt
    if smiles and _is_inorganic(smiles):
        return ("inorganic_salt", "inorganic_check")

    return None


def _role_for_smiles_no_stereo(smiles: str, db) -> Optional[str]:
    """Match SMILES against DB after stripping stereochemistry."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        Chem.RemoveStereochemistry(mol)
        flat_smi = Chem.MolToSmiles(mol)

        for smi_key, entry in db._by_smiles.items():
            mol2 = Chem.MolFromSmiles(smi_key)
            if mol2 is None:
                continue
            Chem.RemoveStereochemistry(mol2)
            if flat_smi == Chem.MolToSmiles(mol2):
                return entry.get("role")
    except ImportError:
        pass
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Tier 2 — RDKit MCS
# ---------------------------------------------------------------------------

def mcs_ratio(reagent_smiles: str, product_smiles: str) -> Optional[float]:
    """Compute MCS heavy-atom ratio: MCS_atoms / reagent_heavy_atoms."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS

        reagent_mol = Chem.MolFromSmiles(reagent_smiles)
        product_mol = Chem.MolFromSmiles(product_smiles)
        if reagent_mol is None or product_mol is None:
            return None

        reagent_heavy = reagent_mol.GetNumHeavyAtoms()
        if reagent_heavy == 0:
            return None

        result = rdFMCS.FindMCS(
            [reagent_mol, product_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            timeout=10,
        )

        if result.canceled or result.numAtoms == 0:
            return 0.0

        return result.numAtoms / reagent_heavy

    except Exception as e:
        print(f"  [warn] MCS failed: {e}", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Main Classification Logic
# ---------------------------------------------------------------------------

def classify_reagents(reagents: List[ReagentInfo],
                      product_smiles: str,
                      mcs_threshold: float = 0.3,
                      use_rxnmapper: bool = True) -> List[ReagentInfo]:
    """Classify each reagent using the tiered strategy.

    Tiers (applied in order):
      1.  Role lookup — reagent_db name/SMILES match, metal check, inorganic check
      1.5 RXNMapper atom mapping (if use_rxnmapper=True and rxn-experiments env available)
      2.  RDKit MCS ratio (fallback for anything still unclassified)
    """
    # --- Tier 1: role lookup (fast, no subprocess) ---
    for r in reagents:
        # Already classified (e.g. by FM type)?
        if r.classification:
            continue

        result = role_lookup(r.smiles, r.name)
        if result:
            role, method = result
            r.classification = "non_contributing"
            r.classification_method = method
            r.role = role

    # --- Tier 1.5: RXNMapper atom mapping ---
    # For reagents that escaped Tier 1, try ML-based classification.
    # One RXNMapper call covers all unclassified reagents at once.
    if use_rxnmapper:
        _try_rxnmapper_classification(reagents, product_smiles)

    # --- Tier 2: RDKit MCS (fallback for remaining unclassified) ---
    for r in reagents:
        if r.classification:
            continue

        if r.smiles and product_smiles:
            ratio = mcs_ratio(r.smiles, product_smiles)
            if ratio is not None:
                r.mcs_ratio = round(ratio, 3)
                r.classification_method = "mcs"
                r.classification = ("atom_contributing"
                                    if ratio > mcs_threshold
                                    else "non_contributing")
                continue

        # Fallback
        r.classification = "unclassified"
        r.classification_method = "none"

    return reagents


def _try_rxnmapper_classification(reagents: List[ReagentInfo],
                                   product_smiles: str) -> None:
    """Tier 1.5: Use RXNMapper atom maps to classify unclassified reagents.

    Builds a reaction SMILES from all unclassified reagent SMILES + product,
    calls RXNMapper via subprocess (rxn-experiments env), and uses the atom
    map results to determine which reagents are atom-contributing.

    Modifies reagents in place. Silently returns if RXNMapper is unavailable.
    """
    if not product_smiles:
        return

    # Collect unclassified reagents that have SMILES
    unclassified = [r for r in reagents if not r.classification and r.smiles]
    if not unclassified:
        return

    # Build reaction SMILES: all unclassified reagent SMILES >> product
    reactant_smiles_list = [r.smiles for r in unclassified]
    rxn_smi = ".".join(reactant_smiles_list) + ">>" + product_smiles

    # Try to call RXNMapper
    try:
        from experiments.atom_mapping.rxn_atom_mapper import classify_roles
    except ImportError:
        # rxn_atom_mapper not available — skip silently
        return

    try:
        result = classify_roles(rxn_smi)
    except Exception as exc:
        print(f"  [info] RXNMapper classification failed: {exc}",
              file=sys.stderr)
        return

    if result is None:
        return

    confidence = result.get("confidence", 0.0)
    components = result.get("components", [])

    if not components:
        return

    print(f"  RXNMapper classification (confidence={confidence:.4f}):",
          file=sys.stderr)

    # Match results back to reagents by canonical SMILES
    try:
        from rdkit import Chem
        def _canon(smi):
            mol = Chem.MolFromSmiles(smi)
            return Chem.MolToSmiles(mol) if mol else smi
    except ImportError:
        def _canon(smi):
            return smi

    # Build lookup: canonical SMILES → RXNMapper component info
    rxnm_by_smi = {}
    for comp in components:
        canon = _canon(comp["smiles"])
        rxnm_by_smi[canon] = comp

    # Apply to unclassified reagents
    for r in unclassified:
        canon = _canon(r.smiles)
        comp = rxnm_by_smi.get(canon)
        if comp is None:
            continue

        is_contributing = comp.get("atom_contributing")
        if is_contributing is None:
            continue

        if is_contributing:
            r.classification = "atom_contributing"
            r.classification_method = "rxnmapper"
            n_atoms = comp.get("n_product_atoms", 0)
            print(f"    {r.smiles[:50]:50s} → atom_contributing "
                  f"({n_atoms} atoms in product)", file=sys.stderr)
        else:
            r.classification = "non_contributing"
            r.classification_method = "rxnmapper"
            print(f"    {r.smiles[:50]:50s} → non_contributing",
                  file=sys.stderr)

        r.rxnmapper_confidence = confidence


# ---------------------------------------------------------------------------
# CDXML Mode Entry Point
# ---------------------------------------------------------------------------

def classify_from_cdxml(cdxml_path: str,
                        mcs_threshold: float = 0.3,
                        use_rxnmapper: bool = True) -> Dict[str, Any]:
    """Parse a CDXML reaction file and classify all reagents."""
    tree = ET.parse(cdxml_path)
    root = tree.getroot()
    page = _get_page(root)

    # --- Parse <step> metadata ---
    scheme = page.find("scheme")
    step = scheme.find("step") if scheme is not None else None
    if step is None:
        raise SystemExit("ERROR: no <scheme><step> found in CDXML")

    reactant_ids = step.get("ReactionStepReactants", "").split()
    product_ids  = step.get("ReactionStepProducts",  "").split()
    above_ids    = step.get("ReactionStepObjectsAboveArrow", "").split()
    below_ids    = step.get("ReactionStepObjectsBelowArrow", "").split()

    # Build id → element map
    id_to_el: Dict[str, ET.Element] = {}
    for el in page:
        eid = el.get("id", "")
        if eid:
            id_to_el[eid] = el

    # --- Extract product SMILES ---
    product_smiles = None
    for pid in product_ids:
        el = id_to_el.get(pid)
        if el is not None and el.tag == "fragment":
            product_smiles = _fragment_to_smiles(el)
            if product_smiles:
                break
    if not product_smiles:
        raise SystemExit("ERROR: could not extract product SMILES")

    print(f"Product SMILES: {product_smiles}", file=sys.stderr)

    # --- Collect reagents ---
    reagents: List[ReagentInfo] = []
    seen_ids: set = set()

    def _process_element(eid: str, position: str):
        """Process a single element (fragment or text) as a potential reagent."""
        if eid in seen_ids:
            return
        seen_ids.add(eid)

        el = id_to_el.get(eid)
        if el is None:
            return

        fm_type = _get_fm_molecule_type(el)

        # Skip products and condition text
        if fm_type == 3:
            return
        if fm_type == 2:
            return

        ri = ReagentInfo(source_id=eid, position=position)

        # FM type = 1 → solvent (immediate classification)
        if fm_type == 1:
            ri.source_type = el.tag
            ri.classification = "non_contributing"
            ri.classification_method = "fm_type"
            ri.role = "solvent"
            if el.tag == "t":
                ri.name = _get_text_content(el)
            reagents.append(ri)
            return

        # Fragment → extract SMILES via ChemScript
        if el.tag == "fragment":
            ri.source_type = "fragment"
            ri.smiles = _fragment_to_smiles(el)
            if ri.smiles:
                print(f"  Fragment {eid}: {ri.smiles}", file=sys.stderr)

        # Text → resolve name to SMILES via PubChem
        elif el.tag == "t":
            ri.source_type = "text"
            text = _get_text_content(el)
            ri.name = text
            ri.smiles = _text_to_smiles(text)
            if ri.smiles:
                print(f"  Text '{text}' → {ri.smiles}", file=sys.stderr)
            else:
                print(f"  Text '{text}' → no SMILES (name-only)", file=sys.stderr)
        else:
            return

        reagents.append(ri)

    # Process reactants first, then above/below arrow
    for rid in reactant_ids:
        _process_element(rid, "reactant")
    for eid in above_ids:
        _process_element(eid, "above_arrow")
    for eid in below_ids:
        _process_element(eid, "below_arrow")

    # --- Classify ---
    classify_reagents(reagents, product_smiles, mcs_threshold,
                      use_rxnmapper=use_rxnmapper)

    return {
        "cdxml_file": os.path.basename(cdxml_path),
        "product_smiles": product_smiles,
        "reagents": [_reagent_to_dict(r) for r in reagents],
    }


# ---------------------------------------------------------------------------
# SMILES Mode Entry Point
# ---------------------------------------------------------------------------

def classify_from_smiles(reagent_smiles: List[str],
                         product_smiles: str,
                         reagent_names: Optional[List[str]] = None,
                         mcs_threshold: float = 0.3,
                         use_rxnmapper: bool = True) -> Dict[str, Any]:
    """Classify reagents given as SMILES strings."""
    reagents: List[ReagentInfo] = []
    for i, smi in enumerate(reagent_smiles):
        name = reagent_names[i] if reagent_names and i < len(reagent_names) else None
        ri = ReagentInfo(source_type="smiles_input", smiles=smi, name=name)
        reagents.append(ri)
    classify_reagents(reagents, product_smiles, mcs_threshold,
                      use_rxnmapper=use_rxnmapper)
    return {
        "product_smiles": product_smiles,
        "reagents": [_reagent_to_dict(r) for r in reagents],
    }


# ---------------------------------------------------------------------------
# Output Helpers
# ---------------------------------------------------------------------------

def _reagent_to_dict(r: ReagentInfo) -> Dict[str, Any]:
    d = asdict(r)
    # Drop empty/None optional fields for cleaner output
    if d["mcs_ratio"] is None:
        del d["mcs_ratio"]
    if d.get("rxnmapper_confidence") is None:
        d.pop("rxnmapper_confidence", None)
    if d["role"] is None:
        del d["role"]
    if d["name"] is None:
        del d["name"]
    return d


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Classify reaction reagents as atom-contributing "
                    "or non-contributing (role lookup + RDKit MCS).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="mode", required=True,
                                help="Input mode")

    # Shared args for both modes
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--output",
                        help="Output JSON file (default: stdout)")
    common.add_argument("--pretty", action="store_true",
                        help="Pretty-print JSON output")
    common.add_argument("--threshold", type=float, default=0.5,
                        help="MCS ratio threshold (default: 0.5)")

    # CDXML mode
    p_cdxml = sub.add_parser("cdxml", parents=[common],
                              help="Classify from a CDXML reaction file")
    p_cdxml.add_argument("-i", "--input", required=True,
                          help="Input CDXML file")

    # SMILES mode
    p_smi = sub.add_parser("smiles", parents=[common],
                            help="Classify from SMILES strings")
    p_smi.add_argument("--reagents", nargs="+", required=True,
                        help="Reagent SMILES strings")
    p_smi.add_argument("--product", required=True,
                        help="Product SMILES")
    p_smi.add_argument("--names", nargs="+", default=None,
                        help="Reagent names (parallel to --reagents)")

    args = parser.parse_args(argv)

    if args.mode == "cdxml":
        result = classify_from_cdxml(args.input, args.threshold)
    elif args.mode == "smiles":
        result = classify_from_smiles(
            args.reagents, args.product, args.names, args.threshold)
    else:
        parser.print_help()
        return 1

    indent = 2 if args.pretty else None
    json_str = json.dumps(result, indent=indent, ensure_ascii=False)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(json_str + "\n")
        print(f"Written to {args.output}", file=sys.stderr)
    else:
        print(json_str)

    return 0


if __name__ == "__main__":
    sys.exit(main())
