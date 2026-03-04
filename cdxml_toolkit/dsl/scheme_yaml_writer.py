"""
scheme_yaml_writer.py — Generate scheme YAML from reaction_parser JSON.

This is the layout-decision layer between perception (reaction_parser) and
rendering (renderer).  It reads a reaction JSON, decides where each species
goes in the scheme, and writes a YAML file that the renderer can consume.

The three decisions made here:
  1. Structure or text?  (atom-contributing → structure; non-contributing → text)
  2. Position?  (substrate → left; other reactant → above arrow; reagents → below)
  3. Priority ordering of below-arrow text (catalyst > base > solvent > conditions)

Usage:
    python experiments/scheme_dsl/scheme_yaml_writer.py reaction.json -o scheme.yaml

    from scheme_dsl.scheme_yaml_writer import write_scheme_yaml
    yaml_path = write_scheme_yaml("reaction.json", "scheme.yaml")
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Role-based priority for below-arrow text ordering
# ---------------------------------------------------------------------------

# Lower number = higher priority (appears first below arrow)
_ROLE_PRIORITY = {
    "catalyst": 10,
    "ligand": 15,
    "base": 20,
    "coupling_reagent": 25,
    "reducing_agent": 30,
    "oxidant": 30,
    "lewis_acid": 30,
    "activating_agent": 30,
    "halogenating_agent": 30,
    "fluorinating_agent": 30,
    "borylating_agent": 30,
    "protecting_group": 35,
    "deprotecting_agent": 35,
    "acid": 35,
    "additive": 40,
    "reductant": 40,
    "reagent": 45,
    "drying_agent": 50,
    "inorganic_salt": 50,
    "solvent": 80,
}

# Roles that should always be shown as text, never as drawn structures
_DEMOTE_ROLES = {
    "base", "catalyst", "ligand", "coupling_reagent", "reducing_agent",
    "oxidant", "protecting_group", "deprotecting_agent", "acid",
    "activating_agent", "lewis_acid", "drying_agent", "halogenating_agent",
    "fluorinating_agent", "borylating_agent", "additive", "reductant",
    "reagent", "solvent", "inorganic_salt",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_scheme_yaml(
    json_path: str,
    output_path: str,
    layout: str = "auto",
    include_run_arrows: bool = True,
) -> str:
    """Read reaction JSON, make layout decisions, write YAML file.

    Parameters
    ----------
    json_path : str
        Path to reaction_parser JSON file.
    output_path : str
        Path where the YAML will be written.
    layout : str
        Layout type: "linear", "sequential", or "auto" (inferred from step count).
    include_run_arrows : bool
        If True and ELN data has SM mass + product yield, include run_arrows.

    Returns
    -------
    str
        The absolute path to the written YAML file.
    """
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    species = data.get("species", [])
    eln_data = data.get("eln_data") or {}
    conditions = data.get("conditions", [])

    # Build the YAML dict
    yaml_dict = _build_yaml_dict(species, conditions, eln_data,
                                  layout=layout,
                                  include_run_arrows=include_run_arrows)

    # Write YAML
    _write_yaml_file(yaml_dict, output_path)

    return os.path.abspath(output_path)


def build_scheme_yaml_dict(
    json_path: str,
    layout: str = "auto",
    include_run_arrows: bool = True,
) -> Dict[str, Any]:
    """Read reaction JSON and return the YAML dict (without writing to disk).

    Useful for programmatic access when you want to inspect or modify
    the dict before writing.
    """
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    species = data.get("species", [])
    eln_data = data.get("eln_data") or {}
    conditions = data.get("conditions", [])

    return _build_yaml_dict(species, conditions, eln_data,
                            layout=layout,
                            include_run_arrows=include_run_arrows)


# ---------------------------------------------------------------------------
# Core layout logic
# ---------------------------------------------------------------------------

def _build_yaml_dict(
    species: List[Dict[str, Any]],
    conditions: List[str],
    eln_data: Dict[str, Any],
    layout: str = "auto",
    include_run_arrows: bool = True,
) -> Dict[str, Any]:
    """Build the complete YAML dict from reaction data.

    This is where all three layout decisions are made:
      1. Structure or text?  (atom_contributing → draw; else → text)
      2. Position?  (substrate left, other reactant above, reagents below)
      3. Priority ordering of below-arrow text
    """
    # --- Classify species into scheme positions ---
    substrates = []        # left of arrow (drawn structures)
    above_structures = []  # above arrow (drawn structures)
    above_text = []        # above arrow (text, e.g. "(1.2 eq)")
    below_text_items = []  # below arrow (text with priority for sorting)
    products = []          # right of arrow (drawn structures)

    # Track species that will be drawn as structures (need StructureRef entries)
    drawn_ids = set()

    for sp in species:
        sp_id = sp.get("id", "")
        role = sp.get("role", "")
        role_detail = (sp.get("role_detail") or "").lower()
        is_sm = sp.get("is_sm", False)
        is_dp = sp.get("is_dp", False)
        is_substrate = sp.get("is_substrate", False)
        is_solvent = sp.get("is_solvent", False)
        atom_contributing = (role == "atom_contributing")
        source = sp.get("source", "")
        # Use raw name for labels (display_text may already have equiv appended)
        name = sp.get("name", "")
        display = sp.get("display_text") or name
        smiles = sp.get("smiles")

        # Check is_dp flag as well as role for product identification
        if role == "product" or is_dp:
            if smiles:
                products.append(sp_id)
                drawn_ids.add(sp_id)

        # Check is_sm/is_substrate flags — but only honor them when the
        # species is atom-contributing or unclassified.  The CSV may mark
        # reagents like n-BuLi as "substrate" even though RXNMapper says
        # they are non-contributing.  Role classification wins.
        elif (is_substrate or is_sm) and (atom_contributing or role in ("", "unclassified")):
            if smiles:
                substrates.append(sp_id)
                drawn_ids.add(sp_id)
            # else: CSV-only SM entry with no structure — skip (the actual
            # structure should be a separate atom_contributing species)

        elif atom_contributing:
            # Atom-contributing species are drawn as structures above arrow
            above_structures.append(sp_id)
            drawn_ids.add(sp_id)
            # Add equiv text below the structure
            equiv = sp.get("csv_equiv")
            if equiv:
                above_text.append(f"({equiv} eq)")

        elif is_solvent:
            # Solvents → below arrow text, low priority
            priority = _ROLE_PRIORITY.get("solvent", 80)
            below_text_items.append((priority, display))

        else:
            # Non-contributing / unclassified species → text below arrow
            # Check if it should be demoted (drawn structure → text)
            should_demote = role_detail in _DEMOTE_ROLES
            has_structure = source in ("fragment", "rxn") and smiles

            if has_structure and not should_demote:
                # Unusual: non-contributing but not a known reagent type.
                # Draw it above the arrow as a structure.
                above_structures.append(sp_id)
                drawn_ids.add(sp_id)
            else:
                # Text label below arrow — use raw name + format equiv here
                priority = _ROLE_PRIORITY.get(role_detail, 50)
                label = _format_label_with_equiv(name, display, sp.get("csv_equiv"))
                if not label:
                    continue
                below_text_items.append((priority, label))

    # If no substrate was identified but there are atom-contributing species
    # above the arrow, promote the most likely one to substrate (left of arrow).
    # Prefer the one with csv_equiv closest to 1.0, else the first one.
    if not substrates and above_structures:
        above_set = set(above_structures)
        best_id = above_structures[0]
        best_diff = float("inf")
        for sp in species:
            sp_id_check = sp.get("id")
            if sp_id_check in above_set:
                eq = sp.get("csv_equiv")
                if eq:
                    try:
                        diff = abs(float(eq) - 1.0)
                        if diff < best_diff:
                            best_diff = diff
                            best_id = sp_id_check
                    except (ValueError, TypeError):
                        pass
        above_structures.remove(best_id)
        substrates.append(best_id)
        # Rebuild above_text from remaining above_structures
        above_text = []
        remaining = set(above_structures)
        for sp in species:
            if sp.get("id") in remaining:
                eq = sp.get("csv_equiv")
                if eq:
                    above_text.append(f"({eq} eq)")

    # Sort below-arrow text by priority
    below_text_items.sort(key=lambda x: x[0])
    below_text = [item[1] for item in below_text_items]

    # Add conditions (temperature, time, atmosphere) at end
    if conditions:
        below_text.extend(_normalize_conditions(conditions))

    # --- Build structures dict ---
    # Assign sequential compound numbers (1, 2, 3, ...) to substrates and
    # products only.  Above-arrow structures (reagents) don't get numbered
    # in chemistry publications — they're shown for clarity, not reference.
    # This is a layout decision: the JSON carries is_sm/is_dp flags but the
    # actual labels ("1", "2", etc.) are determined here.
    label_order = []
    for sid in substrates:
        if sid not in label_order:
            label_order.append(sid)
    for sid in products:
        if sid not in label_order:
            label_order.append(sid)
    label_counter = 1
    label_map: Dict[str, str] = {}
    for sid in label_order:
        label_map[sid] = str(label_counter)
        label_counter += 1

    structures = {}
    for sp in species:
        sp_id = sp.get("id", "")
        if sp_id not in drawn_ids:
            continue
        smiles = sp.get("smiles")
        if not smiles:
            continue
        entry: Dict[str, Any] = {"smiles": smiles}
        # Compound number label below structure (substrates + products only)
        if sp_id in label_map:
            entry["label"] = label_map[sp_id]
        structures[sp_id] = entry

    # --- Build step ---
    step: Dict[str, Any] = {}
    if substrates:
        step["substrates"] = substrates
    if products:
        step["products"] = products

    above: Dict[str, Any] = {}
    if above_structures:
        above["structures"] = above_structures
    if above_text:
        above["text"] = above_text
    if above:
        step["above_arrow"] = above

    if below_text:
        step["below_arrow"] = {"text": below_text}

    # --- Determine layout ---
    if layout == "auto":
        layout = "linear"  # single step for now; multi-step will be "sequential"

    # --- Build top-level YAML dict ---
    yaml_dict: Dict[str, Any] = {
        "structures": structures,
        "steps": [step],
        "layout": layout,
    }

    # --- Run arrows ---
    # Run arrows already display yield in their output label, so only add
    # yield_ to the step when run arrows are NOT present (avoids duplication).
    run_arrows_added = False
    if include_run_arrows and eln_data:
        run_arrows = _build_run_arrows(eln_data)
        if run_arrows:
            yaml_dict["run_arrows"] = run_arrows
            run_arrows_added = True

    if not run_arrows_added:
        yield_pct = eln_data.get("product_yield", "").strip()
        if yield_pct:
            yield_str = yield_pct.rstrip("%").strip()
            step["yield_"] = f"{yield_str}%"

    return yaml_dict


# ---------------------------------------------------------------------------
# Run arrows
# ---------------------------------------------------------------------------

def _build_run_arrows(eln_data: Dict[str, Any]) -> Optional[List[Dict[str, Any]]]:
    """Build run_arrows list from ELN data (SM mass → product yield)."""
    sm_mass = eln_data.get("sm_mass", "").strip()
    product_obtained = eln_data.get("product_obtained", "").strip()
    product_yield = eln_data.get("product_yield", "").strip()

    if not sm_mass or not product_obtained:
        return None

    input_label = sm_mass if _has_unit(sm_mass) else f"{sm_mass} g"
    obtained_str = (product_obtained if _has_unit(product_obtained)
                    else f"{product_obtained} g")
    if product_yield:
        yield_clean = product_yield.rstrip("%").strip()
        output_label = f"{obtained_str}, {yield_clean}% yield"
    else:
        output_label = obtained_str

    return [{
        "step": 1,
        "runs": [{"input": input_label, "output": output_label}],
    }]


def _format_label_with_equiv(
    name: str, display: str, csv_equiv: Optional[str],
) -> str:
    """Build a text label, adding equiv only if not already present.

    Uses raw ``name`` as the base label, falling back to ``display``.
    Avoids duplicating "(X eq)" when ``display_text`` already contains it.
    """
    base = name or display
    if not base:
        return ""
    # If equiv data exists and the label doesn't already mention "eq"
    if csv_equiv and "eq" not in base.lower():
        return f"{base} ({csv_equiv} eq)"
    # If display has equiv but name doesn't, use display as-is
    if "eq" in display.lower():
        return display
    return base


def _normalize_conditions(conditions: List[str]) -> List[str]:
    """Normalize condition strings for display.

    Fixes temperature formatting:
      "80 C"  → "80 °C"
      "105C"  → "105 °C"
      "80°C"  → "80 °C"
      "-78 °C" → "-78 °C" (no change)
    """
    result = []
    for cond in conditions:
        # "80 C" or "-78 C" → "80 °C" (missing degree symbol)
        cond = re.sub(r"(-?\d+\.?\d*)\s+C\b", r"\1 °C", cond)
        # "80C" or "105C" → "80 °C" (no space, no degree)
        cond = re.sub(r"(-?\d+\.?\d*)C\b", r"\1 °C", cond)
        # "80°C" → "80 °C" (no space before degree)
        cond = re.sub(r"(-?\d+\.?\d*)°C", r"\1 °C", cond)
        result.append(cond)
    return result


def _has_unit(value: str) -> bool:
    """Check if a mass string already contains a unit (g, mg, mL, etc.)."""
    return bool(re.search(r"\d\s*(g|mg|kg|mL|µL|L)\b", value))


# ---------------------------------------------------------------------------
# YAML output
# ---------------------------------------------------------------------------

def _write_yaml_file(data: Dict[str, Any], path: str) -> None:
    """Write YAML dict to file.

    Uses PyYAML if available, otherwise writes a simple manual format.
    """
    if yaml is not None:
        with open(path, "w", encoding="utf-8") as f:
            yaml.dump(data, f, default_flow_style=False, allow_unicode=True,
                      sort_keys=False)
    else:
        # Fallback: write JSON with .yaml extension (valid YAML superset)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a scheme YAML from a reaction_parser JSON file.",
    )
    parser.add_argument("json_path", help="Reaction parser JSON file")
    parser.add_argument("-o", "--output", default=None,
                        help="Output YAML path (default: {stem}-scheme.yaml)")
    parser.add_argument("--layout", default="auto",
                        help="Layout: linear, sequential, auto (default: auto)")
    parser.add_argument("--no-run-arrows", action="store_true",
                        help="Suppress run arrows")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.json_path):
        print(f"Error: {args.json_path} not found", file=sys.stderr)
        sys.exit(1)

    output = args.output
    if output is None:
        stem = os.path.splitext(os.path.basename(args.json_path))[0]
        output = os.path.join(
            os.path.dirname(args.json_path) or ".", f"{stem}-scheme.yaml")

    result = write_scheme_yaml(
        args.json_path, output,
        layout=args.layout,
        include_run_arrows=not args.no_run_arrows,
    )

    if args.verbose:
        print(f"Written: {result}", file=sys.stderr)
    print(result)


if __name__ == "__main__":
    main()
