"""cdxml-toolkit MCP server — 15 chemistry tools via Model Context Protocol.

Tools:
    resolve_name          — Chemical name/formula/CAS → rich molecule JSON
    modify_molecule       — Analyze or transform a molecule (6 operations)
    draw_molecule         — Single molecule → CDXML document
    render_scheme         — YAML/compact text/reaction JSON → publication CDXML
    parse_reaction        — CDXML/CDX/CSV/RXN → semantic reaction JSON
    summarize_reaction    — Reaction JSON → compact context-efficient summary
    extract_structures    — Image/PDF → SMILES + bboxes via DECIMER
    parse_scheme          — CDXML scheme → structured species/steps/topology JSON
    convert_cdx_cdxml     — Bidirectional CDX ↔ CDXML file conversion
    parse_analysis_file   — LCMS/NMR PDF → peaks and data
    format_lab_entry      — Entry dicts → formatted lab book text
    extract_cdxml_from_office — PPTX/DOCX → embedded CDXML files
    embed_cdxml_in_office — CDXML → editable OLE object in PPTX/DOCX
    search_compound       — SMILES → exact/similar matches across experiments
    render_to_png         — CDXML → PNG via ChemDraw COM

Run:
    python -m cdxml_toolkit.mcp_server                    # stdio (default)
    python -m cdxml_toolkit.mcp_server --transport http    # streamable-http
    python -m cdxml_toolkit.mcp_server --transport http --port 8080
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from mcp.server.fastmcp import FastMCP

mcp = FastMCP(
    "cdxml-toolkit",
    instructions=(
        "Chemistry toolkit for parsing, rendering, and reasoning about reaction schemes. "
        "Core workflow: parse_reaction → summarize_reaction (inspect) → render_scheme. "
        "For structure work: resolve_name → modify_molecule → draw_molecule → render_scheme. "
        "Convention: atom-contributing species (reactants, products, key reagents) get drawn "
        "structures; text-only labels are for reagents/conditions. Use (IUPAC name + MW) as "
        "compound identifiers, not LLM-generated SMILES — always resolve names first."
    ),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _validate_file(path: str, label: str) -> Path:
    """Resolve and validate that *path* exists and is a file."""
    p = Path(path).resolve()
    if not p.is_file():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p


# ---------------------------------------------------------------------------
# Tool 1: resolve_name
# ---------------------------------------------------------------------------

@mcp.tool()
def resolve_name(query: str, use_network: bool = True) -> dict:
    """Resolve any chemical identifier to a rich molecule descriptor.

    Converts a name, abbreviation, condensed formula, or CAS number into a
    structured molecule dict with SMILES, formula, MW, exact mass, IUPAC name,
    reagent role, and display text. Uses a 4-tier resolution chain: curated
    reagent DB → condensed formula parser → ChemScript IUPAC → PubChem.

    Do NOT hand-construct SMILES — use this tool instead. The returned dict can
    be passed directly to modify_molecule, draw_molecule, or used to build a
    render_scheme input.

    Args:
        query: Chemical identifier — common name, IUPAC name, abbreviation,
               condensed formula (e.g. "PhB(OH)2"), or CAS number. Examples:
               "aspirin", "Cs2CO3", "2-chloropyridine", "534-17-8", "Et3N".
        use_network: Allow PubChem lookup (requires internet). Default True.

    Returns:
        Dict with keys: ok, name, smiles, formula, mw, exact_mass, iupac_name,
        source (which tier resolved it), role (if in reagent DB), display_text,
        prefix_form (IUPAC substituent prefix, if applicable).
        Returns {ok: False, error: "..."} if unresolvable.
    """
    from cdxml_toolkit.naming.mol_builder import resolve_compound

    return resolve_compound(query, use_network=use_network)


# ---------------------------------------------------------------------------
# Tool 2: modify_molecule
# ---------------------------------------------------------------------------

@mcp.tool()
def modify_molecule(
    mol_json: dict,
    operation: str,
    add: Optional[List[dict]] = None,
    remove: Optional[List[str]] = None,
    new_smiles: Optional[str] = None,
    new_name: Optional[str] = None,
    reaction_name: Optional[str] = None,
    reagent: Optional[dict] = None,
    smarts: Optional[str] = None,
    description: Optional[str] = None,
) -> dict:
    """Analyze or transform a molecule with structural verification.

    Takes a molecule dict (at minimum {"smiles": "..."}) and applies one of 6
    operations. Returns the modified molecule with a structural diff so you can
    verify the change was correct. This is the primary molecular editing tool.

    Do NOT set new_smiles to a hallucinated SMILES — only use "set_smiles" if
    you have a validated SMILES from resolve_name or apply_reaction output.

    Operations:
        "analyze"      — Inspect without modifying: functional groups, IUPAC
                         names, formula, MW, prefix form. No extra kwargs needed.
        "name_surgery" — Modify via IUPAC name: add/remove substituents.
                         Pass add=[{"locant": "2", "prefix": "fluoro"}] and/or
                         remove=["methyl"] kwargs.
        "smarts"       — Apply a SMARTS reaction transform. Pass smarts=
                         "reaction SMARTS" (e.g. "[c:1][F]>>[c:1][Cl]") or
                         reaction_name= from list_reactions output.
        "set_smiles"   — Accept a validated SMILES. Pass new_smiles= (validated
                         by RDKit) and optional description= for context.
        "set_name"     — Set the display name. Pass new_name=.
        "reaction"     — Apply a named template from list_reactions. Pass
                         reaction_name= and optionally reagent={"smiles": ...}
                         for binary reactions (coupling, etc.).

    Args:
        mol_json: Source molecule dict with at least {"smiles": "..."}.
        operation: One of: "analyze", "name_surgery", "smarts", "set_smiles",
                   "set_name", "reaction".
        add: For "name_surgery" — list of {"locant": str, "prefix": str} dicts.
        remove: For "name_surgery" — list of prefix strings to remove.
        new_smiles: For "set_smiles" — validated SMILES string.
        new_name: For "set_name" — new display name string.
        reaction_name: For "smarts"/"reaction" — template name.
        reagent: For "reaction" — dict with "smiles" key for the second reagent.
        smarts: For "smarts" — reaction SMARTS string.
        description: For "set_smiles" — optional context note.

    Returns:
        For "analyze": ok, input_smiles, canonical_name, alternative_names,
        functional_groups, prefix_form, bracket_tree, formula, mw.
        For modifications: ok, input_smiles, output_smiles, input_name,
        output_name, aligned_names, diff (atoms_added, atoms_removed,
        atoms_changed, mcs_smarts, delta_formula, delta_mw), formula, mw.
    """
    from cdxml_toolkit.naming.mol_builder import modify_molecule as _modify

    kwargs: Dict[str, Any] = {}
    if add is not None:
        kwargs["add"] = add
    if remove is not None:
        kwargs["remove"] = remove
    if new_smiles is not None:
        kwargs["new_smiles"] = new_smiles
    if new_name is not None:
        kwargs["new_name"] = new_name
    if reaction_name is not None:
        kwargs["reaction_name"] = reaction_name
    if reagent is not None:
        kwargs["reagent"] = reagent
    if smarts is not None:
        kwargs["smarts"] = smarts
    if description is not None:
        kwargs["description"] = description

    return _modify(mol_json, operation, **kwargs)


# ---------------------------------------------------------------------------
# Tool 3: draw_molecule
# ---------------------------------------------------------------------------

@mcp.tool()
def draw_molecule(mol_json: dict, output_path: Optional[str] = None) -> dict:
    """Render a single molecule to a standalone CDXML document.

    Takes a molecule dict (at minimum {"smiles": "..."}) and generates a
    self-contained CDXML string with 2D coordinates in ACS Document 1996 style
    (BondLength=14.40, Arial 10pt). Optionally places a text label below the
    structure using the "label", "name", or "iupac_name" field (in that order).

    Useful for quick structure visualisation before adding to a scheme, or for
    generating a CDXML snippet to open directly in ChemDraw.

    Args:
        mol_json: Molecule dict with at least {"smiles": "..."}. Optional display
                  fields: "label" (used verbatim), "name", "iupac_name".
                  Typically the output of resolve_name or modify_molecule.
        output_path: If given, also write CDXML to this file path.

    Returns:
        Dict with keys: ok, cdxml (CDXML document string), and output_path if
        a path was specified. Returns {ok: False, error: "..."} on failure.
    """
    from cdxml_toolkit.naming.mol_builder import draw_molecule as _draw

    return _draw(mol_json, output_path=output_path)


# ---------------------------------------------------------------------------
# Tool 4: render_scheme
# ---------------------------------------------------------------------------

@mcp.tool()
def render_scheme(
    yaml_text: Optional[str] = None,
    compact_text: Optional[str] = None,
    json_path: Optional[str] = None,
    layout: str = "auto",
) -> str:
    """Render a chemical reaction scheme to publication-ready CDXML.

    Accepts exactly ONE input mode:
    1. yaml_text    — YAML scheme descriptor string (from write_scheme_yaml,
                      possibly hand-edited). Full layout control.
    2. compact_text — Compact DSL string, like "Mermaid for reactions":
                      "ArBr + Amine --> Product (72%)"
                      "  above: Amine"
                      "  below: \"Pd2(dba)3\", \"BINAP\", \"toluene, 110C\""
    3. json_path    — Path to a reaction JSON file from parse_reaction.
                      Auto-generates YAML layout internally.

    Output is ACS Document 1996 style (BondLength=14.40, Arial 10pt).

    Layout conventions (important for multi-step schemes):
    - Atom-contributing species (reactants, key intermediates, products) go
      on the CENTER LINE as substrates/products — drawn as structures.
    - Additional reagents, coupling partners, and small molecules go ABOVE
      or BELOW the arrow as text or small structures.
    - In sequential multi-step schemes, each step should have only ONE main
      substrate on the center line.  This allows the renderer to share
      intermediates between steps (no duplication).  If you put two
      substrates on the center line (e.g. acid + amine both as substrates),
      the intermediate cannot be shared and will be drawn twice.
    - Example: for an amide coupling of intermediate + amine, put the
      intermediate as the sole substrate and the amine in above_arrow.

    Args:
        yaml_text:    YAML scheme descriptor string.
        compact_text: Compact DSL syntax string.
        json_path:    Path to a reaction JSON file.
        layout:       Layout mode for json_path auto-generation:
                      "auto", "landscape", or "portrait".

    Returns:
        Complete CDXML document string (opens in ChemDraw).
    """
    from cdxml_toolkit.render.renderer import render

    modes = sum(x is not None for x in [yaml_text, compact_text, json_path])
    if modes == 0:
        raise ValueError("Provide exactly one of: yaml_text, compact_text, or json_path.")
    if modes > 1:
        raise ValueError("Provide only ONE of: yaml_text, compact_text, or json_path.")

    if yaml_text is not None:
        from cdxml_toolkit.render.parser import parse_yaml

        scheme = parse_yaml(yaml_text)
        return render(scheme)

    if compact_text is not None:
        from cdxml_toolkit.render.compact_parser import parse_compact

        scheme = parse_compact(compact_text)
        return render(scheme)

    # json_path mode — auto-generate YAML then render
    import yaml

    from cdxml_toolkit.render.parser import parse_yaml
    from cdxml_toolkit.render.scheme_yaml_writer import build_scheme_yaml_dict

    p = _validate_file(json_path, "JSON file")
    yaml_dict = build_scheme_yaml_dict(str(p), layout=layout, include_run_arrows=True)
    yaml_str = yaml.dump(yaml_dict, default_flow_style=False, allow_unicode=True, sort_keys=False)
    scheme = parse_yaml(yaml_str)
    return render(scheme)


# ---------------------------------------------------------------------------
# Tool 5: parse_reaction
# ---------------------------------------------------------------------------

@mcp.tool()
def parse_reaction(
    cdxml: Optional[str] = None,
    cdx: Optional[str] = None,
    csv: Optional[str] = None,
    rxn: Optional[str] = None,
    input_dir: Optional[str] = None,
) -> dict:
    """Parse reaction files into a semantic JSON descriptor.

    Extracts every species with canonical SMILES, role classification (using
    Schneider fingerprint scoring for reactant/reagent binary, plus curated
    database for semantic roles like base/solvent/catalyst), display names,
    equivalents, mass data, and adducts. Produces a single JSON source of truth
    suitable for summarize_reaction, render_scheme, or LCMS analysis.

    Provide at least one file path. Multiple may be combined (e.g. cdxml + csv)
    to merge structural data with ELN metadata.

    Args:
        cdxml:      Path to a .cdxml reaction file.
        cdx:        Path to a .cdx reaction file (converted internally).
        csv:        Path to a Findmolecule ELN CSV export.
        rxn:        Path to a .rxn file.
        input_dir:  Directory containing experiment files (auto-discovers
                    cdxml/cdx/csv/rxn by experiment ID).

    Returns:
        Reaction descriptor dict with keys: version, experiment, input_files,
        reaction_smiles, reaction_class, species (list with role, smiles,
        formula, mw, etc.), conditions, and eln_data.
    """
    from cdxml_toolkit.perception.reaction_parser import parse_reaction as _parse

    if not any([cdxml, cdx, csv, rxn, input_dir]):
        raise ValueError(
            "Provide at least one input: cdxml, cdx, csv, rxn, or input_dir."
        )

    kwargs: Dict[str, Any] = {}
    if cdxml:
        kwargs["cdxml"] = str(_validate_file(cdxml, "CDXML file"))
    if cdx:
        kwargs["cdx"] = str(_validate_file(cdx, "CDX file"))
    if csv:
        kwargs["csv"] = str(_validate_file(csv, "CSV file"))
    if rxn:
        kwargs["rxn"] = str(_validate_file(rxn, "RXN file"))
    if input_dir:
        p = Path(input_dir).resolve()
        if not p.is_dir():
            raise FileNotFoundError(f"input_dir not found: {p}")
        kwargs["input_dir"] = str(p)

    descriptor = _parse(**kwargs, verbose=False)
    return descriptor.to_dict()


# ---------------------------------------------------------------------------
# Tool 6: summarize_reaction
# ---------------------------------------------------------------------------

@mcp.tool()
def summarize_reaction(
    json_path: str,
    species_fields: Optional[List[str]] = None,
    top_fields: Optional[List[str]] = None,
    eln_fields: Optional[List[str]] = None,
) -> dict:
    """Return a compact, context-efficient view of a reaction JSON file.

    The full reaction JSON can be 3,000+ tokens with geometry data. This tool
    returns only the fields you need for a given task, making it practical for
    LLM reasoning without burning context.

    Default fields (when no arguments given):
      species:   id, name, role, role_detail, smiles, display_text, formula, mw
      top-level: experiment, conditions
      eln_data:  product_yield, reaction_type

    Pass ["*"] for any field set to get all fields (equivalent to loading the
    full JSON). Request specific fields by name for task-focused summaries.

    Args:
        json_path:      Path to a reaction JSON file from parse_reaction.
        species_fields: Species fields to include. Available: id, name, role,
                        role_detail, smiles, smiles_neutral, is_sm, is_dp,
                        is_substrate, is_solvent, exact_mass, exact_mass_full,
                        mw, formula, adducts, source, source_id, csv_equiv,
                        csv_mass, csv_name, csv_volume, csv_supplier,
                        display_text, original_geometry. Pass ["*"] for all.
        top_fields:     Top-level fields. Available: version, experiment,
                        input_files, reaction_smiles, reaction_class,
                        reaction_name, classification_confidence, warnings,
                        metadata, conditions. Pass ["*"] for all.
        eln_fields:     ELN data fields. Available: sm_mass, product_obtained,
                        product_yield, procedure_text, procedure_plain,
                        reaction_type, start_date, labbook_name, solvents,
                        solvent_details. Pass ["*"] for all.

    Returns:
        Compact dict with requested fields for each species and top-level keys.
    """
    from cdxml_toolkit.perception.reaction_parser import reaction_summary

    p = _validate_file(json_path, "JSON file")
    return reaction_summary(
        str(p),
        species_fields=species_fields,
        top_fields=top_fields,
        eln_fields=eln_fields,
    )


# ---------------------------------------------------------------------------
# Tool 7: extract_structures_from_image
# ---------------------------------------------------------------------------

@mcp.tool()
def extract_structures_from_image(
    image_path: str,
    detect_labels: bool = True,
) -> dict:
    """Extract chemical structures from an image using DECIMER.

    Takes a PNG, JPG, or PDF image and returns SMILES + confidence scores +
    bounding boxes for every detected chemical structure. Segments the image
    into individual structure regions automatically. Optionally detects nearby
    text labels via OCR.

    DECIMER models download on first run (~570 MB to ~/.data/DECIMER-V2/).
    Requires: DECIMER, opencv-python, and optionally pytesseract/easyocr.

    The returned SMILES should be passed through resolve_name or modify_molecule
    to verify and enrich — DECIMER SMILES may not be canonical and can have
    low confidence for complex structures.

    Args:
        image_path:    Path to PNG, JPG, or PDF file.
        detect_labels: Attempt OCR detection of text labels near structures.
                       Requires pytesseract or easyocr; labels are null without
                       an OCR library. Default True.

    Returns:
        Dict with keys: ok, image_path, structures (list of: smiles, confidence
        in [0,1], bbox [x0,y0,x1,y1], label or null). Returns {ok: False,
        error: "..."} if DECIMER is not installed or extraction fails.
    """
    try:
        from cdxml_toolkit.image.structure_from_image import (
            extract_structures_from_image as _extract,
        )
    except ImportError as e:
        return {
            "ok": False,
            "error": (
                f"DECIMER/OpenCV not available: {e}. "
                "Install with: pip install cdxml-toolkit[decimer]"
            ),
        }

    p = _validate_file(image_path, "Image file")
    return _extract(str(p), detect_labels=detect_labels)


# ---------------------------------------------------------------------------
# Tool 8: parse_scheme
# ---------------------------------------------------------------------------

@mcp.tool()
def parse_scheme(cdxml_path: str) -> dict:
    """Parse a CDXML reaction scheme into a structured description.

    Reads a CDXML file containing a reaction scheme (single- or multi-step) and
    returns a structured JSON with a species registry, reaction graph, topology
    classification, and a natural language narrative suitable for LLM reasoning.

    Uses two strategies in order: step-attribute path (reads <scheme><step>
    attributes if present) then geometry-based fallback (spatial arrow detection).
    Text labels near arrows are classified as "chemical", "condition_ref",
    "footnote", "yield", "compound_label", "citation", or "bioactivity".

    Args:
        cdxml_path: Path to a CDXML file containing a reaction scheme.

    Returns:
        Dict with keys: source_file, species (dict of species records with
        smiles, name, formula, mw, role, text_category), steps (list with
        reactant/product/reagent species IDs and conditions), topology
        (linear/parallel/convergent/divergent), content_type, narrative
        (human-readable summary), and optionally sub_schemes for multi-panel
        files.
    """
    from cdxml_toolkit.perception.scheme_reader import read_scheme

    p = _validate_file(cdxml_path, "CDXML file")
    desc = read_scheme(str(p))

    # SchemeDescription is a dataclass — convert to dict
    try:
        from dataclasses import asdict
        result = asdict(desc)
    except Exception:
        # Fallback: use to_dict if available
        if hasattr(desc, "to_dict"):
            result = desc.to_dict()
        else:
            result = {
                "source_file": getattr(desc, "source_file", str(p)),
                "species": {k: v if isinstance(v, dict) else vars(v)
                            for k, v in (desc.species or {}).items()},
                "steps": [s if isinstance(s, dict) else vars(s)
                          for s in (desc.steps or [])],
                "topology": getattr(desc, "topology", None),
                "content_type": getattr(desc, "content_type", None),
                "narrative": getattr(desc, "narrative", None),
            }
    return result


# ---------------------------------------------------------------------------
# Tool 9: convert_cdx_cdxml
# ---------------------------------------------------------------------------

@mcp.tool()
def convert_cdx_cdxml(
    input_path: str,
    output_path: Optional[str] = None,
) -> dict:
    """Convert bidirectionally between CDX (binary) and CDXML (XML) formats.

    Direction is detected from the file extension:
    - .cdx input  → .cdxml output
    - .cdxml input → .cdx output

    Uses available backends in order: ChemDraw COM (best fidelity, Windows) →
    pycdxml (pure Python, partial support) → OpenBabel.

    ChemDraw COM requires ChemDraw to be installed and closed before running.

    Args:
        input_path:  Path to .cdx or .cdxml file.
        output_path: Output file path. If not given, same directory as input
                     with the swapped extension (e.g. foo.cdx → foo.cdxml).

    Returns:
        Dict with keys: ok, input, output (absolute path to written file).
        Returns {ok: False, error: "..."} if conversion fails.
    """
    from cdxml_toolkit.chemdraw.cdx_converter import convert_file

    p = _validate_file(input_path, "Input file")
    ext = p.suffix.lower()
    if ext not in (".cdx", ".cdxml"):
        return {
            "ok": False,
            "error": f"Unsupported extension: {ext}. Use .cdx or .cdxml.",
        }

    try:
        out = convert_file(str(p), output_path=output_path)
        return {"ok": True, "input": str(p), "output": str(Path(out).resolve())}
    except Exception as e:
        return {"ok": False, "error": str(e), "input": str(p)}


# ---------------------------------------------------------------------------
# Tool 10: parse_analysis_file
# ---------------------------------------------------------------------------

@mcp.tool()
def parse_analysis_file(pdf_path: str) -> dict:
    """Parse an LCMS or NMR analysis PDF to extract peaks and data.

    Supports Waters LCMS reports and MestReNova NMR PDFs. Returns structured
    peak data for LCMS species identification or NMR characterisation.

    This module is under active development. If unavailable, the tool returns
    a graceful error rather than crashing.

    Args:
        pdf_path: Path to an LCMS or NMR PDF report.

    Returns:
        For LCMS: dict with retention_times, peak_areas, masses, UV traces.
        For NMR: dict with chemical_shifts, multiplicities, integrations.
        Returns {ok: False, error: "..."} if module unavailable or parse fails.
    """
    try:
        from cdxml_toolkit.analysis.parse_analysis_file import parse_analysis_file as _parse
    except ImportError as e:
        return {
            "ok": False,
            "error": f"parse_analysis_file module not available: {e}",
        }

    p = _validate_file(pdf_path, "PDF file")
    try:
        result = _parse(str(p))
        if isinstance(result, dict):
            return result
        return {"ok": True, "data": result}
    except Exception as e:
        return {"ok": False, "error": str(e), "pdf_path": str(p)}


# ---------------------------------------------------------------------------
# Tool 11: format_lab_entry
# ---------------------------------------------------------------------------

@mcp.tool()
def format_lab_entry(entries_json: Any) -> dict:
    """Format a list of entry dicts into a structured lab book text entry.

    Takes a list of typed entry dicts (or a JSON string) and produces a
    formatted lab book entry.  The tool re-parses LCMS PDFs to fill in
    exact numbers — you only provide peak identifications (name, approximate
    RT, ion) as search keys.

    IMPORTANT: Do NOT write free-form LCMS text.  Use the structured entry
    types below.  The tool will look up the actual RT, area%, m/z, and UV
    from the PDF.

    Entry types and their required fields:

      {"type": "text", "content": "Procedure paragraph or section header..."}

      {"type": "lcms-species",
       "file": "path/to/report.pdf",
       "label": "t = 0 min",
       "peaks": [
         {"name": "Product",  "rt": 1.02, "ion": {"mode": "ES-", "mz": 444.1}},
         {"name": "SM",       "rt": 0.65, "ion": {"mode": "ES+", "mz": 275.1}},
         {"name": "TPPO",     "rt": 1.02, "ion": {"mode": "ES+", "mz": 279.1}}
       ]}

      {"type": "lcms-areas",
       "file": "path/to/report.pdf",
       "label": "t = 10 min",
       "peaks": [
         {"name": "Product",     "rt": 1.03, "compound_related": true},
         {"name": "Byproduct",   "rt": 1.26, "compound_related": false}
       ]}

      {"type": "lcms-species",
       "file": "path/to/report.pdf",
       "label": "Purified product",
       "peaks": [
         {"name": "Product", "rt": 1.01, "ion": {"mode": "ES-", "mz": 444.2},
          "purity": true, "detector": "220nm"}
       ]}

      {"type": "lcms-manual",
       "file": "path/to/manual_integration.pdf",
       "label": "Manual LC",
       "peaks": [
         {"name": "Product", "rt": 1.01, "compound_related": true}
       ]}

      {"type": "nmr", "content": "1H NMR (400 MHz, DMSO-d6): ..."}

    Workflow: First call parse_analysis_file on each PDF to see peaks/masses.
    Then build entries referencing those PDFs with approximate RT and ion as
    search keys.  This tool re-reads the PDF and fills in exact numbers.

    Args:
        entries_json: List of entry dicts, or a JSON string, or {"entries": [...]}.

    Returns:
        Dict with keys: ok, text (formatted lab book entry string).
    """
    from cdxml_toolkit.analysis.format_procedure_entry import process_entries

    # Accept a JSON string, a list, or a dict with "entries" key
    if isinstance(entries_json, str):
        try:
            entries_json = json.loads(entries_json)
        except json.JSONDecodeError as e:
            return {"ok": False, "error": f"Invalid JSON: {e}"}

    if isinstance(entries_json, dict):
        entries = entries_json.get("entries", [])
    elif isinstance(entries_json, list):
        entries = entries_json
    else:
        return {
            "ok": False,
            "error": "entries_json must be a list, a JSON string, or {\"entries\": [...]}",
        }

    try:
        text = process_entries(entries)
        return {"ok": True, "text": text}
    except Exception as e:
        return {"ok": False, "error": str(e)}


# ---------------------------------------------------------------------------
# Tool 12: extract_cdxml_from_office
# ---------------------------------------------------------------------------

@mcp.tool()
def extract_cdxml_from_office(
    file_path: str,
    output_dir: Optional[str] = None,
) -> dict:
    """Extract embedded ChemDraw objects from a PPTX or DOCX file.

    Office files (PPTX/DOCX) are ZIP archives that may contain ChemDraw OLE
    objects as binary blobs. This tool extracts every ChemDraw object, converts
    it to CDXML, and writes the files to output_dir.

    Requires: olefile. CDX→CDXML conversion uses available backends.

    Args:
        file_path:  Path to a .pptx or .docx file.
        output_dir: Directory for extracted CDXML files. Default: a folder
                    named "<basename>_chemdraw/" next to the input file.

    Returns:
        Dict with keys: ok, input, output_dir, objects (list of: source_path,
        cdxml_output, cdx_output, error for each extracted object).
        Returns {ok: False, error: "..."} if extraction fails entirely.
    """
    from cdxml_toolkit.office.ole_extractor import extract_from_office, ExtractedObject

    p = _validate_file(file_path, "Office file")

    try:
        results = extract_from_office(
            str(p),
            output_dir=output_dir,
            output_format="cdxml",
        )
    except Exception as e:
        return {"ok": False, "error": str(e), "input": str(p)}

    objects = []
    for obj in results:
        entry: Dict[str, Any] = {"source_path": obj.source_path}
        if obj.cdxml_output:
            entry["cdxml_output"] = obj.cdxml_output
        if obj.cdx_output:
            entry["cdx_output"] = obj.cdx_output
        if obj.error:
            entry["error"] = obj.error
        objects.append(entry)

    resolved_output_dir = output_dir
    if results:
        # derive from first result's output path
        first_out = results[0].cdxml_output or results[0].cdx_output
        if first_out:
            resolved_output_dir = str(Path(first_out).parent)

    return {
        "ok": True,
        "input": str(p),
        "output_dir": resolved_output_dir,
        "count": len(objects),
        "objects": objects,
    }


# ---------------------------------------------------------------------------
# Tool 13: embed_cdxml_in_office
# ---------------------------------------------------------------------------

@mcp.tool()
def embed_cdxml_in_office(
    cdxml_path: str,
    office_path: str,
    output_path: Optional[str] = None,
) -> dict:
    """Embed a CDXML file as an editable ChemDraw OLE object in PPTX or DOCX.

    Converts CDXML → CDX + EMF preview via ChemDraw COM, builds a CFB OLE
    compound file, and injects it into a PPTX slide or DOCX paragraph as a
    double-clickable, editable ChemDraw object.

    Requires: ChemDraw COM (Windows, ChemDraw 16+), python-pptx or python-docx.
    ChemDraw must be installed and closed before calling this tool.

    The output format (.pptx or .docx) is detected from office_path extension.
    If office_path does not exist, a new file is created.

    Args:
        cdxml_path:   Path to the CDXML file to embed.
        office_path:  Path to the target .pptx or .docx file. Created if it
                      does not exist.
        output_path:  Output file path. If not given, writes to office_path
                      (modifies in place via temp file).

    Returns:
        Dict with keys: ok, input_cdxml, output (absolute path to written
        Office file), format ("pptx" or "docx"), num_objects_embedded.
        Returns {ok: False, error: "..."} if embedding fails.
    """
    from cdxml_toolkit.office.ole_embedder import (
        batch_convert,
        get_cdxml_content_size,
        build_ole_compound_file,
        build_pptx,
        build_docx,
    )

    cdxml_p = _validate_file(cdxml_path, "CDXML file")
    ext = Path(office_path).suffix.lower()
    if ext not in (".pptx", ".docx"):
        return {
            "ok": False,
            "error": f"Unsupported office format: {ext}. Use .pptx or .docx.",
        }

    if output_path is None:
        output_path = office_path

    try:
        # Step 1: convert CDXML → CDX + EMF via ChemDraw COM
        converted = batch_convert([str(cdxml_p)])
        if not converted:
            return {"ok": False, "error": "ChemDraw COM conversion returned no output."}
        item = converted[0]

        # Step 2: compute display dimensions + build OLE compound file
        w_emu, h_emu = get_cdxml_content_size(str(cdxml_p))
        ole_data = build_ole_compound_file(item["cdx_data"])

        items = [{
            "ole_data": ole_data,
            "emf_data": item["emf_data"],
            "width_emu": w_emu,
            "height_emu": h_emu,
            "name": item["name"],
        }]

        # Step 3: build/inject into Office file
        fmt = ext.lstrip(".")
        if ext == ".pptx":
            build_pptx(items, output_path)
        else:
            build_docx(items, output_path)

        return {
            "ok": True,
            "input_cdxml": str(cdxml_p),
            "output": str(Path(output_path).resolve()),
            "format": fmt,
            "num_objects_embedded": 1,
        }

    except ImportError as e:
        return {
            "ok": False,
            "error": (
                f"ChemDraw COM or Office library not available: {e}. "
                "Requires pywin32 (ChemDraw COM) and python-pptx/python-docx."
            ),
        }
    except Exception as e:
        return {"ok": False, "error": str(e)}


# ---------------------------------------------------------------------------
# Tool 14: search_compound
# ---------------------------------------------------------------------------

@mcp.tool()
def search_compound(
    smiles: str,
    experiment_dir: str,
    similarity_threshold: float = 0.85,
) -> dict:
    """Search for a compound across experiment JSON files by SMILES similarity.

    Scans a directory of reaction JSON files (from parse_reaction) and returns
    exact matches and structurally similar compounds above the given Tanimoto
    threshold. Useful for finding related experiments, checking if a compound
    has been made before, or tracing a compound through a multi-step synthesis.

    This module is under active development. If unavailable, the tool returns
    a graceful error rather than crashing.

    Args:
        smiles:               SMILES string of the compound to search for.
                              Use resolve_name to get a validated SMILES first.
        experiment_dir:       Directory containing reaction JSON files to search.
        similarity_threshold: Tanimoto similarity cutoff (0–1). Default 0.85.

    Returns:
        Dict with keys: ok, query_smiles, exact_matches (list), similar_matches
        (list with similarity scores), total_files_searched.
        Returns {ok: False, error: "..."} if module unavailable or search fails.
    """
    try:
        from cdxml_toolkit.perception.compound_search import search_compound as _search
    except ImportError as e:
        return {
            "ok": False,
            "error": f"compound_search module not available: {e}",
        }

    d = Path(experiment_dir).resolve()
    if not d.is_dir():
        return {"ok": False, "error": f"experiment_dir not found: {d}"}

    try:
        result = _search(
            smiles,
            str(d),
            similarity_threshold=similarity_threshold,
        )
        if isinstance(result, dict):
            return result
        return {"ok": True, "data": result}
    except Exception as e:
        return {"ok": False, "error": str(e)}


# ---------------------------------------------------------------------------
# Tool 15: render_to_png
# ---------------------------------------------------------------------------

@mcp.tool()
def render_to_png(
    cdxml_path: str,
    output_path: Optional[str] = None,
) -> dict:
    """Render a CDXML file to PNG using ChemDraw COM.

    Uses ChemDraw's native rendering engine (via COM automation) at 300 DPI
    with a solid white background. ChemDraw must be installed (Professional
    16+) and closed before calling this tool.

    This tool uses ChemDraw COM exclusively — no RDKit fallback. For a
    quick preview without ChemDraw, use draw_molecule which returns CDXML
    that can be opened directly.

    Args:
        cdxml_path:  Path to the CDXML file to render.
        output_path: Output PNG path. If not given, writes to the same
                     directory as the input with a .png extension.

    Returns:
        Dict with keys: ok, input, output (absolute path to PNG file).
        Returns {ok: False, error: "..."} if ChemDraw COM is unavailable
        or rendering fails.
    """
    try:
        from cdxml_toolkit.chemdraw.cdxml_to_image import cdxml_to_image
    except ImportError as e:
        return {
            "ok": False,
            "error": (
                f"ChemDraw COM not available: {e}. "
                "Requires pywin32 and ChemDraw Professional 16+ on Windows."
            ),
        }

    p = _validate_file(cdxml_path, "CDXML file")

    try:
        out = cdxml_to_image(str(p), output_path=output_path, png_dpi=300)
        return {"ok": True, "input": str(p), "output": str(Path(out).resolve())}
    except Exception as e:
        return {"ok": False, "error": str(e), "input": str(p)}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="cdxml-toolkit MCP server",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--transport",
        choices=["stdio", "http"],
        default="stdio",
        help="MCP transport mode (default: stdio)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port for HTTP transport (default: 8000)",
    )
    args = parser.parse_args()

    if args.transport == "stdio":
        mcp.run(transport="stdio")
    else:
        mcp.run(transport="streamable-http", host="0.0.0.0", port=args.port)


if __name__ == "__main__":
    main()
