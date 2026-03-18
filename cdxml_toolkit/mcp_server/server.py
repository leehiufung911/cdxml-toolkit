"""cdxml-toolkit MCP server — exposes core chemistry tools via Model Context Protocol.

Tools:
    convert_cdx      — CDX → CDXML conversion (pycdxml)
    parse_reaction   — CDXML/CDX/CSV/RXN → semantic reaction JSON
    write_scheme_yaml — reaction JSON(s) → YAML scheme descriptor
    render_scheme    — YAML/compact/JSON → publication-ready CDXML
    resolve_cas      — CAS or compound name → {name, cas, mw, formula, smiles}
    lookup_reagent   — abbreviation or SMILES → {display, role, smiles, tier}

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
from typing import Optional

from mcp.server.fastmcp import FastMCP

mcp = FastMCP(
    "cdxml-toolkit",
    instructions=(
        "Chemistry toolkit for reaction scheme parsing, rendering, and lookup. "
        "Typical workflow: parse_reaction → write_scheme_yaml (inspect/edit) → render_scheme. "
        "For quick rendering without inspection, pass json_path directly to render_scheme."
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
# Tool 1: convert_cdx
# ---------------------------------------------------------------------------

@mcp.tool()
def convert_cdx(cdx_path: str) -> str:
    """Convert a binary CDX file to CDXML (XML) format using pycdxml.

    Args:
        cdx_path: Path to a .cdx file on disk.

    Returns:
        CDXML string (complete XML document).
    """
    from pycdxml import cdxml_converter

    from cdxml_toolkit.chemdraw.cdx_converter import sanitise_cdxml

    p = _validate_file(cdx_path, "CDX file")
    doc = cdxml_converter.read_cdx(str(p))
    raw_cdxml = doc.to_cdxml()
    return sanitise_cdxml(raw_cdxml)


# ---------------------------------------------------------------------------
# Tool 2: parse_reaction
# ---------------------------------------------------------------------------

@mcp.tool()
def parse_reaction(
    cdxml: Optional[str] = None,
    cdx: Optional[str] = None,
    csv: Optional[str] = None,
    rxn: Optional[str] = None,
) -> dict:
    """Parse reaction files into a semantic JSON descriptor.

    Provide at least one file path. Multiple may be combined (e.g. cdxml + csv).
    Reagent role classification uses Schneider fingerprints + curated DB.

    Args:
        cdxml: Path to a .cdxml reaction file.
        cdx:   Path to a .cdx reaction file (converted internally).
        csv:   Path to a Findmolecule ELN CSV export.
        rxn:   Path to a .rxn file.

    Returns:
        Reaction descriptor dict with keys: version, experiment, input_files,
        reaction_smiles, reaction_class, reactants, reagents, products, etc.
    """
    from cdxml_toolkit.perception.reaction_parser import parse_reaction as _parse

    if not any([cdxml, cdx, csv, rxn]):
        raise ValueError("Provide at least one input file (cdxml, cdx, csv, or rxn).")

    kwargs = {}
    if cdxml:
        kwargs["cdxml"] = str(_validate_file(cdxml, "CDXML file"))
    if cdx:
        kwargs["cdx"] = str(_validate_file(cdx, "CDX file"))
    if csv:
        kwargs["csv"] = str(_validate_file(csv, "CSV file"))
    if rxn:
        kwargs["rxn"] = str(_validate_file(rxn, "RXN file"))

    descriptor = _parse(**kwargs, verbose=False)
    return descriptor.to_dict()


# ---------------------------------------------------------------------------
# Tool 3: write_scheme_yaml
# ---------------------------------------------------------------------------

@mcp.tool()
def write_scheme_yaml(
    json_paths: list[str],
    layout: str = "auto",
    include_run_arrows: bool = True,
    use_eln_labels: bool = False,
) -> str:
    """Build a YAML scheme descriptor from one or more reaction JSON files.

    For multiple reactions, auto-detects relationships:
    - Parallel (same SM + DP + shared reagent) → stacked run arrows
    - Sequential (product of A = SM of B) → multi-step chain
    - Unrelated → independent stacked rows

    The returned YAML can be inspected and edited before passing to render_scheme.

    Args:
        json_paths:        List of reaction JSON file paths (from parse_reaction output).
        layout:            Layout mode: "auto", "landscape", or "portrait".
        include_run_arrows: Include run arrows showing individual reaction conditions.
        use_eln_labels:    Use ELN experiment names instead of sequential numbers.

    Returns:
        YAML text string (scheme descriptor).
    """
    import yaml

    from cdxml_toolkit.render.scheme_yaml_writer import (
        build_merged_scheme_yaml_dict,
        build_scheme_yaml_dict,
    )

    validated = [str(_validate_file(p, f"JSON file [{i}]")) for i, p in enumerate(json_paths)]

    if len(validated) == 1:
        yaml_dict = build_scheme_yaml_dict(
            validated[0],
            layout=layout,
            include_run_arrows=include_run_arrows,
            use_eln_labels=use_eln_labels,
        )
    else:
        yaml_dict = build_merged_scheme_yaml_dict(
            validated,
            layout=layout,
            include_run_arrows=include_run_arrows,
            use_eln_labels=use_eln_labels,
        )

    return yaml.dump(yaml_dict, default_flow_style=False, allow_unicode=True, sort_keys=False)


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
    1. yaml_text   — YAML string (from write_scheme_yaml, possibly edited).
    2. compact_text — compact DSL syntax string (hand-authored).
    3. json_path   — path to a reaction JSON (auto-generates YAML internally).

    Args:
        yaml_text:    YAML scheme descriptor string.
        compact_text: Compact syntax string (e.g. "ArBr + Amine --> Product (72%)").
        json_path:    Path to a reaction JSON file.
        layout:       Layout mode for json_path auto-generation: "auto", "landscape", "portrait".

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
    from cdxml_toolkit.render.parser import parse_yaml
    from cdxml_toolkit.render.scheme_yaml_writer import build_scheme_yaml_dict

    import yaml

    p = _validate_file(json_path, "JSON file")
    yaml_dict = build_scheme_yaml_dict(str(p), layout=layout, include_run_arrows=True)
    yaml_str = yaml.dump(yaml_dict, default_flow_style=False, allow_unicode=True, sort_keys=False)
    scheme = parse_yaml(yaml_str)
    return render(scheme)


# ---------------------------------------------------------------------------
# Tool 5: resolve_cas
# ---------------------------------------------------------------------------

@mcp.tool()
def resolve_cas(identifier: str) -> dict:
    """Look up a compound by CAS number or name via PubChem.

    Args:
        identifier: CAS number (e.g. "534-17-8") or compound name (e.g. "cesium carbonate").

    Returns:
        Dict with keys: cas, name, mw, formula, smiles, isomeric_smiles.
        Returns {error: "..."} if lookup fails.
    """
    from cdxml_toolkit.resolve.cas_resolver import resolve_cas as _resolve

    result = _resolve(identifier)
    if result is None:
        return {"error": f"Could not resolve: {identifier}"}
    # Drop coords_2d if present (large, not useful for MCP)
    result.pop("coords_2d", None)
    return result


# ---------------------------------------------------------------------------
# Tool 6: lookup_reagent
# ---------------------------------------------------------------------------

@mcp.tool()
def lookup_reagent(query: str) -> dict:
    """Look up a reagent abbreviation or SMILES in the curated reagent database.

    Tries name lookup first, then SMILES lookup if the query contains chemistry
    characters. Covers common reagents like "DIAD", "Cs2CO3", "TEA", "NaBH4", etc.

    Args:
        query: Reagent abbreviation (e.g. "cs2co3") or SMILES string.

    Returns:
        Dict with keys: display, role, smiles, tier.
        Returns {error: "not found"} if no match.
    """
    from cdxml_toolkit.resolve.reagent_db import get_reagent_db

    db = get_reagent_db()

    # Try name-based lookup first
    entry = db.entry_for_name(query)
    if entry is not None:
        return {
            "display": entry.get("display", query),
            "role": entry.get("role"),
            "smiles": entry.get("smiles"),
            "tier": entry.get("tier", 1),
        }

    # Try SMILES-based lookup
    entry = db.entry_for_smiles(query)
    if entry is not None:
        return {
            "display": entry.get("display", query),
            "role": entry.get("role"),
            "smiles": entry.get("smiles"),
            "tier": entry.get("tier", 1),
        }

    return {"error": "not found", "query": query}


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
