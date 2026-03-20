# cdxml-toolkit — Agent Reference

## What this is

A chemistry office automation toolkit with 15 MCP tools. Agents use these tools to draw molecules, render reaction schemes, parse ELN exports, analyze LCMS/NMR data, complete lab books, and manipulate ChemDraw files in PowerPoint/Word.

**All structure output is CDXML** (ChemDraw XML). Uses ACS Document 1996 style.

## Key rules

1. **Never generate SMILES.** Always call `resolve_name` first. LLM-generated SMILES are unreliable.
2. **Never return large output inline.** Tools write files and return `{ok, output_path, size}`.
3. **CDXML is the interchange format.** Binary CDX files must be converted via `convert_cdx_cdxml`.

## MCP server

```bash
python -m cdxml_toolkit.mcp_server                    # stdio (default)
python -m cdxml_toolkit.mcp_server --transport http    # streamable-http
```

## Tool reference

### resolve_name

Resolve any chemical identifier to a rich molecule descriptor.

```
resolve_name(query="aspirin")
resolve_name(query="Cs2CO3")
resolve_name(query="534-17-8")       # CAS number
resolve_name(query="CF3")            # fragment → prefix form
```

Returns: `{ok, name, smiles, formula, mw, exact_mass, iupac_name, source, role, display_text, prefix_form}`

4-tier chain: curated reagent DB (186 entries) → condensed formula parser → ChemScript IUPAC → PubChem.

### modify_molecule

Analyze or transform a molecule. 6 operations, all returning structural diffs for verification.

```
# Analyze functional groups
modify_molecule(mol_json={smiles:"..."}, operation="analyze")

# IUPAC name surgery (add/remove substituents)
modify_molecule(mol_json={smiles:"..."}, operation="name_surgery",
                add=[{locant:"4", prefix:"fluoro"}])

# SMARTS transform
modify_molecule(mol_json={smiles:"..."}, operation="smarts",
                smarts="[c:1][F]>>[c:1][Cl]")

# Named reaction (162 templates: coupling, deprotection, reduction, ...)
modify_molecule(mol_json={smiles:"..."}, operation="reaction",
                reaction_name="amide_coupling", reagent={smiles:"..."})

# Set validated SMILES (from another tool's output, never hand-written)
modify_molecule(mol_json={smiles:"..."}, operation="set_smiles", new_smiles="...")

# Set display name
modify_molecule(mol_json={smiles:"..."}, operation="set_name", new_name="Product A")
```

Input normalization: `mol_json` accepts a bare SMILES string. Operation names are fuzzy-matched (e.g. `"BOC deprotection"` redirects to `operation="reaction", reaction_name="BOC_deprotection"`). Stringified JSON arrays in `add`/`remove` are auto-parsed.

### draw_molecule

Single molecule to CDXML. Optionally labels with name/iupac_name.

```
draw_molecule(mol_json={smiles:"...", label:"Aspirin"}, output_path="aspirin.cdxml")
```

### render_scheme

Reaction scheme to publication-ready CDXML. Accepts YAML, compact text, or reaction JSON.

```
# From YAML (what agents typically write)
render_scheme(yaml_text="...", output_path="scheme.cdxml")

# From reaction JSON (auto-layout)
render_scheme(json_path="reaction.json", output_path="scheme.cdxml")
```

Call with no arguments to see the full YAML schema reference.

**YAML format:**
```yaml
layout: sequential           # linear | sequential | stacked-rows
structures:
  SM:
    smiles: "CCO"            # from resolve_name output
  Product:
    smiles: "CC=O"
steps:
  - substrates: [SM]
    products: [Product]
    above_arrow:
      structures: [Reagent]  # drawn structures above arrow
      text: ["PCC (1.5 eq)"] # text labels above arrow
    below_arrow:
      text: ["DCM, rt, 2 h"]
```

Convention: ONE substrate on center line per step. Additional reagents go in `above_arrow`.

For stacked-rows (independent reactions), use `sections`:
```yaml
sections:
  - label: "(i)"
    steps: [{substrates: [A], products: [B], ...}]
  - label: "(ii)"
    steps: [{substrates: [C], products: [D], ...}]
```

**Forgiving parser** accepts: `species`/`substrates` as alias for `structures`, `reactants` as alias for `substrates` in steps, inline structures in steps, text as string not list, bare SMILES strings, `above_arrow` as list/string/dict, `reagents` key.

### parse_reaction

Parse ELN export files into semantic reaction JSON.

```
parse_reaction(input_dir="experiments/KL-CC-001/")  # auto-discovers .cdxml/.csv/.rxn
parse_reaction(cdxml="experiment.cdxml", csv="experiment.csv")
parse_reaction(cdxml="experiment.cdxml", output_path="reaction.json")
```

Returns species with canonical SMILES, role classification, display names, equivalents, mass data.

### summarize_reaction

Context-efficient view of a reaction JSON. The full JSON can be 3,000+ tokens; this returns only what you need.

```
summarize_reaction(json_path="reaction.json")
summarize_reaction(json_path="reaction.json", species_fields=["name","smiles","adducts"])
summarize_reaction(json_path="reaction.json", species_fields=["*"])  # all fields
```

Default fields: id, name, role, smiles, display_text, formula, mw.

### extract_structures_from_image

Extract chemical structures from an image via DECIMER neural network.

```
extract_structures_from_image(image_path="screenshot.png")
```

Returns: `{ok, structures: [{smiles, confidence, bbox, label}]}`. Pass returned SMILES through `resolve_name` to validate.

Note: DECIMER runs on CPU (~60s/image). First run downloads ~570 MB of models.

### parse_scheme

Read a CDXML reaction scheme into structured JSON (species, steps, topology, narrative).

```
parse_scheme(cdxml_path="scheme.cdxml")
```

### parse_analysis_file

Parse LCMS (Waters/manual) or NMR (MestReNova) PDFs.

```
parse_analysis_file(pdf_path="lcms_report.pdf")
```

### format_lab_entry

Format structured entries into lab book text. Re-reads LCMS PDFs for exact peak numbers.

```
format_lab_entry(entries_json=[
  {"type": "text", "content": "Added 100 mg SM to 5 mL THF..."},
  {"type": "lcms-species", "file": "report.pdf", "label": "t=0",
   "peaks": [{"name": "Product", "rt": 1.02, "ion": {"mode": "ES+", "mz": 445.1}}]},
  {"type": "nmr", "content": "1H NMR (400 MHz, DMSO-d6): ..."}
])
```

Entry types: `text`, `lcms-species`, `lcms-areas`, `lcms-manual`, `nmr`. Workflow: call `parse_analysis_file` first to see peaks, then build entries referencing those PDFs.

Also accepts shorthand: `{"procedure": "text..."}` or a bare string.

### extract_cdxml_from_office / embed_cdxml_in_office

```
extract_cdxml_from_office(file_path="presentation.pptx")
embed_cdxml_in_office(cdxml_path="scheme.cdxml", office_path="output.pptx")
```

### convert_cdx_cdxml

```
convert_cdx_cdxml(input_path="experiment.cdx")  # → .cdxml (auto-detected)
```

### search_compound

```
search_compound(smiles="c1ccncc1", experiment_dir="experiments/", similarity_threshold=0.85)
```

### render_to_png

```
render_to_png(cdxml_path="scheme.cdxml")  # requires ChemDraw COM
```

## Agent workflow recipes

**"Draw a compound"**
1. `resolve_name` → get SMILES
2. `draw_molecule` → CDXML

**"Draw a reaction scheme I describe"**
1. `resolve_name` each compound
2. `modify_molecule` to compute products (reaction templates, SMARTS, deprotection)
3. `render_scheme` with YAML describing the steps

**"Parse and re-render an ELN export"**
1. `parse_reaction(input_dir=...)` → reaction JSON
2. `render_scheme(json_path=...)` → publication CDXML

**"Extract structures from an image and build a scheme"**
1. `extract_structures_from_image` → SMILES + confidence
2. `modify_molecule` to apply reactions to extracted structures
3. `render_scheme` with YAML

**"Complete a lab book entry"**
1. `parse_reaction(input_dir=...)` → reaction JSON
2. `summarize_reaction` → species + metadata
3. `parse_analysis_file` on each LCMS/NMR PDF
4. `format_lab_entry` with structured entries

**"Search for a compound across experiments"**
1. `resolve_name` (or `extract_structures_from_image`) → SMILES
2. `search_compound` across experiment directory

**"Modify a PowerPoint scheme"**
1. `extract_cdxml_from_office` → embedded CDXML files
2. `parse_scheme` → understand the scheme
3. `modify_molecule` → change structures
4. `render_scheme` → new CDXML
5. `embed_cdxml_in_office` → back into PPTX

## Installation

```bash
pip install -e ".[all]"        # Everything
pip install -e ".[rdkit]"      # Core + RDKit (minimum for rendering)
```

**Required:** `lxml>=4.6`. **Recommended:** `rdkit>=2023.03`.

## Package structure

```
cdxml_toolkit/
├── mcp_server/                 # MCP server (15 tools) — primary interface
│   ├── server.py               # Tool definitions + input normalization
│   ├── __main__.py
│   └── __init__.py
│
├── naming/                     # Agent chemistry reasoning
│   ├── mol_builder.py          # resolve_compound, modify_molecule, draw_molecule
│   ├── reactions_datamol.json  # 162 reaction templates
│   ├── name_decomposer.py     # Structural decomposition
│   └── aligned_namer.py       # Multi-step aligned naming
│
├── perception/                 # Read and understand reactions
│   ├── reaction_parser.py      # ELN exports → JSON (parse_reaction)
│   ├── scheme_reader.py        # CDXML → structured JSON (parse_scheme)
│   ├── compound_search.py      # search_compound across experiments
│   ├── reactant_heuristic.py   # Reagent role classification
│   └── ...
│
├── render/                     # Generate CDXML
│   ├── renderer.py             # Layout engine (5 patterns)
│   ├── parser.py               # YAML → SchemeDescriptor (forgiving parser)
│   ├── compact_parser.py       # Compact text → SchemeDescriptor
│   ├── schema.py               # Dataclass definitions
│   └── ...
│
├── resolve/                    # Chemical name → SMILES (4-tier chain)
│   ├── reagent_db.py           # Tier 1: curated DB (~186 entries)
│   ├── condensed_formula.py    # Tier 2: generative parser
│   ├── cas_resolver.py         # Tier 4: PubChem
│   └── superatom_table.py      # 2,854 fragment abbreviations
│
├── analysis/                   # LCMS/NMR parsing, lab book formatting
│   ├── parse_analysis_file.py  # Unified LCMS + NMR parser
│   ├── format_procedure_entry.py
│   ├── lcms_analyzer.py
│   ├── extract_nmr.py
│   └── deterministic/          # Batch pipeline
│
├── image/                      # Image → SMILES (DECIMER)
│   ├── structure_from_image.py
│   └── reaction_from_image.py
│
├── office/                     # PPTX/DOCX integration (OLE)
│   ├── ole_embedder.py
│   ├── ole_extractor.py
│   └── doc_from_template.py
│
├── chemdraw/                   # ChemDraw COM automation (Windows)
│   ├── cdx_converter.py
│   ├── cdxml_to_image.py
│   └── ...
│
├── layout/                     # Layout tools
│   ├── reaction_cleanup.py     # Pure-Python reaction layout
│   ├── alignment.py            # Structure alignment
│   └── scheme_merger.py        # Multi-scheme merging
│
├── constants.py                # ACS 1996 style values
├── cdxml_utils.py              # CDXML geometry utilities
├── rdkit_utils.py              # CDXML ↔ RDKit bridge
├── text_formatting.py          # Chemical text (subscripts, italics)
└── cdxml_builder.py            # Build CDXML from atom/bond data
```

## YAML parser normalization rules

The forgiving YAML parser in `render/parser.py` accepts these LLM-common mistakes:

1. `species` or `substrates` as alias for `structures` (top-level)
2. `reactants` as alias for `substrates` (in steps)
3. `structures` as list of dicts → converted to keyed mapping
4. Inline structure dicts in `substrates`/`products` → auto-registered
5. Bare SMILES strings in `substrates`/`products` → auto-registered
6. `text` as string not list → wrapped
7. `reagents` key in steps → distributed to above/below arrow
8. `above_arrow` as list of dicts, bare list, or string → normalized to dict
9. Redundant `id` field inside structure defs → silently accepted

## MCP tool input normalization

Applied in `mcp_server/server.py`:

- `mol_json` accepts bare SMILES string → wrapped as `{"smiles": "..."}`
- `add`/`remove`/`species_fields` accept stringified JSON → auto-parsed
- `operation` fuzzy-matches reaction template names (e.g. `"BOC deprotection"` → `reaction` + `BOC_deprotection`)
- `format_lab_entry` accepts `{"procedure": "text"}` shorthand
- `parse_reaction(input_dir=...)` auto-discovers .cdxml/.cdx/.csv/.rxn files

## Layout patterns

| Pattern | Keyword | Description |
|---------|---------|-------------|
| Single-step | `linear` (default) | One arrow, substrates left, products right |
| Multi-step | `sequential` | Linear chain of steps, auto-wraps |
| Fan-out | `divergent` | One SM giving multiple products |
| Independent rows | `stacked-rows` | Multiple unrelated sequences stacked |
| Wrap (serpentine) | `wrap: serpentine` | Zigzag rows with vertical arrows |

## CDXML conventions

- Carbon atoms: no `Element` attribute. Heteroatoms: `Element` number + `NumHydrogens`.
- Coordinates in points (1/72 inch), y-axis down.
- ACS Document 1996 style: BondLength=14.40, ChainAngle=120, Arial 10pt Bold.
- Abbreviation groups report expanded BoundingBox — use `fragment_bbox()` (atom-only) instead.

## Reagent database

Edit `resolve/reagent_abbreviations.json` to add reagents:

```json
{
  "cs2co3": { "display": "Cs2CO3", "role": "base", "smiles": "O=C([O-])[O-].[Cs+].[Cs+]" }
}
```

Roles: catalyst, ligand, base, solvent, coupling_reagent, reducing_agent, oxidant, protecting_group, deprotecting_agent, acid, activating_agent, lewis_acid.

## Running tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## Dependencies

| Dependency | Purpose |
|-----------|---------|
| `lxml` | CDXML parsing (required) |
| `rdkit` | SMILES, 2D coords, MW, MCS |
| `pywin32` | ChemDraw COM (Windows) |
| `pyyaml` | YAML parsing |
| `python-pptx` / `python-docx` / `olefile` | Office files |
| `opencv-python` / `Pillow` | Image processing |
| `mcp` | MCP server |

ChemDraw COM tools require ChemDraw to be installed and **closed** before running.
