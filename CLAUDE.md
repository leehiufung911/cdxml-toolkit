# cdxml-toolkit — LLM Reference

## What this toolkit does

A Python package for processing ChemDraw CDXML files and generating publication-ready reaction schemes. Designed to be orchestrated by an LLM for chemistry office work: parsing ELN exports, classifying reagents, laying out reactions, and producing CDXML output that opens correctly in ChemDraw.

## Key rules

- **All structure output is CDXML.** Must open in ChemDraw. Uses ACS Document 1996 style (BondLength=14.40, ChainAngle=120, Arial 10pt).
- **Never trust LLM-generated SMILES.** Resolve compound names via PubChem API, CAS lookup, or reagent database. Use existing structure files when available.
- **CDXML is the interchange format.** Tools consume and produce CDXML. Binary CDX files must be converted first via `cdxml-convert`.

## Installation

```bash
pip install -e .                    # Core (lxml only)
pip install -e ".[rdkit]"           # + RDKit (SMILES, 2D coords, MW)
pip install -e ".[chemdraw]"        # + ChemDraw COM automation (Windows)
pip install -e ".[all]"             # Everything
```

**Required:** `lxml>=4.6` (CDXML is XML).
**Strongly recommended:** `rdkit>=2023.03` (needed for scheme rendering, structure alignment, MW calculation).

## CLI tools

| Command | Module | Description |
|---------|--------|-------------|
| `cdxml-parse` | `reaction_parser` | Parse ELN exports (CDXML/CDX/RXN/CSV) into a JSON reaction descriptor |
| `cdxml-render` | `dsl.render_scheme` | Render YAML, compact text, or JSON into a CDXML reaction scheme |
| `cdxml-polish` | `scheme_polisher_v2` | Clean up and enrich an ELN-exported CDXML scheme |
| `cdxml-layout` | `reaction_cleanup` | Re-layout a CDXML reaction (pure Python, no COM) |
| `cdxml-merge` | `scheme_merger` | Merge multiple schemes (parallel, sequential, or auto-detect) |
| `cdxml-convert` | `cdx_converter` | Convert between CDX and CDXML (ChemDraw COM) |
| `cdxml-image` | `cdxml_to_image` | Render CDXML to PNG/SVG (ChemDraw COM) |
| `cdxml-build` | `cdxml_builder` | Build CDXML from atom/bond data |
| `cdxml-ole` | `ole_embedder` | Embed CDXML as editable OLE objects in PPTX/DOCX |

All CLIs use `python -m cdxml_toolkit.<module>` or the installed console script name.

## The JSON-first pipeline

The primary workflow builds publication-ready reaction schemes from ELN export files:

```bash
# 1. Convert CDX to CDXML (needs ChemDraw COM — close ChemDraw first)
cdxml-convert experiment.cdx -o experiment.cdxml

# 2. Parse reaction into semantic JSON
cdxml-parse experiment.cdxml --csv experiment.csv -o reaction.json

# 3. Render JSON into publication-ready CDXML
cdxml-render --from-json reaction.json -o scheme.cdxml

# 4. Render to image for visual inspection (needs ChemDraw COM)
cdxml-image scheme.cdxml
```

**What each step does:**

1. **Reaction parsing** (`reaction_parser`) — extracts every species with canonical SMILES, role classification (3-tier: reagent_db lookup, RXNMapper atom mapping, RDKit MCS), display names, equivalents, mass data, and adducts. Produces a single JSON source of truth.

2. **Layout decisions** (`scheme_yaml_writer`) — determines which species are drawn structures vs text labels, positions them above/below the arrow, sorts by role priority.

3. **CDXML rendering** (`renderer`) — generates 2D coordinates from SMILES via RDKit, computes bounding boxes, sizes the arrow to fit content, formats chemical text (subscripts, italics), and outputs ACS-styled CDXML.

Steps 2-3 require only RDKit. No ChemDraw COM needed for rendering.

## Scheme DSL: three input modes

### From reaction JSON (automated)
```bash
cdxml-render --from-json reaction.json -o scheme.cdxml
```

### From YAML (hand-authored)
```yaml
layout: sequential
structures:
  ArBr:
    smiles: "Brc1ncnc2sccc12"
  Morph:
    smiles: "C1COCCN1"
  Product:
    smiles: "c1nc(N2CCOCC2)c2ccsc2n1"
steps:
  - substrates: [ArBr]
    products: [Product]
    above_arrow:
      structures: [Morph]
    below_arrow:
      text:
        - "Pd2(dba)3 (0.05 eq.)"
        - "rac-BINAP (0.1 eq.)"
        - "Cs2CO3 (2 eq.)"
        - "Dioxane, 105 °C, 24 h"
```
```bash
cdxml-render scheme.yaml -o scheme.cdxml
```

### From compact text ("Mermaid for reactions")
```
ArBr: {Brc1ncnc2sccc12}
Morph: {C1COCCN1}

ArBr + Morph --> Product{c1nc(N2CCOCC2)c2ccsc2n1} (72%)
  above: Morph
  below: "Pd2(dba)3", "BINAP", "Cs2CO3", "toluene, 110 °C, 16 h"
```
```bash
cdxml-render scheme.txt -o scheme.cdxml
```

Full syntax references: `experiments/scheme_dsl/YAML_SYNTAX.md` and `COMPACT_SYNTAX.md`.

## Layout patterns

| Pattern | Keyword | Description |
|---------|---------|-------------|
| Single-step | `linear` (default) | One arrow, substrates left, products right |
| Multi-step | `sequential` | Linear chain of steps, auto-wraps |
| Wrap (repeat) | `wrap: repeat` | Multi-row L→R with repeated structures |
| Wrap (serpentine) | `wrap: serpentine` | Zigzag L→R, R→L with vertical arrows |
| Fan-out | `divergent` | One SM giving multiple products |
| Independent rows | `stacked-rows` | Multiple unrelated sequences stacked |

## Annotations (orthogonal to layout)

- **Run arrows** — SM mass → product yield annotations below the scheme
- **Arrow styles** — `solid` (default), `dashed`, `failed` (X overlay)
- **Compound labels** — numbered labels below structures (e.g. "KL-CC-001")
- **Letter conditions** — `a`, `b`, `c` labels on arrows with legend

## Polisher pipeline (alternative to JSON-first)

Surgically modifies an ELN-exported CDXML rather than building from scratch. Preserves original ChemDraw drawings while cleaning layout, classifying reagents, and injecting ELN data.

```bash
cdxml-polish experiment.cdxml -o polished.cdxml \
    --approach chemdraw_mimic \
    --align-mode rxnmapper \
    --classify-method rxnmapper \
    --eln-csv experiment.csv
```

## Package structure

```
cdxml_toolkit/
├── __init__.py                 # Version, core exports
├── constants.py                # ACS Document 1996 style, layout gaps
├── cdxml_utils.py              # CDXML geometry (bbox, IO, id map)
├── rdkit_utils.py              # CDXML fragment → RDKit Mol/SMILES/MW
├── text_formatting.py          # Chemical subscripts + italic prefixes
├── reagent_db.py               # Two-tier reagent database (172 + 5,837 entries)
├── superatom_table.py          # Abbreviation label → SMILES (~2,850 entries)
├── eln_csv_parser.py           # Findmolecule ELN CSV parser
├── alignment.py                # Structure alignment (Kabsch, MCS, RXNMapper)
├── reaction_cleanup.py         # Pure Python reaction layout (6 approaches)
├── reaction_parser.py          # Reaction → JSON semantic descriptor
├── reactant_heuristic.py       # Reagent role classification
├── scheme_polisher.py          # Scheme polishing (classify, swap, align)
├── scheme_polisher_v2.py       # Full polishing pipeline
├── scheme_merger.py            # Multi-scheme merging (parallel/sequential/auto)
├── eln_enrichment.py           # ELN CSV → scheme annotation
├── cdxml_builder.py            # Build CDXML from atom/bond data
├── coord_normalizer.py         # Coordinate normalization to ACS 1996
├── cas_resolver.py             # CAS/name → SMILES via PubChem
├── cdx_converter.py            # CDX ↔ CDXML (ChemDraw COM)
├── cdxml_to_image.py           # CDXML → PNG/SVG (ChemDraw COM)
├── cdxml_to_image_rdkit.py     # CDXML → image (backup, RDKit only)
├── eln_cdx_cleanup.py          # ELN CDX cleanup
├── chemscript_bridge.py        # ChemScript .NET bridge
├── ole_embedder.py             # CDXML → editable OLE in PPTX/DOCX
├── ole_extractor.py            # Extract ChemDraw from Office files
├── doc_from_template.py        # Fill PPTX/DOCX templates
├── rdf_parser.py               # SciFinder RDF parser
├── reaction_from_image.py      # Image → reaction scheme CDXML
├── structure_from_image.py     # Image → structure CDXML (DECIMER)
├── scheme_aligner.py           # MCS-based orientation alignment
├── scheme_maker.py             # Experimental: JSON → scheme
├── reagent_abbreviations.json  # Tier-1 reagent DB (172 curated entries)
├── chemscanner_abbreviations.json  # Tier-2 ChemScanner DB (5,837 entries)
├── superatom_data.json         # Superatom abbreviation data
└── dsl/                        # Scheme DSL subpackage
    ├── schema.py               # Dataclass definitions
    ├── parser.py               # YAML → SchemeDescriptor
    ├── compact_parser.py       # Compact text → SchemeDescriptor
    ├── renderer.py             # CDXML rendering engine (5 layouts)
    ├── render_scheme.py        # CLI entry point
    ├── auto_layout.py          # Zero-effort JSON → CDXML
    └── scheme_yaml_writer.py   # JSON → YAML layout decisions
```

## Shared modules — what they do

### constants.py
All ACS Document 1996 style values: `ACS_BOND_LENGTH` (14.40), `ACS_CHAIN_ANGLE` (120), font constants, `ACS_STYLE` dict, `CDXML_HEADER`/`CDXML_FOOTER` templates. Also layout gap constants.

### cdxml_utils.py
CDXML geometry utilities. Key functions:
- `fragment_bbox(frag)` — atom-only bounding box (NOT XML BoundingBox — those are unreliable for abbreviation groups)
- `fragment_bbox_with_label_extension(frag)` — extends bbox for hanging N-H/P-H labels
- `parse_cdxml(path)` / `write_cdxml(tree, path)` — IO with DOCTYPE preservation
- `build_id_map(parent)` — `{id: element}` dict

### rdkit_utils.py
CDXML fragment ↔ RDKit bridge. Key functions:
- `frag_to_mol(frag)` / `frag_to_smiles(frag)` / `frag_to_mw(frag)` — convert CDXML fragments
- `cleanup_fragment_rdkit(frag)` — 2D coordinate cleanup with Kabsch orientation preservation

### text_formatting.py
Chemistry-specific text formatting for CDXML `<s>` elements:
- Subscript digits in formulas: "CH3OH" → "CH₃OH"
- Italic IUPAC prefixes: "n-BuLi" → "*n*-BuLi"

### reagent_db.py
Two-tier reagent database. Tier-1 (172 curated entries with roles) always wins. Tier-2 (5,837 ChemScanner entries, no roles) is fallback. Key methods: `display_for_name()`, `role_for_name()`, `display_for_smiles()`, `smiles_role_display()`.

### eln_csv_parser.py
Parses Findmolecule ELN CSV exports (semicolon-delimited, @TYPE sections) into `ExperimentData` dataclass with reactants, solvents, product, and metadata. Pure stdlib, no external dependencies.

## CDXML conventions

- Carbon atoms: no `Element` attribute. Heteroatoms: `Element` number + `NumHydrogens`.
- Coordinates in points (1/72 inch), y-axis down.
- Reactions: `<fragment>` for molecules, `<arrow>` for arrows, `<t>` for text labels.
- `<scheme>` / `<step>` elements hold reaction metadata (`ReactionStepObjectsAboveArrow`, etc.).
- Abbreviation groups: `<n NodeType="Fragment">` — report expanded BoundingBox, not visible size. Always use `fragment_bbox()` (atom-only) instead.

## Reagent database

**To add a reagent:** Edit `reagent_abbreviations.json`. Each entry needs at minimum `display`. Add `role` and `smiles` if known. Alternate spellings go in `aliases`.

```json
{
    "cs2co3": {
        "display": "Cs2CO3",
        "role": "base",
        "smiles": "O=C([O-])[O-].[Cs+].[Cs+]"
    }
}
```

Role categories: catalyst, ligand, base, solvent, coupling_reagent, reducing_agent, oxidant, protecting_group, deprotecting_agent, acid, activating_agent, lewis_acid, and others.

## Dependencies

| Dependency | Required by | Purpose |
|-----------|-------------|---------|
| `lxml` | Core (required) | CDXML XML parsing and writing |
| `rdkit` | DSL renderer, rdkit_utils, alignment, reaction_parser | SMILES, 2D coords, MW, MCS |
| `pywin32` | cdx_converter, cdxml_to_image, ole_embedder | ChemDraw COM automation (Windows) |
| `pyyaml` | dsl.parser, dsl.scheme_yaml_writer | YAML input parsing |
| `python-pptx` / `python-docx` / `olefile` | ole_embedder, ole_extractor, doc_from_template | Office file manipulation |
| `opencv-python` / `Pillow` | structure_from_image, cdxml_to_image | Image processing |

ChemDraw COM tools require ChemDraw to be installed and **closed** before running.

## Running tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

327 tests. All pure Python — no ChemDraw COM required for the test suite.

## Bundled samples

`samples/` contains two real ELN exports with full pipeline output:
- `KL-CC-001/` — Buchwald coupling (source CDX/CSV/RXN + reaction.json + scheme.yaml + scheme.cdxml)
- `KL-CC-002/` — SNAr alkylation (same structure)
- `consolidated/` — Hand-authored two-step sequential YAML combining both reactions
