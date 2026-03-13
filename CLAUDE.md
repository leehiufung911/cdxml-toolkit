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
| `cdxml-parse` | `perception.reaction_parser` | Parse ELN exports (CDXML/CDX/RXN/CSV) into a JSON reaction descriptor |
| `cdxml-render` | `render.render_scheme` | Render YAML, compact text, or JSON into a CDXML reaction scheme |
| `cdxml-polish` | `layout.scheme_polisher_v2` | Clean up and enrich an ELN-exported CDXML scheme |
| `cdxml-layout` | `layout.reaction_cleanup` | Re-layout a CDXML reaction (pure Python, no COM) |
| `cdxml-merge` | `layout.scheme_merger` | Merge multiple schemes (parallel, sequential, or auto-detect) |
| `cdxml-convert` | `chemdraw.cdx_converter` | Convert between CDX and CDXML (ChemDraw COM) |
| `cdxml-image` | `chemdraw.cdxml_to_image` | Render CDXML to PNG/SVG (ChemDraw COM) |
| `cdxml-build` | `cdxml_builder` | Build CDXML from atom/bond data |
| `cdxml-ole` | `office.ole_embedder` | Embed CDXML as editable OLE objects in PPTX/DOCX |

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

1. **Reaction parsing** (`reaction_parser`) — extracts every species with canonical SMILES, role classification (2-tier: Schneider FP scoring for binary reactant/reagent, curated DB for semantic roles), display names, equivalents, mass data, and adducts. Produces a single JSON source of truth.

2. **Layout decisions** (`scheme_yaml_writer`) — determines which species are drawn structures vs text labels, positions them above/below the arrow, sorts by role priority. Also handles **multi-reaction merging**: classifies reaction pairs as parallel/sequential/unrelated via canonical SMILES, clusters with Union-Find, and produces combined scheme YAML with stacked run arrows and range equiv notation.

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
    --align-mode mcs \
    --eln-csv experiment.csv
```

## Package structure

```
cdxml_toolkit/
│
│  # ── Foundation (top level) ─────────────────────────────
├── __init__.py                 # Version, core exports
├── constants.py                # ACS Document 1996 style, layout gaps
├── cdxml_utils.py              # CDXML geometry (bbox, IO, id map)
├── rdkit_utils.py              # CDXML fragment → RDKit Mol/SMILES/MW
├── text_formatting.py          # Chemical subscripts + italic prefixes
├── cdxml_builder.py            # Build CDXML from atom/bond data
├── coord_normalizer.py         # Coordinate normalization to ACS 1996
│
│  # ── Subpackages ────────────────────────────────────────
├── perception/                 # Read and understand reaction schemes
│   ├── scheme_reader.py        # CDXML → SchemeDescription JSON
│   ├── scheme_segmenter.py     # Multi-panel CDXML auto-segmentation
│   ├── spatial_assignment.py   # Geometry-based element → arrow assignment
│   ├── reaction_parser.py      # ELN exports → JSON reaction descriptor
│   ├── reactant_heuristic.py   # Reagent role classification (FP + MCS)
│   ├── eln_csv_parser.py       # Findmolecule ELN CSV parser
│   ├── rdf_parser.py           # SciFinder .rdf reaction parser
│   ├── scheme_reader_audit.py  # Quality audit tool
│   ├── scheme_reader_verify.py # Visual HTML verification report
│   └── scheme_refine.py        # LLM-based refinement layer
│
├── resolve/                    # Chemical name/formula → SMILES
│   ├── reagent_db.py           # Tier 1: curated reagent DB (~186 entries)
│   ├── condensed_formula.py    # Tier 2: generative parser (PhB(OH)2 → SMILES)
│   ├── cas_resolver.py         # Tier 4: PubChem name/CAS → SMILES
│   ├── superatom_table.py      # Fragment vocabulary (2,854 entries)
│   ├── reagent_abbreviations.json
│   ├── chemscanner_abbreviations.json
│   └── superatom_data.json
│
├── render/                     # Generate CDXML from descriptions (was dsl/)
│   ├── schema.py               # Dataclass definitions
│   ├── parser.py               # YAML → SchemeDescriptor
│   ├── compact_parser.py       # Compact text → SchemeDescriptor
│   ├── renderer.py             # CDXML rendering engine (5 layouts)
│   ├── render_scheme.py        # CLI entry point (cdxml-render)
│   ├── auto_layout.py          # Zero-effort JSON → CDXML
│   ├── scheme_yaml_writer.py   # JSON → YAML layout decisions
│   └── scheme_maker.py         # JSON → CDXML scheme (experimental)
│
├── layout/                     # Polish, align, merge existing schemes
│   ├── scheme_polisher.py      # Surgical CDXML modification
│   ├── scheme_polisher_v2.py   # Full polishing pipeline
│   ├── reaction_cleanup.py     # Pure-Python layout (6 approaches)
│   ├── alignment.py            # Structure alignment (Kabsch, MCS, RXNMapper)
│   ├── scheme_aligner.py       # Product-relative MCS alignment
│   ├── scheme_merger.py        # Multi-scheme merging (parallel/seq/auto)
│   └── eln_enrichment.py       # ELN CSV → scheme annotation
│
├── naming/                     # IUPAC name analysis and alignment
│   ├── name_decomposer.py      # IUPAC name decomposition (ChemScript)
│   └── aligned_namer.py        # Multi-step aligned naming (Viterbi DP)
│
├── chemdraw/                   # ChemDraw-specific integrations (COM, ChemScript)
│   ├── cdx_converter.py        # CDX ↔ CDXML (COM / pycdxml / obabel)
│   ├── cdxml_to_image.py       # CDXML → PNG/SVG (ChemDraw COM)
│   ├── cdxml_to_image_rdkit.py # CDXML → image (RDKit-only backup)
│   ├── chemscript_bridge.py    # ChemScript .NET bridge
│   ├── _chemscript_server.py   # 32-bit subprocess server
│   └── eln_cdx_cleanup.py      # ELN CDX file cleanup
│
├── office/                     # Office document integration (OLE, PPTX, DOCX)
│   ├── ole_embedder.py         # CDXML → editable OLE in PPTX/DOCX
│   ├── ole_extractor.py        # Extract ChemDraw from Office files
│   └── doc_from_template.py    # Fill PPTX/DOCX templates
│
├── image/                      # Image-based structure extraction
│   ├── structure_from_image.py # Image → SMILES + 2D coords (DECIMER)
│   └── reaction_from_image.py  # Screenshot → reaction scheme CDXML
│
└── mcp_server/                 # Model Context Protocol server (unchanged)
    ├── server.py
    ├── __main__.py
    └── __init__.py
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

### resolve/reagent_db.py
Two-tier reagent database. Tier-1 (~186 curated entries with roles) always wins. Tier-2 (5,837 ChemScanner entries, no roles) is fallback. Key methods: `display_for_name()`, `role_for_name()`, `display_for_smiles()`, `smiles_role_display()`.

**Name normalization** (progressive, in `_lookup_name_entry()`): Unicode subscript digits → ASCII (Pd₂(dba)₃ → pd2(dba)3), solvate suffix stripping (·CHCl3, .DCM, ·HCl, etc.), rac-/(±)- prefix stripping with fallback to base ligand name. All three name-lookup methods use this cascade.

### resolve/condensed_formula.py
Generative parser for condensed structural formulae (PhB(OH)₂, Et₃N, Me₃SiCl, etc.). Uses the superatom table (2,854 entries) as fragment vocabulary. Two phases:

1. **Tokenizer** (`tokenize()`): greedy longest-match against superatom table + periodic table elements. Token types: `abbrev`, `element`, `count`, `paren_open`, `paren_close`. Two-letter elements (Na, Ag, Si) always match before abbreviations. Single-letter element collisions (n→N, s→S, etc. in superatom table) are excluded from abbreviation matching.

2. **Assembler** (`resolve_condensed_formula()`): RDKit-based molecule assembly from token stream. Handles five patterns: linear (MeI), multiplied groups (Et₃N), parenthesised branches (PhB(OH)₂), hydrogen subscripts (PhCH₂Br), and ionic/metal (Ag₂O). Validates with `Chem.SanitizeMol()`. Returns canonical SMILES or None.

**Resolution chain** in `perception/reaction_parser._resolve_text_label()`:
1. reagent_db (curated dictionary)
2. condensed formula parser (generative, offline)
3. OPSIN (offline, IUPAC/systematic names)
4. PubChem (online, if `use_network=True`)

### perception/scheme_reader.py
Reads CDXML reaction schemes into structured `SchemeDescription` JSON. Dual-strategy parsing: step-attribute (ChemDraw native `<scheme><step>`) or geometry-based (spatial arrow detection). Key concepts:
- **Text classification** (`_classify_text_species`): categorises `<t>` elements near arrows as `"chemical"`, `"condition_ref"`, `"footnote"`, `"yield"`, `"compound_label"`, `"citation"`, or `"bioactivity"`. Only `"chemical"` species undergo SMILES resolution; others are metadata.
- **Footnote linking**: condition_ref letters (a, b, c above arrows) are linked to `(a) conditions...` footnote text blocks elsewhere on the page.
- **Multi-line text splitting**: `"chemical"` text blocks (e.g. "Pd2(dba)3 (5 mol%)\nrac-BINAP (10 mol%)\nCs2CO3 (2.0 eq)\n1,4-dioxane\nreflux, 24 h") are split on newlines, then on comma+space or semicolons (protecting names like "1,4-dioxane") into individual `SpeciesRecord`s. Condition tokens (temperature, time, atmosphere) are filtered via `_is_condition_token()` and go to `step.conditions` only.
- **`SpeciesRecord` fields**: `is_solvent` (bool, set when `reagent_db.role_for_name()` returns `"solvent"`), `equiv_text` (e.g. "1.2 eq", "5 mol%", extracted from parenthetical annotations), `text_category`.
- **`elem_to_species`** maps CDXML element IDs to `List[str]` of species IDs (one-to-many for split text blocks).
- Returns `SchemeDescription` with `species` dict, `steps` list, `topology`, `content_type`, `narrative`, and optional `scope_entries`.

### perception/eln_csv_parser.py
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
| `rdkit` | render/, rdkit_utils, layout/alignment, perception/reaction_parser | SMILES, 2D coords, MW, MCS |
| `pywin32` | chemdraw/ (cdx_converter, cdxml_to_image), office/ole_embedder | ChemDraw COM automation (Windows) |
| `pyyaml` | render/parser, render/scheme_yaml_writer | YAML input parsing |
| `python-pptx` / `python-docx` / `olefile` | office/ (ole_embedder, ole_extractor, doc_from_template) | Office file manipulation |
| `opencv-python` / `Pillow` | image/, chemdraw/cdxml_to_image | Image processing |

ChemDraw COM tools require ChemDraw to be installed and **closed** before running.

## Running tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

573 tests. All pure Python — no ChemDraw COM required for the test suite.

## Bundled samples

`samples/` contains two real ELN exports with full pipeline output:
- `KL-CC-001/` — Buchwald coupling (source CDX/CSV/RXN + reaction.json + scheme.yaml + scheme.cdxml)
- `KL-CC-002/` — SNAr alkylation (same structure)
- `consolidated/` — Hand-authored two-step sequential YAML combining both reactions
