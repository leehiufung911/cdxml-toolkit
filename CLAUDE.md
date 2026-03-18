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
| `cdxml-polish` | `deterministic_pipeline.legacy.scheme_polisher_v2` | Clean up and enrich an ELN-exported CDXML scheme (legacy pipeline) |
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

## Agent playbook

This section maps common agent intents to the right tool chains. The toolkit is designed for LLM orchestration — the agent reasons about chemistry and calls tools; it never generates SMILES directly.

### Core principle: never generate SMILES

LLMs are unreliable at producing valid SMILES strings. Instead, use grounded tools:

| Need SMILES for... | Use this |
|---------------------|----------|
| A known compound name/abbreviation | `mol_builder.resolve_to_smiles("Cs2CO3")` |
| A structure in an image | `image/structure_from_image.py` (DECIMER extraction) |
| A reaction product | `mol_builder.apply_reaction()` or `mol_builder.deprotect()` |
| A modified molecule | `mol_builder.modify_name()` or `mol_builder.assemble_name()` |

### Key Python APIs for agents

**Name resolution** (public API — use this, not the private `_resolve_text_label`):
```python
from cdxml_toolkit.naming.mol_builder import resolve_to_smiles
result = resolve_to_smiles("2-fluoronicotinic acid")
# {'ok': True, 'smiles': '...', 'source': 'opsin'}
```

**Molecule construction** (reaction templates, deprotection, name surgery):
```python
from cdxml_toolkit.naming.mol_builder import (
    apply_reaction, deprotect, list_reactions,
    assemble_name, modify_name, get_prefix_form, validate_name,
    get_tool_definitions,   # export schemas for LLM function calling
)

# Deprotect a Boc-protected amine
result = deprotect(boc_amine_smiles)
# {'ok': True, 'product_smiles': '...', 'removed': ['Boc']}

# Amide coupling
result = apply_reaction("amide_coupling", amine_smiles, acid_smiles)
# {'ok': True, 'products': [{'smiles': '...', 'name': '...'}]}

# See available reaction templates
templates = list_reactions()
```

**Structure understanding** (decompose a molecule to reason about reactivity):
```python
from cdxml_toolkit.naming.name_decomposer import decompose_name
result = decompose_name(smiles)
# DecompositionResult with canonical_name, bracket_tree, alternatives
# Shows parent scaffold, substituents, functional groups — lets the agent
# "see" the molecule from multiple angles before choosing a reaction.
```

**Scheme rendering** (direct Python API — no YAML files needed):
```python
from cdxml_toolkit.render.schema import (
    SchemeDescriptor, StepDescriptor, StructureRef, ArrowContent,
)
from cdxml_toolkit.render.renderer import render, render_to_file

scheme = SchemeDescriptor(
    layout="sequential",
    structures={
        "SM": StructureRef(id="SM", smiles="..."),
        "Product": StructureRef(id="Product", smiles="..."),
    },
    steps=[
        StepDescriptor(
            substrates=["SM"], products=["Product"],
            below_arrow=ArrowContent(text=["TFA", "DCM"]),
        ),
    ],
)
render_to_file(scheme, "scheme.cdxml")
```

### Workflow recipes

**"Draw a scheme from an image"**
1. `image/reaction_from_image.py` — extract structures from screenshot via DECIMER
2. `perception/reaction_parser.py` — parse into semantic JSON (species, roles, SMILES)
3. `render/render_scheme.py --from-json` — render JSON into CDXML scheme

**"Draw a scheme I describe in words"**
1. `mol_builder.resolve_to_smiles()` — resolve each named compound to SMILES
2. `naming/name_decomposer.decompose_name()` — understand structure/reactivity if needed
3. `mol_builder.apply_reaction()` / `deprotect()` — compute products
4. Build `SchemeDescriptor` with steps → `render()` → CDXML

**"What would this reaction produce?"**
1. `mol_builder.resolve_to_smiles()` — get SMILES for all reactants
2. `naming/name_decomposer.decompose_name()` — analyze functional groups present
3. `mol_builder.list_reactions()` — find applicable reaction templates
4. `mol_builder.apply_reaction()` — get product SMILES

**"Resolve a compound name"**
- `mol_builder.resolve_to_smiles(name)` — 4-tier chain: reagent DB → condensed formula → OPSIN → PubChem

**"Modify an existing structure"**
- `mol_builder.modify_name(smiles, add=[...], remove=[...])` — add/swap/remove substituents via IUPAC name manipulation
- `mol_builder.assemble_name(parent, substituents)` — build from scratch

### LLM function-calling integration

`mol_builder.get_tool_definitions()` exports all tool schemas in Anthropic tool-use format. Register these as available tools in your LLM orchestrator:

```python
from cdxml_toolkit.naming.mol_builder import get_tool_definitions
tools = get_tool_definitions()
# [{"name": "resolve_to_smiles", "description": "...", "input_schema": {...}}, ...]
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
│  # ── Agent-friendly subpackages ─────────────────────────
├── perception/                 # Read and understand reaction schemes
│   ├── scheme_reader.py        # CDXML → SchemeDescription JSON
│   ├── scheme_segmenter.py     # Multi-panel CDXML auto-segmentation
│   ├── spatial_assignment.py   # Geometry-based element → arrow assignment
│   ├── reaction_parser.py      # ELN exports → JSON reaction descriptor
│   ├── reactant_heuristic.py   # Reagent role classification (FP + MCS)
│   ├── eln_csv_parser.py       # Findmolecule ELN CSV parser
│   ├── rdf_parser.py           # SciFinder .rdf reaction parser
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
├── render/                     # Generate CDXML from descriptions
│   ├── schema.py               # Dataclass definitions
│   ├── parser.py               # YAML → SchemeDescriptor
│   ├── compact_parser.py       # Compact text → SchemeDescriptor
│   ├── renderer.py             # CDXML rendering engine (5 layouts)
│   ├── render_scheme.py        # CLI entry point (cdxml-render)
│   ├── auto_layout.py          # Zero-effort JSON → CDXML
│   ├── scheme_yaml_writer.py   # JSON → YAML layout decisions
│   └── scheme_maker.py         # JSON → CDXML scheme (experimental)
│
├── layout/                     # Composable layout tools
│   ├── reaction_cleanup.py     # Pure-Python layout (6 approaches)
│   ├── alignment.py            # Structure alignment (Kabsch, MCS, RXNMapper)
│   └── scheme_merger.py        # Multi-scheme merging (parallel/seq/auto)
│
├── naming/                     # Agent chemistry reasoning tools
│   ├── mol_builder.py          # Agent API: resolve names, apply reactions, deprotect, name surgery
│   ├── name_decomposer.py      # Structural decomposition: see a molecule's functional groups
│   └── aligned_namer.py        # Multi-step aligned naming (Viterbi DP)
│
├── chemdraw/                   # ChemDraw-specific integrations (COM, ChemScript)
│   ├── cdx_converter.py        # CDX ↔ CDXML (COM / pycdxml / obabel)
│   ├── cdxml_to_image.py       # CDXML → PNG/SVG (ChemDraw COM)
│   ├── cdxml_to_image_rdkit.py # CDXML → image (RDKit-only backup)
│   ├── chemscript_bridge.py    # ChemScript .NET bridge
│   └── _chemscript_server.py   # 32-bit subprocess server
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
├── analysis/                   # LCMS parsing, species ID, lab book generation
│   ├── lcms_analyzer.py        # Single-file LCMS PDF parser (Waters + manual)
│   ├── format_procedure_entry.py # Agent-driven entries-based lab book formatter
│   ├── extract_nmr.py          # NMR data extraction from MestReNova PDFs
│   └── deterministic/          # Original deterministic batch pipeline
│       ├── multi_lcms_analyzer.py  # Cross-file LCMS collation
│       ├── procedure_writer.py     # Lab book entry assembler
│       ├── lcms_identifier.py      # LCMS species identification by mass
│       ├── lab_book_formatter.py   # Output section formatting
│       ├── mass_resolver.py        # Structure-based mass determination
│       ├── lcms_file_categorizer.py # LCMS filename classification
│       └── discover_experiment_files.py # Experiment file discovery
│
├── deterministic_pipeline/     # Rigid multi-step workflows (not for agent composition)
│   ├── scheme_reader_audit.py  # Quality audit for scheme_reader output
│   ├── scheme_reader_verify.py # Visual HTML verification report
│   └── legacy/                 # Old polisher pipeline (non-JSON-first)
│       ├── scheme_polisher.py      # Surgical CDXML modification (9-step chain)
│       ├── scheme_polisher_v2.py   # COM-free polishing pipeline (cdxml-polish CLI)
│       ├── eln_enrichment.py       # ELN CSV → scheme annotation (two-phase)
│       ├── scheme_aligner.py       # Product-relative MCS alignment
│       └── eln_cdx_cleanup.py      # ELN CDX file cleanup (ChemDraw COM)
│
└── mcp_server/                 # Model Context Protocol server
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

## Legacy polisher pipeline

The older non-JSON-first pipeline lives in `deterministic_pipeline/legacy/`. It surgically modifies ELN-exported CDXML rather than building from scratch. Preserved for backward compatibility but not recommended for new agent workflows.

```bash
cdxml-polish experiment.cdxml -o polished.cdxml \
    --approach chemdraw_mimic \
    --align-mode mcs \
    --eln-csv experiment.csv
```

Pipeline: `eln_cdx_cleanup` (CDX→CDXML + COM cleanup) → `scheme_polisher_v2` (ACS normalization + reagent classification + structure/text swaps + alignment + formatting) → `eln_enrichment` (equiv injection + run arrows).

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

646 tests (+ 1 pre-existing failure). All pure Python — no ChemDraw COM required for the test suite.

## Bundled samples

`samples/` contains two real ELN exports with full pipeline output:
- `KL-CC-001/` — Buchwald coupling (source CDX/CSV/RXN + reaction.json + scheme.yaml + scheme.cdxml)
- `KL-CC-002/` — SNAr alkylation (same structure)
- `consolidated/` — Hand-authored two-step sequential YAML combining both reactions
