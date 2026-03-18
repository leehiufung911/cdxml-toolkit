# cdxml-toolkit

> **WORK IN PROGRESS** — This project is under active development. Expect breaking changes, missing features, and rough edges. Do NOT be surprised if things don't work.

Python toolkit for ChemDraw CDXML reaction scheme processing, layout, and rendering.

Built for organic/medicinal chemists who work with reaction schemes in ChemDraw. Reads, writes, and manipulates CDXML files programmatically — build publication-ready reaction schemes from SMILES, clean up ELN exports, render schemes from declarative YAML, and more.

TLDR: Best feature currently is being able to specify a reaction scheme in a text format (YAML or mermaid like) and get a .cdxml chemdraw scheme. Lots of other little utilities also, like OLE embedding into docx/pptx

PROJECT INTENT:
Broadly/in the long run, to build out a toolkit to make organic chemistry office work more automatable (A LOT of which involves chemdraw)
Such tools may be run on their own, or as I envision, used/called by an LLM agent

(Prime example: DSL which allows LLM to write a YAML or mermaid description of a scheme, which is then rendered into a .cdxml using renderer.py)

There are a bunch of dependencies that are easy to install (RDKit), some that may be more finnicky (Chemscript, Chemdraw if it doesn't autodetect. **Chemdraw (COM) is an IMPORTANT DEPENDENCY. There are fallbacks, but I haven't validated how well they work**)
This was made on a machine with Chemdraw 16.
I will validate and fix things gradually.

Professionally, I do not have much programming background. I am a just a PhD organic chemist trying to make life easier hopefully for myself and my fellow chemists.

**I leaned heavily on Claude Code in the making of this project.** I directed, refined, debugged, and figured out the broad design of things, but Claude Code (Opus 4.6) did basically all of the actual coding.

Note that LLMs (even frontier models) are terrible at judging whether a chemical structure is correct, as well as whether a scheme is correctly spaced... they're even bad at writing SMILES. This is why the design decision was made to have a source of chemical truth (JSON, parsed from CDX files/other ELN files), have a script or an LLM declare the layout of a scheme with a YAML file (but not actually place objects), and have a layout engine (renderer.py) deterministically place objects at the right places. This affords flexibility for the LLM to choose the scheme layout it wants, but avoids its weaknesses of being bad at chemistry and layout/design.

I would like to express immense gratitude to Anthropic for making such a revolutionary product. Without it, these ideas would only be random musings in my brain, since I do not have the programming skill to implement any of this in any reasonable amount of time otherwise.

## Installation

```bash
# Core package (CDXML manipulation, no external dependencies beyond lxml)
pip install cdxml-toolkit

# With RDKit support (structure alignment, MW calculation, 2D coord generation)
pip install cdxml-toolkit[rdkit]

# Everything (RDKit + ChemDraw COM + PDF parsing + Office document support)
pip install cdxml-toolkit[all]

# Development install from source
git clone https://github.com/leehiufung911/cdxml-toolkit.git
cd cdxml-toolkit
pip install -e ".[dev]"
```

## Quick(?) Start

Honestly? Maybe try using this repository with Claude Code/Opencode/Windsurf/another agent. The endgame is meant to be like "Hey, help me make the scheme for these reactions", and the agent just goes and does the things you ask.
I uploaded a CLAUDE.MD. This is meant to be read as reference by both humans and LLMs.

## Quick Start — JSON-First Pipeline

The recommended workflow builds publication-ready schemes from ELN export files in two steps:

```bash
# Step 1: Parse reaction files -> semantic JSON
cdxml-parse experiment.cdxml --csv experiment.csv -o reaction.json

# Step 2: Render JSON -> publication-ready CDXML
cdxml-render --from-json reaction.json -o scheme.cdxml
```

The reaction parser extracts every species (structures, reagents, solvents, products) with canonical SMILES, roles, equivalents, and mass data. The renderer turns that into a properly laid-out scheme — no manual positioning needed.

If your input is CDX (binary ChemDraw), convert first:

```bash
cdxml-convert experiment.cdx -o experiment.cdxml
```

### Python API

```python
from cdxml_toolkit.perception import parse_reaction
from cdxml_toolkit.render import render_to_file

# Parse
desc = parse_reaction(cdxml="experiment.cdxml", csv="experiment.csv")
desc.to_json("reaction.json")

# Render
from cdxml_toolkit.render.auto_layout import auto_layout_to_cdxml
auto_layout_to_cdxml("reaction.json", "scheme.cdxml")
```

### Other input formats

The renderer also accepts hand-authored YAML or compact text syntax:

```bash
# YAML input
cdxml-render scheme.yaml -o scheme.cdxml

# Compact text ("Mermaid for reactions")
cdxml-render scheme.txt -o scheme.cdxml
```

### Reagent database

```python
from cdxml_toolkit import get_reagent_db

db = get_reagent_db()
db.display_for_name("cs2co3")    # "Cs2CO3"
db.role_for_name("cs2co3")       # "base"
db.display_for_name("hatu")      # "HATU" (from ChemScanner tier-2)
```

## Features

- **Reaction parser** — Unified semantic layer: CDX/CDXML/RXN/CSV -> JSON with canonical SMILES, roles, equivalents, masses, and adducts. The single source of truth for all downstream tools.
- **Scheme DSL renderer** — JSON/YAML/compact text -> publication-ready CDXML. 6 layout engines (linear, sequential, serpentine, divergent, stacked-rows, numbered-parallel). Run arrows, dashed/failed arrows, compound labels. No ChemDraw COM needed — uses RDKit for 2D coordinate generation.
- **Mol builder** — Agent API for compound name resolution (reagent DB -> condensed formula -> OPSIN -> PubChem), reaction template application, deprotection, and IUPAC name surgery. Exports LLM function-calling schemas.
- **LCMS analysis** — Parse Waters and manual LCMS PDFs, identify species by mass, cross-file collation.
- **Lab book formatting** — Agent-driven procedure entry formatting from parsed reaction data.
- **Scheme merger** — Combine multiple reaction schemes with auto-detected relationships (parallel, sequential, unrelated).
- **CDXML reading/writing** — Parse and write CDXML with DOCTYPE preservation. Atom-only bounding boxes (reliable for abbreviation groups). Hanging label detection.
- **Chemical text formatting** — Subscript digits in formulas (CH3OH), italic IUPAC prefixes (*n*-BuLi, *tert*-BuOH).
- **Reagent database** — Two-tier lookup: ~186 curated entries with roles + ~5,837 ChemScanner entries. SMILES-based and name-based matching.
- **Superatom table** — ~2,850 abbreviation->SMILES mappings for MW calculation of ChemDraw abbreviation groups (OTs, Boc, Me, etc.).
- **RDKit bridge** — CDXML fragment -> RDKit Mol, SMILES, MW. 2D coordinate cleanup with Kabsch orientation preservation.
- **Reaction layout** — Pure Python replacement for ChemDraw COM "Clean Up Reaction". 6 layout approaches including `chemdraw_mimic` (default).
- **Structure alignment** — Kabsch, RDKit MCS, and RXNMapper-based orientation alignment.
- **OLE embedding** — Embed CDXML as editable ChemDraw objects in PPTX/DOCX.
- **MCP server** — Model Context Protocol server for tool integration (skeleton).

## Requirements

- **Python >= 3.9**
- **lxml** — core dependency (always required)

| Optional dependency | What it enables | Install extra |
|---|---|---|
| RDKit | Structure alignment, MW calculation, SMILES, 2D coords | `[rdkit]` |
| PyYAML | YAML scheme parsing | `[yaml]` |
| pywin32 | ChemDraw COM automation (Windows only) | `[chemdraw]` |
| python-pptx, python-docx, olefile | Office document support | `[office]` |
| OpenCV, Pillow | Image-based structure extraction | `[image]` |
| pdfplumber | LCMS PDF parsing | `[analysis]` |
| mcp | MCP server | `[mcp]` |

## CLI Tools

Installed as console scripts:

| Command | Description |
|---|---|
| `cdxml-parse` | Parse reaction files -> JSON descriptor (**start here**) |
| `cdxml-render` | Render JSON/YAML/compact text -> CDXML scheme |
| `cdxml-convert` | CDX <-> CDXML conversion |
| `cdxml-image` | Render CDXML to PNG/SVG (ChemDraw COM) |
| `cdxml-merge` | Merge multiple reaction schemes |
| `cdxml-layout` | Clean up reaction layout (pure Python, no COM) |
| `cdxml-build` | Build CDXML from atom/bond data |
| `cdxml-polish` | Full scheme polishing pipeline (legacy) |
| `cdxml-ole` | Embed CDXML as editable OLE in PPTX/DOCX |
| `cdxml-mcp` | MCP server |
| `cdxml-lcms` | Parse single LCMS PDF report |
| `cdxml-multi-lcms` | Cross-file LCMS collation |
| `cdxml-procedure` | Assemble lab book entries |
| `cdxml-discover` | Discover experiment files in a directory |
| `cdxml-format-entry` | Agent-driven lab book entry formatting |
| `cdxml-nmr` | Extract NMR data from MestReNova PDFs |

All tools also work as Python modules:
```bash
python -m cdxml_toolkit.perception.reaction_parser experiment.cdxml --csv exp.csv -o reaction.json
python -m cdxml_toolkit.render.render_scheme --from-json reaction.json -o scheme.cdxml
```

## Layout Patterns

The DSL renderer supports 5 layout patterns. The layout is auto-selected based on the reaction JSON, or you can override it:

| Layout | Description |
|---|---|
| `linear` | Default. Single-step reaction |
| `sequential` | Multi-step with shared intermediates drawn once |
| `divergent` | One starting material -> multiple products |
| `stacked-rows` | Independent rows with optional section labels |
| `serpentine` | Alternating direction rows (left->right, right->left) |

No ChemDraw COM needed — structures are generated from SMILES via RDKit.

## Bundled Samples

The `samples/` directory contains two real ELN exports (Buchwald coupling, SNAr alkylation) with full pipeline output:

```
samples/
├── KL-CC-001/         # Buchwald coupling (single-step)
│   ├── KL-CC-001.cdx  # ELN export (binary ChemDraw)
│   ├── KL-CC-001.csv  # ELN export (reagent table)
│   ├── KL-CC-001.rxn  # ELN export (RXN file)
│   ├── reaction.json  # Parsed reaction (cdxml-parse output)
│   ├── scheme.yaml    # Auto-generated YAML
│   └── scheme.cdxml   # Rendered scheme (cdxml-render output)
│
├── KL-CC-002/         # SNAr alkylation (single-step)
│   └── ...            # Same structure as KL-CC-001
│
└── consolidated/      # Two-step sequential scheme
    ├── two-step-scheme.yaml   # Hand-authored sequential YAML
    └── two-step-scheme.cdxml  # Rendered two-step scheme
```

## Package Structure

```
cdxml_toolkit/
│
│  # ── Foundation ─────────────────────────────────────────
├── __init__.py                 # Version, core exports
├── constants.py                # ACS Document 1996 style, layout gaps
├── cdxml_utils.py              # CDXML geometry (bbox, IO, id map)
├── rdkit_utils.py              # CDXML fragment <-> RDKit Mol/SMILES/MW
├── text_formatting.py          # Chemical subscripts + italic prefixes
├── cdxml_builder.py            # Build CDXML from atom/bond data
├── coord_normalizer.py         # Coordinate normalization to ACS 1996
│
│  # ── Subpackages ────────────────────────────────────────
├── perception/                 # Read and understand reaction schemes
│   ├── reaction_parser.py      # ELN exports -> JSON reaction descriptor
│   ├── scheme_reader.py        # CDXML -> SchemeDescription JSON
│   ├── scheme_segmenter.py     # Multi-panel CDXML auto-segmentation
│   ├── reactant_heuristic.py   # Reagent role classification (FP + MCS)
│   └── ...
│
├── resolve/                    # Chemical name/formula -> SMILES
│   ├── reagent_db.py           # Curated reagent DB (~186 entries)
│   ├── condensed_formula.py    # Generative parser (PhB(OH)2 -> SMILES)
│   ├── cas_resolver.py         # PubChem name/CAS -> SMILES
│   └── superatom_table.py      # Fragment vocabulary (2,854 entries)
│
├── render/                     # Generate CDXML from descriptions
│   ├── renderer.py             # CDXML rendering engine (5 layouts)
│   ├── schema.py               # Dataclass definitions
│   ├── parser.py               # YAML parser
│   ├── compact_parser.py       # Compact text parser
│   ├── auto_layout.py          # Zero-effort JSON -> CDXML
│   └── scheme_yaml_writer.py   # JSON -> YAML layout decisions
│
├── layout/                     # Composable layout tools
│   ├── reaction_cleanup.py     # Pure-Python layout (6 approaches)
│   ├── alignment.py            # Structure alignment (Kabsch, MCS)
│   └── scheme_merger.py        # Multi-scheme merging
│
├── naming/                     # Agent chemistry reasoning tools
│   ├── mol_builder.py          # Resolve names, apply reactions, deprotect
│   ├── name_decomposer.py      # Structural decomposition
│   └── aligned_namer.py        # Multi-step aligned naming (Viterbi DP)
│
├── analysis/                   # LCMS parsing, lab book generation
│   ├── lcms_analyzer.py        # Single-file LCMS PDF parser
│   ├── format_procedure_entry.py # Agent-driven lab book formatting
│   ├── extract_nmr.py          # NMR data extraction
│   └── deterministic/          # Batch pipeline (cross-file LCMS, etc.)
│
├── chemdraw/                   # ChemDraw COM automation (Windows)
├── office/                     # OLE embedding in PPTX/DOCX
├── image/                      # Image-based structure extraction (DECIMER)
├── deterministic_pipeline/     # Legacy polishing workflows
└── mcp_server/                 # Model Context Protocol server
```

## ChemDraw COM Warning

Tools using ChemDraw COM automation (`cdx_converter`, `cdxml_to_image`, `scheme_polisher`, `ole_embedder`) launch and control ChemDraw programmatically. **Close ChemDraw before running these tools.**

## ACS Document 1996 Style

All CDXML output conforms to ACS Document 1996 style:
- Bond length: 14.40 pt
- Chain angle: 120 degrees
- Font: Arial 10pt Bold
- Line width: 0.6 pt

## License

[MIT](LICENSE)

## Attribution

See [NOTICE.md](NOTICE.md) for third-party data attribution (ChemScanner, RDKit).

## Author

Hiu Fung Kevin Lee ([@leehiufung911](https://github.com/leehiufung911))
