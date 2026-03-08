# cdxml-toolkit

Python toolkit for ChemDraw CDXML reaction scheme processing, layout, and rendering.

Built for organic/medicinal chemists who work with reaction schemes in ChemDraw. Reads, writes, and manipulates CDXML files programmatically — build publication-ready reaction schemes from SMILES, clean up ELN exports, render schemes from declarative YAML, and more.

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

## Quick Start — JSON-First Pipeline

The recommended workflow builds publication-ready schemes from ELN export files in two steps:

```bash
# Step 1: Parse reaction files → semantic JSON
cdxml-parse experiment.cdxml --csv experiment.csv -o reaction.json

# Step 2: Render JSON → publication-ready CDXML
cdxml-render --from-json reaction.json -o scheme.cdxml
```

The reaction parser extracts every species (structures, reagents, solvents, products) with canonical SMILES, roles, equivalents, and mass data. The renderer turns that into a properly laid-out scheme — no manual positioning needed.

If your input is CDX (binary ChemDraw), convert first:

```bash
cdxml-convert experiment.cdx -o experiment.cdxml
```

### Python API

```python
from cdxml_toolkit.reaction_parser import parse_reaction
from cdxml_toolkit.dsl import render_to_file

# Parse
desc = parse_reaction(cdxml="experiment.cdxml", csv="experiment.csv")
desc.to_json("reaction.json")

# Render
from cdxml_toolkit.dsl.auto_layout import auto_layout_to_cdxml
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

- **Reaction parser** — Unified semantic layer: CDX/CDXML/RXN/CSV → JSON with canonical SMILES, roles, equivalents, masses, and adducts. The single source of truth for all downstream tools.
- **Scheme DSL renderer** — JSON/YAML/compact text → publication-ready CDXML. 6 layout engines (linear, sequential, serpentine, divergent, stacked-rows, numbered-parallel). Run arrows, dashed/failed arrows, compound labels. No ChemDraw COM needed — uses RDKit for 2D coordinate generation.
- **Scheme merger** — Combine multiple reaction schemes with auto-detected relationships (parallel, sequential, unrelated).
- **CDXML reading/writing** — Parse and write CDXML with DOCTYPE preservation. Atom-only bounding boxes (reliable for abbreviation groups). Hanging label detection.
- **Chemical text formatting** — Subscript digits in formulas (CH₃OH), italic IUPAC prefixes (*n*-BuLi, *tert*-BuOH).
- **Reagent database** — Two-tier lookup: ~172 curated entries with roles + ~5,837 ChemScanner entries. SMILES-based and name-based matching.
- **Superatom table** — ~2,850 abbreviation→SMILES mappings for MW calculation of ChemDraw abbreviation groups (OTs, Boc, Me, etc.).
- **RDKit bridge** — CDXML fragment → RDKit Mol, SMILES, MW. 2D coordinate cleanup with Kabsch orientation preservation.
- **Reaction layout** — Pure Python replacement for ChemDraw COM "Clean Up Reaction". 6 layout approaches including `chemdraw_mimic` (default).
- **Structure alignment** — Kabsch, RDKit MCS, and RXNMapper-based orientation alignment.
- **ELN enrichment** — Inject equivalents, run arrows, and mass/yield annotations from Findmolecule CSV exports.
- **Scheme polishing (legacy)** — Surgically modifies ELN-exported CDXML. Fragment cleanup → normalization → reagent classification → text↔structure swap → alignment → layout. The JSON-first pipeline above is recommended for new work.

## Requirements

- **Python ≥ 3.9**
- **lxml** — core dependency (always required)

| Optional dependency | What it enables | Install extra |
|---|---|---|
| RDKit | Structure alignment, MW calculation, SMILES, 2D coords | `[rdkit]` |
| PyYAML | YAML scheme parsing | `[yaml]` |
| pywin32 | ChemDraw COM automation (Windows only) | `[chemdraw]` |
| pdfplumber | LCMS PDF parsing | `[pdf]` |
| python-pptx, python-docx, olefile | Office document support | `[office]` |
| OpenCV, Pillow | Image-based structure extraction | `[image]` |

## CLI Tools

Installed as console scripts:

| Command | Description |
|---|---|
| `cdxml-parse` | Parse reaction files → JSON descriptor (**start here**) |
| `cdxml-render` | Render JSON/YAML/compact text → CDXML scheme |
| `cdxml-convert` | CDX ↔ CDXML conversion |
| `cdxml-image` | Render CDXML to PNG/SVG (ChemDraw COM) |
| `cdxml-merge` | Merge multiple reaction schemes |
| `cdxml-layout` | Clean up reaction layout (pure Python, no COM) |
| `cdxml-build` | Build CDXML from atom/bond data |
| `cdxml-polish` | Full scheme polishing pipeline (legacy) |
| `cdxml-ole` | Embed CDXML as editable OLE in PPTX/DOCX |

All tools also work as Python modules:
```bash
python -m cdxml_toolkit.reaction_parser experiment.cdxml --csv exp.csv -o reaction.json
python -m cdxml_toolkit.dsl.render_scheme --from-json reaction.json -o scheme.cdxml
```

## Layout Patterns

The DSL renderer supports 5 layout patterns. The layout is auto-selected based on the reaction JSON, or you can override it:

| Layout | Description |
|---|---|
| `linear` | Default. Single-step reaction |
| `sequential` | Multi-step with shared intermediates drawn once |
| `divergent` | One starting material → multiple products |
| `stacked-rows` | Independent rows with optional section labels |
| `serpentine` | Alternating direction rows (left→right, right→left) |

No ChemDraw COM needed — structures are generated from SMILES via RDKit.

## Bundled Samples

The `samples/` directory contains two real ELN exports (Buchwald coupling and SNAr alkylation) with full pipeline output:

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

## Project Structure

```
cdxml-toolkit/
├── cdxml_toolkit/              # The pip-installable package
│   ├── __init__.py
│   ├── constants.py            # ACS Document 1996 style constants
│   ├── cdxml_utils.py          # CDXML geometry and IO
│   ├── rdkit_utils.py          # RDKit-CDXML bridge
│   ├── text_formatting.py      # Chemical text formatting
│   ├── reagent_db.py           # Two-tier reagent database
│   ├── superatom_table.py      # Abbreviation→SMILES lookup
│   ├── reaction_cleanup.py     # Reaction layout engine
│   ├── reaction_parser.py      # Reaction → JSON semantic layer
│   ├── scheme_polisher_v2.py   # Full polishing pipeline
│   ├── scheme_merger.py        # Scheme merging
│   ├── alignment.py            # Structure alignment strategies
│   ├── ...                     # 30+ modules total
│   └── dsl/                    # Scheme DSL subpackage
│       ├── renderer.py         # CDXML rendering engine
│       ├── parser.py           # YAML parser
│       ├── compact_parser.py   # Compact text parser
│       ├── schema.py           # Dataclass definitions
│       ├── render_scheme.py    # CLI entry point
│       ├── auto_layout.py      # Zero-effort JSON→CDXML
│       └── scheme_yaml_writer.py  # JSON→YAML layout decisions
│
├── tests/                      # pytest test suite
├── experiments/                # ML tools (RXNMapper, FlowER, RXN Insight)
├── pyproject.toml
├── LICENSE                     # MIT
└── NOTICE.md                   # Third-party attribution
```

## ChemDraw COM Warning

Tools using ChemDraw COM automation (`cdx_converter`, `cdxml_to_image`, `scheme_polisher`, `ole_embedder`) launch and control ChemDraw programmatically. **Close ChemDraw before running these tools.**

## ACS Document 1996 Style

All CDXML output conforms to ACS Document 1996 style:
- Bond length: 14.40 pt
- Chain angle: 120°
- Font: Arial 10pt Bold
- Line width: 0.6 pt

## License

[MIT](LICENSE)

## Attribution

See [NOTICE.md](NOTICE.md) for third-party data attribution (ChemScanner, RDKit).

## Author

Hiu Fung Kevin Lee ([@leehiufung911](https://github.com/leehiufung911))
