# Chemistry Toolkit — Project Instructions

## Who I am
I'm an organic medicinal chemist in drug discovery. I work with reactions like Buchwald, SNAr, Suzuki, amide couplings, Boc deprotections, reductive aminations, Grignard additions, etc.

## What this project is
A modular toolkit of Python scripts that let an LLM assist with chemistry office work. Each tool does one thing well. Claude orchestrates them as needed. Version 0.4 — unified codebase (portable/ merged into root). 29 working tools + 1 experimental + 3 ML experimental (atom mapping, byproduct prediction, role classification) + 1 backup renderer + 7 shared modules + 3 support modules + Scheme DSL (7 modules: schema, YAML parser, compact parser, renderer, render_scheme CLI, auto_layout, scheme_yaml_writer).

## Key rules
- **ChemDraw output (CDXML) is non-negotiable.** All structure output must open in ChemDraw 16. Use ACS Document 1996 style (BondLength=14.40, ChainAngle=120, Arial 10pt).
- **Never trust LLM-generated SMILES.** Use compound name resolution (PubChem API, CAS lookup, ChemScript) or existing structure files.
- **Don't generate stoichiometry.** My ELN (Findmolecule) handles mass/mmol/equiv. Procedure text uses reagent names as placeholders.
- **Each tool = standalone Python script with CLI interface.** Use argparse. CDXML is the interchange format. Library modules live in `cdxml_toolkit/` (the pip-installable package); private workflow scripts stay at the project root.
- **Python environment:** Use `/c/Users/mic23/miniconda3/envs/LLMChem` for Python. RDKit is available in this environment.
- **Test with real files.** Test data lives outside the project tree in `C:\Users\mic23\chem-test-data\` (configurable via `CHEM_TEST_DATA` env var). Always validate output by opening in ChemDraw when possible.
- **`chem-test-data/additional_test/`** contains files for final user validation only. Do not use those files while building or testing tools.
- **Always render CDXML to image for visual inspection.** Whenever a CDXML file is produced or modified, run `cdxml_to_image.py` on it (ChemDraw backend — this is the primary, faithful backend) and visually inspect the result before reporting back. Do not rely solely on XML structure checks. The rendered image is the ground truth for correctness.

## Canonical workflow: JSON-first pipeline (CDX → JSON → CDXML)

The primary workflow builds publication-ready reaction schemes from ELN export files. The reaction parser extracts all species into a canonical JSON descriptor; the DSL renderer turns that into a properly laid-out CDXML scheme. No ChemDraw COM needed for rendering — only for CDX→CDXML conversion and image export.

**Standard invocation:**
```bash
# Convert CDX to CDXML (ChemDraw COM — close ChemDraw first)
python -m cdxml_toolkit.cdx_converter experiment.cdx -o experiment.cdxml

# Parse reaction into semantic JSON
python -m cdxml_toolkit.reaction_parser experiment.cdxml --csv experiment.csv -o reaction.json

# Render JSON → publication-ready CDXML
python -m cdxml_toolkit.dsl.render_scheme --from-json reaction.json -o scheme.cdxml

# Render to image for visual inspection
python -m cdxml_toolkit.cdxml_to_image scheme.cdxml
```

**What happens under the hood:**
1. **Reaction parsing** (`reaction_parser.py`) — extracts every species with canonical SMILES, role classification, display names, equivalents, mass data, and adducts. Produces a single JSON source of truth.
2. **Layout decisions** (`scheme_yaml_writer.py`) — determines which species are drawn structures vs text labels, positions them above/below the arrow, sorts by role priority.
3. **CDXML rendering** (`renderer.py`) — generates 2D coordinates from SMILES via RDKit, computes bounding boxes, sizes the arrow to fit content, formats chemical text (subscripts, italics), and outputs ACS Document 1996 styled CDXML.

**Graceful degradation:** The renderer requires only RDKit (for SMILES → 2D coordinates). No ChemDraw COM, no ChemScript. Falls back to text labels when SMILES resolution fails.

**Batch orchestrator:** `run_pipeline.py` wraps this pipeline with batch CDX→CDXML conversion (single COM session), YAML config, experiment discovery, LCMS analysis, and procedure writing. Set `scheme.renderer: "dsl"` in config.yaml to use the JSON-first pipeline. Setup infrastructure lives in `deploy/`.

**Architecture — three layers:**
```
reaction_parser JSON  ──→  scheme_yaml_writer  ──→  YAML  ──→  renderer  ──→  CDXML
                           (layout decisions)       (or compact text)      (spatial engine)
```

**Other input modes (for hand-authored schemes):**
```bash
# From YAML
python -m cdxml_toolkit.dsl.render_scheme scheme.yaml -o scheme.cdxml

# From compact text ("Mermaid for reactions")
python -m cdxml_toolkit.dsl.render_scheme scheme.txt -o scheme.cdxml

# Zero-effort (bypasses YAML, builds SchemeDescriptor in memory)
python -m cdxml_toolkit.dsl.auto_layout reaction.json -o scheme.cdxml
```

## Legacy workflow: CDX → polished enriched CDXML (scheme_polisher_v2)

The polisher pipeline surgically modifies an ELN-exported CDXML rather than building from scratch. It preserves the original ChemDraw drawing while cleaning up layout, classifying reagents, and injecting ELN data. Use this when you need to keep the exact original structures from the ELN export.

**Standard invocation:**
```bash
python -m cdxml_toolkit.scheme_polisher_v2 experiment.cdxml -o polished.cdxml \
    --approach chemdraw_mimic \
    --align-mode rxnmapper \
    --classify-method rxnmapper \
    --eln-csv experiment.csv
```

**What `scheme_polisher_v2.py` does, step by step:**

1. **Fragment geometry cleanup** — RDKit `cleanup_fragment_rdkit()` fixes bond angles while preserving orientation via Kabsch rotation. Falls back to ChemScript.
2. **Bond length normalization** — scales all fragments to ACS bond length (14.40 pt).
3. **ACS Document 1996 settings** — sets document-level style properties (ChainAngle=120, fonts, etc.).
4. **Font normalization** — standardizes all text to Arial 10pt Bold, fixes narrow ELN text, anchors orphan labels.
5. **Scheme polishing** (`scheme_polisher.py`) — reagent classification (RXNMapper or heuristic), structure-text swaps, deduplication, subscript formatting, condition text merging.
6. **Reactant/reagent alignment** (`alignment.py`) — orients all structures to match the product using RXNMapper atom maps or RDKit MCS + `GenerateDepictionMatching2DStructure`.
7. **ELN enrichment Phase A** (`eln_enrichment.py`, before layout) — injects equivalents into text labels (e.g. "Cs2CO3" → "Cs2CO3 (2 eq.)") so arrow length accounts for wider text. Optionally repositions non-substrate reactant above arrow.
8. **Compact toward arrow** — moves above/below-arrow objects closer to the arrow.
9. **Reaction layout** (`reaction_cleanup.py`) — final spatial positioning using `chemdraw_mimic` approach with atom-only bounding boxes.
10. **ELN enrichment Phase B** (`eln_enrichment.py`, after layout) — adds run arrow with substrate mass → product yield, fragment equivalents labels.

**Graceful degradation:** If `rxn-experiments` conda env is unavailable, classification and alignment fall back to heuristic role lookup + RDKit MCS. If RDKit is unavailable, alignment falls back to Kabsch. The pipeline always produces output.

**What the LLM specifies (semantic layer):**
- Structure references (SMILES, compound names, CDXML file paths — never coordinates)
- Roles (substrate, reagent, catalyst, product)
- Conditions text (reagents, solvents, temperature, time)
- Layout pattern hint (single keyword: `linear`, `sequential`, `divergent`, `stacked-rows`)
- Wrapping mode (`repeat`, `serpentine`, `none`) and `steps_per_row`
- Annotations (yields, mass amounts, run arrows, compound labels)

**What the renderer handles (spatial layer):**
- 2D coordinate generation (RDKit, no ChemDraw COM needed)
- Bounding box measurement via `cdxml_utils.fragment_bbox()`
- Arrow sizing from content width (min 5× ACS_BOND_LENGTH)
- Text formatting (subscripts, italics via `text_formatting.py`)
- Gap computation (ACS style from `constants.py`)
- MCS-based structure alignment to product orientation

**Layout patterns (fixed vocabulary):**
| Pattern | Description | Status |
|---------|-------------|--------|
| `linear` | Single-step reaction (default) | ✅ Implemented |
| `sequential` | Multi-step linear synthesis (2-10+ steps) | ✅ Implemented |
| `wrap: repeat` | Multi-row L→R with repeated structures (default wrapping) | ✅ Implemented |
| `wrap: serpentine` | Zigzag L→R, R→L with vertical arrows | ✅ Implemented |
| `divergent` | One SM giving multiple products (Y-shape fan-out) | ✅ Implemented |
| `stacked-rows` | Multiple independent sequences stacked vertically | ✅ Implemented |
| `numbered-parallel` | One-pot sequential additions | 🔜 Planned |
| `convergent` | Two streams merging | 🔜 Planned |

**Annotations (orthogonal to layout):**
- Run arrows (SM mass → product yield) — ✅ implemented
- Dashed arrows (`arrow_style: dashed`) — ✅ implemented
- Failed arrows with X overlay (`arrow_style: failed`) — ✅ implemented
- Compound number labels (`label:` on structures) — ✅ implemented
- Letter conditions (`condition_key:` block) — ✅ implemented

**Compact text syntax ("Mermaid for reaction schemes"):**
```
# Simple Buchwald coupling — 3 lines
ArBr: {BrC1=CC=CC=C1}
Morph: {C1COCCN1}

ArBr + Morph --> ArMorph{C1=CC=C(N2CCCCC2)C=C1} (72%)
  above: Morph
  below: "Pd2(dba)3", "BINAP", "Cs2CO3", "toluene, 110 °C, 16 h"
```

**Pipeline integration:** Activated in `run_pipeline.py` by setting `scheme.renderer: "dsl"` in config.yaml. When enabled, reaction_parser runs first (Phase 3.05), then `render_scheme.py --from-json` consumes the JSON (Phase 3a-dsl). Per-experiment override: `experiment_overrides: { KL-7001-008: { renderer: dsl } }`.

**Graceful degradation:** The renderer requires only RDKit (for SMILES → 2D coordinates). No ChemDraw COM, no ChemScript. Falls back to text labels when SMILES resolution fails.

**Internal documentation:** Detailed design docs live in `experiments/scheme_dsl/`: `DESIGN.md` (architecture, YAML format v0.2, real-world examples), `LAYOUT_PATTERNS.md` (formal rendering specification for all patterns), `COMPACT_SYNTAX.md` (compact text syntax with EBNF grammar), plus 30 showcase YAML examples in `experiments/scheme_dsl/showcase/`.

## Project layout
```
chem-tools-0.4/
├── cdxml_toolkit/                  # THE PACKAGE (pip-installable, import name: cdxml_toolkit)
│   ├── __init__.py                 # Version, top-level exports
│   ├── constants.py                # ACS Document 1996 style, LCMS defaults, layout gaps
│   ├── cdxml_utils.py              # CDXML geometry utilities (bbox, IO, id map)
│   ├── rdkit_utils.py              # RDKit-based CDXML fragment utilities (mol, SMILES, MW, cleanup)
│   ├── text_formatting.py          # Chemical text formatting (subscripts, italics)
│   ├── reagent_db.py               # Two-tier reagent database (singleton)
│   ├── superatom_table.py          # Superatom label→SMILES lookup (~2,850 entries)
│   ├── alignment.py                # Alignment strategies (Kabsch, RDKit MCS, RXNMapper)
│   ├── reaction_cleanup.py         # Pure Python reaction layout (6 approaches)
│   ├── reaction_parser.py          # Reaction → JSON semantic layer
│   ├── scheme_polisher.py          # Scheme polishing (classification, swaps, alignment)
│   ├── scheme_polisher_v2.py       # Full polishing pipeline
│   ├── scheme_merger.py            # Multi-scheme merging
│   ├── cdxml_builder.py            # CDXML construction from atom/bond data
│   ├── coord_normalizer.py         # Coordinate normalization to ACS 1996
│   ├── reactant_heuristic.py       # Reagent role classification
│   ├── chemscript_bridge.py        # ChemScript .NET bridge
│   ├── _chemscript_server.py       # Internal subprocess helper
│   ├── cas_resolver.py             # CAS/name → SMILES via PubChem
│   ├── cdx_converter.py            # CDX ↔ CDXML conversion
│   ├── cdxml_to_image.py           # CDXML → PNG/SVG (ChemDraw COM)
│   ├── cdxml_to_image_rdkit.py     # CDXML → image (backup, RDKit only)
│   ├── eln_enrichment.py           # ELN CSV → scheme annotation
│   ├── eln_cdx_cleanup.py          # ELN CDX cleanup
│   ├── ole_embedder.py             # CDXML → editable OLE in PPTX/DOCX
│   ├── ole_extractor.py            # Extract embedded ChemDraw from Office
│   ├── doc_from_template.py        # Fill PPTX/DOCX templates
│   ├── rdf_parser.py               # SciFinder RDF parser
│   ├── reaction_from_image.py      # Image → reaction scheme CDXML
│   ├── structure_from_image.py     # Image → structure CDXML (DECIMER)
│   ├── scheme_maker.py             # Experimental: JSON → complete scheme
│   ├── scheme_aligner.py           # Experimental: MCS-based alignment
│   ├── reagent_abbreviations.json  # Tier-1 reagent DB (172 entries)
│   ├── chemscanner_abbreviations.json  # Tier-2 ChemScanner DB (5,837 entries)
│   ├── superatom_data.json         # Superatom abbreviation data
│   └── dsl/                        # Scheme DSL subpackage
│       ├── __init__.py
│       ├── schema.py               # Dataclass definitions
│       ├── parser.py               # YAML → SchemeDescriptor
│       ├── compact_parser.py       # Compact text → SchemeDescriptor
│       ├── renderer.py             # CDXML rendering engine (6 layouts)
│       ├── render_scheme.py        # CLI entry point
│       ├── auto_layout.py          # Zero-effort JSON → CDXML
│       └── scheme_yaml_writer.py   # JSON → YAML layout decisions
│
├── # --- PRIVATE (not in package, stay at root) ---
├── run_pipeline.py             # Batch orchestrator
├── procedure_writer.py         # Lab book entry assembler
├── mass_resolver.py            # Structure-based mass determination
├── lcms_identifier.py          # LCMS species identification
├── lab_book_formatter.py       # Lab book output formatting
├── lcms_analyzer.py            # Single-file LCMS PDF parser
├── multi_lcms_analyzer.py      # Cross-file LCMS collation
├── lcms_file_categorizer.py    # LCMS file classification
├── discover_experiment_files.py # Experiment file discovery
├── query_status.py             # Pipeline status query
├── deploy_check.py             # Deployment verification
├── setup_wizard.py             # Interactive config wizard
├── prepare_reaction_scheme.py  # Deprecated (use scheme_polisher_v2)
│
├── # --- PROJECT FILES ---
├── pyproject.toml              # Package build config
├── LICENSE                     # MIT
├── NOTICE.md                   # Third-party attribution
├── README.md                   # GitHub-facing README
├── CLAUDE.md                   # Claude project instructions (this file)
├── pytest.ini                  # Pytest configuration
├── .gitignore
│
├── tests/                      # Test suite (pytest)
│   ├── conftest.py
│   ├── test_constants.py       # Unit tests for constants (41 tests)
│   ├── test_text_formatting.py # Unit tests for text formatting (29 tests)
│   ├── test_cdxml_utils.py     # Unit tests for CDXML utils (27 tests)
│   ├── test_rdkit_utils.py     # Unit tests for RDKit utils (23 tests)
│   ├── test_reagent_db.py      # Unit tests for reagent DB (28 tests)
│   ├── test_superatom_table.py # Unit tests for superatom table (20 tests)
│   ├── test_reaction_parser.py # Unit tests for reaction parser (60 tests)
│   ├── test_scheme_maker.py    # Unit tests for scheme maker
│   ├── test_smoke.py           # CLI smoke tests for core tools
│   ├── test_smoke_extended.py  # CLI smoke tests for extended tools
│   └── test_builder.py         # Integration tests (ChemDraw COM)
├── deploy/                     # Deployment infrastructure
├── experiments/                # ML tools and data extraction
│   ├── atom_mapping/           # RXNMapper (rxn-experiments env)
│   ├── byproduct_prediction/   # FlowER (flower env)
│   ├── role_classification/    # RXN Insight (rxn-experiments env)
│   ├── scheme_dsl/             # DSL docs + showcase (code moved to cdxml_toolkit/dsl/)
│   └── build_superatom_json.py # One-time data extraction scripts
├── exploration/                # Reverse-engineering, POC work
├── docs/                       # Reference material
└── local_llm/                  # Local LLM integration

# Test data lives OUTSIDE the project tree:
C:\Users\mic23\chem-test-data\
├── KL-CC-001/                  # Buchwald ELN export (.cdx, .csv, .pdf)
├── KL-CC-002/                  # SNAr ELN export (.cdx, .csv, .pdf)
├── polished-scheme/            # Test inputs for scheme_polisher
├── procedurefilltest/          # 14 experiments (Mitsunobu + Boc deprotection)
├── additional_test/            # User-only validation files
└── *.cdxml, *.png, *.pdf       # Standalone test files (Buchwald examples, etc.)
```

## File formats I work with
- `.cdxml` / `.cdx` — ChemDraw (CDXML is XML, CDX is binary)
- `.rdf` — SciFinder reaction export (V3000 MOL blocks + CAS + metadata)
- `.csv` — Findmolecule ELN export (semicolon-delimited, @TYPE rows)
- `.pdf` — LCMS reports (Waters MassLynx), ELN pages
- `.pptx` — Weekly reports with embedded ChemDraw OLE objects

## Tool inventory
See `README.md` for the full tool table with status and descriptions. Key points:
- **29 working tools + 1 experimental + 3 ML experimental + Scheme DSL (7 modules)** covering perception (LCMS, RDF, images, OLE), manipulation (building, normalizing, classifying, polishing, reaction layout), assembly (procedure writing, document templates), enrichment (ELN CSV → scheme annotation), output (rendering, converting), embedding (OLE into PPTX/DOCX), ELN integration (`findmolecule_export.py`), ML-assisted analysis (atom mapping, byproduct prediction, role classification), and DSL-based scheme rendering.
- **Two scheme processing pipelines exist:**
  - **Polisher pipeline** (`scheme_polisher_v2.py`): Surgically modifies ELN-exported CDXML. COM-free polishing: per-fragment RDKit/ChemScript cleanup + bond normalization + classification + alignment + ELN enrichment + reaction_cleanup layout. Default alignment and classification use RXNMapper; falls back to RDKit MCS. `--eln-csv` enables enrichment. See "Canonical workflow" above.
  - **DSL pipeline** (`experiments/scheme_dsl/`): Builds schemes from scratch via `reaction_parser JSON → scheme_yaml_writer → renderer → CDXML`. No ChemDraw COM needed at any stage. Supports 6 layout patterns, 3 arrow styles, run arrows, compound labels, letter conditions. See "Alternative workflow" above.
- **`scheme_aligner.py` is experimental.** Standalone MCS-based orientation alignment (separate from `alignment.py`). Uses RDKit MCS for structure alignment without the full polishing pipeline.
- **ML experimental tools** live in `experiments/` and require dedicated conda environments. `rxn_atom_mapper.py` wraps IBM RXNMapper for atom mapping and reagent role classification (requires `rxn-experiments` env). `flower_predictor.py` wraps FlowER for byproduct prediction via beam search (requires `flower` env). `rxn_role_classifier.py` wraps RXN Insight for reaction role classification (requires `rxn-experiments` env). All communicate via JSON-over-subprocess and degrade gracefully when their envs are missing.
- **`procedure_writer.py` was split** into three support modules: `mass_resolver.py` (structure-based mass determination), `lcms_identifier.py` (species identification from ion m/z), `lab_book_formatter.py` (output section generation). The CLI interface is unchanged — procedure_writer.py orchestrates the three modules.
- **`cdxml_combiner.py` was removed (deprecated).** Multi-step scheme combination had too many positioning edge-cases. Re-engineered as `scheme_merger.py`.
- **`cdxml_to_image_rdkit.py` is backup only.** Single molecules, no reactions, limited quality. Use `cdxml_to_image.py` (ChemDraw COM) as primary.

## Shared modules (v0.4)

Seven shared modules provide common functionality to all tools:

### constants.py — Centralized constants
All hardcoded magic numbers previously scattered across tool scripts. Organized into sections:
1. **ACS Document 1996 style** — `ACS_BOND_LENGTH` (14.40), `ACS_CHAIN_ANGLE` (120), font/size/face constants, `ACS_STYLE` dict, `CDXML_HEADER` template, `CDXML_MINIMAL_HEADER`, `CDXML_FOOTER`
2. **LCMS analysis** — `LCMS_RT_TOLERANCE` (0.02), `LCMS_MZ_TOLERANCE` (0.5), `LCMS_TREND_THRESHOLD` (0.2), `LCMS_MIN_SUMMARY_AREA` (2.0), column boundary, axis ticks, UV wavelength range
3. **Mass matching** — `MW_MATCH_TOLERANCE` (2.0 Da), `MW_MATCH_TOLERANCE_LOOSE` (5.0 Da), `MASS_TOLERANCE` (1.5 Da), `MIN_REPORT_AREA_PCT` (20.0), `MIN_SIGNIFICANT_AREA` (2.0%)
4. **Layout** — `LAYOUT_ABOVE_GAP` (8pt), `LAYOUT_BELOW_GAP` (4pt), `LAYOUT_HANGING_LABEL_GAP` (16pt), fragment/inter-fragment gap constants
5. **Image/structure** — `EXPAND_SCALE_BOND` (10.0)

**To add a new constant:** Add it to the appropriate section in `constants.py` with a comment noting its purpose and which file(s) it was originally in. Import it in the consuming tool.

### text_formatting.py — Chemical text formatting
Handles two chemistry-specific typographic conventions for ChemDraw CDXML `<s>` elements:
- **Subscript digits in chemical formulas** — "CH3OH" → "CH₃OH", "Pd2(dba)3" → "Pd₂(dba)₃". Plain numbers (temperatures, durations, percentages) are left as normal text.
- **Italic prefixes in IUPAC nomenclature** — "n-BuLi" → "*n*-BuLi", "tert-BuOH" → "*tert*-BuOH", "N-Boc" → "*N*-Boc".

Key functions: `needs_subscript()`, `split_italic_prefix()`, `build_formatted_s_xml()`. Backward-compatible aliases (`_build_formatted_s_xml`, `_build_subscripted_s_xml`) are provided.

Previously duplicated in `scheme_polisher.py` and `reaction_from_image.py`.

### cdxml_utils.py — CDXML geometry and IO utilities
Extracted from duplicated code across `reaction_cleanup.py`, `eln_enrichment.py`, and `scheme_polisher.py`.

Key functions:
- `fragment_bbox(frag)` — Atom-only bounding box using direct-child `<n>` atom `p` positions (NOT recursive, NOT XML BoundingBox). Falls back to XML BoundingBox only if no atoms have `p`.
- `fragment_bbox_with_label_extension(frag)` — Like `fragment_bbox` but extends `max_y` by `LAYOUT_HANGING_LABEL_GAP` (16pt) when the fragment has a hanging N-H/P-H label. Use this for layout spacing calculations.
- `fragment_centroid(frag)` — Center of `fragment_bbox`.
- `fragment_bottom_has_hanging_label(frag)` — Detects N/P at bottom with ≤2 bonds (hanging H label, needs extra gap).
- `recompute_text_bbox(t_elem)` — Estimates BoundingBox for `<t>` elements from anchor position and text content.
- `build_id_map(parent)` — Builds `{id: element}` dict for all descendants (recursive). Note: `reaction_cleanup.py` keeps its own `_build_id_map()` which is intentionally shallow (direct children only).
- `parse_cdxml(path)` / `write_cdxml(tree, path)` — IO with DOCTYPE preservation.

**Key design decision:** Atom-only bounding boxes are used because XML BoundingBox attributes are unreliable for `NodeType="Fragment"` abbreviation groups (OTs, Boc) which report the expanded inner structure, not the visible abbreviation.

### rdkit_utils.py — RDKit-based CDXML fragment utilities
Complements `cdxml_utils.py` (pure XML geometry) with RDKit-powered chemical operations on CDXML fragments:

Key functions:
- `frag_to_mol(frag_elem)` — CDXML `<fragment>` → RDKit Mol with atom metadata. Abbreviation groups (NodeType="Fragment") become dummy atoms (element 0).
- `frag_to_smiles(frag_elem)` — CDXML `<fragment>` → canonical SMILES string.
- `frag_to_mw(frag_elem)` — CDXML `<fragment>` → molecular weight (average MW via `Descriptors.MolWt`). Handles abbreviation groups (OTs, Boc, etc.) via `superatom_table.py` lookup — core MW + abbreviation MW - 1.008 Da per attachment bond.
- `frag_to_molblock(frag_elem)` — CDXML `<fragment>` → MOL block with CDXML coordinates preserved.
- `cleanup_fragment_rdkit(frag_elem)` — 2D coordinate cleanup with Kabsch orientation preservation. Generates new 2D coords via RDKit, then uses Kabsch rotation to align them back to the original orientation. For multi-component fragments (salt products like amine+HCl), disconnected counterions are repositioned to preserve their original offset from the main structure (RDKit's `Compute2DCoords` places disconnected components arbitrarily).
- `set_cdxml_conformer(mol, atoms_data)` — Sets RDKit conformer from CDXML atom coordinates (with y-axis flip).
- `rdkit_default_bond_length()` — Returns RDKit's default 2D depiction bond length.
- `avg_bond_length_from_atoms(atoms_data)` — Average bond length from CDXML atom coordinate data.

All RDKit imports are lazy — the module can be imported even if RDKit is not installed (functions raise ImportError at call time).

Uses `ACS_BOND_LENGTH` from `constants.py` for scale conversion.

### alignment.py — Shared alignment strategies
Centralized alignment geometry and orientation-matching strategies for reaction scheme polishing. Three independent layers, all consumed by `scheme_polisher.py` and `scheme_polisher_v2.py`:

**Layer 1 — Geometry primitives** (stdlib only):
- `fragment_centroid(frag)` — Centroid from direct-child node positions.
- `get_visible_carbon_positions(frag)` — Carbon backbone positions for Kabsch alignment (excludes heteroatoms, abbreviation groups).
- `compute_rigid_rotation_2d(source, target)` — Kabsch 2D rotation angle from matched point sets.
- `rotate_fragment_in_place(frag, angle, cx, cy)` — Rotates all coordinates in a fragment.
- `translate_subtree(elem, dx, dy)` — Translates all `p` and `BoundingBox` attributes.
- `make_abbrev_dummy_copy(frag)` — Deep copy with abbreviation groups replaced by dummy atoms.

**Layer 2 — Kabsch alignment** (ChemScript + RDKit, lazy imports):
- `kabsch_align_fragment_to_product(frag, product_frag, ...)` — Single fragment alignment via MOL export + MCS + rigid rotation.
- `kabsch_align_to_product(tree, ...)` — Aligns all reactant/reagent fragments in a scheme to the product.

**Layer 3 — RDKit MCS alignment** (RDKit only, lazy imports):
- `rdkit_align_to_product(tree, ...)` — Per-fragment alignment via `GenerateDepictionMatching2DStructure` with scale conversion (CDXML pts ↔ RDKit units).
- `align_product_to_reference(tree, ref_tree, ...)` — Aligns a product to a known-good reference structure from another CDXML file.
- `rxnmapper_align_to_product(tree, mapped_rxn, ...)` — Uses RXNMapper atom maps (from `rxn_atom_mapper.py`) for alignment instead of MCS. Falls back to RDKit MCS, then Kabsch.

**Scale conversion:** CDXML bond lengths are ~34.4 pts; RDKit uses ~1.5 units. All Layer 3 functions handle this conversion internally. Y-axis is flipped (CDXML down, RDKit up).

**Dependencies:** `constants.py` (`ACS_BOND_LENGTH`, `CDXML_MINIMAL_HEADER`), `cdxml_utils.py` (`write_cdxml`). RDKit and ChemScript are lazy-imported.

### superatom_table.py — Superatom fragment abbreviation lookup
Maps ChemDraw abbreviation group labels (OTs, Boc, Me, Et, etc.) to SMILES strings for MW calculation of fragments containing `NodeType="Fragment"` abbreviation nodes. Data backed by `superatom_data.json` (~2,850 entries generated from ChemScanner's `superatom.txt`) supplemented by RDKit built-in abbreviations (~40 entries).

Key functions:
- `get_superatom_table()` — Returns the label→SMILES dict (singleton, built on first call). Keys are lowercase.
- `lookup_smiles(label)` — Look up a superatom label, return its SMILES or None.
- `lookup_mw(label)` — Look up a superatom label, return its standalone MW or None. Requires RDKit. The returned MW is for the standalone fragment; callers subtract 1.008 per attachment bond.
- `get_abbrev_label(node)` — Extract visible abbreviation label text from a CDXML `<n NodeType="Fragment">` element (from `<t><s>...</s></t>` child).

Data sources (in priority order):
1. `superatom_data.json` — ~2,850 entries from ChemScanner's superatom.txt (forward + reverse forms). Generated by `experiments/build_superatom_json.py`.
2. RDKit `rdAbbreviations.GetDefaultAbbreviations()` — ~40 entries (BSD). Only adds entries not already in JSON.

## Inter-tool dependencies
Several tools import from each other. The seven shared modules are imported by most tools:
- `constants.py` is imported by essentially everything — any tool that references ACS style, LCMS defaults, layout gaps, or mass matching tolerances
- `text_formatting.py` is imported by `scheme_polisher.py`, `reaction_from_image.py`, `eln_enrichment.py`
- `cdxml_utils.py` is imported by `reaction_cleanup.py`, `eln_enrichment.py`, `scheme_polisher.py`, `scheme_aligner.py`, `alignment.py`
- `rdkit_utils.py` is imported by tools needing RDKit-CDXML bridge (fragment→mol, SMILES, MW, cleanup)
- `alignment.py` is imported by `scheme_polisher.py` and `scheme_polisher_v2.py` (Kabsch, RDKit MCS, and RXNMapper alignment strategies)
- `superatom_table.py` is imported by `rdkit_utils.py` (for abbreviation group MW calculation)
- `reagent_db.py` is the shared reagent database (two-tier: curated + ChemScanner) — imported by `reactant_heuristic.py`, `scheme_polisher.py`, `reaction_from_image.py`, `eln_enrichment.py`, `reaction_parser.py`
- `scheme_polisher.py` imports from `reactant_heuristic.py`, `reagent_db.py`, `chemscript_bridge.py`, `reaction_from_image.py`, `eln_cdx_cleanup.py`, `alignment.py`
- `scheme_polisher_v2.py` imports from `scheme_polisher.py`, `alignment.py`, `reaction_cleanup.py`, `eln_enrichment.py` (when `--eln-csv` is provided)
- `eln_enrichment.py` imports from `procedure_writer.py` (CSV parsing), `reagent_db.py`, `text_formatting.py`, `cdxml_utils.py`, `reaction_cleanup.py`, `chemscript_bridge.py` (MW via SMILES), RDKit
- `reaction_from_image.py` imports from `reagent_db.py`, `structure_from_image.py`, `cdxml_builder.py`, `chemscript_bridge.py`, `cas_resolver.py`
- `reactant_heuristic.py` imports from `reagent_db.py`, `chemscript_bridge.py`, `cas_resolver.py`, `experiments.atom_mapping.rxn_atom_mapper` (when `use_rxnmapper=True`)
- `multi_lcms_analyzer.py` imports from `lcms_analyzer.py`
- `procedure_writer.py` orchestrates three support modules: `mass_resolver.py` → `lcms_identifier.py` → `lab_book_formatter.py`. Also imports `discover_experiment_files.py`, `lcms_analyzer.py`, `constants.py`
- `mass_resolver.py` imports from `constants.py`, `chemscript_bridge.py` (optional), RDKit (optional)
- `lcms_identifier.py` imports from `mass_resolver.py`, `lcms_analyzer.py`, `multi_lcms_analyzer.py`, `constants.py`
- `lab_book_formatter.py` imports from `mass_resolver.py`, `lcms_identifier.py`, `constants.py`
- `reaction_cleanup.py` imports from `cdxml_utils.py` and `constants.py` only (no ChemDraw COM needed)
- `scheme_merger.py` imports from `cdxml_utils.py`, `rdkit_utils.py`, `alignment.py` (lazy), `constants.py` (no ChemDraw COM needed)
- `reaction_parser.py` imports from `constants.py`, `reagent_db.py`, `cdxml_utils.py`, `rdkit_utils.py`, `reactant_heuristic.py` (lazy), `procedure_writer.py` (CSV parsing), `cas_resolver.py` (lazy), `experiments.atom_mapping.rxn_atom_mapper` (lazy), `experiments.role_classification.rxn_role_classifier` (lazy), RDKit (lazy)
- **Scheme DSL modules** (`experiments/scheme_dsl/`):
  - `renderer.py` imports from `constants.py` (ACS style values), `text_formatting.py` (`build_formatted_s_xml`), `.schema` (dataclasses). Uses RDKit for SMILES → 2D coords.
  - `render_scheme.py` imports from `.parser`, `.compact_parser`, `.scheme_yaml_writer`, `.renderer`
  - `auto_layout.py` imports from `.schema`, `.renderer`. Adds project root to sys.path to access `constants.py`.
  - `scheme_yaml_writer.py` imports from `.schema` only (reads JSON files directly)
  - `parser.py` and `compact_parser.py` import from `.schema` only

Library modules live in the `cdxml_toolkit` package (using relative imports). Private root-level scripts import from the package via `from cdxml_toolkit.X import ...`. The package is installed in development mode (`pip install -e .`).

## LCMS specifics
- Instrument: Waters Acquity UPLC with QDa detector, MassLynx software
- Instrument names: `PPIMSA05`, `PPIMSA13`, etc. — extracted from line after "Paraza Pharma Inc" in PDF header
- PDF extraction: use `pdfplumber` — works reliably on MassLynx reports
- Peak tables (pages 1-2): parsed from full-page extracted text. Three detectors: TAC (all-wavelength), 220nm (Ch1), 254nm (Ch2). Each peak has raw area + area%.
- Mass spectra + UV (pages 3+): two-column layout parsed with spatial word-level extraction (`extract_words()`) to prevent column interleaving. Column boundary at x=306 (half of 612pt page width).
- Duplicate peak numbers: when a wavelength table has two rows with the same peak number at different RTs, they get suffixed as `"2a"`, `"2b"`. All `peak_num` values are strings.
- m/z extraction: MassLynx labels use 1 decimal place. `_extract_mz_values()` splits joined labels (e.g. `569.1814.6874.9` → `[569.1, 814.6, 874.9]`) using `\d+\.\d` regex, which also filters out integer y-axis labels.
- UV lambda-max: wavelength values from `3:UV Detector` panels, filtered to 150-400nm range excluding axis ticks (multiples of 50).
- Common adducts: [M+H]+, [M-H]-, [M+Na]+, [M+formate]-
- Watch for in-source fragmentation (e.g. Boc loss gives -99 Da in ESI+)
- Run time: `Time:HH:MM:SS` extracted from PDF header (line 4 of text). Used by `multi_lcms_analyzer.py` for chronological sorting.

## multi_lcms_analyzer.py — Cross-File LCMS Collation
Collates peaks across multiple LCMS PDF files from the same reaction. Imports `lcms_analyzer.py` for single-file parsing. Also used by `procedure_writer.py` for tracking LCMS analysis.

**CLI:**
```bash
python multi_lcms_analyzer.py file1.pdf file2.pdf file3.pdf ... \
    --rt-tolerance 0.02 --mz-tolerance 0.5 --trend-threshold 0.2
python multi_lcms_analyzer.py *.pdf --out-of-order --hide-other-ions -o report.txt
python multi_lcms_analyzer.py *.pdf --json -o report.json
```

**Peak matching algorithm:**
- Peaks matched across files by RT proximity (default ±0.02 min) with UV ratio (area_220nm / area_254nm) as discriminant
- UV ratio compatibility: only rejects when both ratios are finite and outside 2× of each other; missing data = inconclusive (not rejected)
- Within each file, largest peaks matched first for stability
- Unmatched peaks become new compounds
- Canonical RT: majority vote (mode) of rounded values across files

**Ion merging:**
- m/z values clustered within ±0.5 Da (running mean)
- Split into recurring (seen in ≥2 files) vs single-observation ("other ions")
- Rank preserved: position in `top_ions` list (0 = base peak)
- Sorted by: occurrences (desc) → best rank → m/z

**UV lambda-max consensus:**
- Collected from matched `ChromPeak.uv_lambda_max` across all files
- Deduplicated by grouping within 10 nm; mean of each cluster reported
- Half-up rounding for display (222.5 → 223, not banker's rounding)

**Trend analysis:**
- Computes area% trend as (last − first) / max over full timeline
- Files where compound is absent treated as 0% (between first and last observation)
- Area% fallback: TAC → 220nm → 254nm
- Threshold: ±20% change = increasing/decreasing (default)
- Outlier/blank files excluded from trend (but still shown in timeline)
- Ambiguous-timing files (e.g. "beforeadd", sort_key=500) excluded from trend unless run-time sorting resolves their position

**File sorting:**
- **Default:** Sort by actual LCMS acquisition timestamp extracted from PDF (`Date:` + `Time:` in header). Resolves ambiguous filenames automatically.
- **`--out-of-order`:** Sort by filename heuristics (sort_key from `categorize_lcms_file`). Files with ambiguous timing flagged and excluded from trend.
- Files grouped by instrument (unless `--ignore-instrument`)

**Outlier detection:**
- Peak count < 40% of median → likely blank
- Single peak > 95% TAC area with ≤5 total peaks → likely blank
- Excluded files marked `**EXCLUDED**` in output

**Output sections:**
1. **Reaction summary** (top) — increasing/decreasing/stable compounds sorted by max area, with recurring ions. Minor compounds below `--min-summary-area` (default 2%) hidden with count.
2. **Compound details** — per-compound: trend, UV ratio, UV λmax, recurring ions, other ions, area% timeline across all files.

**Key flags:**
- `--rt-tolerance` (default 0.02): RT matching window in minutes
- `--mz-tolerance` (default 0.5): ion clustering window in Da
- `--trend-threshold` (default 0.2): fraction change for increasing/decreasing
- `--min-summary-area` (default 2.0): hide compounds below this max area% from summary
- `--hide-other-ions`: suppress single-observation ions in compound details
- `--out-of-order`: use filename heuristic sorting instead of PDF timestamps
- `--ignore-instrument`: analyze all files together regardless of instrument
- `--json`: structured JSON output

**Test data:** `test_data/procedurefilltest/KL-7001-incomplete/LCMS files/KL-7001-004-*.pdf` (6 files: beforeadd, 0min, 10min, 30min, add50mgDEAD-10min, add50mgDEAD-20min). Reference output in `KL-7001-004-multi-lcms-output.txt`.

## findmolecule_export.py — Automated ELN Export

Downloads lab book exports (CDX, CSV, RXN files) from Findmolecule ELN by replicating browser HTTP requests. Reverse-engineered from HAR capture + `export.js` analysis (Feb 2026). Lives in `exploration/findmolecule exploration/` (not at project root — this tool runs on the **work machine**, not the development machine).

**Dependencies:** `requests` only (stdlib + requests). No ChemDraw, no RDKit.

**Credentials:** Stored in `fm_credentials.json` (gitignored) next to the script. Created via `--save-credentials` interactive prompt. Resolution order: `--credentials FILE` → `fm_credentials.json` → `FM_USER`/`FM_PASSWORD` env vars → `--user`/`--password` flags → interactive prompt.

**CLI:**
```bash
python findmolecule_export.py --save-credentials              # one-time setup
python findmolecule_export.py --list                          # show available lab books
python findmolecule_export.py --lab-book KL-7001              # export one lab book
python findmolecule_export.py --lab-book KL-7001 KL-1001     # export multiple
python findmolecule_export.py --all -o ./exports/             # export all to directory
python findmolecule_export.py --lab-book KL-7001 --options cdx csv rxn separate-csv
```

**How it works:**
1. `POST /welcome/loginAuthentication` — form login with username, password, FingerprintJS hash, timezone
2. Extract userId from `<option value="XXXXX" selected>` in the post-login page's user dropdown (position ~90K in a 130K HTML page — must search the full page, not a truncated snippet)
3. `GET /labBook/export/openExportLabBooks?userId=XXXXX` — returns HTML with lab book checkboxes (`data-id` attributes)
4. `GET /labBook/export/exportLabBooks?ids=ID&delimiter=;&options=1,2,7,4` — returns ZIP download
5. `POST /labBook/export/cleanUpExportSession` — server cleanup

**Export option codes** (from HTML checkbox `value` attributes):
- 1=cdx, 2=csv, 3=attachments, 4=separate-csv, 5=pdf, 6=updated-only, 7=rxn, 8=all-users

**Known quirks:**
- The userId is buried ~90K into a 130K-char HTML page (in a `<select class="chosen-project-users">` with 145 options). The `selected` attribute marks the logged-in user — no manual userId input needed.
- FingerprintJS hash is sent at login. We generate a stable MD5 from the username. Findmolecule logs but does not enforce it (as of Feb 2026).
- Large exports may be split into multiple ZIP files via a polling mechanism (`/labBook/export/isExportReady?fileNo=N`).
- `fm_credentials.json` stores password in **plain text** — keep it private. Added to `.gitignore`.

**Exploration artifacts:** `exploration/findmolecule exploration/` also contains `fm_probe.py` (endpoint discovery), `analyze_har.py` (HAR file parser), `README_capture_guide.md` (HAR capture instructions), and debug HTML snapshots. These were used during reverse-engineering and are kept for reference.

## Findmolecule ELN export format
- CSV: semicolon-delimited, @TYPE rows (REACTANT, SOLVENT, EVENT, PRODUCT), columns for mass, mmol, equiv, MW, supplier
- CDX: process into CDXML via `eln_cdx_cleanup.py` (scales + styles + cleans reaction layout)
- Procedure HTML: `<a href="#reagent-row0">` links to reagent table rows

## ChemDraw CDXML conventions
- Carbon atoms: no Element attribute. Heteroatoms: Element number + NumHydrogens + NeedsClean="yes" + `<t>` label
- Coordinates in points (1/72 inch), y-axis down
- Reactions: `<fragment>` for molecules, `<arrow>`, `<t>` for conditions, `<scheme>`/`<step>` for metadata
- `ReactionStepObjectsAboveArrow` / `ReactionStepObjectsBelowArrow` reference element IDs

## procedure_writer.py — Lab Book Entry Assembler
Assembles a polished, copy-paste-ready lab book entry from raw experiment files. Reaction-agnostic — no hardcoded species or reaction-type-specific logic. Split into four files for maintainability:
- `procedure_writer.py` — CLI orchestrator, file discovery, CSV parsing (~520 LOC)
- `mass_resolver.py` — Structure-based mass determination via ChemScript + RDKit (~570 LOC)
- `lcms_identifier.py` — LCMS species identification from ion m/z values (~420 LOC)
- `lab_book_formatter.py` — Output section generation (PROCEDURE/CHARACTERIZATION/NOTES) (~690 LOC)

**Dependencies:** `discover_experiment_files.py` (file discovery), `lcms_analyzer.py`, `multi_lcms_analyzer.py`, `mass_resolver.py`, `lcms_identifier.py`, `lab_book_formatter.py`, `constants.py`. Optional: `chemscript_bridge.py`, RDKit — enable structure-based mass determination from CDX/RXN files; without them, falls back to CSV MW values.

**CLI:**
```bash
python procedure_writer.py --input-dir path/ --experiment KL-7001-004 --output result.txt
python procedure_writer.py --input-dir path/ --experiment KL-7001-004 --sm-mass 274 --product-mass 459
python procedure_writer.py --input-dir path/ --experiment KL-7001-004 --flower-json predictions.json
```

**Inputs:**
- Findmolecule ELN CSV (semicolon-delimited, @TYPE sections for REACTANT/SOLVENT/PRODUCT/ANALYSIS)
- LCMS PDFs from MassLynx (tracking, workup, purification, final product)
- FlowER predictions JSON (optional, from `run_pipeline.py` Phase 3.15 or `flower_predictor.py --output-json`)
- NMR PDFs from MestReNova (extracts `1H NMR (400 MHz, solvent) δ ...` data strings)
- CDX/RXN structure files (loaded via ChemScript to extract SMILES → RDKit exact masses)

**Output: three sections:**
1. **PROCEDURE** — setup text from CSV + auto-generated tracking narrative + inferred workup/purification
2. **CHARACTERIZATION** — LCMS annotation lines (tracking + purified product) + NMR data strings
3. **NOTES** — per-file SM/DP area% timeline, unidentified compounds ≥20%, missing data flags

**File discovery and categorization:**
- LCMS PDFs: searches experiment subdirs, `LCMS files/`, parent dir. NOT `DATA/` (reserved for NMR).
- NMR PDFs: scans `DATA/` directories for PDFs matching experiment name that contain NMR data strings (content-based detection via `\d+[A-Z]\s+NMR` regex).
- File categorization delegated to `multi_lcms_analyzer.categorize_lcms_file()`.

**Mass determination pipeline:**
1. Try CDX file via ChemScript `load_reaction()` → if fails, try RXN file
2. For each reactant/product: extract SMILES → `compute_masses()` returns `(neutral_mass, full_mass)` where neutral = largest fragment after salt splitting (for ESI adduct matching), full = entire molecule including counterions (for CSV MW matching). Note: `compute_masses()` and `build_adducts()` are public functions (backward-compatible aliases `_compute_masses`/`_build_adducts` still work)
3. Match each species to CSV reagent rows by MW (tries both neutral and full mass within 2 Da tolerance) to get the chemist's reagent name
4. SM identified by matching against CSV substrate MW; product labeled "DP" if single product
5. Fallback: if no CDX/RXN or no ChemScript/RDKit, uses CSV MW values directly

**LCMS analysis:**
- **Tracking:** Delegates to `multi_lcms_analyzer.analyze()` for cross-file compound matching, trend detection, and area% timelines. Compounds matched to expected species by ion m/z.
- **Purified product:** Single-file analysis via `lcms_analyzer.parse_report()`. Peaks matched to expected species by ion m/z.
- **Adduct matching priority:** [M+H]+/[M-H]- preferred over [M+Na]+/[M+formate]-. Among same priority, higher-intensity ion preferred. Among same intensity, ESI+ preferred.
- **Reporting threshold:** Only compounds with ≥20% area% (in any detector) appear in CHARACTERIZATION. Must also match an expected species ion.
- **UV lambda max:** Reported for every peak in the characterization section (from multi-LCMS consensus or single-report data).
- **Purity:** All three detectors reported (TAC, 220nm, 254nm).
- **Instrument name:** Kept as-is from PDF header (e.g. `PPIMSA05`, not reformatted).

**NMR extraction:**
- Regex matches `1H NMR (400 MHz, solvent) δ ...` patterns from pdfplumber text
- Handles 1H, 13C, 19F, 31P (any `\d+[A-Z] NMR` pattern)
- Two-level deduplication: within-file (MestReNova repeats data on each page) and cross-file (same compound analyzed in multiple PDFs)

**Known limitations:**
- LCMS conversion % can be inflated when SM/DP peaks co-elute with other compounds
- NMR extraction fails on spectrum-image-only PDFs (no data string in text layer)
- Workup narrative is basic (keyword-based); doesn't capture solvents or detailed procedures
- Byproducts not drawn in the reaction scheme (e.g. TPPO from Mitsunobu) appear as "unidentified" in notes

**Test data:** `test_data/procedurefilltest/KL-7001-incomplete/` (14 experiments, Mitsunobu + Boc deprotections).

## reaction_parser.py — Unified Reaction Semantic Layer
Parses ELN export files (any combination of CDX, CDXML, RXN, CSV) into a single persisted JSON descriptor listing every chemical species with canonical SMILES, role classification, display names, and mass data. The JSON output serves as a **single source of truth** for downstream tools (procedure_writer, scheme_merger, flower_predictor) instead of each tool re-parsing raw files independently.

**Key design decisions:**
- **Does NOT rely on `<step>` attributes** for role assignment. Arrow position-based detection: fragments right of arrow head = products, others = reactants/reagents.
- **Dual mass paradigm:** `exact_mass` (neutral fragment, for LCMS adducts) vs `exact_mass_full` (salt form, for CSV matching).
- **Text label recovery:** Resolves text labels from polished schemes (where structures were demoted) back to SMILES via reagent_db chain.
- **Condition text splitting:** `split_condition_text()` handles merged condition blocks (e.g. "DEAD\nPPh3\nTHF, rt, 16 h"), filtering out temperatures, times, and atmosphere tokens.

**CLI:**
```bash
python reaction_parser.py experiment.cdxml -o reaction.json
python reaction_parser.py experiment.cdxml --csv exp.csv --rxn exp.rxn --pretty
python reaction_parser.py --input-dir path/ --experiment KL-7001-004
```

**Python API:**
```python
from reaction_parser import parse_reaction, ReactionDescriptor
desc = parse_reaction(cdxml="scheme.cdxml", csv="exp.csv")
desc.to_json("reaction.json")
sm = desc.get_sm()    # SpeciesDescriptor for starting material
dp = desc.get_dp()    # SpeciesDescriptor for desired product
expected = desc.get_expected_species()  # ExpectedSpecies-compatible dicts for LCMS
```

**Extraction pipeline (10 steps):**
1. Input resolution (CDX, CDXML, RXN, CSV; CDX auto-converted)
2. Extract species from structural source (CDXML > RXN)
3. Resolve text labels (reagent_db + OPSIN + PubChem, `--no-network` disables online)
4. Match CSV data (3-pass: name match, MW match fragments, MW match text via DB SMILES)
5. Classify roles (3-tier: role_lookup + RXNMapper + MCS, then optional RXN Insight)
6. Identify SM and DP (CSV substrate flag > largest atom-contributing reactant; single/largest product)
7. Apply display names (SM/DP > reagent_db > CSV name > formula > SMILES)
8. Compute masses (salt splitting, adducts: [M+H]+, [M-H]-, [M+Na]+, [M+formate]-)
9. Deduplicate species by canonical SMILES
10. Build reaction SMILES and assemble JSON

**Dataclasses:** `SpeciesDescriptor` (per-species: SMILES, neutral SMILES, role, name, masses, adducts, CSV match data) and `ReactionDescriptor` (complete reaction: experiment, input files, species list, reaction SMILES, reaction class, warnings, metadata).

**Dependencies:** `constants.py`, `reagent_db.py`, `cdxml_utils.py`, `rdkit_utils.py`, `reactant_heuristic.py` (lazy), `procedure_writer.py` (CSV parsing), `cas_resolver.py` (lazy), RDKit (lazy). ML: `rxn_atom_mapper.py` (lazy), `rxn_role_classifier.py` (lazy).

**Flags:** `--no-rxnmapper`, `--no-rxn-insight`, `--no-network`, `--pretty`, `--json-errors`, `-v`/`--verbose`.

**Pipeline integration:** Phase 3.05 in `run_pipeline.py` (between scheme polishing and FlowER). Configured via `reaction_parser:` section in config.yaml. FlowER Phase 3.15 reads `reaction_smiles` from the JSON when available instead of re-extracting from CDXML.

**Test data:** `chem-test-data/KL-CC-001/` (Buchwald), `chem-test-data/KL-CC-002/` (SNAr).

## reaction_cleanup.py — Pure Python reaction layout
Pure Python replacement for ChemDraw COM's "Clean Up Reaction". No external dependencies.

**6 layout approaches:** `bbox_center`, `arrow_driven`, `proportional`, `compact`, `golden_ratio`, `chemdraw_mimic` (default). All use the same core helpers; they differ in gap ratios and arrow sizing. `--all` runs all approaches for side-by-side comparison.

**Key design decisions:**
- **Atom-only bounding boxes** for fragments via `cdxml_utils.fragment_bbox()` — XML BoundingBox attributes are unreliable (especially for `NodeType="Fragment"` abbreviation groups like OTs, Boc, which report the expanded inner structure, not the visible abbreviation).
- **Hanging label edge case:** Detected by `cdxml_utils.fragment_bottom_has_hanging_label()`. If N (Element=7) or P (Element=15) is the bottommost atom and has ≤2 explicit bonds, the H label renders as a vertical stack below the atom symbol (e.g. morpholine NH). These fragments get `LAYOUT_HANGING_LABEL_GAP` (16pt) above-arrow gap instead of `LAYOUT_ABOVE_GAP` (8pt). All gap values are defined in `constants.py`.
- **Text always below arrow:** All `<t>` elements (conditions text) are placed below the arrow with `LAYOUT_BELOW_GAP` (4pt) gap, regardless of their original `ReactionStepObjectsAboveArrow` / `ReactionStepObjectsBelowArrow` assignment.
- **Arrow length from content:** Arrow length is computed to fit the widest above/below-arrow content, with a minimum of 5× `ACS_BOND_LENGTH`.

## eln_enrichment.py — ELN CSV → Scheme Annotation
Enriches a polished CDXML reaction scheme with data from a Findmolecule ELN CSV. Used by `scheme_polisher_v2.py` when `--eln-csv` is provided.

**Two-phase design:**
- **Phase A** (before layout): Inject equivalents into text labels (e.g. "Cs2CO3" → "Cs2CO3 (2 eq.)") so arrow length computation accounts for the wider text.
- **Phase B** (after layout): Add above-arrow/side eq labels for fragment structures + run arrow with SM mass → product yield.

**Reactant repositioning (Step 4.5, before Phase A):**
When ≥2 fragment reactants sit left of the arrow and nothing is above it, the non-substrate is moved above the arrow (modifies `<step>` metadata only; physical repositioning is handled downstream by `reaction_cleanup`). Substrate identified by `is_substrate=True` or `equiv=1.0` in CSV, matched to a fragment by MW.

**CSV-to-scheme matching (3 passes):**
1. **Name match:** CSV reagent name resolved via `reagent_db.resolve_display()`, compared to scheme `<t>` text.
2. **MW match (fragments):** Fragment MW computed via ChemScript SMILES export + RDKit `Descriptors.MolWt()`, matched to CSV MW within `MW_MATCH_TOLERANCE` (2.0 Da, from `constants.py`). Closest match wins.
3. **MW match (text via reagent_db SMILES):** Scheme text looked up in reagent_db → SMILES → RDKit MW, compared to CSV MW. Closest match wins.

**MW computation:**
- **Primary:** Wrap `<fragment>` in minimal CDXML → ChemScript `write_data()` → SMILES → RDKit `Descriptors.MolWt()` (average MW, matches CSV values exactly).
- **Fallback:** Manual atom counting from CDXML Element attributes. Only estimates implicit H for carbon; heteroatom H requires `NumHydrogens` attribute. Typically ~2 Da error on molecules with NH/OH groups.

**Fragment bounding boxes:**
Uses `cdxml_utils.fragment_bbox_with_label_extension()` (atom-only bbox + hanging N-H/P-H label extension). No artificial padding.

**Dependencies:** `procedure_writer.py` (CSV parsing), `reagent_db.py`, `text_formatting.py` (`build_formatted_s_xml`), `cdxml_utils.py` (`recompute_text_bbox`, `fragment_bbox_with_label_extension`), `reaction_cleanup.py`, `chemscript_bridge.py` (SMILES export), RDKit.

**Test data:** `test_data/KL-CC-001/` (Buchwald), `test_data/KL-CC-002/` (SNAr), `test_data/procedurefilltest/KL-7001-incomplete/` (14 Mitsunobu + Boc deprotection reactions).

## scheme_merger.py — Merge Multiple Reaction Schemes
Merges two or more ELN-enriched CDXML reaction schemes (from `scheme_polisher_v2.py`) into a single scheme. Three relationship categories, auto-detected by default:

**Parallel:** Same reactants + products → stacked run arrows (e.g. 50 mg + 2 g of the same Mitsunobu).
**Sequential:** Product of step A = SM of step B → multi-step linear scheme. Shared intermediate drawn once.
**Unrelated:** Neither → placed side by side (`--adjacent`, default) or rejected (`--no-adjacent`).

**Auto-detect algorithm:** Builds pairwise classification matrix, clusters parallel reactions via Union-Find, chains sequential groups via topological sort (Kahn's algorithm), places unrelated groups adjacent.

**CLI:**
```bash
# Auto-detect (default — recommended)
python scheme_merger.py s1.cdxml s2.cdxml s3.cdxml s4.cdxml -o merged.cdxml

# Explicit mode override
python scheme_merger.py --mode parallel s1.cdxml s2.cdxml
python scheme_merger.py --mode sequential s1.cdxml s2.cdxml

# Equiv handling
python scheme_merger.py s1.cdxml s2.cdxml --no-equiv
python scheme_merger.py s1.cdxml s2.cdxml --equiv-range

# Alignment reference
python scheme_merger.py s1.cdxml s2.cdxml --ref-cdxml ref.cdxml

# Unrelated handling
python scheme_merger.py s1.cdxml s2.cdxml --no-adjacent  # error if any are unrelated

# Deprecated (still works with warning)
python scheme_merger.py --parallel s1.cdxml s2.cdxml
python scheme_merger.py --sequential s1.cdxml s2.cdxml
```

**Key functions:**
- `classify_pair(ps_a, ps_b)` — Returns `"parallel"`, `"sequential_ab"`, `"sequential_ba"`, or `"unrelated"`. Product-centric: compares SMILES sets, MW fallback for abbreviation groups.
- `auto_detect(schemes)` → `MergePlan` — Union-Find for parallel groups, Kahn's topological sort for sequential chain.
- `execute_merge_plan(schemes, plan)` — Parallel merge within groups, sequential merge across chain, adjacent placement for unrelated.
- `adjacent_place(trees)` — Side-by-side placement with 3x bond-length gap.

**Parallel validation:** `parallel_merge()` has `strict=True` by default — rejects mismatched products with `ValueError`. `strict=False` preserves old warn-only behavior.

**Sequential layout fixes:** Multi-component bbox centering (uses largest connected component for salt products), text stacking accounts for multi-line element heights, standalone equiv labels stripped (redundant with conditions text), SupersededBy cross-references correctly remapped across element subtrees.

**Alignment cascade (sequential mode):** Works backwards from the final product. Last step's product is aligned to `--ref-cdxml` (if provided) via MCS, then within-step structures aligned via RXNMapper. Earlier steps chain alignment through the shared intermediate.

**Pipeline integration:** Configured via `merge:` section in config.yaml. Called by `run_pipeline.py` in Phase 3.25 (after scheme polishing, before batch render). When `mode` is omitted in config, auto-detect is used. Supports `adjacent: false` to reject unrelated inputs in automated pipelines.

**Dependencies:** `cdxml_utils.py`, `rdkit_utils.py`, `alignment.py` (lazy), `constants.py`. No ChemDraw COM.

**Test data:** Reference parallel merge: `C:\Users\mic23\chem-test-data\KL-7001-004-with-009-scheme.cdxml`. Sequential merge reference: `C:\Users\mic23\chem-test-data\Report-scheme.cdx` (KL-CC-001 + KL-CC-002). Auto-detect test: KL-7001-004+009+011+013 (2 Mitsunobu + 2 Boc deprotection → hybrid parallel+sequential).

## scheme_aligner.py — MCS-Based Structure Orientation Alignment (Experimental)
Aligns all drawn structures in a CDXML reaction scheme to match the product's orientation using RDKit MCS + `GenerateDepictionMatching2DStructure`. Standalone — no inter-tool imports.

**Motivation:** In reaction schemes, the same substructure (e.g. a thienopyrimidine core) may be drawn at different orientations in the reactant vs product. This tool makes all structures visually consistent by aligning shared scaffolds.

**CLI:**
```bash
python scheme_aligner.py reaction.cdxml
python scheme_aligner.py reaction.cdxml -o aligned.cdxml --svg
python scheme_aligner.py reaction.cdxml --timeout 60
```

**How it works:**
1. Parse CDXML `<scheme><step>` to identify product, reactants, and above/below-arrow objects
2. Convert each `<fragment>` to an RDKit Mol (abbreviation groups like OTs become dummy atoms, element 0)
3. Product is the reference — its coordinates define the "truth" orientation
4. For each other fragment: find MCS with product (`FindMCS` with `ringMatchesRingOnly`, `completeRingsOnly`), skip if < 3 atoms
5. Scale CDXML coords to RDKit scale (`scale = rdkit_bond_length / cdxml_bond_length`), set reference conformer
6. Call `GenerateDepictionMatching2DStructure(target, ref, atom_map)` to align target to reference
7. Scale aligned coords back to CDXML space (`1/scale`), translate to preserve original fragment center
8. Write updated coordinates back to XML, recompute fragment BoundingBox

**Key technical details:**
- **Scale conversion is mandatory:** CDXML bond lengths are ~34.4 pts in coordinate space; RDKit uses ~1.5 units. Without conversion, `GenerateDepictionMatching2DStructure` compresses MCS atoms into a tiny area.
- **Atom map format:** `(ref_idx, target_idx)` pairs, passed as positional argument (not keyword). Despite some docs suggesting otherwise.
- **Y-axis flip:** CDXML y-axis points down; RDKit y-axis points up. Applied during scale conversion and writeback.
- **Center preservation:** After alignment, each fragment is translated back to its original centroid so the overall scheme layout is preserved.
- **Abbreviation groups:** `NodeType="Fragment"` nodes (e.g. OTs, Boc) are treated as dummy atoms (element 0) — they participate in connectivity but not MCS element matching.

**Limitations:**
- **No collision avoidance:** Reoriented fragments may overlap the arrow or other elements. Layout cleanup (via `reaction_cleanup.py`) is a separate step.
- **Single-step only:** Processes first product per step. Multi-product steps use only the first product as reference.
- **No stereochemistry preservation:** `GenerateDepictionMatching2DStructure` may alter wedge/dash bond directions.

**Dependencies:** RDKit only. No ChemDraw COM, no ChemScript, no inter-tool imports.

**Test data:** `test_data/KL-CC-002/KL-CC-002.cdxml` (SNAr, 4 fragments), `test_data/KL-CC-002/extrarottest.cdx` (SNAr with rotated structures, 3 fragments — convert to CDXML first).

## rxn_atom_mapper.py — RXNMapper Atom Mapping Wrapper (ML Experimental)

Wraps IBM's RXNMapper transformer model for reaction atom mapping and reagent role classification. Lives in `experiments/atom_mapping/`. Requires the `rxn-experiments` conda environment (RXNMapper + dependencies). Communicates via JSON-over-subprocess — the main LLMChem env calls `_rxnmapper_worker.py` as a subprocess in the `rxn-experiments` env.

**Public API (importable from LLMChem env):**
```python
from experiments.atom_mapping.rxn_atom_mapper import map_reaction, classify_roles

# Atom mapping — returns mapped SMILES with atom indices
result = map_reaction("CC(=O)O.OCC>>CC(=O)OCC")
# {"mapped_rxn": "[CH3:1]...", "confidence": 0.98}

# Role classification — determines which reactants contribute atoms to product
result = classify_roles("BrC1=CC=CC=C1.C1CCNCC1.O=C([O-])[O-].[Cs+].[Cs+]>>C1=CC=C(N2CCCCC2)C=C1")
# {"mapped_rxn": "...", "confidence": 0.98, "components": [
#   {"smiles": "BrC1=CC=CC=C1", "atom_contributing": true, "role": "reactant"},
#   {"smiles": "C1CCNCC1", "atom_contributing": true, "role": "reactant"},
#   {"smiles": "O=C([O-])[O-].[Cs+].[Cs+]", "atom_contributing": false, "role": "reagent"}
# ]}
```

**How classification works:**
1. RXNMapper produces atom-mapped SMILES where atom indices show which atoms move from reactants to products
2. For each reactant component, count how many of its mapped atoms appear in the product
3. If any mapped atoms contribute → `atom_contributing=True` (reactant); otherwise → `atom_contributing=False` (reagent/catalyst)

**Used by:**
- `reactant_heuristic.py` — Tier 1.5 classification (between role_lookup and MCS fallback)
- `alignment.py` — `rxnmapper_align_to_product()` uses atom maps for structure orientation alignment
- `scheme_polisher_v2.py` — default `--classify-method rxnmapper` and `--align-mode rxnmapper`

**Graceful degradation:** If the `rxn-experiments` conda env is not installed, all callers silently fall back to heuristic methods (role_lookup + RDKit MCS for classification, RDKit MCS for alignment). No errors raised.

**Dependencies:** Requires `rxn-experiments` conda env with `rxnmapper`, `torch`, `rdkit`. The worker script (`_rxnmapper_worker.py`) runs in this env. The public API module (`rxn_atom_mapper.py`) runs in LLMChem and uses subprocess to call the worker.

## flower_predictor.py — FlowER Byproduct Prediction Wrapper (ML Experimental)

Wraps the FlowER (Flow-based Enumerator of Reactions) model for predicting reaction byproducts via beam search. Lives in `experiments/byproduct_prediction/`. Requires the `flower` conda environment. Returns predicted byproducts as `ExpectedSpecies`-compatible dicts with adducts for LCMS matching.

**Public API:**
```python
from experiments.byproduct_prediction.flower_predictor import (
    predict_byproducts,              # from RXN file path
    predict_byproducts_from_smiles,  # from reaction SMILES string
)

# From SMILES (used by run_pipeline.py Phase 3.15)
species = predict_byproducts_from_smiles(
    "CC(=O)O.OCC>>CC(=O)OCC",
    max_depth=10, beam_size=3, device="auto"
)
# Returns: [{"name": "BP-C2H4O", "role": "byproduct", "exact_mass": 44.026,
#            "smiles": "CCO", "mw": 46.07, "formula": "C2H4O",
#            "adducts": {"[M+H]+": 45.033, "[M-H]-": 43.018}}]

# From RXN file (used by flower_predictor CLI)
species = predict_byproducts(rxn_path="reaction.rxn", csv_path="eln.csv")
```

**Features:**
- **Result caching** — SHA256 of reaction SMILES + beam params → cached JSON in `.cache/flower/`
- **Graceful failure** — returns empty list if flower env unavailable
- **MW filtering** — skips tiny fragments (MW < 50), polymeric artifacts (MW > 2× product), known species
- **Adduct computation** — ESI+/ESI- adducts for LCMS matching
- **Device auto-detection** — `device="auto"` uses CUDA if available, CPU otherwise
- **JSON output** — `--output-json` saves predictions as a JSON sidecar file for downstream consumption

**Used by:**
- `run_pipeline.py` — Phase 3.15 calls `flower_predictor.py --smiles` after scheme polishing, saves JSON sidecar
- `procedure_writer.py` — loads FlowER JSON via `--flower-json` for LCMS species matching
- `mass_resolver.py` — `--predict-byproducts` flag triggers inline FlowER prediction (standalone use)

**CLI:**
```bash
python flower_predictor.py --rxn path/to/reaction.rxn --csv path/to/eln.csv --pretty
python flower_predictor.py --smiles "BrC1=CC=CC=C1.C1CCNCC1>>C1=CC=C(N2CCCCC2)C=C1" -o predictions.json
python flower_predictor.py --rxn reaction.rxn --device cpu --no-cache
```

**Dependencies:** Requires `flower` conda env with FlowER model + dependencies. The worker script (`_flower_worker.py`) runs in this env. Subprocess timeout: 30 minutes (FlowER beam search can be slow).

## run_pipeline.py — Batch Pipeline Orchestrator

Orchestrates the full experiment processing pipeline across multiple experiments. Wraps `scheme_polisher_v2.py`, `multi_lcms_analyzer.py`, and `procedure_writer.py` with batch COM sessions, YAML configuration, experiment discovery, and ELN export integration.

**CLI:**
```bash
python run_pipeline.py --config deploy/config.yaml
python run_pipeline.py --config deploy/config.yaml --experiments KL-7001-008 KL-7001-009
python run_pipeline.py --config deploy/config.yaml --range KL-7001-008 KL-7001-014
python run_pipeline.py --config deploy/config.yaml --skip-export --input-dir ./exports/extracted
python run_pipeline.py --config deploy/config.yaml --resume   # skip completed phases
python run_pipeline.py --validate
```

**Pipeline phases:**
1. **Phase 1: ELN Export** — calls `findmolecule_export.py` (if at root) to download CDX/CSV/RXN from Findmolecule. Skippable with `--skip-export`.
2. **Phase 2: Experiment Discovery** — scans directories for experiment files (CDX, CSV, LCMS PDFs, NMR PDFs).
3. **Phase 2.5: Batch CDX → CDXML** — converts all CDX files in one ChemDraw COM session via `cdx_converter.py --batch`.
4. **Phase 3: Per-experiment processing** — for each experiment:
   - **3a.** Scheme polishing via `scheme_polisher_v2.py` (subprocess) — **OR** DSL rendering (see below)
   - **3.05** Reaction parsing via `reaction_parser.py` (subprocess). Parses polished CDXML + CSV + RXN into `{experiment}-reaction.json` — canonical species list with SMILES, roles, names, masses. Configured via `reaction_parser:` config section.
   - **3.15** FlowER byproduct prediction (optional, `flower.enabled: true`). Reads `reaction_smiles` from reaction parser JSON (Phase 3.05) when available, falls back to CDXML extraction. Saves predictions as `{experiment}-flower.json` sidecar file.
   - **3b.** Multi-LCMS analysis via `multi_lcms_analyzer.py` (subprocess). Tracking files only; purification/final files excluded from trend.
   - **3c.** Procedure writing via `procedure_writer.py` (subprocess). Receives `--flower-json` if Phase 3.15 produced predictions.
5. **Phase 3.5: Batch Render** — renders all polished schemes to PNG in one COM session via `cdxml_to_image.py --batch` (optional, controlled by `scheme.render` config).

**Scheme renderer selection (`scheme.renderer` config):**
- `"polisher"` (default) — uses `scheme_polisher_v2.py`. Execution order: scheme polishing (3a) → reaction parsing (3.05). The polished CDXML is passed to reaction_parser as input.
- `"dsl"` — uses the Scheme DSL renderer. Execution order is **reversed**: reaction parsing (3.05) runs **first** (with `scheme_cdxml=None`), then `render_scheme.py --from-json` consumes the JSON (Phase 3a-dsl). Config option `scheme.run_arrows: false` suppresses run arrows.

```yaml
# In config.yaml:
scheme:
  renderer: "dsl"          # "polisher" (default) or "dsl"
  run_arrows: true          # include SM mass → product yield arrows (DSL only)
```

**Key design decisions:**
- **`script_dir`** is always `os.path.dirname(os.path.abspath(__file__))` — the project root. All subprocess calls use `os.path.join(script_dir, "tool.py")`.
- **No inter-tool imports** — all tools are called as subprocesses. Only `lcms_file_categorizer` and `lcms_analyzer` are imported directly for lightweight file categorization and validation.
- **ML config propagation** — `RXN_EXPERIMENTS_PYTHON`, `FLOWER_PYTHON`, `FLOWER_DIR` environment variables passed to subprocesses.
- **LCMS file filtering** — non-standard PDFs (not Waters MassLynx reports) are filtered via `lcms_analyzer.is_waters_report()`.

**Experiment selection:** Supports `--experiments` (full names), `--range` (start/end), or interactive prompt. Flexible parsing: `"KL-7001-008 to KL-7001-014"`, `"KL-7001-008-014"`, `"5, 8, 9, 11"` (bare numbers with prefix from config).

**Validation mode:** `--validate` checks all required scripts exist, imports resolve, ChemScript is configured, and ML environments are available. Returns exit code 0 (pass) or 1 (fail).

**Deployment:** Invoked via `deploy/run.bat` which initializes conda, activates the `chem-pipeline` environment, and cds to the project root before calling `python run_pipeline.py %*`. Config file lives at `deploy/config.yaml` (generated from `deploy/config_template.yaml` by `setup_wizard.py`). Use `deploy_check.py` to verify what will work on a given machine (see "Deployment verification" section).

**Checkpoint/resume (`--resume`):** Skips phases whose output files already exist and are newer than input files. Checks:
- CDX→CDXML: skip if CDXML exists (handled by `batch_convert_cdx` inherently)
- Scheme polishing: skip if `{exp}-scheme.cdxml` exists and is newer than input CDX
- LCMS analysis: skip if `{exp}-lcms.json` or `{exp}-lcms.txt` exists
- Procedure: skip if `{exp}-procedure.txt` exists

Resumed phases log `RESUMED` status instead of `OK`. Re-running with `--resume` after a failure only re-processes the failed experiments/phases.

**Per-experiment config overrides:** The `experiment_overrides:` section in config.yaml allows overriding global settings per experiment:
```yaml
experiment_overrides:
  KL-7001-005:
    align_mode: rdkit        # override global scheme.align_mode
    skip_lcms: true          # skip LCMS analysis for this experiment
  KL-7001-008:
    ref_cdxml: custom-ref.cdxml  # experiment-specific reference
```
Supported keys: `align_mode`, `classify_method`, `layout_approach`, `ref_cdxml`, `chemscript_cleanup`, `merge_conditions`, `eln_enrichment`, `flower`, `renderer` (`"polisher"` or `"dsl"`), `skip_scheme`, `skip_reaction`, `skip_lcms`, `skip_procedure`, `skip_flower`. Overrides are applied via `_apply_experiment_overrides()` at the start of `process_experiment()`.

**Structured error output (`--json-errors`):** Key tools (`scheme_polisher_v2.py`, `procedure_writer.py`, `multi_lcms_analyzer.py`) support `--json-errors` flag that outputs structured JSON error objects to stderr on failure:
```json
{"error": "smiles_parse_failed", "detail": "Fragment 74: RDKit could not parse CDXML", "file": "input.cdxml"}
```
Error codes are machine-readable (e.g. `file_not_found`, `cdxml_parse_failed`, `smiles_parse_failed`, `enrichment_failed`, `pipeline_failed`, `no_parseable_files`, `analysis_empty`, `csv_parse_failed`). This lets agent orchestrators programmatically determine what went wrong without parsing free-text error messages.

**Dependencies:** `pyyaml` (config), subprocess calls to most tool scripts. `findmolecule_export.py` is optional (only needed for ELN export phase; lives in `exploration/` on dev machine, deployed to root on work machine).

## Unified Reagent Database (`reagent_db.py` — two-tier)

Two-tier reagent database with cascade lookup. Tier-1 (`reagent_abbreviations.json`, ~172 curated entries with roles) always wins. Tier-2 (`chemscanner_abbreviations.json`, ~5,837 entries from ChemScanner, no roles) is consulted only when tier-1 returns None. Replaced three previously separate data sources (a hardcoded SMILES→role dict, a regex→role pattern list, and a flat name→display JSON file).

**Tier-1: `reagent_abbreviations.json`** — Curated entries with roles, display names, SMILES, aliases. Med-chem-specific entries (Buchwald ligands, protecting group reagents, etc.) not found in any public database.

**Tier-2: `chemscanner_abbreviations.json`** — 5,837 entries (8,384 findable names via keys + aliases) extracted from ChemScanner's abbreviations.yaml, solvents.yaml, and superatom.txt. No role field. Sources: TCI/Sigma catalog reagents, chiral ligands, phosphines, organometallics, ionic liquids, fragment abbreviations. Generated by `experiments/build_chemscanner_json.py`.

**JSON format:** Both files use the same format. Each entry is keyed by lowercase primary name:
```json
{
    "cs2co3": {
        "display": "Cs2CO3",
        "role": "base",
        "smiles": "O=C([O-])[O-].[Cs+].[Cs+]"
    },
    "n-buli": {
        "display": "n-BuLi",
        "role": "base",
        "smiles": ["[Li]CCCC", "[Li][CH2]CCC"],
        "aliases": ["nbuli"]
    }
}
```
Fields: `display` (required), `role` (optional, tier-1 only), `smiles` (optional string or list), `aliases` (optional list of alternate lowercase keys).

**Cascade behavior:**
- `display_for_name()`, `display_for_smiles()`, `resolve_display()`, `entry_for_name()`, `entry_for_smiles()` — check tier-1 first, then tier-2
- `role_for_name()`, `role_for_smiles()`, `smiles_role_display()` — tier-1 only (ChemScanner has no roles)
- If `chemscanner_abbreviations.json` is missing, degrades gracefully to tier-1 only (warning to stderr)

**`reagent_db.py` API:**
```python
from reagent_db import get_reagent_db
db = get_reagent_db()                          # singleton, loads both JSONs once
db.display_for_name("cs2co3")                  # "Cs2CO3"  (tier-1)
db.display_for_name("hatu")                    # "HATU"    (tier-2 fallback)
db.role_for_name("tea")                        # "base" (alias of et3n, tier-1 only)
db.display_for_smiles("CCN(CC)CC")             # "Et3N" (RDKit canonicalization)
db.role_for_smiles("CCN(CC)CC")                # "base"
db.resolve_display("cs2co3")                   # "Cs2CO3"  (or original if unknown)
db.resolve_display("hatu")                     # "HATU"    (tier-2, not "hatu")
db.smiles_role_display("O=C([O-])[O-].[Cs+].[Cs+]")  # ("base", "Cs2CO3")
```

**Role categories:** catalyst, ligand, base, solvent, coupling_reagent, reducing_agent, oxidant, protecting_group, deprotecting_agent, acid, activating_agent, lewis_acid, drying_agent, halogenating_agent, fluorinating_agent, borylating_agent, additive, reductant, reagent.

**To add a new reagent:** Edit `reagent_abbreviations.json` directly. Add a new entry with at minimum `display`. Add `role` and `smiles` if known. To make an alternate spelling resolve to an existing entry, add it to that entry's `aliases` list. ChemScanner entries should not be edited — they are auto-generated.

## ChemScript (.NET library)
ChemScript is a 32-bit .NET library (`CambridgeSoft.ChemScript16.dll`) from ChemDraw 16. Access it only through `chemscript_bridge.py` (never import `_chemscript_server.py` directly).

**When to use ChemScript vs other tools:**
- Format conversion → prefer ChemScript over COM or openbabel
- Structure generation from names/SMILES → prefer ChemScript (PubChem for CAS lookup; web search as last resort)
- Coordinate cleanup → prefer ChemScript `.CleanupStructure()`
- Structure editing → ChemScript (atom/bond manipulation + substructure search)
- CAS → compound name → still use PubChem API

**Known limitations:**
- Name resolution fails for abbreviations (BINAP, Pd2dba3, Cs2CO3) — use CAS resolver or manual SMILES
- `ChemicalName()` throws on unusual charges (e.g. cesium salts) — bridge returns None
- `LoadFile()` fails on reaction schemes — use `ReactionData.LoadFile()` (bridge tries both)

## ole_embedder.py — ChemDraw OLE Embedding
Embeds one or more CDXML files as editable ChemDraw OLE objects into PowerPoint (.pptx) or Word (.docx). Double-clicking the embedded object opens ChemDraw for in-place editing. Uses binary OLE construction (no Office COM required — only ChemDraw COM for format conversion).

**CLI:**
```bash
# Single structure
python ole_embedder.py scheme.cdxml --pptx -o report.pptx
python ole_embedder.py scheme.cdxml --docx -o report.docx

# Multiple structures (one per slide / one per paragraph)
python ole_embedder.py s1.cdxml s2.cdxml s3.cdxml --pptx -o report.pptx
```

**How it works:**
1. Opens ChemDraw COM once, batch-converts all CDXML files to CDX binary + EMF preview
2. Builds OLE compound files (CFB format) from scratch — 10KB per structure
3. Computes display size from CDXML root BoundingBox (+ 2% scale factor)
4. Creates PPTX/DOCX scaffolding via python-pptx/python-docx, injects OLE + EMF into the ZIP

**OLE compound file structure (CFB):**
- Root Entry with ChemDraw CLSID (`41BA6D21-A02E-11CE-8FD9-0020AFD1F20C`)
- `CONTENTS` stream: raw CDX binary data (regular sectors)
- `\x01CompObj`, `\x01Ole`, `\x02OlePres000`: constant metadata (mini-stream)
- Layout matches Office COM's known-good output

**PPTX XML:** Each slide uses `mc:AlternateContent` wrapper with `mc:Choice` + `mc:Fallback` containing `p:pic` with EMF blip for the preview image. OLE objects are centered on each slide.

**DOCX XML:** Uses VML `v:shape` with `type="#_x0000_t75"` and `o:ole=""` attribute, plus `o:OLEObject` with matching ShapeID. Each object is an inline `w:object` in its own paragraph.

**Sizing:**
- Uses CDXML root `BoundingBox` attribute (tight content box in points) with a 2% scale factor
- This matches ChemDraw's native OLE display size (the size you get from manual copy-paste)
- EMF preview dimensions are NOT used for sizing (they're ~2x too large)
- Conversion: 1 pt = 12700 EMU

**Dependencies:** ChemDraw COM (CDX + EMF conversion), python-pptx, python-docx, lxml

**Known limitations:**
- Input must be CDXML (not CDX) — CDX files lack a parseable BoundingBox for sizing
- CDX files >64KB would exceed single-FAT CFB capacity (not an issue for typical structures)
- ChemDraw must not be open when running (COM automation conflict)

**Research & exploration:** Full technical findings in `exploration/FINDINGS.md`. POC scripts and diagnostic tools in `exploration/`.

## Scheme DSL — Text-Based Reaction Scheme Renderer (`experiments/scheme_dsl/`)

A declarative system for generating publication-ready CDXML reaction schemes from YAML or compact text input. The LLM specifies semantic content; the deterministic renderer handles all spatial layout. No ChemDraw COM needed — uses RDKit for 2D coordinate generation.

**Motivation:** LLMs cannot handle spatial layout but are excellent at structured text. The DSL separates semantic content (what to draw) from spatial rendering (where to draw it).

### schema.py — YAML Schema as Dataclasses

Defines the data model consumed by the renderer:

**`SchemeDescriptor`** (top-level):
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `source` | `Optional[str]` | `None` | Path to reaction_parser JSON (for structure resolution) |
| `structures` | `dict[str, StructureRef]` | `{}` | Structure definitions (SMILES, names, files) |
| `steps` | `list[StepDescriptor]` | `[]` | Reaction steps |
| `layout` | `str` | `"linear"` | Layout pattern keyword |
| `wrap` | `str` | `"repeat"` | `"repeat"`, `"serpentine"`, `"none"` |
| `steps_per_row` | `Optional[int]` | `None` | Auto-computed if omitted |
| `run_arrows` | `list[StepRunArrows]` | `[]` | Run arrows (SM mass → product yield) |
| `condition_key` | `Optional[dict]` | `None` | Letter conditions mapping |
| `sections` | `list[SectionDescriptor]` | `[]` | For `stacked-rows` layout |

**`StructureRef`**: `id` (required), `smiles`, `name`, `file`, `cdxml_id`, `label` (compound number displayed below).

**`StepDescriptor`**: `substrates`, `products` (lists of structure ref IDs), `above_arrow` / `below_arrow` (`ArrowContent` with `structures` and `text` lists), `yield_`, `arrow_style` (`"solid"`, `"dashed"`, `"failed"`), `number`, `id`.

**`ArrowContent`**: `structures` (list of ref IDs for drawn structures above arrow), `text` (list of condition text lines below arrow).

**`RunArrowEntry`**: `input_label` (e.g. `"2.15 g"`), `output_label` (e.g. `"1.60 g, 72% yield"`).

**`StepRunArrows`**: `step` (1-indexed), `runs` (list of `RunArrowEntry`).

**Valid layout keywords:** `linear`, `sequential`, `divergent`, `stacked-rows`, `numbered-parallel`, `convergent`.

### scheme_yaml_writer.py — Layout Decision Layer (JSON → YAML)

Reads a `reaction_parser` JSON file, makes automated layout decisions, and writes a YAML file the renderer can consume.

**Three decisions:**
1. **Structure or text?** — `atom_contributing` role = drawn structure; non-contributing (catalyst, ligand, base, solvent, etc.) = text label.
2. **Position?** — SM goes left of arrow; other atom-contributing reactants go above arrow as structures; reagents go below arrow as text.
3. **Ordering** — below-arrow text sorted by role priority (catalyst=10, ligand=15, base=20, ..., solvent=80).

**CLI:**
```bash
python experiments/scheme_dsl/scheme_yaml_writer.py reaction.json -o scheme.yaml
python experiments/scheme_dsl/scheme_yaml_writer.py reaction.json --layout sequential --no-run-arrows
```

**Python API:**
```python
from scheme_dsl.scheme_yaml_writer import write_scheme_yaml, build_scheme_yaml_dict
write_scheme_yaml("reaction.json", "scheme.yaml", layout="auto", include_run_arrows=True)
d = build_scheme_yaml_dict("reaction.json")  # returns dict, no file write
```

### render_scheme.py — CLI Entry Point

Single CLI for all input formats: YAML, compact text, or direct from reaction_parser JSON.

**CLI:**
```bash
# From YAML
python experiments/scheme_dsl/render_scheme.py scheme.yaml -o scheme.cdxml

# From compact text
python experiments/scheme_dsl/render_scheme.py scheme.rxn -o scheme.cdxml

# Direct from reaction JSON (two-step: JSON → YAML → CDXML)
python experiments/scheme_dsl/render_scheme.py --from-json reaction.json -o scheme.cdxml
```

**Format auto-detection:** `.yaml`/`.yml` → YAML; `.rxn`/`.scheme`/`.txt` → compact; content sniffing (`-->`, `..>`, `X>`) as fallback.

### auto_layout.py — Zero-Effort JSON → CDXML

Bypasses YAML entirely — builds a `SchemeDescriptor` directly in memory from a reaction_parser JSON, then renders to CDXML.

**CLI:**
```bash
python experiments/scheme_dsl/auto_layout.py reaction.json -o scheme.cdxml
```

**Python API:**
```python
from scheme_dsl.auto_layout import auto_layout, auto_layout_to_cdxml
scheme = auto_layout("reaction.json")  # returns SchemeDescriptor
path = auto_layout_to_cdxml("reaction.json", "scheme.cdxml")  # renders directly
```

### renderer.py — CDXML Rendering Engine

The core spatial engine. Takes a `SchemeDescriptor` and produces a complete CDXML document string.

**Public API:**
```python
from scheme_dsl.renderer import render, render_to_file
cdxml_str = render(scheme, yaml_dir="/path/to/yaml/dir")
render_to_file(scheme, "output.cdxml", yaml_dir="/path/to/yaml/dir")
```

**Supported layouts:**
| Layout | Internal function | Description |
|--------|-------------------|-------------|
| `linear` | `_layout_linear()` | Single-step reaction |
| `sequential` | `_layout_sequential()` | Multi-step in single row or wrapped |
| `serpentine` | `_layout_serpentine()` | Zigzag L→R, R→L with vertical arrows |
| `divergent` | `_layout_divergent()` | Fan-out to multiple products |
| `stacked-rows` | `_layout_stacked_rows()` | Independent rows stacked vertically |

**Key rendering constants:**
```python
_CHAR_WIDTH = 4.7       # average character width (proportional Arial)
_LINE_ADVANCE = 11.5    # line-to-line distance
_CAP_HEIGHT = 9.1       # baseline to top of uppercase
_FRAG_PAD_W = 3.3       # horizontal bbox padding beyond atom centers
_FRAG_PAD_H = 1.45      # vertical bbox padding beyond atom centers
ROW_GAP = 55.0           # vertical gap between wrapped rows
```

**Dependencies:** `constants.py` (ACS style values), `text_formatting.py` (`build_formatted_s_xml`), `schema.py` (dataclasses). RDKit for SMILES → 2D coordinates. No ChemDraw COM.

### parser.py and compact_parser.py — Input Parsers

**`parser.py`** — Parses YAML files into `SchemeDescriptor`. Validates layout, wrap, and arrow_style values. Handles `structures:`, `steps:`, `run_arrows:`, `condition_key:`, and `sections:` (for stacked-rows).

**`compact_parser.py`** — Parses compact text syntax into `SchemeDescriptor`. Syntax elements:
- Structure definitions: `id: {SMILES}`, `id: name "compound name"`, `id: file "path"`
- Reaction chain: `ArBr + Morph --> Product (72%)`
- Arrow types: `-->` (solid), `..>` (dashed), `X>` (failed), `-->|a|` (letter conditions)
- Conditions: indented `above:` / `below:` / `step N:` / `run:` blocks
- Directives: `@layout`, `@wrap`, `@steps_per_row`, `@conditions`

See `experiments/scheme_dsl/COMPACT_SYNTAX.md` for the full specification with EBNF grammar.

### Inter-module data flow

```
reaction_parser JSON
        │
        ├──→ scheme_yaml_writer.py ──→ YAML file ──→ parser.py ──┐
        │    (layout decisions)                                    │
        └──→ auto_layout.py ─────────────────────────────────────┤
             (zero-effort, no YAML)                                │
                                                                   ▼
compact text ──→ compact_parser.py ──→ SchemeDescriptor (schema.py)
                                              │
                                              ▼
                                       renderer.py ──→ CDXML string
```

### Test data and examples

- `experiments/scheme_dsl/showcase/` — 30+ example YAML and compact text files demonstrating all layout patterns
- `experiments/scheme_dsl/test_renderer.py` — integration tests for all layout engines
- `experiments/scheme_dsl/test_compact_syntax.py` — compact parser tests
- Reference reactions for testing: `chem-test-data/KL-CC-001/` (Buchwald), `chem-test-data/KL-CC-002/` (SNAr)

## Running tests
```bash
cd /c/Users/mic23/chem-tools-0.4
python -m pytest tests/ -v
```

288 tests total. Custom pytest marks registered in `pytest.ini`: `network` (tests requiring internet).

**Pure Python tests (no ChemDraw needed):**
- `test_constants.py` (41 tests) — verifies constant values and types, all 39 constants
- `test_text_formatting.py` (29 tests) — subscript detection, italic prefix splitting, XML output
- `test_cdxml_utils.py` (27 tests) — fragment bbox, hanging label detection, bbox with label extension, text bbox, id map, IO
- `test_rdkit_utils.py` (23 tests) — frag_to_mol, frag_to_smiles, frag_to_mw (including superatom-assisted MW), cleanup, scale helpers. Skipped if RDKit unavailable
- `test_reagent_db.py` (28 tests) — singleton behavior, name/SMILES lookups, JSON integrity, two-tier cascade (tier-2 fallback, tier-1 wins, graceful degradation)
- `test_superatom_table.py` (20 tests) — table building, case-insensitivity, L/R forms, MW lookup, label extraction, SMILES validity, ChemScanner-exclusive entries
- `test_reaction_parser.py` (60 tests) — SpeciesDescriptor/ReactionDescriptor roundtrip, split_condition_text (14 tests), _is_condition_token (9 tests), text label resolution, display name precedence, SM/DP identification, deduplication, build_reaction_smiles, mass computation (with salt splitting + adducts)
- `test_smoke.py` (14 tests) — CLI smoke tests for `lcms_analyzer`, `multi_lcms_analyzer`, `reaction_cleanup`, `coord_normalizer`, `rdf_parser`, `cdxml_builder`, `reagent_db`
- `test_smoke_extended.py` (42 tests) — CLI smoke tests for `discover_experiment_files`, `cas_resolver`, `scheme_aligner`, `prepare_reaction_scheme`, `ole_extractor`, `procedure_writer` (6 tests), `eln_enrichment`, `reactant_heuristic` (5 tests), `scheme_merger` (5 tests), `reaction_parser` (8 tests). Tests marked `@pytest.mark.network` require PubChem API

**ChemDraw COM tests (skipped in CI):**
- `test_builder.py` — integration tests for `cdxml_builder` with ChemDraw rendering
- COM-dependent tools (`scheme_polisher`, `cdx_converter`, `cdxml_to_image`, `eln_cdx_cleanup`, `ole_embedder`, `chemscript_bridge`) are not covered by smoke tests.

## Deployment verification
```bash
# Quick check: what works on this machine?
python deploy_check.py

# Full check with config validation:
python deploy_check.py --config deploy/config.yaml

# With PDF parsing test:
python deploy_check.py --config deploy/config.yaml --test-pdf path/to/lcms.pdf
```

Reports `[PASS]`/`[CRITICAL]`/`[WARNING]`/`[INFO]` for 14 check categories (Python version, packages, imports, ChemDraw COM, ChemScript, Java/OPSIN, ML envs, config parse, write permissions, etc.). Exit code 1 if any CRITICAL, 0 otherwise.

## Pipeline integration test
A self-contained test lives **outside** the project tree at `C:\Users\mic23\pipeline-test\`:
```
pipeline-test/
├── run_test.bat      # Runs pipeline on test data (double-click or CLI)
├── config.yaml       # Test config (output → pipeline-test\output\)
├── clean.bat         # Delete output dir for a fresh re-run
└── output/           # Created by run_test.bat (gitignored in project)
```

**Test data** (in `C:\Users\mic23\chem-test-data\`):
- Export ZIP: `lab-books.zip` (19 experiments, KL-7001-001 through KL-7001-019)
- LCMS PDFs: `testing-LCMS-NMR\LCMS\` (KL-7001-004 has 11 files, KL-7001-019 has 3)
- NMR PDFs: `testing-LCMS-NMR\DATA\` (KL-7001-019 has NMR data)

**Usage:**
```batch
REM Full test (KL-7001-004 + KL-7001-019):
C:\Users\mic23\pipeline-test\run_test.bat

REM Single experiment:
run_test.bat --experiments KL-7001-004

REM Scheme polishing only:
run_test.bat --step scheme

REM Clean up and re-run:
C:\Users\mic23\pipeline-test\clean.bat
C:\Users\mic23\pipeline-test\run_test.bat
```

The test uses `--export-zip` (no Findmolecule connection needed) and writes all output to `pipeline-test\output\`, keeping the project directory clean.

## ChemDraw COM Warning
Tools using ChemDraw COM automation (`eln_cdx_cleanup.py`, `cdx_converter.py`, `cdxml_to_image.py`, `scheme_polisher.py`, `ole_embedder.py`) launch and control ChemDraw programmatically. **Close ChemDraw before running these tools.** If ChemDraw is already open, unpredictable behavior may occur.
