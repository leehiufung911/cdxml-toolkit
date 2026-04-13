# cdxml-toolkit: What It Does and What Problem Does It Solve

This document goes beyond the README to explain what problems cdxml-toolkit solves, the design decisions behind it, and what it looks like in practice across 60,000 lines of code.

---

## The Problem

Chemistry office work is tedious and repetitive. Every day, chemists:

- Draw and redraw reaction schemes in ChemDraw
- Parse ELN (Electronic Lab Notebook) exports into usable formats
- Cross-reference LCMS and NMR characterization data with their reactions
- Write lab book entries with exact mass spec peak numbers
- Build weekly report slides with publication-quality schemes
- Search old lab books: "Have I ever made this compound before?"

LLMs can reason about chemistry — they know what a Buchwald coupling is, what reagents to use, what products to expect. But they can't:

- **Generate valid SMILES reliably.** SMILES are a character-level molecular encoding. One wrong character = wrong molecule. LLMs hallucinate SMILES constantly.
- **Read chemical structures from images.** Vision models can tell you "this is a reaction scheme" but can't determine exact atom connectivity.
- **Produce publication-quality ChemDraw output.** CDXML is a complex XML format with 2D coordinates, bond angles, text formatting, and ACS style rules.

cdxml-toolkit bridges this gap. It provides grounded, validated chemistry tools that LLMs call via MCP (Model Context Protocol). The agent reasons; the tools handle the chemistry.

---

## The Origin Story

The vision started simple: could I ask an agent to draw stuff in ChemDraw for me?

> "Hey, draw 2-chloropyridine, with a CF2Me on the 3-position!"
> "Take this structure but replace the pyridine with a pyrazine!"

Initially this seemed to require decomposing the problem into API calls to ChemDraw or MarvinSketch. It seemed hard, so it was shelved.

What came next was a series of deterministic scripts: clean up schemes from ELN exports, parse species identities, cross-reference with LCMS files, help write lab book entries. During this process, I used Claude Code to help evaluate ~8 different reaction product/byproduct prediction tools — installing each one, troubleshooting, passing test reactions, reading outputs. This was the moment I realized the potential of agentic AI: it would be much better to make **composable tools** that an agent can chain together to adapt to any situation, than hard-coded deterministic scripts.

The key insight was that I needed two things:
1. A **DSL** (domain-specific language) to let the LLM describe the *shape* of a scheme in text, which a deterministic renderer then translates into an actual chemical scheme
2. A **scheme parser** that does the reverse — takes CDXML schemes and extracts where the arrows are, what the species are relative to the arrows, what their roles are

From there, the toolkit grew into 15 tools covering the full chemistry office workflow.

---

## Design Decisions That Matter

### 1. Never Trust LLM-Generated SMILES

This is the single most important design decision. Every molecule in cdxml-toolkit comes from a grounded source: a curated reagent database (186 entries), a condensed formula parser, ChemScript's IUPAC engine, OPSIN (bundled offline), or PubChem. The agent is **never allowed** to write a SMILES string from its built-in chemistry knowledge, from reading an image with vision, or by hand-editing a string.

Why? Because SMILES hallucination is the #1 source of chemistry errors in LLM outputs. `c1ccccc1` is benzene. `c1cccc1` is... nothing valid, but an LLM might generate it. `CC(=O)Oc1ccccc1C(=O)O` is aspirin. `CC(=O)Oc1ccccc1C(O)=O` is also aspirin. `CC(=O)Oc1cccc1C(=O)O` is not. The difference is one character.

### 2. Verify Every Transformation (Molecular Diffs)

When a chemist edits a molecule in ChemDraw, they click to add a bond, see the methyl group appear, and conclude the edit worked. Agents need the same feedback loop.

Every `modify_molecule` operation returns:
- **Aligned IUPAC name diffs**: `"5-aminopentyl" → "3-aminopropoxy"` — the agent can read exactly what changed in human-chemistry terms
- **MCS-based molecular diffs**: which atoms changed, what they changed from/to
- **Delta formula and MW**: `+C6H4, -CH3, delta_mw: +52.03`

This means the agent can confirm that the transformation it requested is actually what happened — not just trust blindly.

### 3. IUPAC Names as an LLM-Native Molecular Representation

SMILES are machine-friendly but LLM-hostile. `c1nc(N2CCOCC2)c2ccsc2n1` tells an LLM almost nothing. But "4-(morpholin-4-yl)thieno[3,2-d]pyrimidine" — now the LLM can reason. It knows there's a morpholine. It knows there's a thienopyrimidine.

This led to the idea of **aligned non-canonical IUPAC names**. When a scheme shows:

> 2-bromoquinoline → 2-morpholinoquinoline → 2-morpholino-3-bromoquinoline

...it's immediately clear that step 1 swaps bromine for morpholine, and step 2 is a bromination. This works even on larger molecules — at every step, usually only one thing changes, and the names make that change visible.

### 4. Name Surgery: Solving the Original Problem in Reverse

Building the IUPAC name decomposer accidentally solved the original vision — "draw me 2-chloropyridine with a CF2Me on the 3-position." Because if the agent understands it needs CF2Me in prefix form → "1,1-difluoroethyl" → then constructs "2-chloro-3-(1,1-difluoroethyl)pyridine" → we're done. IUPAC names are a **lossless molecular representation**. For substitutions (which is most of medicinal chemistry), simple "name surgery" just works.

Building the aligned naming function was non-trivial. IUPAC rules define one canonical name per molecule, and no software exists to enumerate *non-canonical* forms. (In fact, no open-source software exists to do structure-to-IUPAC-name at all.) The workaround: use ChemDraw's naming engine as an oracle. Parse the canonical IUPAC name's bracket hierarchy to identify substituent boundaries, surgically cleave the molecule at each boundary, attach each fragment to a "naming probe" (a long-chain acid) to force IUPAC to treat it as a substituent, and read off the -yl name ChemDraw produces. Non-canonical forms can then be enumerated by string manipulation.

### 5. Forgiving Inputs, Actionable Errors

LLMs make predictable mistakes when writing YAML. The scheme parser accepts 9+ common errors:
- `species` or `substrates` as alias for `structures`
- `reactants` as alias for `substrates` in steps
- Inline structure dicts in step arrays
- Bare SMILES strings where structure refs are expected
- `text` as a string instead of a list
- `reagents` key (auto-distributed to above/below arrow)
- `above_arrow` as a list, string, or dict

This forgiving parser reduced `render_scheme` retries from 17 (on Ministral 8B) to 1-3. Every error message tells the agent what to do: "Did you mean: BOC_deprotection?" instead of "KeyError: 'BOC deprotection'".

### 6. Never Flood the Agent's Context

CDXML files are 8-37 KB of XML. Reaction JSON can be 34 KB. Dumping that into the agent's context window is fatal — it wastes tokens and confuses the agent. All large outputs write to files and return `{ok: true, output_path: "...", size: 23456}`. The agent gets a handle, not a wall of XML. A separate `summarize_reaction` tool lets the agent request only the fields it needs from a parsed reaction.

---

## What It Looks Like in Practice

cdxml-toolkit was validated across 10 chemistry tasks of increasing complexity, tested with Sonnet 4.6 (Phase 2) and Haiku 4.5 (Phase 3). All results below are from real MCP tool calls — no manual intervention.

### Task 1: Draw a Single Molecule

**Prompt**: "Draw deucravacitinib"

The agent calls `resolve_name("deucravacitinib")` → gets SMILES with isotopic deuterium labels → `draw_molecule()` → 8.8 KB publication-ready CDXML. Two tool calls, ~27 seconds.

This is the simplest workflow, but it validates the entire resolution chain. Deucravacitinib contains CD3 groups (deuterated methyl), which the resolution pipeline handles correctly.

### Task 2: Molecular Modification with Verification

**Prompt**: "Replace the CD3 amide with a benzyl amide in deucravacitinib"

The agent resolves deucravacitinib, then uses `modify_molecule(operation="smarts")` with a SMARTS pattern that swaps the CD3 amide for benzyl. The tool returns a verification diff:

```
delta_formula: +C6H4
atoms_changed: [{position: "amide N-substituent", from: "CD3", to: "benzyl"}]
```

The agent reads this, confirms the swap is correct, and draws the modified structure. The diff is the agent's "eyes" — without it, the agent would have to trust blindly that the SMARTS worked.

### Task 3: Multi-Step Reaction Scheme

**Prompt**: "Show the hydrolysis of deucravacitinib's CD3 amide, then HATU coupling with 2-aminoaniline"

The agent:
1. Resolves deucravacitinib, HATU, DIPEA, 2-aminoaniline
2. Uses `modify_molecule(operation="smarts")` for targeted hydrolysis — must cleave only the CD3 amide while preserving the cyclopropane amide
3. Uses `modify_molecule(operation="reaction", reaction_name="amide_coupling")` for step 2
4. Writes YAML and calls `render_scheme()` → 24 KB sequential CDXML

This task tests selective chemistry reasoning: the molecule has two amide bonds, and only one should be hydrolyzed.

### Task 4: Image Extraction + Reaction Building

**Input**: A PNG image of a deprotection reaction

The agent:
1. Calls `extract_structures_from_image()` → DECIMER neural network extracts 3 structures with confidence scores (0.9957, 0.9959, 0.7544)
2. Applies `modify_molecule(operation="reaction", reaction_name="BOC_deprotection")` to the extracted starting material
3. Chains with an amide coupling
4. Renders the 2-step scheme

This demonstrates the full pipeline: LLM reasoning → ML perception (DECIMER) → grounded chemistry (reaction templates) → publication output (CDXML). The LLM never sees the image's molecules directly — DECIMER handles perception, and the SMILES are validated.

### Task 5: Complex 3-Step Synthesis from an Image

**Input**: Image of a Boc-protected amine

The agent builds a 3-step route: Boc deprotection → amide coupling with 2-bromonicotinic acid → Buchwald-Hartwig amination with morpholine. Seven tool calls, 27 KB CDXML.

The chemical subtlety: the Buchwald reaction uses an aryl bromide, so the aryl-Br bond from step 2 must survive through the coupling and be selectively activated in step 3. The agent reasons about this correctly through the tool's reaction template system.

### Task 6: Lab Book Entry Assembly

**Input**: ELN CDX file + 4 LCMS PDFs (premix, 0 min, 10 min, final purified)

The agent:
1. Converts CDX → CDXML, parses the reaction (species, roles, equivalents, masses)
2. Parses all 4 LCMS PDFs — extracts retention times, m/z values, UV areas
3. Assembles a structured lab book entry with exact peak numbers

The reaction is a Mitsunobu etherification at 50 mg scale. The tool correctly identifies the product ([M-H]- 444.1), the triphenylphosphine oxide byproduct ([M+H]+ 279.1), and a double-addition artifact ([M+H]+ 517.2). Final C18 purity: 100% by 220 nm UV.

### Task 7: Multi-Experiment Weekly Report

**Input**: 6 ELN experiments (KL-7001-009 through KL-7001-014)

The agent parses all 6 experiments, generates individual CDXML schemes, and groups them into two distinct routes by comparing product SMILES and molecular weights:

- **Route A** (5-carbon linker, MW 459→396): experiments 009, 011, 013
- **Route B** (4-carbon linker, MW 445→382): experiments 010, 012, 014

**The agent caught an ELN mislabeling error**: experiment KL-7001-012 had its starting material labeled "KL-7001-009" in the ELN, but the SMILES and MW showed it was a 4-carbon linker (Route B), not 5-carbon (Route A). This is the kind of error that can propagate through months of lab work if undetected.

### Task 8: Compound Search Across 19 Experiments

**Input**: Image of a target compound (L5-C1)

The agent extracts the structure via DECIMER, then searches across 19 lab book experiments:

- **2 exact matches**: KL-7001-011 (80.3 mg, 81% yield) and KL-7001-013 (1998 mg, 88% yield, scale-up)
- **3 similar matches**: L4-C1 homologues (one CH2 shorter, 97.96% Tanimoto similarity)
- **Total L5-C1 produced**: 2.08 g across both batches

### Task 9: Patent Route Extraction from PowerPoint

**Input**: PPTX file where slide 2 contains a PNG image of a Chinese patent route (CN118027001A)

The agent extracts slide images from the PPTX, runs DECIMER on the patent route image (5 structures detected, confidence 0.68-1.00), and reconstructs the 3-step synthesis as a 22 KB CDXML scheme:
1. Boc protection (Boc2O, NaHCO3)
2. Mitsunobu coupling (PPh3, DEAD)
3. Boc deprotection (TFA)

### Task 10: Modify an Embedded PPTX Scheme

**Input**: PPTX with 2 stacked reaction schemes containing ChemDraw OLE objects

**Task**: Shorten the linker by 2 carbons in both reactions (5C→3C top, 4C→2C bottom)

The agent:
1. Extracts embedded CDXML from the PPTX OLE objects
2. Decodes the fragment SMILES
3. Uses `modify_molecule(operation="set_smiles")` for each structure, with verification diffs confirming the carbon count change
4. Renders new CDXML with stacked-rows layout
5. Round-trips the modified CDXML back into the PPTX as OLE objects

This is the most ambitious task — a full PPTX → CDXML → modify → CDXML → CDX → OLE → PPTX round-trip. The OLE embedder is still barebones, so this workflow is more proof-of-concept than production-ready.

---

## Test Results

### Phase 2 (Sonnet 4.6) and Phase 3 (Haiku 4.5) — Full Suite

| Task | Description | Sonnet 4.6 | Haiku 4.5 |
|------|-------------|:----------:|:---------:|
| 1 | Draw molecule | Pass | Pass |
| 2 | Modify + verify | Pass | Pass |
| 3 | 2-step scheme | Pass | Pass |
| 4 | Image + coupling | Pass | Pass |
| 5 | Image + 3-step | Pass | Pass |
| 6 | Lab book entry | Pass | Pass |
| 7 | Merge 6 experiments | Pass | Pass |
| 8 | Compound search | Pass | Pass |
| 9 | Patent route from PPTX | Pass | Pass |
| 10 | Modify PPTX scheme | Pass | Pass |

**Phase 2**: 444k tokens, 118 output files (22 CDXML, 19 JSON).
**Phase 3**: 342k tokens (23% fewer), 82 output files.

Both phases: **10/10, zero manual intervention.**

For comparison, Phase 1 (Sonnet 4.6 with raw Python APIs, no MCP tools) also scored 10/10 but required 551k tokens and 4 code patches during execution. MCP tools made the same work **41% faster** and eliminated all patching.

---

## Sample Output

The scheme below was rendered entirely by cdxml-toolkit from a YAML description — 2D coordinates computed by RDKit, layout by the toolkit's rendering engine, formatted in ACS Document 1996 style. It opens directly in ChemDraw as an editable, publication-ready file.

The YAML that produces it:

```yaml
layout: sequential

structures:
  ArBr:
    smiles: "Brc1ncnc2sccc12"
  Morph:
    smiles: "C1COCCN1"
  KL_CC_001:
    smiles: "c1nc(N2CCOCC2)c2ccsc2n1"
    label: "KL-CC-001"
  TsOPyrrole:
    smiles: "Cc1ccc(S(=O)(=O)OCc2ccc[nH]2)cc1"
  KL_CC_002:
    smiles: "c1c[nH]c(Cc2cc3c(N4CCOCC4)ncnc3s2)c1"
    label: "KL-CC-002"

steps:
  - substrates: [ArBr]
    products: [KL_CC_001]
    above_arrow:
      structures: [Morph]
    below_arrow:
      text:
        - "Pd2(dba)3 (0.05 eq.)"
        - "rac-BINAP (0.1 eq.)"
        - "Cs2CO3 (2 eq.)"
        - "Dioxane, 105 C, 24 h"
  - substrates: [KL_CC_001]
    products: [KL_CC_002]
    above_arrow:
      structures: [TsOPyrrole]
    below_arrow:
      text:
        - "n-BuLi (1.05 eq.)"
        - "THF, -78 C, 1 h"

run_arrows:
  - step: 1
    runs:
      - input: "2.15 g"
        output: "1.60 g, 72% yield"
  - step: 2
    runs:
      - input: "1.25 g"
        output: "600 mg, 44% yield"
```

(See `samples/consolidated/two-step-scheme.png` for the rendered output.)

---

## Maturity

Not all parts of the toolkit are equally battle-tested:

- **Most mature**: Scheme drawing and rendering (`resolve_name`, `modify_molecule`, `draw_molecule`, `render_scheme`). These are the core workflow and have been tested extensively across multiple models and task types.
- **Well tested**: Reaction parsing (`parse_reaction`, `summarize_reaction`, `parse_scheme`) and the YAML forgiving parser.
- **Less tested**: Lab book entry assembly (`format_lab_entry`, `parse_analysis_file`) and compound search (`search_compound`). These work in the tested scenarios but haven't seen as much real-world use.
- **Proof of concept**: OLE embedding into PPTX/DOCX (`embed_cdxml_in_office`). Extraction works well; re-embedding is barebones.

---

## Who It's For

Any chemist who wants an AI assistant that can actually *do* chemistry office work, not just talk about it. The MCP interface means it works with Claude Desktop, Claude Code, opencode, or any MCP-compatible agent. The goal: a chemist with a consumer GPU runs a local LLM that handles the tedious parts of their day — drawing schemes, parsing data, writing lab books — while the tools guarantee chemical correctness.
