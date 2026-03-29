# cdxml-toolkit

Chemistry office automation toolkit with MCP (Model Context Protocol) server. Lets LLM agents draw reaction schemes, parse ELN exports, analyze LCMS data, and produce publication-ready ChemDraw (CDXML) output.

The goal: any chemist with a consumer GPU can run a local LLM agent that helps with routine chemistry office tasks. The toolkit provides 15 grounded, validated chemistry tools that LLMs call via MCP — the agent reasons about chemistry while the tools handle SMILES resolution, 2D coordinate generation, and CDXML layout.

> Built and tested with Claude Code (Opus 4.6). I directed the design and architecture; Claude did the implementation. I'm a PhD organic chemist, not a programmer — this project wouldn't exist without Claude Code, and I thank Anthropic. 

## Installation

**Prerequisites:** Windows with ChemDraw (ChemOffice 2015+) installed. Python 3.10–3.13 (3.14 is not yet supported by TensorFlow/DECIMER).

```bash
# 1. Create a conda environment and install
conda create -n cdxml python=3.12 pip -y
conda activate cdxml
pip install cdxml-toolkit

# 2. Run the doctor to check your setup
cdxml-doctor --no-tests
```

Everything is included by default: RDKit, MCP server, ChemDraw COM, Office support, PDF analysis, image processing, DECIMER neural image extraction, OPSIN, and OCR.

On first run, `cdxml-doctor` will extract the bundled JRE for OPSIN (~45 MB, one-time) and download DECIMER neural models (~570 MB). Subsequent runs are fast.

If ChemScript is not configured, `cdxml-doctor` will detect your ChemDraw installation, show what it found, and offer to set everything up automatically:

```
=== ChemScript setup ===

  Found ChemScript DLLs:
    Managed:  C:\...\CambridgeSoft.ChemScript16.dll (32-bit)
    Native:   C:\...\ChemScript160.dll (32-bit)

  32-bit ChemScript requires a 32-bit Python environment.
  The doctor will run the following commands:

    set CONDA_SUBDIR=win-32 && conda create -n chemscript32 python=3.10 pip -y
    C:\Users\YOU\miniconda3\envs\chemscript32\python.exe -m pip install pythonnet

  Proceed? [y/N] y

  Creating chemscript32 conda env...
  chemscript32 env created.
  Installing pythonnet in chemscript32...
  pythonnet installed.
  Saving config...
  ChemScript configured. Run cdxml-doctor again to verify.
```

Run `cdxml-doctor --no-tests` again to confirm ChemScript shows OK.

ChemScript is optional — without it, OPSIN handles IUPAC name resolution as an offline fallback. ChemScript adds bidirectional name-to-structure conversion and aligned naming.

Alternatively, install from GitHub for the latest development version:

```bash
pip install "cdxml-toolkit @ git+https://github.com/leehiufung911/cdxml-toolkit.git@main"
```

## MCP server

The primary interface is the MCP server. Connect it to any MCP-compatible agent (Claude Desktop, Claude Code, opencode, qwen-agent, etc.) and chat naturally: "Draw deucravacitinib", "Help me complete my lab book", "Extract structures from this image".

Edit your MCP config to point to the Python in your conda environment (replace `YOUR_USERNAME`):

```json
{
  "mcpServers": {
    "cdxml-toolkit": {
      "command": "C:\\Users\\YOUR_USERNAME\\miniconda3\\envs\\cdxml\\python.exe",
      "args": ["-m", "cdxml_toolkit.mcp_server"]
    }
  }
}
```

For Claude Desktop, this file is at `%APPDATA%\Claude\claude_desktop_config.json`.

### Verify it works

```
> Resolve "aspirin", then draw it.
```

Expected: 2 tool calls (resolve_name, draw_molecule), produces an aspirin CDXML file.

## MCP tools (15)

### Chemistry resolution
| Tool | Description |
|------|-------------|
| `resolve_name` | Name/abbreviation/CAS/formula to rich molecule JSON (5-tier: reagent DB, condensed formula, ChemScript, OPSIN, PubChem) |
| `modify_molecule` | 6 operations: analyze, name_surgery, smarts, set_smiles, set_name, reaction. 162 named reaction templates. Returns MCS-based structural diffs. |

### Structure rendering
| Tool | Description |
|------|-------------|
| `draw_molecule` | Single molecule to CDXML |
| `render_scheme` | YAML/compact text/reaction JSON to publication-ready CDXML. Forgiving parser handles common LLM YAML mistakes. |

### Perception (reading existing chemistry)
| Tool | Description |
|------|-------------|
| `parse_reaction` | ELN exports (CDXML/CDX/CSV/RXN) to semantic JSON with species, roles, SMILES, equivalents |
| `summarize_reaction` | Context-efficient view of reaction JSON (select only the fields you need) |
| `extract_structures_from_image` | Image to SMILES + confidence scores via DECIMER neural network |
| `parse_scheme` | CDXML scheme to structured species/steps/topology JSON |

### Analysis
| Tool | Description |
|------|-------------|
| `parse_analysis_file` | LCMS (Waters/manual) or NMR (MestReNova) PDF to structured peak data |
| `format_lab_entry` | Structured entry dicts to formatted lab book text. Re-reads LCMS PDFs for exact numbers. |

### Office integration
| Tool | Description |
|------|-------------|
| `extract_cdxml_from_office` | Pull embedded ChemDraw OLE objects from PPTX/DOCX |
| `embed_cdxml_in_office` | Inject CDXML as editable ChemDraw OLE into PPTX/DOCX |
| `convert_cdx_cdxml` | Bidirectional CDX/CDXML conversion |
| `search_compound` | Find a molecule across experiment directories by SMILES similarity |
| `render_to_png` | CDXML to PNG via ChemDraw COM |

## Design principles

**Never trust LLM-generated SMILES.** The agent always goes through `resolve_name` to get grounded SMILES from databases. Direct SMILES generation is the #1 source of chemistry hallucination.

**Verify every transformation.** `modify_molecule` returns aligned IUPAC name diffs and MCS-based molecular diffs after every edit. The agent can confirm the transformation is correct.

**Never flood the agent.** Large outputs (CDXML, JSON) always write to files and return `{ok: true, output_path: "...", size: 23456}`. The agent never gets 30KB of XML in its context window.

**Forgiving inputs.** The YAML parser accepts 9+ common LLM mistakes (inline structures, `substrates` as alias for `structures`, text as string not list, bare SMILES, `above_arrow` as list/string). Input parameters accept bare SMILES strings, stringified JSON arrays, and fuzzy operation names.

**Actionable errors.** Every error tells the agent what to do instead: "Did you mean: BOC_deprotection?", not "KeyError".

**Progressive discovery.** Call any tool with no arguments to get usage examples and schema reference.

## CLI tools

All tools are also available as command-line scripts:

| Command | Description |
|---------|-------------|
| `cdxml-mcp` | MCP server (primary interface) |
| `cdxml-parse` | Parse reaction files to JSON |
| `cdxml-render` | Render JSON/YAML/compact text to CDXML |
| `cdxml-convert` | CDX/CDXML bidirectional conversion |
| `cdxml-image` | CDXML to PNG/SVG (ChemDraw COM) |
| `cdxml-merge` | Merge multiple reaction schemes |
| `cdxml-layout` | Clean up reaction layout (pure Python) |
| `cdxml-ole` | Embed CDXML as editable OLE in PPTX/DOCX |
| `cdxml-lcms` | Parse LCMS PDF reports |
| `cdxml-nmr` | Extract NMR data from MestReNova PDFs |
| `cdxml-format-entry` | Format lab book entries |
| `cdxml-discover` | Discover experiment files in a directory |
| `cdxml-doctor` | Diagnostics, test runner, and ChemScript setup guide |

## Scheme DSL

The renderer accepts three input formats:

**YAML** (what agents typically write):
```yaml
layout: sequential
structures:
  SM:
    smiles: "Brc1ncnc2sccc12"
  Product:
    smiles: "c1nc(N2CCOCC2)c2ccsc2n1"
steps:
  - substrates: [SM]
    products: [Product]
    above_arrow:
      structures: [Morph]
    below_arrow:
      text: ["Pd2(dba)3", "BINAP", "Cs2CO3", "Dioxane, 105 C"]
```

**Compact text** ("Mermaid for reactions"):
```
SM: {Brc1ncnc2sccc12}
SM --> Product{c1nc(N2CCOCC2)c2ccsc2n1}
  above: Morph{C1COCCN1}
  below: "Pd2(dba)3", "BINAP", "Cs2CO3"
```

**Reaction JSON** (from parse_reaction):
```bash
cdxml-render --from-json reaction.json -o scheme.cdxml
```

## Running tests

```bash
# Using cdxml-doctor (recommended — also prints diagnostics)
cdxml-doctor

# Or directly with pytest
pytest tests/ -v
```

## License

[MIT](LICENSE)

## Attribution

See [NOTICE.md](NOTICE.md) for third-party data attribution (ChemScanner, RDKit).

## Author

Hiu Fung Kevin Lee ([@leehiufung911](https://github.com/leehiufung911))
