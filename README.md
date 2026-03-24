# cdxml-toolkit

Chemistry office automation toolkit with MCP (Model Context Protocol) server. Lets LLM agents draw reaction schemes, parse ELN exports, analyze LCMS data, and produce publication-ready ChemDraw (CDXML) output.

The goal: any chemist with a consumer GPU can run a local LLM agent that helps with routine chemistry office tasks. The toolkit provides 15 grounded, validated chemistry tools that LLMs call via MCP — the agent reasons about chemistry while the tools handle SMILES resolution, 2D coordinate generation, and CDXML layout.

> Built and tested with Claude Code (Opus 4.6). I directed the design and architecture; Claude did the implementation. I'm a PhD organic chemist, not a programmer — this project wouldn't exist without Claude Code, and I thank Anthropic. 

## Quick start: MCP server

The primary interface is the MCP server. Connect it to any MCP-compatible agent (Claude Desktop, opencode, qwen-agent, etc.) and just chat naturally: "Draw deucravacitinib", "Help me complete my lab book", "Extract structures from this image".

### Claude Desktop

Edit `%APPDATA%\Claude\claude_desktop_config.json` (Windows) or `~/Library/Application Support/Claude/claude_desktop_config.json` (Mac):

```json
{
  "mcpServers": {
    "cdxml-toolkit": {
      "command": "python",
      "args": ["-m", "cdxml_toolkit.mcp_server"]
    }
  }
}
```

### opencode (for OpenRouter / local models)

Create `opencode.json`:

```json
{
  "provider": {
    "openrouter": {
      "models": { "qwen/qwen3.5-27b": {} }
    }
  },
  "mcp": {
    "cdxml-toolkit": {
      "type": "local",
      "command": ["python", "-m", "cdxml_toolkit.mcp_server"],
      "enabled": true,
      "timeout": 120000
    }
  }
}
```

### Verify it works

```
> Use cdxml-toolkit. Resolve "aspirin", then draw it.
```

Expected: 2 tool calls (resolve_name, draw_molecule), produces an aspirin CDXML file.

## MCP tools (15)

### Chemistry resolution
| Tool | Description |
|------|-------------|
| `resolve_name` | Name/abbreviation/CAS/formula to rich molecule JSON (4-tier: reagent DB, condensed formula, ChemScript, PubChem) |
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

## Installation

**Prerequisites:** Windows with ChemDraw (ChemOffice 2015+) and ChemScript installed.

```bash
# From GitHub — recommended install (includes RDKit, MCP server, ChemDraw COM,
# Office support, PDF parsing, image processing, and ChemScript bridge)
pip install "cdxml-toolkit[all] @ git+https://github.com/leehiufung911/cdxml-toolkit.git@main"

# With DECIMER neural image extraction (extract_structures_from_image)
pip install "cdxml-toolkit[all,decimer] @ git+https://github.com/leehiufung911/cdxml-toolkit.git@main"

# Development (editable install)
git clone https://github.com/leehiufung911/cdxml-toolkit.git
cd cdxml-toolkit
pip install -e ".[dev]"
```

**Required:** `lxml>=4.6`. **Recommended:** `rdkit>=2023.03` (needed for scheme rendering).

### Extras

| Extra | What it includes | Notes |
|-------|-----------------|-------|
| `[all]` | RDKit, pywin32, image, Office, YAML, PDF analysis, MCP server, pythonnet | **Use this.** Everything most users need. |
| `[decimer]` | TensorFlow, DECIMER, PyMuPDF | Neural image-to-SMILES. Adds ~1 GB. |
| `[full]` | `[all]` + `[decimer]` + `[ocr]` | Everything pip-installable. |
| `[rdkit]` | RDKit only | Minimal install for scripting. |
| `[mcp]` | MCP server + PyYAML | MCP server only (no RDKit/Office). |
| `[dev]` | `[all]` + pytest | For running the test suite. |

### System dependencies (not pip-installable)

| Dependency | Required for | Setup |
|-----------|-------------|-------|
| **ChemDraw** (ChemOffice 2015+) | CDX conversion, PNG rendering | Assumed installed. Must be **closed** before running. |
| **ChemScript .NET** | Name resolution (Tier 3: IUPAC to structure) | Comes with ChemOffice. Needs a 32-bit Python env (`chemscript32`) with pythonnet and a `~/.chemscript_config.json`. |
| **Microsoft Office** | OLE embedding into PPTX/DOCX | Optional. Only needed for `embed_cdxml_in_office`. |

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
pip install -e ".[dev]"
pytest tests/ -v
```

## License

[MIT](LICENSE)

## Attribution

See [NOTICE.md](NOTICE.md) for third-party data attribution (ChemScanner, RDKit).

## Author

Hiu Fung Kevin Lee ([@leehiufung911](https://github.com/leehiufung911))
