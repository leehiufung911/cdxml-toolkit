"""
parser.py — Parse YAML scheme files into SchemeDescriptor dataclasses.

Validates structure references, layout keywords, and produces clear error
messages for malformed input.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml

from .schema import (
    ArrowContent,
    RunArrowEntry,
    SchemeDescriptor,
    SectionDescriptor,
    StepDescriptor,
    StepRunArrows,
    StructureRef,
)

# Valid layout and wrap keywords
VALID_LAYOUTS = {
    "linear", "sequential", "divergent", "stacked-rows",
    "numbered-parallel", "convergent",
}
VALID_WRAPS = {"repeat", "serpentine", "none"}
VALID_ARROW_STYLES = {"solid", "dashed", "failed"}


class SchemeParseError(Exception):
    """Raised when YAML content is invalid or violates schema rules."""
    pass


def parse_yaml(source: Union[str, Path]) -> SchemeDescriptor:
    """
    Parse a YAML file or string into a SchemeDescriptor.

    Parameters
    ----------
    source : str or Path
        If a Path or string ending in .yaml/.yml, read as file.
        Otherwise, interpret as raw YAML text.

    Returns
    -------
    SchemeDescriptor

    Raises
    ------
    SchemeParseError
        On any validation failure.
    """
    text = _load_yaml_text(source)
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as e:
        raise SchemeParseError(f"Invalid YAML syntax: {e}") from e

    if not isinstance(data, dict):
        raise SchemeParseError("Top-level YAML must be a mapping (dict)")

    # Allow either top-level keys directly or wrapped in 'scheme:'
    if "scheme" in data and isinstance(data["scheme"], dict):
        data = data["scheme"]

    return _parse_scheme(data)


def _load_yaml_text(source: Union[str, Path]) -> str:
    """Load YAML text from a file path or return raw string."""
    if isinstance(source, Path):
        return source.read_text(encoding="utf-8")
    if isinstance(source, str) and (
        source.endswith(".yaml") or source.endswith(".yml")
    ):
        return Path(source).read_text(encoding="utf-8")
    return source


def _parse_scheme(data: Dict[str, Any]) -> SchemeDescriptor:
    """Parse the scheme-level dict."""
    # --- Source (reaction_parser JSON) ---
    source = data.get("source")
    if source is not None:
        source = str(source)

    # --- Structures ---
    raw_structs = data.get("structures", {})
    if raw_structs is None:
        raw_structs = {}
    if not isinstance(raw_structs, dict):
        raise SchemeParseError("'structures' must be a mapping")
    structures = {}
    for key, val in raw_structs.items():
        key = str(key)
        structures[key] = _parse_structure(key, val)

    # --- Sections (for stacked-rows layout) ---
    raw_sections = data.get("sections", [])
    if not isinstance(raw_sections, list):
        raise SchemeParseError("'sections' must be a list")
    sections = [_parse_section(i, s) for i, s in enumerate(raw_sections)]

    # --- Steps ---
    raw_steps = data.get("steps", [])
    if not isinstance(raw_steps, list):
        raise SchemeParseError("'steps' must be a list")
    if not raw_steps and not sections:
        raise SchemeParseError("At least one step is required (or use 'sections' for stacked-rows)")
    steps = [_parse_step(i, s) for i, s in enumerate(raw_steps)]

    # --- Validate structure refs in steps and sections ---
    # When source is present, refs may be resolved from JSON at render time,
    # so we only validate refs that are NOT in the declared structures block
    # (the renderer will handle resolution failures for source-backed refs).
    def _validate_step_refs(step_list, context_prefix=""):
        for i, step in enumerate(step_list):
            prefix = f"{context_prefix}step {i+1}"
            _validate_refs(f"{prefix} substrates", step.substrates, structures)
            _validate_refs(f"{prefix} products", step.products, structures)
            if step.above_arrow:
                _validate_refs(
                    f"{prefix} above_arrow.structures",
                    step.above_arrow.structures,
                    structures,
                )
            if step.below_arrow:
                _validate_refs(
                    f"{prefix} below_arrow.structures",
                    step.below_arrow.structures,
                    structures,
                )

    if not source:
        _validate_step_refs(steps)
        for sec_idx, sec in enumerate(sections):
            _validate_step_refs(sec.steps, f"section {sec_idx+1} ")

    # --- Layout ---
    layout = str(data.get("layout", "linear"))
    if layout not in VALID_LAYOUTS:
        raise SchemeParseError(
            f"Invalid layout '{layout}'. Must be one of: {sorted(VALID_LAYOUTS)}"
        )

    # --- Wrap ---
    wrap = str(data.get("wrap", "repeat"))
    if wrap not in VALID_WRAPS:
        raise SchemeParseError(
            f"Invalid wrap '{wrap}'. Must be one of: {sorted(VALID_WRAPS)}"
        )

    # --- Steps per row ---
    steps_per_row = data.get("steps_per_row")
    if steps_per_row is not None:
        steps_per_row = int(steps_per_row)
        if steps_per_row < 1:
            raise SchemeParseError("steps_per_row must be >= 1")

    # --- Title ---
    title = data.get("title")
    if title is not None:
        title = str(title)

    # --- Run arrows ---
    raw_runs = data.get("run_arrows", [])
    if not isinstance(raw_runs, list):
        raise SchemeParseError("'run_arrows' must be a list")
    run_arrows = [_parse_run_arrows(r) for r in raw_runs]

    # --- Condition key ---
    condition_key = data.get("condition_key")
    if condition_key is not None and not isinstance(condition_key, dict):
        raise SchemeParseError("'condition_key' must be a mapping")

    return SchemeDescriptor(
        source=source,
        structures=structures,
        steps=steps,
        layout=layout,
        wrap=wrap,
        steps_per_row=steps_per_row,
        title=title,
        run_arrows=run_arrows,
        condition_key=condition_key,
        sections=sections,
    )


def _parse_section(index: int, data: Any) -> SectionDescriptor:
    """Parse a single section entry (for stacked-rows layout)."""
    if not isinstance(data, dict):
        raise SchemeParseError(
            f"Section {index+1} must be a mapping, got {type(data).__name__}"
        )
    label = data.get("label")
    if label is not None:
        label = str(label)

    raw_steps = data.get("steps", [])
    if not isinstance(raw_steps, list):
        raise SchemeParseError(f"Section {index+1} 'steps' must be a list")
    if not raw_steps:
        raise SchemeParseError(f"Section {index+1} must have at least one step")
    steps = [_parse_step(i, s) for i, s in enumerate(raw_steps)]

    layout = str(data.get("layout", "linear"))

    return SectionDescriptor(label=label, steps=steps, layout=layout)


def _parse_structure(key: str, val: Any) -> StructureRef:
    """Parse a single structure entry."""
    if isinstance(val, str):
        # Shorthand: just a SMILES string
        return StructureRef(id=key, smiles=val)
    if not isinstance(val, dict):
        raise SchemeParseError(
            f"Structure '{key}' must be a mapping or SMILES string, got {type(val).__name__}"
        )
    return StructureRef(
        id=key,
        smiles=val.get("smiles"),
        name=val.get("name"),
        file=val.get("file"),
        cdxml_id=val.get("cdxml_id"),
        label=val.get("label"),
    )


def _parse_step(index: int, data: Any) -> StepDescriptor:
    """Parse a single step entry."""
    if not isinstance(data, dict):
        raise SchemeParseError(f"Step {index+1} must be a mapping, got {type(data).__name__}")

    substrates = _as_str_list(data.get("substrates", []), f"step {index+1} substrates")
    products = _as_str_list(data.get("products", []), f"step {index+1} products")

    if not substrates:
        raise SchemeParseError(f"Step {index+1} must have at least one substrate")
    if not products:
        raise SchemeParseError(f"Step {index+1} must have at least one product")

    above = _parse_arrow_content(data.get("above_arrow"), f"step {index+1} above_arrow")
    below = _parse_arrow_content(data.get("below_arrow"), f"step {index+1} below_arrow")

    return StepDescriptor(
        substrates=substrates,
        products=products,
        above_arrow=above,
        below_arrow=below,
        yield_=data.get("yield"),
        number=data.get("number"),
        id=data.get("id"),
        arrow_style=_validate_arrow_style(data.get("arrow_style", "solid"), index),
    )


def _parse_arrow_content(data: Any, context: str) -> Optional[ArrowContent]:
    """Parse above_arrow or below_arrow content."""
    if data is None:
        return None
    if not isinstance(data, dict):
        raise SchemeParseError(f"'{context}' must be a mapping, got {type(data).__name__}")
    structures = _as_str_list(data.get("structures", []), f"{context}.structures")
    text = _as_str_list(data.get("text", []), f"{context}.text")
    if not structures and not text:
        return None
    return ArrowContent(structures=structures, text=text)


def _parse_run_arrows(data: Any) -> StepRunArrows:
    """Parse a run_arrows entry."""
    if not isinstance(data, dict):
        raise SchemeParseError(f"run_arrows entry must be a mapping")
    step = data.get("step")
    if step is None:
        raise SchemeParseError("run_arrows entry must have a 'step' field")
    step = int(step)
    raw_runs = data.get("runs", [])
    if not isinstance(raw_runs, list):
        raise SchemeParseError(f"run_arrows step {step} 'runs' must be a list")
    runs = []
    for r in raw_runs:
        if not isinstance(r, dict):
            raise SchemeParseError(f"run entry must be a mapping")
        inp = r.get("input", "")
        out = r.get("output", "")
        if not inp:
            raise SchemeParseError(
                f"run entry must have an 'input' field"
            )
        note = r.get("note")
        if note is not None:
            note = str(note)
        runs.append(RunArrowEntry(input_label=str(inp), output_label=str(out),
                                  note=note))
    return StepRunArrows(step=step, runs=runs)


def _validate_arrow_style(style: str, step_idx: int) -> str:
    """Validate arrow_style value."""
    style = str(style)
    if style not in VALID_ARROW_STYLES:
        raise SchemeParseError(
            f"Step {step_idx+1}: invalid arrow_style '{style}'. "
            f"Must be one of: {sorted(VALID_ARROW_STYLES)}"
        )
    return style


def _validate_refs(
    context: str,
    refs: List[str],
    structures: Dict[str, StructureRef],
) -> None:
    """Check that all refs point to defined structures.

    Undeclared refs are auto-registered as bare StructureRefs so the renderer
    can attempt resolution via reagent_db (name → SMILES lookup).
    """
    for ref in refs:
        if ref not in structures:
            # Auto-register as a bare ref — renderer will try reagent_db
            structures[ref] = StructureRef(id=ref)


def _as_str_list(val: Any, context: str) -> List[str]:
    """Coerce a value to a list of strings."""
    if val is None:
        return []
    if isinstance(val, str):
        return [val]
    if isinstance(val, list):
        return [str(v) for v in val]
    raise SchemeParseError(f"'{context}' must be a list or string, got {type(val).__name__}")
