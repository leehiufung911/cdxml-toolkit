#!/usr/bin/env python3
"""
render_scheme.py — CLI entry point for scheme → CDXML rendering.

Accepts both YAML and compact syntax input files.

Usage:
    python experiments/scheme_dsl/render_scheme.py examples/simple_linear.yaml -o output.cdxml
    python experiments/scheme_dsl/render_scheme.py examples/buchwald.yaml
    python experiments/scheme_dsl/render_scheme.py scheme.txt --format compact
"""

import argparse
import os
import sys

from .parser import SchemeParseError, parse_yaml
from .compact_parser import parse_compact_file, ParseError
from .scheme_yaml_writer import write_scheme_yaml
from .renderer import render, render_to_file

# File extensions that trigger compact syntax parsing
_COMPACT_EXTENSIONS = {".rxn", ".scheme", ".txt"}
_YAML_EXTENSIONS = {".yaml", ".yml"}


def _detect_format(input_path: str) -> str:
    """Detect input format from file extension or content.

    Returns "yaml" or "compact".
    """
    ext = os.path.splitext(input_path)[1].lower()
    if ext in _YAML_EXTENSIONS:
        return "yaml"
    if ext in _COMPACT_EXTENSIONS:
        return "compact"

    # Fallback: sniff content — if first non-comment line contains -->, it's compact
    try:
        with open(input_path, "r", encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                if "-->" in stripped or "..>" in stripped or "X>" in stripped or "x>" in stripped:
                    return "compact"
                break
    except OSError:
        pass

    return "yaml"


def main():
    parser = argparse.ArgumentParser(
        description="Render a reaction scheme (YAML, compact syntax, or JSON) to CDXML.",
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="Input file (YAML or compact syntax)",
    )
    parser.add_argument(
        "--from-json",
        default=None,
        metavar="JSON",
        help="Render directly from a reaction_parser JSON file",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output CDXML file (default: input stem + .cdxml)",
    )
    parser.add_argument(
        "--format",
        choices=["yaml", "compact"],
        default=None,
        help="Input format (auto-detected from extension if omitted)",
    )
    parser.add_argument(
        "--layout",
        default="auto",
        help="Layout override for --from-json (default: auto)",
    )
    parser.add_argument(
        "--no-run-arrows",
        action="store_true",
        help="Suppress run arrows (SM mass → product yield)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print progress to stderr",
    )
    args = parser.parse_args()

    # --- Handle --from-json mode ---
    if args.from_json:
        json_path = args.from_json
        if not os.path.exists(json_path):
            print(f"Error: JSON file not found: {json_path}", file=sys.stderr)
            sys.exit(1)

        # Output CDXML path
        output_path = args.output
        if output_path is None:
            stem = os.path.splitext(os.path.basename(json_path))[0]
            output_path = os.path.join(
                os.path.dirname(json_path) or ".", f"{stem}-scheme.cdxml")

        # Step 1: JSON → YAML (layout decisions happen here)
        yaml_path = os.path.splitext(output_path)[0] + ".yaml"
        if args.verbose:
            print(f"Writing scheme YAML: {json_path} → {yaml_path}",
                  file=sys.stderr)

        try:
            write_scheme_yaml(
                json_path, yaml_path,
                layout=args.layout,
                include_run_arrows=not args.no_run_arrows,
            )
        except Exception as e:
            print(f"YAML writer error: {e}", file=sys.stderr)
            sys.exit(1)

        # Step 2: YAML → SchemeDescriptor → CDXML (pure rendering)
        if args.verbose:
            print(f"Rendering: {yaml_path} → {output_path}", file=sys.stderr)

        try:
            scheme = parse_yaml(yaml_path)
        except SchemeParseError as e:
            print(f"Parse error in generated YAML: {e}", file=sys.stderr)
            sys.exit(1)

        if args.verbose:
            print(
                f"  {len(scheme.structures)} structures, "
                f"{len(scheme.steps)} steps, "
                f"layout={scheme.layout}",
                file=sys.stderr,
            )

        try:
            yaml_dir = os.path.dirname(os.path.abspath(yaml_path))
            render_to_file(scheme, output_path, yaml_dir=yaml_dir)
        except Exception as e:
            print(f"Render error: {e}", file=sys.stderr)
            sys.exit(1)

        print(f"Written: {output_path}")
        return

    # --- Standard YAML/compact mode ---
    if not args.input:
        print("Error: input file required (or use --from-json)", file=sys.stderr)
        sys.exit(1)

    # Resolve input path
    input_path = args.input
    if not os.path.isabs(input_path):
        # Try relative to CWD first, then relative to script dir
        if not os.path.exists(input_path):
            alt = os.path.join(os.path.dirname(__file__), input_path)
            if os.path.exists(alt):
                input_path = alt

    if not os.path.exists(input_path):
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Resolve output path
    output_path = args.output
    if output_path is None:
        stem = os.path.splitext(os.path.basename(input_path))[0]
        output_path = os.path.join(os.path.dirname(input_path), f"{stem}.cdxml")

    # Detect format
    fmt = args.format or _detect_format(input_path)

    # Parse
    if args.verbose:
        print(f"Parsing {input_path} (format: {fmt})...", file=sys.stderr)
    try:
        if fmt == "compact":
            scheme = parse_compact_file(input_path)
        else:
            scheme = parse_yaml(input_path)
    except (SchemeParseError, ParseError) as e:
        print(f"Parse error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        source_info = f", source={scheme.source}" if scheme.source else ""
        print(
            f"  {len(scheme.structures)} structures, "
            f"{len(scheme.steps)} steps, "
            f"layout={scheme.layout}{source_info}",
            file=sys.stderr,
        )

    # Render
    if args.verbose:
        print(f"Rendering to {output_path}...", file=sys.stderr)
    try:
        yaml_dir = os.path.dirname(os.path.abspath(input_path))
        render_to_file(scheme, output_path, yaml_dir=yaml_dir)
    except Exception as e:
        print(f"Render error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Written: {output_path}")


if __name__ == "__main__":
    main()
