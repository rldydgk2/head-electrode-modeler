"""Command-line entry points for inspection, previews, and the research editor."""

import argparse
import json
from pathlib import Path

import numpy as np

from .poly import read_poly


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser("inspect", help="Summarize a TetGen .poly surface")
    inspect.add_argument("input", type=Path)
    render = commands.add_parser("render", help="Save a two-view PNG preview")
    render.add_argument("input", type=Path)
    render.add_argument("--output", type=Path, default=Path("outputs/preview.png"))
    render.add_argument("--title", default="Head model")
    edit = commands.add_parser("edit", help="Launch the experimental research editor")
    edit.add_argument("input", type=Path)
    edit.add_argument("--output-dir", type=Path, default=Path("outputs"))
    edit.add_argument(
        "--single", action="store_true", help="Place one electrode instead of arrays"
    )
    args = parser.parse_args()
    if args.command == "inspect":
        surface = read_poly(args.input)
        markers, counts = np.unique(surface.markers, return_counts=True)
        print(
            json.dumps(
                {
                    "points": len(surface.points),
                    "facets": surface.facet_count,
                    "preview_triangles": len(surface.triangles),
                    "triangles_by_marker": dict(
                        zip(map(str, markers), map(int, counts))
                    ),
                    "regions": len(surface.regions),
                    "bounds_mm": [
                        surface.points.min(axis=0).tolist(),
                        surface.points.max(axis=0).tolist(),
                    ],
                },
                indent=2,
            )
        )
    elif args.command == "render":
        from .render import render as render_surface

        render_surface(args.input, args.output, title=args.title)
        print(args.output)
    else:
        from .editor import run

        run(args.input, args.output_dir, single=args.single)


if __name__ == "__main__":
    main()
