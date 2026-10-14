"""Read the ASCII TetGen PLC subset used by this project.

Point IDs need not be contiguous. Polygon facets are triangulated as fans for
display only; the original files are never rewritten by this reader.
"""

from dataclasses import dataclass
from pathlib import Path
import gzip

import numpy as np


@dataclass
class Surface:
    points: np.ndarray
    triangles: np.ndarray
    markers: np.ndarray
    point_markers: np.ndarray
    facet_count: int
    regions: list[list[float]]


def read_poly(path: str | Path) -> Surface:
    """Read inline 3-D points and facets, rejecting unsupported facet holes."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as stream:
        rows = (line.split("#", 1)[0].split() for line in stream)
        rows = (row for row in rows if row)

        def take():
            try:
                return next(rows)
            except StopIteration as exc:
                raise ValueError("Unexpected end of .poly file") from exc

        count, dimension, attributes, marked = map(int, take())
        if dimension != 3 or count == 0:
            raise ValueError(
                "Expected inline 3-D points; external .node files are unsupported"
            )
        points = np.empty((count, 3))
        point_markers = np.zeros(count, dtype=np.int32)
        lookup = {}
        for index in range(count):
            row = take()
            point_id = int(row[0])
            if point_id in lookup:
                raise ValueError(f"Duplicate point ID: {point_id}")
            lookup[point_id] = index
            points[index] = row[1:4]
            if marked:
                point_markers[index] = int(row[4 + attributes])
        facet_count, marked = map(int, take())
        triangles, markers = [], []
        for _ in range(facet_count):
            header = take()
            polygons, holes = map(int, header[:2])
            marker = int(header[2]) if marked else 0
            if holes:
                raise ValueError(
                    "Facet holes are unsupported by the preview triangulator"
                )
            for _ in range(polygons):
                row = take()
                size = int(row[0])
                if size < 3 or len(row) != size + 1:
                    raise ValueError("Expected a polygon with at least three point IDs")
                try:
                    ids = [lookup[int(value)] for value in row[1:]]
                except KeyError as exc:
                    raise ValueError(
                        f"Facet references missing point {exc.args[0]}"
                    ) from exc
                for corner in range(1, size - 1):
                    triangles.append((ids[0], ids[corner], ids[corner + 1]))
                    markers.append(marker)
        for _ in range(int(take()[0])):
            take()
        region_header = next(rows, ["0"])
        regions = [list(map(float, take())) for _ in range(int(region_header[0]))]
    if not np.isfinite(points).all():
        raise ValueError("Non-finite point coordinates")
    return Surface(
        points,
        np.asarray(triangles, dtype=np.int32).reshape(-1, 3),
        np.asarray(markers, dtype=np.int32),
        point_markers,
        facet_count,
        regions,
    )
