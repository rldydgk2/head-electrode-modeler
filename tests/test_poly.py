"""Format-level checks for previews and point-ID handling."""

import gzip
from pathlib import Path
import tempfile
import unittest

from head_model.poly import read_poly

POLY = """# Non-contiguous IDs, an attribute and boundary markers
4 3 1 1
10 0 0 0 0.5 1005
20 1 0 0 0.5 1005
30 1 1 0 0.5 1005
40 0 1 0 0.5 1005
1 1
1 0 1005
4 10 20 30 40
0
1
1 0.5 0.5 -1 7
"""


class PolyTests(unittest.TestCase):
    def test_polygon_ids_attributes_and_gzip(self):
        with tempfile.TemporaryDirectory() as directory:
            for filename in ("sample.poly", "sample.poly.gz"):
                path = Path(directory) / filename
                opener = gzip.open if filename.endswith(".gz") else open
                with opener(path, "wt") as stream:
                    stream.write(POLY)
                surface = read_poly(path)
                self.assertEqual(surface.triangles.tolist(), [[0, 1, 2], [0, 2, 3]])
                self.assertEqual(surface.markers.tolist(), [1005, 1005])
                self.assertEqual(surface.point_markers.tolist(), [1005] * 4)
                self.assertEqual(surface.regions, [[1, 0.5, 0.5, -1, 7]])

    def test_missing_reference_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.poly"
            path.write_text(POLY.replace("4 10 20 30 40", "4 10 20 30 99"))
            with self.assertRaisesRegex(ValueError, "missing point 99"):
                read_poly(path)

    def test_truncated_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.poly"
            path.write_text("4 3 0 0\n")
            with self.assertRaisesRegex(ValueError, "Unexpected end"):
                read_poly(path)


if __name__ == "__main__":
    unittest.main()
