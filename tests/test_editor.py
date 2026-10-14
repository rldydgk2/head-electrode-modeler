"""Regression checks for the small corrections to the historical editor."""

import unittest
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import numpy as np

from head_model import editor


class EditorTests(unittest.TestCase):
    def test_spacing_remains_defined_after_four_placements(self):
        self.assertEqual(
            [editor.array_spacing(i) for i in range(6)],
            [(36, 25), (36, 25), (25, 36), (25, 36), (36, 25), (36, 25)],
        )

    def test_orientation_uses_facet_size_not_dictionary_size(self):
        skin = {
            1: np.array([0.0, 0.0, 1.0]),
            2: np.array([0.0, 1.0, 1.0]),
            3: np.array([1.0, 0.0, 1.0]),
        }
        facets = {9: [1005, 1, 2, 3]}
        skin_facets = {9: [1, 2, 3]}
        with patch.multiple(
            editor,
            create=True,
            skin=skin,
            facet=facets,
            skin_facet=skin_facets,
            mass_p_dic={9: np.array([1 / 3, 1 / 3, 1.0])},
        ):
            editor.orient_facet(9)
        self.assertEqual(facets[9], [1005, 3, 2, 1])


if __name__ == "__main__":
    unittest.main()
