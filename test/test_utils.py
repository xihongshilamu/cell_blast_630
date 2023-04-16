import os
import sys
import unittest

import numpy as np

if os.environ.get("TEST_MODE", "INSTALL") == "DEV":
    sys.path.insert(0, "..")
import Cell_BLAST as cb


class TestJsonCellTypeTree(unittest.TestCase):

    def setUp(self):
        self.dag = cb.utils.CellTypeDAG.load("./tree.json")

    def test_load(self):
        with self.assertRaises(ValueError):
            dag = cb.utils.CellTypeDAG.load("./asd.asd")

    def test_alias(self):
        self.assertIn("root", self.dag.vdict)
        l1 = ["a", "b", "c", "d", "e", "f", "g", "h"]
        l2 = ["A", "B", "C", "D", "E", "F", "G", "H"]
        for n1, n2 in zip(l1, l2):
            self.assertEqual(self.dag.get_vertex(n1), self.dag.get_vertex(n2))

    def test_relation(self):
        self.assertFalse(self.dag.is_related("b", "c"))
        self.assertFalse(self.dag.is_related("d", "i"))
        self.assertTrue(self.dag.is_related("g", "e"))
        self.assertTrue(self.dag.is_descendant_of("g", "e"))
        self.assertFalse(self.dag.is_descendant_of("e", "g"))
        self.assertTrue(self.dag.is_ancestor_of("e", "g"))
        self.assertFalse(self.dag.is_ancestor_of("g", "e"))

    def test_similarity(self):
        self.assertAlmostEqual(self.dag.similarity("g", "h"), 0.5)
        self.assertAlmostEqual(self.dag.similarity("g", "f"), 0.75)
        self.assertAlmostEqual(self.dag.similarity("g", "e"), 0.75)
        self.assertAlmostEqual(self.dag.similarity("b", "d"), 0.25)
        self.assertAlmostEqual(self.dag.similarity("h", "c"), 1 / 6)
        self.assertAlmostEqual(self.dag.similarity("h", "i"), 0)

    def test_value(self):
        self.dag.value_reset()
        for item in self.dag.vdict:
            self.assertEqual(self.dag.get_vertex(item)["value"], 0)
        self.dag.value_set("g", 5)
        self.dag.value_set("f", 2)
        self.dag.value_set("c", 3)
        self.dag.value_update()
        self.assertEqual(self.dag.get_vertex("root")["value"], 10)
        self.assertEqual(self.dag.get_vertex("a")["value"], 3)
        self.assertEqual(self.dag.get_vertex("b")["value"], 0)
        self.assertEqual(self.dag.get_vertex("c")["value"], 3)
        self.assertEqual(self.dag.get_vertex("d")["value"], 0)
        self.assertEqual(self.dag.get_vertex("e")["value"], 7)
        self.assertEqual(self.dag.get_vertex("f")["value"], 7)
        self.assertEqual(self.dag.get_vertex("g")["value"], 5)
        self.assertEqual(self.dag.get_vertex("h")["value"], 0)
        result = self.dag.best_leaves(4, min_path=0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "g")
        result = self.dag.best_leaves(6, min_path=0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "f")


class TestOBOCellTypeTree(TestJsonCellTypeTree):

    def setUp(self):
        self.dag = cb.utils.CellTypeDAG.load("./tree.obo")

    def test_value(self):
        self.dag.value_reset()
        for item in self.dag.vdict:
            self.assertEqual(self.dag.get_vertex(item)["value"], 0)
        self.dag.value_set("g", 5)
        self.dag.value_set("f", 2)
        self.dag.value_set("c", 3)
        self.dag.value_update()
        self.assertEqual(self.dag.get_vertex("root")["value"], 10)
        self.assertEqual(self.dag.get_vertex("a")["value"], 3)
        self.assertEqual(self.dag.get_vertex("b")["value"], 0)
        self.assertEqual(self.dag.get_vertex("c")["value"], 3)
        self.assertEqual(self.dag.get_vertex("d")["value"], 0)
        self.assertEqual(self.dag.get_vertex("e")["value"], 7)
        self.assertEqual(self.dag.get_vertex("f")["value"], 7)
        self.assertEqual(self.dag.get_vertex("g")["value"], 5)
        self.assertEqual(self.dag.get_vertex("h")["value"], 0)
        result = self.dag.best_leaves(4, retrieve="cell_ontology_class", min_path=0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "g")
        result = self.dag.best_leaves(6, retrieve="cell_ontology_class", min_path=0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "f")


class TestOBOCellTypeDAG(TestJsonCellTypeTree):

    def setUp(self):
        self.dag = cb.utils.CellTypeDAG.load("./dag.obo")

    def test_similarity(self):
        self.assertAlmostEqual(self.dag.similarity("g", "h"), 7 / 24)
        self.assertAlmostEqual(self.dag.similarity("g", "f"), 0.75)
        self.assertAlmostEqual(self.dag.similarity("g", "e"), 0.75)
        self.assertAlmostEqual(self.dag.similarity("b", "d"), 0.25)
        self.assertAlmostEqual(self.dag.similarity("h", "c"), 1 / 12)
        self.assertAlmostEqual(self.dag.similarity("h", "i"), 0)

    def test_value(self):
        self.dag.value_reset()
        for item in self.dag.vdict:
            self.assertEqual(self.dag.get_vertex(item)["value"], 0)
        self.dag.value_set("g", 5)
        self.dag.value_set("f", 2)
        self.dag.value_set("c", 3)
        self.dag.value_update()
        self.assertEqual(self.dag.get_vertex("root")["value"], 10)
        self.assertEqual(self.dag.get_vertex("a")["value"], 3)
        self.assertEqual(self.dag.get_vertex("b")["value"], 0)
        self.assertEqual(self.dag.get_vertex("c")["value"], 3)
        self.assertEqual(self.dag.get_vertex("d")["value"], 3)
        self.assertEqual(self.dag.get_vertex("e")["value"], 7)
        self.assertEqual(self.dag.get_vertex("f")["value"], 7)
        self.assertEqual(self.dag.get_vertex("g")["value"], 5)
        self.assertEqual(self.dag.get_vertex("h")["value"], 0)
        result = self.dag.best_leaves(4, retrieve="cell_ontology_class", min_path=0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "g")
        result = self.dag.best_leaves(6, retrieve="cell_ontology_class", min_path=0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], "f")


class TestHybridPath(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.m1 = np.array([
            [2, 1, 0],
            [3, 0, 0],
            [0, 4, 5]
        ])
        cls.v2 = np.array(["a", "s", "d"])

    def test_hybrid_path(self):
        cb.utils.write_hybrid_path(self.m1, "./test.h5//a")
        cb.utils.write_hybrid_path(self.v2, "./test.h5//b/c")
        m1 = cb.utils.read_hybrid_path("./test.h5//a")
        v2 = cb.utils.read_hybrid_path("./test.h5//b/c")
        self.assertTrue(np.all(self.m1 == m1))
        self.assertTrue(np.all(self.v2 == v2))

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("./test.h5"):
            os.remove("./test.h5")


if __name__ == "__main__":
    unittest.main()
