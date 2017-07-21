import unittest
from main import *
import pandas as pd


class BestOneTester(unittest.TestCase):

    def test_match(self):

        row = {"ccode": "="}

        self.assertEqual(bestone(row), 10)

    def test_various_class_codes(self):

        vals = [("=", 10),
                ("_", 9),
                ("f, =", 8),
                ("f, _", 7),
                ("n", 6),
                ('g', 5),
                ("m", 4),
                ('rI', 3),
                ('x', 2)]

        for ccode, val in vals:
            row = {"ccode": ccode}
            self.assertEqual(bestone(row), val)

    def test_wrong_input(self):

        row = []
        with self.assertRaises(TypeError):
            bestone(row)

        row = dict()
        with self.assertRaises(KeyError):
            bestone(row)

    def test_match2(self):

        row = {"ccode": "f, _"}

        self.assertEqual(bestone(row), 7)

class CategTester(unittest.TestCase):

    def test_match(self):

        row = {'rank': 9}

        self.assertEqual(categ(row), 'Match')

    def test_various_class_codes(self):

        vals = [(9, 'Match'),
                (7, 'Fusion-Match'),
                (6, 'Extension'),
                (5, 'Alternative Splicing'),
                (4, 'Overlap'),
                (3, 'Intronic'),
                (2, 'Fragment'),
                (1, 'Unknown')]

        for category, val in vals:
            row = {"rank": category}
            self.assertEqual(categ(row), val)


class refloadTester(unittest.TestCase):

    # @unittest.skip
    def test_wrong_input(self):

        row = ''
        with self.assertRaises(FileNotFoundError):
            load_ref_stats(row)

        row = 0
        with self.assertRaises(ValueError):
            load_ref_stats(row)

    def test_name_change(self):
        refcol = ['TID_x', '# coding exons_x']
        with open("sample.tsv", "r") as a:
            df = load_ref_stats(a)
            self.assertEqual(list(df.columns), refcol)





if __name__ == "__main__":
    unittest.main()
