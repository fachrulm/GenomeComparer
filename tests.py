import unittest
from main import *
import pandas as pd


class BestOneTester(unittest.TestCase):

    def test_match(self):

        row = {"CCODE": "="}

        self.assertEqual(bestone(row), 10)

    def test_various_class_codes(self):

        vals = [("=", 10),
                ("_", 9),
                ("f, =", 8),
                ("f, _", 7),
                ("n", 6),
                ('g', 5),
                ('f', 4),
                ("m", 3),
                ('rI', 2),
                ('x', 1)]

        for ccode, val in vals:
            row = {"CCODE": ccode}
            self.assertEqual(bestone(row), val)

    def test_wrong_input(self):

        row = []
        with self.assertRaises(TypeError):
            bestone(row)

        row = dict()
        with self.assertRaises(KeyError):
            bestone(row)

    def test_match2(self):

        row = {"CCODE": "f, _"}

        self.assertEqual(bestone(row), 7)

    #def test_dummy(self):
     #   with open("az.tsv", "r") as ax:
      #      df = bestone(ax)
       #     self.assertEqual(df[0], 10)

class CategTester(unittest.TestCase):

    def test_match(self):

        row = {'rank': 9}

        self.assertEqual(categ(row), 'Match')

    def test_various_class_codes(self):

        vals = [(9, 'Match'),
                (7, 'Fusion-Match'),
                (6, 'Extension'),
                (5, 'Alternative Splicing'),
                (4, 'Fusion'),
                (3, 'Overlap'),
                (2, 'Intronic'),
                (1, 'Fragment'),
                (0, 'Unknown')]

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

    def test_dummy_ref(self):
        refcol = ['TID_x', '# coding exons_x']
        with open("chr1A.reference.stats.tsv", "r") as testref:
            df = load_ref_stats(testref)
            self.assertEqual(list(df.columns), refcol)

class statsloadTester (unittest.TestCase):

    def test_name_change(self):
        refcol = ['TID_y', 'GID_y', '# Exon number_y']
        with open("1A_on_1B.stats.tsv", "r") as a:
            df = load_aligned_stats(a)
            self.assertEqual(list(df.columns), refcol)

class comparloadTester (unittest.TestCase):

    def test_name_change(self):
        refcol = ['TID', 'GID', 'CCODE', 'REF_ID', 'REF_GENE', 'NF1', 'EF1', 'JF1', "CONFIDENCE", "rank", "category"]
        with open("1D_on_1A.compare.refmap", "r") as a:
            df = load_comparisons(a)
            self.assertIsInstance(df, pd.DataFrame)
            self.assertEqual(list(df.columns), refcol)


class F1FileTester (unittest.TestCase):

    def test_F1(self):
        TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), '1D_on_1A.compare.refmap')
        self.testdata = open(TESTDATA_FILENAME).read()

        comparison = load_comparisons(TESTDATA_FILENAME)
        f1 = getf1(comparison, '1D', '1A')
        self.assertFalse(isinstance(f1, str))


class InitMergeTester(unittest.TestCase):

    def test_init(self):
        def init_merge(ref, aligned, comparison):
            merge = pd.merge(pd.merge(ref, aligned, left_on='TID_x', right_on='TID_y'),
                             comparison, left_on='TID_x', right_on='TID')
            return merge
        ref = pd.DataFrame(columns=['TID_x', '# coding exons_x'])
        aligned = pd.DataFrame(columns=['TID_y', "GID_y", '# Exon number_y'])
        compa = pd.DataFrame(columns=['TID', "GID", 'CCODE', "REF_ID", 'REF_GENE', 'rank', 'category'])
        merge = init_merge(ref, aligned, compa)

        cols = list(merge.columns)
        self.assertEqual(cols,
                         ['TID_x', '# coding exons_x', 'TID_y', 'GID_y', '# Exon number_y', 'TID', "GID", 'CCODE',
                          "REF_ID", 'REF_GENE', 'rank', 'category'])




    def test_one(self):

        df = pd.DataFrame(columns=["TID_y", 'GID', "REF_ID", 'REF_GENE', "# coding exons_x", 'CCODE', "extra"])
        df.loc[0] = ["tid_1", "gid_1", "ref_1", "ref_gid_1", 2, "=", "foo"]
        df.loc[1] = ["tid_2", "gid_2", "-", "-", 2, "-", "bar"]
        self.assertEqual(len(df), 2)
        new_df = pre_ref_merge(df, "A", "B")
        self.assertEqual(len(new_df), 1, new_df)
        cols = list(new_df.columns)
        self.assertEqual(cols,
                         ["A id", "A gene", "ref (B) id", "ref (B) gene", "A Exon(s)", "A-B ccode"])
        row = new_df[new_df["A id"] == "tid_2"]
        self.assertEqual(len(row), 0, row)
        row = new_df[new_df["A id"] == "tid_1"]
        self.assertEqual(len(row), 1, new_df)

    # def ref_merge(mergexz, mergeyz, x, y, z):
    #     mergexz = pre_ref_merge(mergexz, x, z)
    #     mergeyz = pre_ref_merge(mergeyz, y, z)
    #     zmerge = pd.merge(mergexz, mergeyz, how='inner', on=["ref ({}) id".format(z), "ref ({}) gene".format(z)])
    #     zmerge = zmerge[["{} id".format(x),
    #                      "{} gene".format(x),
    #                      "{} Exon(s)".format(x),
    #                      '{}-{} ccode'.format(x, z),
    #                      "{} id".format(y),
    #                      "{} gene".format(y),
    #                      "{} Exon(s)".format(y),
    #                      '{}-{} ccode'.format(y, z),
    #                      "ref ({}) id".format(z),
    #                      "ref ({}) gene".format(z)]]
    #     return zmerge

    def test_ref_merge_one(self):
        df = pd.DataFrame(columns=["TID_y", 'GID', "REF_ID", 'REF_GENE', "# coding exons_x", 'CCODE', "extra"])
        df_one = pd.DataFrame(columns=["TID_y", 'GID', "REF_ID", 'REF_GENE', "# coding exons_x", 'CCODE', "extra"])

        df.loc[0] = ["tidA_1", "gidA_1", "refB_1", "refB_gid_1", 2, "=", "foo"]
        df.loc[1] = ["tidA_2", "gidA_2", "-", "-", 2, "-", "bar"]

        df_one.loc[0] = ["tidD_1", "gidD_1", "refB_1", "refB_gid_1", 2, "=", "foo"]
        df_one.loc[1] = ["tidD_2", "gidD_2", "-", "-", 2, "-", "bar"]

        merged = ref_merge(df, df_one, "A", "D", "B")
        cols = list(merged.columns)
        self.assertEqual(cols,
                         ["A id", "A gene", "A Exon(s)", "A-B ccode", "D id", "D gene", "D Exon(s)", "D-B ccode",
                          "ref (B) id", "ref (B) gene"])
        self.assertEqual(len(merged), 1, merged)


    def test_sixway(self):
        df = pd.DataFrame(columns=["A id", "A gene", "A Exon(s)", "A-D ccode", "B id", "B gene", "B Exon(s)", "B-D ccode",
                     "ref (D) id", "ref (D) gene"])

        df_two = pd.DataFrame(columns=["A id", "A gene", "A Exon(s)", "A-B ccode", "D id", "D gene", "D Exon(s)", "D-B ccode",
                          "ref (B) id", "ref (B) gene"])

        df_three = pd.DataFrame(columns=["B id", "B gene", "B Exon(s)", "B-A ccode", "D id", "D gene", "D Exon(s)", "D-A ccode",
                          "ref (A) id", "ref (A) gene"])


        df.loc[0] = ["as_1", "gidA_1", 2, "=", "bs_1", "gidB_1",  2, "=", "ds_1", "gidD_1"]
        df.loc[1] = ["tidA_2", "gidA_2", 2, "-", "-", "gidD_2",  2, "-", "tidB_2", "refB_gid_2"]

        df_two.loc[0] = ["as_1", "gidA_1", 2, "=", "ds_1", "gidD_1",  3, "=", "bs_1", "gidB_1"]
        df_two.loc[1] = ["tidA_2", "gidA_2", 3, "-", "-", "-",  '-', "=", "tidD_2", "refD_gid_2"]

        df_three.loc[0] = ["bs_1", "gidB_1", 2, "=", "ds_1", "gidD_1", 2, "=", "as_1", "as_1"]
        df_three.loc[1] = ["tidB_2", "gidB_2", 3, "-", "-", "-", 3, "=", "tidA_2", "refA_gid_2"]

        sixwayy = sixway(df, df_two, df_three, "A", "B", "D")
        coln = list(sixwayy.columns)
        self.assertEqual(coln, ["A", "B", "D", 'A Exon(s)', 'B Exon(s)', 'D Exon(s)', 'A-B ccode', 'A-D ccode',
                                'B-A ccode', 'B-D ccode', 'D-A ccode', 'D-B ccode'])
        self.assertEqual(len(sixwayy), 1, sixwayy)

    def test_allmatch(self):

        df = pd.DataFrame(columns=["A", "B", "D", 'A Exon(s)', 'B Exon(s)', 'D Exon(s)', 'A-B ccode', 'A-D ccode',
                                'B-A ccode', 'B-D ccode', 'D-A ccode', 'D-B ccode'])
        df.loc[0] = ["tidA_1", "tidB_1", "tidD_1", 2, 2, 2, "=","=","=","=","=","="]
        df.loc[1] = ["tidA_1", "tidB_1", "tidD_1", 2, 3, 2, "=", "=", "=", "=", "=", "="]
        df.loc[2] = ["tidA_1", "tidB_1", "tidD_1", 2, 2, 2, "=", "c", "=", "=", "j", "="]
        df.loc[3] = ["tidA_1", "tidB_1", "tidD_1", 10, 2, 5, "=", "x", "=", "j", "=", "c"]

        math = allmatches(df, "A", "B", "D")




if __name__ == "__main__":
    unittest.main()
