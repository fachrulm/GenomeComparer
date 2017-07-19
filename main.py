#!/usr/bin/env python3

import sys  # System functions and streams, e.g. sys.stdout for STandard OUTput or sys.stderr for STandard ERRor
import os  # Library to perform OS-level functions, such as location of files, their sizes, etc.
import argparse  # Library to create the script interface on the shell
import pandas as pd
import numpy as np
import itertools
import re


__doc__ = """This script will analyse triplets of alignments to find gene triplets"""


def bestone(row):
    """Function to create rank column"""
    # Between the triple " is a doctstring, i.e. the documentation for the funtion


    if row['ccode'] == '=':
        val = 10
    elif row['ccode'] == '_':
        val = 9
    elif 'f' and '=' in row['ccode']:
        val = 8
    elif 'f' and '_' in row['ccode']:
        val = 7
    elif row['ccode'] in ('n', 'J', 'C', 'c'):
        val = 6
    elif row['ccode'] in ('j', 'h', 'g', 'G'):
        val = 5
    elif row['ccode'] in ('o', 'e', 'm'):
        val = 4
    elif row['ccode'] in ('i', 'I', 'ri', 'rI'):
        val = 3
    elif row['ccode'] in ('x', 'X', 'p', 'P'):
        val = 2
    else:
        val = 1
    return val


def strip_path(df):

    """Snippet to remove the .mrna and .path suffices from the dataframe"""

    df['TID'] = df['TID'].str.replace('.mrna[0-9]*', '')
    df['GID'] = df['GID'].str.replace('\.[0-9]*\.path[0-9]*', '')

    return df


def load_ref_stats(path):
    ref = pd.read_csv(path, sep='\t')
    ref = ref.loc[:, ['TID', '# coding exons']]
    ref.columns = ['TID_x', '# coding exons_x']
    return ref

def load_aligned_stats(path):
    ali = pd.read_csv(path, sep='\t')
    ali = ali.loc[:, ['TID', "GID", 'Exon number']]
    strip_path(ali)
    ali.columns = ['TID_y', "GID_y", '# Exon number_y']
    return ali

def load_comparisons(path):
    com = pd.read_csv(path, sep='\t')
    com = com.loc[:, ['tid', "gid", 'ccode', "ref_id", 'ref_gene']]
    com.columns = com.columns.str.upper()
    strip_path(com)
    com['RANK'] = com.apply(bestone, axis=1)
    com = com.loc[(com['RANK'] >= 9) | ~com['GID'].duplicated()]
    com = com.loc[(com['RANK'] >= 9) | ~com['REF_GENE'].duplicated()]
    return com

def init_merge(ref, aligned, comparison):
    merge = pd.merge(pd.merge(ref,aligned,left_on='TID_x',right_on='TID_y'),
                     comparison,left_on='TID_x',right_on='TID')
    return merge

def pre_ref_merge(merge, x, z):
    merge = merge[["TID_y", 'GID', "REF_ID", 'REF_GENE', "# coding exons_x", 'CCODE']]
    merge.columns = ["{} id".format(x),
                     '{} gene'.format(x),
                     "ref ({}) id".format(z),
                     "ref ({}) gene".format(z),
                     "{} Exon(s)".format(x),
                     '{}-{} ccode'.format(x, z)]
    merge.replace("-", np.nan, inplace=True)
    merge.dropna(inplace=True)
    return merge

def ref_merge(mergexz, mergeyz, x, y, z):
    mergexz = pre_ref_merge(mergexz, x, z)
    mergeyz = pre_ref_merge(mergeyz, y, z)
    zmerge = pd.merge(mergexz, mergeyz, how='inner', on=["ref ({}) id".format(z),"ref ({}) gene".format(z)])
    zmerge = zmerge[["{} id".format(x),
                     "{} gene".format(x),
                     "{} Exon(s)".format(x),
                     '{}-{} ccode'.format(x, z),
                     "{} id".format(y),
                     "{} gene".format(y),
                     "{} Exon(s)".format(y),
                     '{}-{} ccode'.format(y, z),
                     "ref ({}) id".format(z),
                     "ref ({}) gene".format(z)]]
    return zmerge

def sixway(a_merge,b_merge,d_merge):
    pretrip = pd.merge(a_merge,b_merge,left_on=['1D id','ref (1A) id','1B id'],right_on=['1D id','1A id','ref (1B) id'])
    pretrip = pretrip.drop_duplicates(subset='1D id', keep=False)
    triplets = pd.merge(pretrip,d_merge,left_on=['1A id','1B id','1D id'],right_on=['1A id','1B id','ref (1D) id'])
    triplets = triplets[['1A id', '1B id', '1D id', 'A Exon(s)_x', 'B Exon(s)_x', 'D Exon(s)_x', 'A-B ccode',
                              'A-D ccode', 'B-A ccode', 'B-D ccode', 'D-A ccode', 'D-B ccode']]
    triplets.columns = ['1A', '1B', '1D', 'A Exon(s)', 'B Exon(s)', 'D Exon(s)', 'A-B ccode', 'A-D ccode',
                          'B-A ccode', 'B-D ccode', 'D-A ccode', 'D-B ccode']
    return triplets


def crossch(ref2, triplets):
    """Crosscheck with consortium list"""


def main():

    """Main function for the utility"""

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("A", help="\"A\" genome name (eg. 1A)")
    parser.add_argument("B", help="\"B\" genome name (eg. 1B)")
    parser.add_argument("D", help="\"D\" genome name (eg. 1D)")
    parser.add_argument("--tmap",
                        help="Template for the TMAP files, eg ?_on_?.compare.tmap. It must contain two ?",
                        default="?_on_?.compare.tmap")
    parser.add_argument("--ref",
                        help="Template for the reference statistics, eg chr?.reference.stats.tsv. It must contain one ?",
                        default="chr?.reference.stats.tsv")
    parser.add_argument("--aligned",
                        help="Template for the aligned statistics, eg ?_on_?.stats.tsv. It must contain two ?",
                        default="?_on_?.stats.tsv")
    parser.add_argument("--inmerge",
                        help="Initial merge of the reference statistics, aligned statistics, and TMAP files. It must contain two ?",
                        default="nt_??")
    args = parser.parse_args()


    comparisons = dict()
    aligned_stats = dict()
    ref_stats = dict()

    genomes = [args.A, args.B, args.D]
    # Load the reference statistics into the dictionary
    for genome in genomes:
        ref_stats[genome] = load_ref_stats(re.sub("\?", genome, args.ref))

        aligned_stats[genome] = load_aligned_stats(re.sub("\?", genome, args.aligned))

        comparisons[genome] = load_comparisons(re.sub("\?", genome, args.tmap))

        ref = ref_stats[genome]
        aligned = aligned_stats[genome]
        comparison = comparisons[genome]

        init_merge[genome] = init_merge(re.sub("\?", genome, args.inmerge))

        merge = init_merge[genome]

        pre_ref_merge[genome] = pre_ref_merge()







    ref1A = pd.read_csv('chr1A.reference.stats.tsv', sep='\t')
    st_1A1B = pd.read_csv('1A_on_1B.stats.tsv', sep='\t')
    comp_1A1B = pd.read_csv('1A_on_1B.compare.tmap', sep='\t')







# If the script is called as a script (instead of being imported as a library), execute main()
if __name__ == "__main__":
    main()