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

    # assert isinstance(row, pd.DataFrame), (row, row.keys())

    # val = 1
    if row["CCODE"] is np.nan:
        val = 1
    elif row["CCODE"] == '=':
        val = 10
    elif row["CCODE"] == '_':
        val = 9
    elif 'f' in row["CCODE"]:
        if "=" in row["CCODE"]:
            val = 8
        elif '_' in row["CCODE"]:
            val = 7
        else:
            val = 4
    elif row["CCODE"] in ('n', 'J', 'C', 'c'):
        val = 6
    elif row["CCODE"] in ('j', 'h', 'g', 'G'):
        val = 5
    elif row["CCODE"] in ('o', 'e', 'm'):
        val = 3
    elif row["CCODE"] in ('i', 'I', 'ri', 'rI'):
        val = 2
    elif row["CCODE"] in ('x', 'X', 'p', 'P'):
        val = 1
    else:
        val = 0

    return val


def categ(row):
    if row['rank'] in (10, 9):
        val = 'Match'
    elif row['rank'] in (8, 7):
        val = 'Fusion-Match'
    elif row['rank'] == 6:
        val = 'Extension'
    elif row['rank'] == 5:
        val = 'Alternative Splicing'
    elif row['rank'] == 4:
        val = 'Fusion'
    elif row['rank'] == 3:
        val = 'Overlap'
    elif row['rank'] == 2:
        val = 'Intronic'
    elif row['rank'] == 1:
        val = 'Fragment'
    else:
        val = 'Unknown'
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

    #print(ref.columns)
    return ref


def load_aligned_stats(path):
    ali = pd.read_csv(path, sep='\t')
    ali = ali.loc[:, ['TID', "GID", 'Exon number']]
    strip_path(ali)
    ali.columns = ['TID_y', "GID_y", '# Exon number_y']

    #print(ali.columns)
    return ali


def load_comparisons(path):
    """Input is refmap file"""
    # refmap = dict()
    com = pd.read_csv(path, sep='\t')
    com = com.loc[:, ['tid', "gid", 'ccode', "ref_id", 'ref_gene', "nF1", "eF1", "jF1"]]
    com.columns = com.columns.str.upper()

    strip_path(com)
    a = ~(com['REF_ID'].str.contains('LC'))

    com["CONFIDENCE"] = a
    #com["CONFIDENCE"] = str()

    #for index, row in com.iterrows():
    #    if row["TID"] is np.nan:
    #        #com.set_value(row, "CONFIDENCE", 'None')
    #        row["CONFIDENCE"] = 'None'
    #        print(row)
    #    elif row['REF_ID'].str.contains('LC'):
    #        row["CONFIDENCE"] = False
    #    else:
    #        row["CONFIDENCE"] = True

    # refmap['LC'] = com.loc[a]
    # refmap['HC'] = com.loc[~a]
    com.dropna()
    com["rank"] = com.apply(bestone, axis=1)
    com["category"] = com.apply(categ, axis=1)
    com = com.replace(np.nan, "-")

    # for i in ['LC','HC']:
    #     refmap[i] = refmap[i].dropna()
    #     refmap[i]['rank'] = refmap[i].apply()
    #     refmap[i]['category'] = refmap[i]

    #print(com.columns)
    with open("compa.tsv", "wt") as out:
        com.to_csv(out, sep="\t")

    #comparHC = pd.DataFrame
    #comparLC = pd.DataFrame

    for index, row in com.iterrows():
        # print(row['CONFIDENCE'])
        if row['CONFIDENCE'] is 'False':
            comparHC = com.drop(com.index[row])
        elif row['CONFIDENCE'] is 'True':
            comparLC = com.drop(com.index[row])
    #print(com)
    #with open("comparHC.tsv", "wt") as out:
        #comparHC.to_csv(out, sep="\t")


    return com

def getf1(comparison):

    print("Category", *["{} {}".format(*_) for _ in itertools.product(["NF1", "EF1", "jF1"], ["mean", "StDEV"])],
          sep="\t")
    for category in comparison["category"].unique():
        row = [category]
        means = comparison[comparison["category"] == category][["NF1", "EF1", "JF1"]].mean().astype(list)
        stdev = comparison[comparison["category"] == category][["NF1", "EF1", "JF1"]].std().astype(list)
        for md, std in zip(means, stdev):
            row.extend([round(md, 2), round(std, 2)])
        print(*row, sep="\t")



def init_merge(ref, aligned, comparison):
    merge = pd.merge(pd.merge(ref,aligned,left_on='TID_x',right_on='TID_y'),
                     comparison,left_on='TID_x',right_on='TID')

    #with open("merge.tsv", "wt") as out:
    #    merge.to_csv(out, sep="\t")
    #print(merge)
    return merge

def pre_ref_merge(merge, x, z):
    new_df = merge[["TID_y", 'GID', "REF_ID", 'REF_GENE', "# coding exons_x", 'CCODE']]
    new_df.columns = ["{} id".format(x),
                     '{} gene'.format(x),
                     "ref ({}) id".format(z),
                     "ref ({}) gene".format(z),
                     "{} Exon(s)".format(x),
                     '{}-{} ccode'.format(x, z)]
    new_df = new_df.replace("-", np.nan)
    new_df = new_df.dropna()

    #print(new_df)
    return new_df

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
    with open("zmerge.tsv", "wt") as out:
        zmerge.to_csv(out, sep="\t")
    #print(zmerge)
    return zmerge


def sixway(xy_z_merge, xz_y_merge, yz_x_merge, x, y, z):

    """

    :param xy_z_merge:
    :type xy_z_merge: pd.DataFrame
    :param xz_y_merge:
    :type xz_y_merge: pd.DataFrame
    :param yz_x_merge:
    :type yz_x_merge: pd.DataFrame
    :param x:
    :param y:
    :param z:
    :return:
    """

    pre_trip = pd.merge(xy_z_merge, xz_y_merge,
                        left_on=['{} id'.format(x), 'ref ({}) id'.format(z),
                                 '{} id'.format(y)],
                        right_on=["{} id".format(x), "{} id".format(z), "ref ({}) id".format(y)])

    pre_trip = pre_trip.drop_duplicates(subset='{} id'.format(x), keep=False)

    triplets = pd.merge(pre_trip,
                        yz_x_merge,
                        left_on=['{} id'.format(x), '{} id'.format(y), '{} id'.format(z)],
                        right_on=['ref ({}) id'.format(x), '{} id'.format(y), '{} id'.format(z)])
    triplets = triplets[['{} id'.format(z), '{} id'.format(y), '{} id'.format(x),
                         '{} Exon(s)_x'.format(z), '{} Exon(s)_x'.format(y), '{} Exon(s)_x'.format(x),
                         '{}-{} ccode'.format(z, y), '{}-{} ccode'.format(z, x), '{}-{} ccode'.format(y, z),
                         '{}-{} ccode'.format(y, x), '{}-{} ccode'.format(x, z), '{}-{} ccode'.format(x, y)]]

    triplets.columns = [x, y, z, '{} Exon(s)'.format(x), '{} Exon(s)'.format(y), '{} Exon(s)'.format(z),
                        '{}-{} ccode'.format(x, y), '{}-{} ccode'.format(x, z), '{}-{} ccode'.format(y, x),
                        '{}-{} ccode'.format(y, z), '{}-{} ccode'.format(z, x), '{}-{} ccode'.format(z, y)]

    with open("triplets.tsv", "wt") as out:
        triplets.to_csv(out, sep="\t")

    #print(triplets)
    return triplets


def crossch(existl, triplets, x, y, z):
    """Crosscheck with consortium list"""
    for i in [x, y ,z]:
        triplets[i] = triplets[i].str.replace('\.[0-9]*','')
    bothm = pd.merge(triplets, existl, left_on=[x, y, z], right_on=['A', 'B', 'D'],
                        indicator=True, how='inner')
    return bothm


def excode(xlist, x, y, z):
    """Get list of exons and ccode"""
    refs = dict()
    bothm_ex = dict()
    for i in [x, y ,z]:

        refs[i] = pd.read_csv('chr{}.reference.stats.tsv'.format(i),sep='\t')
        refs[i] = refs[i][['GID', '# coding exons']]
        refs[i] = refs[i].sort_values(['GID','# coding exons'], ascending=False)
        refs[i].columns = ['GID', '{} coding exons'.format(i)]
        refs[i] = refs[i].drop_duplicates(subset='GID', keep='first')


        bothm_ex = pd.merge(xlist,refs[i],left_on=['A'],right_on=['GID'])
    bothm_ex = bothm_ex[['1A', '1B', '1D', 'A Exon(s)', 'B Exon(s)', 'D Exon(s)', 'A-B ccode', 'A-D ccode', 'B-A ccode',
         'B-D ccode', 'D-A ccode', 'D-B ccode']]
    return bothm_ex

def allmatches(df,x,y,z):
    xx = (df['{} Exon(s)'.format(x)] == df['{} Exon(s)'.format(y)]) & \
        (df['{} Exon(s)'.format(y)] == df['{} Exon(s)'.format(z)])
    a = df['{}-{} ccode'.format(x, y)].isin(['=', '_']) & df['{}-{} ccode'.format(x, z)].isin(['=', '_']) &\
        df['{}-{} ccode'.format(y, x)].isin(['=', '_']) & df['{}-{} ccode'.format(y, z)].isin(['=', '_']) & \
        df['{}-{} ccode'.format(z, x)].isin(['=', '_']) & df['{}-{} ccode'.format(z, y)].isin(['=', '_'])
    triplets_eq = df[xx & a]

    with open("triplets_eq.tsv", "wt") as out:
        triplets_eq.to_csv(out, sep="\t")

    #print(triplets_eq)

    return triplets_eq



def main():

    """Main function for the utility"""

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("A", help="\"A\" genome name (eg. 1A)")
    parser.add_argument("B", help="\"B\" genome name (eg. 1B)")
    parser.add_argument("D", help="\"D\" genome name (eg. 1D)")
    parser.add_argument("--chroms", help="Chromosome numbers of wheat; 1 to 7")
    parser.add_argument("--refmap",
                        help="Template for the REFMAP files, eg ?_on_?.compare.refmap. It must contain two ?",
                        default="?_on_?.compare.refmap")
    parser.add_argument("--ref",
                        help="Template for the reference statistics, eg chr?.reference.stats.tsv. It must contain one ?",
                        default="chr?.reference.stats.tsv")
    parser.add_argument("--aligned",
                        help="Template for the aligned statistics, eg ?_on_?.stats.tsv. It must contain two ?",
                        default="?_on_?.stats.tsv")
    parser.add_argument("--existl",
                        help="Existing df to compare to, eg wheat.homeolog_groups.TRIADS.chr?. The ? is the chromosome number",
                        default="wheat.homeolog_groups.TRIADS.chr?")
    parser.add_argument("--triplets",
                        help="Df of generated triplets.")
    parser.add_argument("--comparcon", help="Output file")
    parser.add_argument("--out", help="Output file")

    # parser.add_argument("--inmerge",
    #                     help="Initial merge of the reference statistics, aligned statistics, and TMAP files. It must contain two ?",
    #                     default="nt_??")
    args = parser.parse_args()


    comparisons = dict()
    aligned_stats = dict()
    ref_stats = dict()
    initial_merge = dict()
    refr_merge = dict()
    triplets = dict()
    getF1 = dict()
    genomes = [args.A, args.B, args.D]
    comparHC = pd.DataFrame(columns=list(itertools.permutations(genomes, 2)))
    comparLC = pd.DataFrame(columns=list(itertools.permutations(genomes, 2)))

    # Load the reference statistics into the dictionary
    for genome in genomes:
        ref_stats[genome] = load_ref_stats(re.sub("\?", genome, args.ref))

    template = re.sub("\?", "{}", args.aligned)
    comp_template = re.sub("\?", "{}", args.refmap)

    for x, y in itertools.permutations(genomes, 2):
        aligned_stats[(x, y)] = load_aligned_stats(template.format(x, y))

    for x, y in itertools.permutations(genomes, 2):
        comparisons[(x, y)] = load_comparisons(comp_template.format(x, y))
        #print(comparisons[(x, y)].columns)
        for index, row in comparisons[(x, y)].iterrows():
            #print(row)
            if row['CONFIDENCE'] == 'False':
                comparHC[(x, y)].append(row)
                #comparHC[(x, y)] = comparisons[(x, y)].drop(comparisons[(x, y)].index[row])
            elif row['CONFIDENCE'] == 'True':
                comparLC[(x, y)].append(row)
                #comparLC[(x, y)] = comparisons[(x, y)].drop(comparisons[(x, y)].index[row])
        #print(comparHC[(x, y)])

        getF1[(x, y)] = getf1(comparisons[(x, y)])


    for x, y in itertools.permutations(genomes, 2):
        initial_merge[(x, y)] = init_merge(ref_stats[x], aligned_stats[(x, y)], comparisons[(x, y)])
        pre_ref_merge(initial_merge[(x, y)], x, y)

    for x, y, z in itertools.permutations(genomes, 3):
        refr_merge[(x, y)] = ref_merge(initial_merge[(x, z)], initial_merge[(y, z)], x, y, z)

    for z, y, x in itertools.combinations(genomes, 3):
        triplets = sixway(refr_merge[(x, y)], refr_merge[(x, z)], refr_merge[(y, z)], x, y, z)

    for x, y, z in itertools.combinations(genomes, 3):
        triplets_eq = allmatches(triplets, x, y, z)


    #for x, y, z in genomes:
    #    for n in list(range(1, 8)):
    #        crosschecked = crossch(re.sub("\?", n, args.existl), triplets, x, y, z)

    # for



    # Let's say the result is a pandas DataFrame
    df = pd.DataFrame()
    df.to_csv(args.out, sep="\t")





# If the script is called as a script (instead of being imported as a library), execute main()
if __name__ == "__main__":
    main()