#!/usr/bin/env python3

import sys  # System functions and streams, e.g. sys.stdout for STandard OUTput or sys.stderr for STandard ERRor
import os  # Library to perform OS-level functions, such as location of files, their sizes, etc.
import argparse  # Library to create the script interface on the shell
import pandas as pd
import numpy as np
import itertools
import re
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
plt.style.use('bmh')

__doc__ = """This script will analyse triplets of alignments to find gene triplets"""


def bestone(row):
    """Function to create rank column"""
    # Between the triple " is a doctstring, i.e. the documentation for the funtion

    # assert isinstance(row, pd.DataFrame), (row, row.keys())

    # val = 1
    if row["CCODE"] == '=':
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
    elif row['CCODE'] is 'blank':
        val = 0
    else:
        val = 66

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
    elif row['rank'] == 0:
        val = 'Not aligned'
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

    com = com.replace(np.nan, "blank")
    #com.dropna()
    com["rank"] = com.apply(bestone, axis=1)
    com["category"] = com.apply(categ, axis=1)
    #com = com.replace(np.nan, "-")

    # for i in ['LC','HC']:
    #     refmap[i] = refmap[i].dropna()
    #     refmap[i]['rank'] = refmap[i].apply()
    #     refmap[i]['category'] = refmap[i]

    #print(com.columns)
    with open("compa.tsv", "wt") as out:
        com.to_csv(out, sep="\t")

    #comparHC = pd.DataFrame
    #comparLC = pd.DataFrame

    #for index, row in com.iterrows():
        # print(row['CONFIDENCE'])
    #    if row['CONFIDENCE'] is 'False':
    #        comparHC = com.drop(com.index[row])
    #    elif row['CONFIDENCE'] is 'True':
    #        comparLC = com.drop(com.index[row])
    #print(com)
    #with open("comparHC.tsv", "wt") as out:
        #comparHC.to_csv(out, sep="\t")


    return com

def pieplot(comparison,x,y,confidence):
    leng = len(comparison)
    a = comparison.groupby(["category"]).size().reset_index().rename(columns={0: 'count'})
    #print("N =", sum(a['count']))
    a.plot(labels=a["category"], y='count', kind='pie', autopct='%1.1f%%', startangle=90, figsize=(9, 9), fontsize=11)
    plt.axis('equal')
    plt.axis('off')
    plt.title('{} on {} {} (n = {})'.format(x, y, confidence, leng), y=1.15, fontsize=24)
    plt.legend(loc="best", labels=a["category"], fontsize=12)
    plt.tight_layout()

    plotsv = "%s.png" % '{}{}_{}'.format(x, y, confidence)
    plt.savefig(plotsv)
    plt.close()

def getf1(comparison,x,y):
    #F1 = open('{}_on_{}_F1.text'.format(x,y), 'w')
    #sys.stdout = F1
    print("For {} on {}".format(x,y))
    print("Category", *["{} {}".format(*_) for _ in itertools.product(["NF1", "EF1", "jF1"], ["mean", "StDEV"])],
          sep="\t")
    for category in comparison["category"].unique():
        row = [category]
        means = comparison[comparison["category"] == category][["NF1", "EF1", "JF1"]].mean().astype(list)
        stdev = comparison[comparison["category"] == category][["NF1", "EF1", "JF1"]].std().astype(list)
        for md, std in zip(means, stdev):
            row.extend([round(md, 2), round(std, 2)])
        print(*row, sep="\t")
    print("###########\n")
    #F1.close()


def init_merge(ref, aligned, comparison):
    merge = pd.merge(pd.merge(ref,aligned,left_on='TID_x',right_on='TID_y'),
                     comparison,left_on='TID_x',right_on='TID')
    #with open("merge.tsv", "wt") as out:
    #    merge.to_csv(out, sep="\t")
    #print(merge)
    return merge

def pre_ref_merge(merge, x, z):
    new_df = merge[["TID_y", 'GID', "REF_ID", 'REF_GENE', "# coding exons_x", 'CCODE',"NF1", "EF1", "JF1"]]
    new_df.columns = ["{} id".format(x),
                     '{} gene'.format(x),
                     "ref ({}) id".format(z),
                     "ref ({}) gene".format(z),
                     "{} Exon(s)".format(x),
                     '{}-{} ccode'.format(x, z),
                      '{}-{} nF1'.format(x, z),
                      '{}-{} eF1'.format(x, z),
                      '{}-{} jF1'.format(x, z)]
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
                     '{}-{} nF1'.format(x, z),
                     '{}-{} eF1'.format(x, z),
                     '{}-{} jF1'.format(x, z),
                     "{} id".format(y),
                     "{} gene".format(y),
                     "{} Exon(s)".format(y),
                     '{}-{} ccode'.format(y, z),
                     '{}-{} nF1'.format(y, z),
                     '{}-{} eF1'.format(y, z),
                     '{}-{} jF1'.format(y, z),
                     "ref ({}) id".format(z),
                     "ref ({}) gene".format(z)]]
    #with open("zmerge.tsv", "wt") as out:
    #    zmerge.to_csv(out, sep="\t")
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

    triplets['nF1_mean'] = triplets[['{}-{} nF1'.format(x, y), '{}-{} nF1'.format(x, z),
                                     '{}-{} nF1'.format(y, x), '{}-{} nF1'.format(y, z), '{}-{} nF1'.format(z, x),
                            '{}-{} nF1'.format(z, y)]].mean(axis=1)

    triplets['eF1_mean'] = triplets[['{}-{} eF1'.format(x, y), '{}-{} eF1'.format(x, z),
                                     '{}-{} eF1'.format(y, x), '{}-{} eF1'.format(y, z), '{}-{} eF1'.format(z, x),
                                     '{}-{} eF1'.format(z, y)]].mean(axis=1)

    triplets['jF1_mean'] = triplets[['{}-{} jF1'.format(x, y), '{}-{} jF1'.format(x, z),
                                     '{}-{} jF1'.format(y, x), '{}-{} jF1'.format(y, z), '{}-{} jF1'.format(z, x),
                                     '{}-{} jF1'.format(z, y)]].mean(axis=1)

    triplets = triplets[['{} id'.format(x), '{} id'.format(y), '{} id'.format(z),
                         '{} Exon(s)_x'.format(x), '{} Exon(s)_x'.format(y), '{} Exon(s)_x'.format(z),
                         '{}-{} ccode'.format(x, y), '{}-{} ccode'.format(x, z), '{}-{} ccode'.format(y, z),
                         '{}-{} ccode'.format(y, x), '{}-{} ccode'.format(z, x), '{}-{} ccode'.format(z, y), 'nF1_mean',
                         'eF1_mean', 'jF1_mean']]

    triplets.columns = [x, y, z, '{} Exon(s)'.format(x), '{} Exon(s)'.format(y), '{} Exon(s)'.format(z),
                        '{}-{} ccode'.format(x, y), '{}-{} ccode'.format(x, z), '{}-{} ccode'.format(y, x),
                        '{}-{} ccode'.format(y, z), '{}-{} ccode'.format(z, x), '{}-{} ccode'.format(z, y), 'nF1_mean',
                        'eF1_mean', 'jF1_mean']

    #print(triplets)
    return triplets

def histogram(triplets,confidence,categ,x,y,z):
    a = triplets.groupby(["nF1_mean","eF1_mean","jF1_mean"]).size().reset_index().rename(columns={0: 'count'})
    #aq = triplets_eq.groupby(["nF1_mean", "eF1_mean", "jF1_mean"]).size().reset_index().rename(columns={0: 'count'})


    a.hist(column=["nF1_mean", "eF1_mean", "jF1_mean"], bins=10, figsize=[8, 8], range=(0, 101))
    plt.suptitle('F1 Stats of {} Triplets ({})'.format(confidence,categ), fontsize=22, fontweight='bold')
    #plt.tight_layout()
    plotsv = "%s.png" % 'F1 Statistics of {} {}{}{} Triplets ({})'.format(confidence,x,y,z,categ)
    plt.savefig(plotsv)

    #aq.hist(column=["nF1_mean", "eF1_mean", "jF1_mean"], bins=10, figsize=[8, 8], range=(0, 101))
    #plt.suptitle("F Stats (All Categories)", fontsize=22, fontweight='bold')
    # plt.tight_layout()
    #plotsv = "%s.png" % 'F1 Statstics of {} Triplets (Exact Matches)'.format(confidence)
    #plt.savefig(plotsv)

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

    #with open("triplets_eq.tsv", "wt") as out:
    #    triplets_eq.to_csv(out, sep="\t")

    #print(triplets_eq)

    return triplets_eq

def nomatches(df,x,y,z):
    xx = (df['{} Exon(s)'.format(x)] == df['{} Exon(s)'.format(y)]) & \
        (df['{} Exon(s)'.format(y)] == df['{} Exon(s)'.format(z)])
    a = df['{}-{} ccode'.format(x, y)].isin(['=', '_']) & df['{}-{} ccode'.format(x, z)].isin(['=', '_']) &\
        df['{}-{} ccode'.format(y, x)].isin(['=', '_']) & df['{}-{} ccode'.format(y, z)].isin(['=', '_']) & \
        df['{}-{} ccode'.format(z, x)].isin(['=', '_']) & df['{}-{} ccode'.format(z, y)].isin(['=', '_'])
    triplets_nomatch = df[~xx & ~a]

    #with open("triplets_eq.tsv", "wt") as out:
    #    triplets_eq.to_csv(out, sep="\t")

    #print(triplets_eq)

    return triplets_nomatch

def percent(cross,lst):
    trpcrs = (len(cross)/len(lst))*100
    trpcrs = str(round(trpcrs, 2))

    return trpcrs



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
    parser.add_argument("--F1", help="Output file for nF1, eF1, and jF1 statistics")
    parser.add_argument("--out", help="Output file")

    # parser.add_argument("--inmerge",
    #                     help="Initial merge of the reference statistics, aligned statistics, and TMAP files. It must contain two ?",
    #                     default="nt_??")
    args = parser.parse_args()


    comparisons = dict()
    aligned_stats = dict()
    ref_stats = dict()
    initial_merge_HC = dict()
    initial_merge_LC = dict()
    refr_merge_HC = dict()
    refr_merge_LC = dict()
    triplets_all = dict()
    triplets_all_eq = dict()
    triplets_all_noneq = dict()
    getF1 = dict()
    genomes = [args.A, args.B, args.D]
    comparHC = dict()  # pd.DataFrame(columns=list(itertools.permutations(genomes, 2)))
    comparLC = dict()  # pd.DataFrame(columns=list(itertools.permutations(genomes, 2)))
    triplets_HC = dict()
    triplets_LC = dict()
    triplets_HC_eq = dict()
    triplets_LC_eq = dict()
    triplets_HC_noneq = dict()
    triplets_LC_noneq = dict()

    # Load the reference statistics into the dictionary
    for genome in genomes:
        ref_stats[genome] = load_ref_stats(re.sub("\?", genome, args.ref))

    template = re.sub("\?", "{}", args.aligned)
    comp_template = re.sub("\?", "{}", args.refmap)

    for x, y in itertools.product(genomes, repeat=2):
        aligned_stats[(x, y)] = load_aligned_stats(template.format(x, y))

    for x, y in itertools.product(genomes, repeat=2):
        comparisons[(x, y)] = load_comparisons(comp_template.format(x, y))
        #print(comparisons[(x, y)].columns)
        #comparLC[(x, y)] = pd.DataFrame()
        #comparHC[(x, y)] = pd.DataFrame()
        comparLC[(x, y)] = comparisons[(x, y)][comparisons[(x, y)]["CONFIDENCE"] == False]
        #print(comparLC[(x, y)])

        comparHC[(x, y)] = comparisons[(x, y)][comparisons[(x, y)]["CONFIDENCE"] == True]

        # for index, row in comparisons[(x, y)].iterrows():
        #     #print(row)
        #     if row['CONFIDENCE'] is 'False':
        #         comparLC[(x, y)].append(row, ignore_index=True)
        #         #comparHC[(x, y)] = comparisons[(x, y)].drop(comparisons[(x, y)].index[row])
        #     elif row['CONFIDENCE'] is'True':
        #         comparHC[(x, y)].append(row, ignore_index=True)
        #         #comparLC[(x, y)] = comparisons[(x, y)].drop(comparisons[(x, y)].index[row])
        pieplot(comparisons[(x, y)], x, y, 'Overall')
        pieplot(comparHC[(x, y)], x, y, 'HC')
        pieplot(comparLC[(x, y)], x, y, 'LC')


    with open('All_F1.text', 'w') as F1:
        for x, y in itertools.permutations(genomes, 2):
            sys.stdout = F1
            getF1[(x, y)] = getf1(comparisons[(x, y)],x,y)
    F1.close()


    for x, y in itertools.permutations(genomes, 2):
        initial_merge_HC[(x, y)] = init_merge(ref_stats[x], aligned_stats[(x, y)], comparHC[(x, y)])
        pre_ref_merge(initial_merge_HC[(x, y)], x, y)

        initial_merge_LC[(x, y)] = init_merge(ref_stats[x], aligned_stats[(x, y)], comparLC[(x, y)])
        pre_ref_merge(initial_merge_LC[(x, y)], x, y)

    for x, y, z in itertools.permutations(genomes, 3):
        refr_merge_HC[(x, y)] = ref_merge(initial_merge_HC[(x, z)], initial_merge_HC[(y, z)], x, y, z)


        refr_merge_LC[(x, y)] = ref_merge(initial_merge_LC[(x, z)], initial_merge_LC[(y, z)], x, y, z)


    for x, y, z in itertools.combinations(genomes, 3):
        triplets_HC = sixway(refr_merge_HC[(x, y)], refr_merge_HC[(x, z)], refr_merge_HC[(y, z)], x, y, z)
        with open("triplets_HC_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_HC.to_csv(out, sep="\t")
        histogram(triplets_HC, 'HC','All Categories',x,y,z)

        triplets_LC = sixway(refr_merge_LC[(x, y)], refr_merge_LC[(x, z)], refr_merge_LC[(y, z)], x, y, z)
        with open("triplets_LC_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_LC.to_csv(out, sep="\t")
        histogram(triplets_LC, 'LC', 'All Categories',x,y,z)

        triplets_all = pd.concat([triplets_HC, triplets_LC],join='outer')
        with open("triplets_all_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_all.to_csv(out, sep="\t")
        histogram(triplets_all, 'all', 'All Categories',x,y,z)

    for x, y, z in itertools.combinations(genomes, 3):
        triplets_HC_eq = allmatches(triplets_HC, x, y, z)
        with open("triplets_HC_eq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_HC_eq.to_csv(out, sep="\t")
        histogram(triplets_HC_eq, 'HC', 'Exact Matches',x,y,z)

        triplets_LC_eq = allmatches(triplets_LC, x, y, z)
        with open("triplets_LC_eq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_LC_eq.to_csv(out, sep="\t")
        histogram(triplets_LC_eq, 'LC', 'Exact Matches',x,y,z)

        triplets_all_eq = allmatches(triplets_all, x, y, z)
        with open("triplets_all_eq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_all_eq.to_csv(out, sep="\t")
        histogram(triplets_all_eq, 'all', 'Exact Matches',x,y,z)

    for x, y, z in itertools.combinations(genomes, 3):
        triplets_HC_noneq = nomatches(triplets_HC, x, y, z)
        with open("triplets_HC_noneq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_HC_noneq.to_csv(out, sep="\t")
        histogram(triplets_HC_noneq, 'HC', 'Non-Matches',x,y,z)

        triplets_LC_noneq = nomatches(triplets_LC, x, y, z)
        with open("triplets_LC_noneq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_LC_noneq.to_csv(out, sep="\t")
        histogram(triplets_LC_noneq, 'LC', 'Non-Matches',x,y,z)

        triplets_all_noneq = nomatches(triplets_all, x, y, z)
        with open("triplets_all_noneq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            triplets_all_noneq.to_csv(out, sep="\t")
        histogram(triplets_all_noneq, 'all', 'Non-Matches',x,y,z)

    for x, y, z in itertools.combinations(genomes, 3):
        rawexistl = pd.read_csv('wheat.homeolog_groups.release.nonTE.TRIADS.tsv', sep='\t')
        existl = rawexistl[rawexistl['chrs'].astype(str).str.contains('{}'.format(args.A[:1]))]
        crosscheck_HC = crossch(existl, triplets_HC, x, y, z)
        with open("crosscheck_HC_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            crosscheck_HC.to_csv(out, sep="\t")
        with open("percent_HC_{}-{}-{}.txt".format(x, y, z), "wt") as out:
            for i in [existl, triplets_HC]:
                percent_HC = percent(crosscheck_HC,i)
                if i is existl:
                    out.write('{}-{}-{} proportion in consortium list ('.format(x,y,z) + str(len(existl)) + '): ')
                    out.write(str(percent_HC) + '% (' + str(len(crosscheck_HC)) + ')\n')
                else:
                    out.write('{}-{}-{} proportion in generated triplets ('.format(x,y,z) + str(len(triplets_HC))+'): ')
                    out.write(str(percent_HC) + '% (' + str(len(crosscheck_HC)) + ')\n')

        crosscheck_HC_eq = crossch(existl, triplets_HC_eq, x, y, z)
        with open("crosscheck_HC_eq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
            crosscheck_HC_eq.to_csv(out, sep="\t")
        with open("percent_HCeq_{}-{}-{}.txt".format(x, y, z), "wt") as out:
            for i in [existl, triplets_HC_eq]:
                percent_HC_eq = percent(crosscheck_HC_eq,i)
                if i is existl:
                    out.write('{}-{}-{} proportion in consortium list ('.format(x,y,z) + str(len(existl)) + '): ')
                    out.write(str(percent_HC_eq) + '% (' + str(len(crosscheck_HC_eq)) + ')\n')
                else:
                    out.write('{}-{}-{} proportion in generated triplets ('.format(x,y,z) + str(len(triplets_HC_eq))+'): ')
                    out.write(str(percent_HC_eq) + '% (' + str(len(crosscheck_HC_eq)) + ')\n')

        #crosscheck_LC = crossch(existl, triplets_LC, x, y, z)
        #with open("crosscheck_LC_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
        #    crosscheck_LC.to_csv(out, sep="\t")
        #with open("percent_LC_{}-{}-{}.txt".format(x, y, z), "wt") as out:
        #    for i in [existl, triplets_LC]:
        #        percent_LC = percent(crosscheck_LC,i)
        #        if i is existl:
        #            out.write('{}-{}-{} Against existing list: '.format(x,y,z))
        #            out.write(str(percent_LC)+'%\n')
        #        else:
        #            out.write('{}-{}-{} Against generated triplets: '.format(x,y,z))
        #            out.write(str(percent_LC) + '%\n')

        #crosscheck_all = crossch(existl, triplets_all, x, y, z)
        #with open("crosscheck_all_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
        #    crosscheck_all.to_csv(out, sep="\t")
        #with open("percent_all_{}-{}-{}.txt".format(x, y, z), "wt") as out:
        #    for i in [existl, triplets_all]:
        #        percent_all = percent(crosscheck_all,i)
        #        if i is existl:
        #            out.write('{}-{}-{} Against existing list: '.format(x,y,z))
        #            out.write(str(percent_all)+'%\n')
        #        else:
        #            out.write('{}-{}-{} Against generated triplets: '.format(x,y,z))
        #            out.write(str(percent_all) + '%\n')

        #crosscheck_all_eq = crossch(existl, triplets_all_eq, x, y, z)
        #with open("crosscheck_all_eq_{}-{}-{}.tsv".format(x,y,z), "wt") as out:
        #    crosscheck_all_eq.to_csv(out, sep="\t")
        #with open("percent_alleq_{}-{}-{}.txt".format(x, y, z), "wt") as out:
        #    for i in [existl, triplets_all_eq]:
        #        percent_all_eq = percent(crosscheck_all_eq,i)
        #        if i is existl:
        #            out.write('{}-{}-{} Against existing list: '.format(x,y,z))
        #            out.write(str(percent_all_eq)+'%\n')
        #        else:
        #            out.write('{}-{}-{} Against generated triplets: '.format(x,y,z))
        #            out.write(str(percent_all_eq) + '%\n')






# If the script is called as a script (instead of being imported as a library), execute main()
if __name__ == "__main__":
    main()
