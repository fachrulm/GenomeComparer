import sys
import glob
import re
import itertools
import pandas as pd
import numpy as np

CHROMOSOMES = list(range(1, 8))
SUBGENOMES = ['A', 'B', 'D']

chroms = []
chromset = []
stringy = []
for chrom in CHROMOSOMES:
    alist = list(map(''.join, itertools.chain(itertools.product([str(chrom)], SUBGENOMES))))
    chromset.append(alist)
    for _ in chromset:
        stringy.append(" ".join(str(elm) for elm in _))
    for _ in alist:
        chroms.append(_)

comb = {}
for chrom in CHROMOSOMES:
    comb[chrom] = []
    chromss = list(map(''.join, itertools.chain(itertools.product([str(chrom)], SUBGENOMES))))
    #print(chromss)
    for x, y in itertools.permutations(chromss, 2):
        comb[chrom].append((x, y))

allcomb = [comb[x] for x in comb]

TRIPLETS = None

rule all:
    input:
        #chroms = chroms,
        #chromset = chromset,
        #stringy = stringy,
        #CHROMOSOMES = CHROMOSOMES,
        #SUBGENOMES = SUBGENOMES,
        fasta = expand("chr{chrom}{sub}.fasta", chrom=CHROMOSOMES, sub=SUBGENOMES),
        index = expand("Index/chr{chrom}{sub}/chr{chrom}{sub}.salcpchilddc", chrom=CHROMOSOMES, sub=SUBGENOMES),
        gff3 = expand("chr{chrom}{sub}.reference.gff3", chrom=CHROMOSOMES, sub=SUBGENOMES),
        cdna = expand("chr{chrom}{sub}.cdnas.fa", chrom=CHROMOSOMES, sub=SUBGENOMES),
        cds = expand("chr{chrom}{sub}.cds.fa", chrom=CHROMOSOMES, sub=SUBGENOMES),
        prot = expand("chr{chrom}{sub}.proteins.fa", chrom=CHROMOSOMES, sub=SUBGENOMES),
        ali = expand("{chrom}{sub_query}_on_{chrom}{sub_target}.gff3", chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        cpare = expand("{chrom}{sub_query}_on_{chrom}{sub_target}.compare.stats", chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        stats_cm=expand("{chrom}{sub_query}_on_{chrom}{sub_target}.stats.tsv",chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        stats_cx=expand("{chrom}{sub_query}_on_{chrom}{sub_target}.stats.txt",chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        stats_rf=expand("chr{chrom}{sub}.reference.stats.tsv",chrom=CHROMOSOMES, sub=SUBGENOMES),
        stats_rx=expand("chr{chrom}{sub}.reference.stats.txt",chrom=CHROMOSOMES, sub=SUBGENOMES),
        triplets_hc = expand("triplets_HC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_lc = expand("triplets_LC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_all = expand("triplets_all_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_hc_eq = expand("triplets_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_lc_eq = expand("triplets_LC_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_all_eq = expand("triplets_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        crossc_hc = expand("crosscheck_HC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        #crossc_lc = expand("crosscheck_LC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        #crossc_all = expand("crosscheck_all_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        crossc_hc_eq = expand("crosscheck_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        perc = expand("percent_HC_{chrom}A-{chrom}B-{chrom}D.txt", chrom=CHROMOSOMES),
        #crossc_all_eq = expand("crosscheck_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        #tripin = expand("{chrom}A {chrom}B {chrom}D", chrom=CHROMOSOMES)
    output: touch("all.done")

##Retrieve & index the chromosomes and extract GFF3S
rule start:
    input: genome="../references/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
    output:
        fasta="chr{chrom}{sub}.fasta"
    params:
        chroms=lambda wildcards: "chr%s" % wildcards.chrom,
	subg=lambda wildcards: wildcards.sub
    log: "{chrom}{sub}.log"         
    shell: """set +u && source samtools-1.3_dm && samtools faidx {input.genome} {params.chroms}{params.subg} > {output.fasta} && set -u"""

rule infas:
    input: fasta = rules.start.output.fasta
    output:
        index="Index/chr{chrom}{sub}/chr{chrom}{sub}.salcpchilddc"
    params:
        chroms=lambda wildcards: wildcards.chrom,
        subg=lambda wildcards: wildcards.sub
    log: "{chrom}{sub}.log"
    shell: """set +u && source samtools-1.3_dm && source gmap-20170508 && mkdir -p Index && gmap_build -D $(pwd)/Index -d chr{params.chroms}{params.subg} {input.fasta} > {log} 2> {log} && set -u"""

rule ref_gff3:
    input:   
        annot="../references/iwgsc_refseqv1.0_UTR_2017May05.gff3"
    output:
        gff3="chr{chrom}{sub}.reference.gff3"
    params:
        chroms=lambda wildcards: "chr%s" % wildcards.chrom,
	subg=lambda wildcards: wildcards.sub
    log:
        "{chrom}{sub}.ref_gff3.log"
    shell:"egrep \"^(#|{params.chroms}{params.subg})\" {input.annot} | uniq > {output.gff3}"

##Retrieve FASTAs
rule fastas:
    input:
        fasta = rules.start.output.fasta,
        gff3 = rules.ref_gff3.output.gff3
    output:
        cdna="chr{chrom}{sub}.cdnas.fa",
        cds="chr{chrom}{sub}.cds.fa",
        prot="chr{chrom}{sub}.proteins.fa"
    shell:
       "" "set +u && source cufflinks-2.2.2_beta_20150904 && gffread -g {input.fasta} -w {output.cdna} -x {output.cds} -y {output.prot} {input.gff3} && set -u"""


##Perform Alignments
rule gmap:
    input:
        cds="chr{chrom}{sub_query}.cds.fa",
        index="Index/chr{chrom}{sub_target}/chr{chrom}{sub_target}.salcpchilddc"
    output:
        ali="{chrom}{sub_query}_on_{chrom}{sub_target}.gff3"
    params:
        chroms=lambda wildcards: wildcards.chrom,
	subg=lambda wildcards: wildcards.sub_target
    log: "logs/{chrom}{sub_query}_on_{chrom}{sub_target}.gmap.log"
    shell: """set +u && source gmap-20170508 && gmap -F -D $(pwd)/Index -d chr{params.chroms}{params.subg} --max-intronlength-middle=71000 --max-intronlength-ends=71000 --no-chimeras -t 10 -f gff3_gene -n 1 {input.cds} > {output.ali} 2> {log} && set -u"""


###Get initial pairwise comparison files (*.compare.*)
rule filter_gmap:
    input:
      ali=rules.gmap.output.ali
    output:
      ali="{chrom}{sub_query}_on_{chrom}{sub_target}.no_cds.gff3"
    log: "logs/{chrom}{sub_query}_on_{chrom}{sub_target}.filter_gmap.log"
    shell: """egrep -v "(CDS|UTR)\W" {input.ali} > {output.ali} 2> {log}"""

rule index_mikado:
    input:
      gff3="chr{chrom}{sub}.reference.gff3"
    output:
      midx="chr{chrom}{sub}.reference.gff3.midx"
    log: "logs/chr{chrom}{sub}.mikado_index.log"
    shell: """set +u && source mikado-1.0.1 && mikado compare -r {input.gff3} --index -l {log} && set -u"""

rule mikado:
    input:
        ali=rules.filter_gmap.output.ali,
        gff3="chr{chrom}{sub_target}.reference.gff3",
        midx="chr{chrom}{sub_target}.reference.gff3.midx"
    output:
        cpare="{chrom}{sub_query}_on_{chrom}{sub_target}.compare.stats",
        tmap="{chrom}{sub_query}_on_{chrom}{sub_target}.compare.tmap",
        refmap="{chrom}{sub_query}_on_{chrom}{sub_target}.compare.refmap"
    params:
        cpare=lambda wildcards: "{chrom}{sub_query}_on_{chrom}{sub_target}.compare".format(chrom=wildcards.chrom, sub_query=wildcards.sub_query, sub_target=wildcards.sub_target)
    log: "{chrom}{sub_query}_on_{chrom}{sub_target}.compare.log"
    shell: """set +u && source mikado-1.0.1 && mikado compare -eu -r {input.gff3} -p {input.ali} --log {log} -o {params.cpare} && set -u"""

rule statsdo:
    input:
        ali = rules.filter_gmap.output.ali
    output:
        stats_cm="{chrom}{sub_query}_on_{chrom}{sub_target}.stats.tsv",
        stats_cx="{chrom}{sub_query}_on_{chrom}{sub_target}.stats.txt"
    shell: "set +u && source mikado-1.0.1 && mikado util stats --tab-stats {output.stats_cm} {input.ali} {output.stats_cx} && set -u"

rule ref_stats:
    input:
       gff3 = rules.ref_gff3.output.gff3
    output:
       stats_rf="chr{chrom}{sub}.reference.stats.tsv",
       stats_rx="chr{chrom}{sub}.reference.stats.txt"
    shell: "set +u && source mikado-1.0.1 && mikado util stats --tab-stats {output.stats_rf} {input.gff3} {output.stats_rx} && set -u"


rule triplets:
    input:
        #tripin="{chrom}A {chrom}B {chrom}D"
    output:
        triplets_hc = "triplets_HC_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_lc = "triplets_LC_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_all = "triplets_all_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_hc_eq = "triplets_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_lc_eq = "triplets_LC_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_all_eq = "triplets_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        crossc_hc = "crosscheck_HC_{chrom}A-{chrom}B-{chrom}D.tsv",
        #crossc_lc = "crosscheck_LC_{chrom}A-{chrom}B-{chrom}D.tsv",
        #crossc_all = "crosscheck_all_{chrom}A-{chrom}B-{chrom}D.tsv",
        crossc_hc_eq = "crosscheck_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        #crossc_all_eq = "crosscheck_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        perc = "percent_HC_{chrom}A-{chrom}B-{chrom}D.txt",
    params: chroms=lambda wildcards: "%s" % wildcards.chrom
    shell:
        "set +u && source python-3.5.1 && source pandas-0.18.0 && python3 main.py {params.chroms}A {params.chroms}B {params.chroms}D && set -u"

rule calper:
    output:
        finperchc = "percentage_perchrom_HC.txt",
        finperchceq = "percentage_perchrom_HCequal.txt",
    shell: "set +u && cat percent_HC_* > {output.finperchc} && cat percent_HCeq_* > {output.finperchceq} && set -u"

rule totper:
    input:
        finperchc = "percentage_perchrom_HC.txt",
        finperchceq = "percentage_perchrom_HCequal.txt"
    #output:
    shell: """
        set +u &&
        contrip=$(tail -n +2 wheat.homeolog_groups.release.nonTE.TRIADS.tsv| wc -l) &&
        crosshc=$(tail -n +2 crosscheck_HC_[1-7]A-[1-7]B-[1-7]D.tsv|wc -l) &&
        crosshc_eq=$(tail -n +2 crosscheck_HC_eq_[1-7]A-[1-7]B-[1-7]D.tsv|wc -l) &&
        triplhc=$(tail -n +2 triplets_HC_[1-7]A-[1-7]B-[1-7]D.tsv|wc -l) &&
        triplhc_eq=$(tail -n +2 triplets_HC_eq_[1-7]A-[1-7]B-[1-7]D.tsv|wc -l) &&
        echo ""Overal percentage against consortium list: $(echo "scale=2; $crosshc_eq*100/$contrip" | bc)%"" >> {input.finperchceq} &&
        echo ""$crosshc_eq out of $contrip"" >> {input.finperchceq} &&
        echo ""Overal percentage against generated triplets: $(echo "scale=2; $crosshc_eq*100/$triplhc_eq" | bc)%"" >> {input.finperchceq} &&
        echo ""$crosshc_eq out of $triplhc_eq"" >> {input.finperchceq} &&
        echo ""Overal percentage against consoritum list: $(echo "scale=2; $crosshc*100/$contrip" | bc)%"" >> {input.finperchc} &&
        echo ""$crosshc out of $contrip"" >> {input.finperchc} &&
        echo ""Overal percentage against generated triplets: $(echo "scale=2; $crosshc*100/$triplhc" | bc)%"" >> {input.finperchc} &&
        echo ""$crosshc out of $triplhc"" >> {input.finperchc} && set -u"""
