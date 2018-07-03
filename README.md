# GenomeComparer
This tool was designed to compare the subgenomes of the hexaploid wheat in terms of the gene structure; in other words, it performs a six-way comparison based on pairwise comparisons between the three subgenomes and reports which genes are conserved between the three.


The main Python-based parsing script __(main.py)__ was mainly built based on [Pandas](https://pandas.pydata.org/) and [argparse](https://github.com/python/cpython/blob/3.7/Lib/argparse.py).
The overall pipeline was built using [SnakeMake](https://snakemake.readthedocs.io/en/stable/), and utilized [samtools](http://www.htslib.org/), [GMAP](http://research-pub.gene.com/gmap/), [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks), and [Mikado](https://github.com/lucventurini/mikado).

