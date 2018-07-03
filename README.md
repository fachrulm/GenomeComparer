# GenomeComparer
This tool was designed to compare the subgenomes of the hexaploid wheat in terms of the gene structure; in other words, it performs a six-way comparison based on pairwise comparisons between the three subgenomes and reports which genes are conserved between the three.


The main Python-based parsing script __(main.py)__ was mainly built based on [Pandas](https://pandas.pydata.org/) and [argparse](https://github.com/python/cpython/blob/3.7/Lib/argparse.py).
The overall pipeline __(Pipeline.py)__ was built using [SnakeMake](https://snakemake.readthedocs.io/en/stable/), and utilized [samtools](http://www.htslib.org/), [GMAP](http://research-pub.gene.com/gmap/), [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks), and [Mikado](https://github.com/lucventurini/mikado).


Only __Pipeline.py__ and __main.py__ are necessary, _tests.py_ was created for testing __main.py__.

To use the pipeline: `snakemake -s Pipeline.py`, given that both __Pipeline.py__ and __main.py__ are in the same directory.
