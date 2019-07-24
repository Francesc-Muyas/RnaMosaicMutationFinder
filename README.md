# RnaMosaicMutationFinder
Somatic variant calling tools for RNAseq data

RnaMosaicMutationFinder provides a workflow and individual custom scripts to perform somatic variant calling in RNA-seq data. 
As observed in the workflow figure, the workflow is divided in 4 parts:

* Step 1. Single sample somatic variant calling. 
* Step 2. Re-genotyping variant sites in all samples.
* Step 3. 3D-Matrix: multi-tissue, multi-individual
* Step 4. 3D-Matrix: Filtering

![Workflow](https://github.com/Francesc-Muyas/RnaMosaicMutationFinder/blob/master/pictures/Workflow_github.png)

## Get RnaMosaicMutationFinder tools  
You will need to run `git clone ` to get RnaMosaicMutationFinder tools. 

```
git clone https://github.com/Francesc-Muyas/RnaMosaicMutationFinder
```

## Install dependencies
Firstly, you will need to install some software

- [Python 2.7](https://www.python.org/download/releases/2.7/), with next python modules:
    * argparse
    * numpy
    * pysam
    * time   
    * sys
    * subprocess
    * os
    * bgzf (from Bio)
    * re
    * ntpath
    * struct
    * hashlib
    * collections
    * itertools

- R version 3.3.2, with next packages and respective dependencies:   
    * argparse
    * bbmle
    * data.table
    * GenomicRanges (Bioconductor)
    * ggplot2
    * randomForest
    * survcomp (Bioconductor)
    * VGAM

## Documentation and usage

* Step 1. Single sample somatic calling. 

Raw BAM files are post-processed in order to remove alignment artefacts. PCR duplicates are marked using Picard (version 2.10.1) and reads mapping to different exons are split using SplitNCigar (part of GATK 3.7 package). Furthermore, low quality reads are removed with the python tool `bam_quality_filter.py`. `python bam_quality_filter.py` requires:

```
usage: bam_quality_filter.py [-h] --infile INFILE --outfile OUTFILE
                             [--minMQ MINMQ] [--maxMM MAXMM] [--maxGAP MAXGAP]

optional arguments:
  -h, --help         show this help message and exit
  --infile INFILE    Input BAM file.
  --outfile OUTFILE  Output BAM file.
  --minMQ MINMQ      Minimum required mapping quality of aligned read. Default
                     = 30
  --maxMM MAXMM      Maximum number of mismatches of aligned read . Default =
                     4
  --maxGAP MAXGAP    Maximum number of GAPs . Default = 1
```

Finally, in order to avoid systematic alignment errors at the extremes of the reads, we suggest trimming the first and last 4 bases from each read-end or read-breakpoint (BamUtil version 1.0.14).

Mpileup file is obtained using Samtools mpileup (version 1.3.1). Pileup files are transformed to tsv files using `pileup2tsv.py` python script, which requires:

```
usage: pileup2tsv.py [-h] --pileup PILEUP -o VARIANTF [--minBQ MINBQ]
                     [-m MINDEPTH]


optional arguments:
  -h, --help            show this help message and exit
  --pileup PILEUP       Pileup file.
  -o VARIANTF, --variantf VARIANTF
                        Output file for predicted SNPs and indels.
  --minBQ MINBQ         Minimum base quality to consider nucleotide in SNP
                        analysis. Default = 10
  -m MINDEPTH, --mindepth MINDEPTH
                        Minimum minimum total depth
```

Once pileup files are transformed to read counts in tsv files, we model errors using a Beta Binomial distribution where:

    Error alternative counts ~ Bin(Coverage, error rate)

    Error rate ~ Beta(alpha, beta)

As the error rate differs depending on the nucleotide change, we modeled error distributions independently for each possible nucleotide change (A>C, A>T, A>G, C>A, C>T, C>G). Finally, we identified all sites showing alternative allele counts significantly deviating from the ER distribution after FDR correction.

This step is performed by the R script `tsv_betabinomial.r`. `Rscript tsv_betabinomial.r` requires:

```
usage: tsv_betabinomial.r [-h] -t file [-n NORMAL_FILE] [-tr TRAIN_FILE] -o
                          file [-b BED]

optional arguments:
  -h, --help            show this help message and exit
  -t file, --tumor_file file
                        Tumor - Read count input file
  -n NORMAL_FILE, --normal_file NORMAL_FILE
                        Normal - Read count input file (not mandatory)
  -tr TRAIN_FILE, --train_file TRAIN_FILE
                        File to trian error rate distribution [optional]
  -o file, --out_file file
                        output_file
  -b BED, --bed BED     Bed file with positions to ignore for modelling (not mandatory)
```

