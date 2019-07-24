# RnaMosaicMutationFinder
Somatic variant calling tools for RNAseq data

RnaMosaicMutationFinder provides a workflow and individual custom scripts to perform somatic variant calling in RNA-seq data. 
As observed in the workflow figure, the workflow is divided in 4 parts:

* Step 1. Single sample somatic variant calling. 
* Step 2. Re-genotyping variant sites in all samples.
* Step 3. 3D-Matrix: multi-tissue, multi-individual
* Step 4. 3D-Matrix: Filtering

![Workflow](https://github.com/Francesc-Muyas/RnaMosaicMutationFinder/blob/master/pictures/Workflow_github.pdf)

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
Raw BAM files are post-processed in order to remove alignment artefacts. PCR duplicates are marked using Picard (version 2.10.1) and reads mapping to different exons are split using SplitNCigar (part of GATK 3.7 package). Furthermore, low quality reads are removed with the python tool `bam_quality_filter.py`. `bam_quality_filter.py` requires:

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