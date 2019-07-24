# RnaMosaicMutationFinder
Somatic variant calling tools for RNAseq data

RnaMosaicMutationFinder provides a workflow and individual custom scripts to perform somatic variant calling in RNA-seq data. 
As observed in the workflow figure, the workflow is divided in 4 parts:

* Step 1. Single sample somatic variant calling. 
* Step 2. Re-genotyping variant sites in all samples.
* Step 3. 3D-Matrix: multi-tissue, multi-individual
* Step 4. 3D-Matrix: Filtering

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
