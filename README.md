# RnaMosaicMutationFinder
Somatic variant calling tools for RNAseq data

RnaMosaicMutationFinder provides a workflow and individual custom scripts to perform somatic variant calling in RNA-seq data. 
As observed in the workflow figure, the workflow is divided in 4 parts:

![Workflow](https://github.com/Francesc-Muyas/RnaMosaicMutationFinder/blob/master/pictures/Workflow_github.png)

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

Then, the resulting file is transformed to a Variant Calling File (vcf) using python script `tsv2vcf.py`:

```
usage: tsv2vcf.py [-h] -i INFILE [-tID TUMORID] [-nID NORMALID] -ref REFERENCE
                  -o OUTFILE [-cov MIN_COV] [-ac MIN_AC]
                  [-variant_dist MIN_DIST] [-dup1 DUPLICATE1_FILTER]
                  [-dup2 DUPLICATE2_FILTER] [-str {0,1}] [-af MIN_AF]
                  [-end {0,1}] [-som {single,paired}]

Getting barcodes in fastq file and labels

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Tsv table
  -tID TUMORID, --tumorid TUMORID
                        Tumor sample id
  -nID NORMALID, --normalid NORMALID
                        Normal sample id
  -ref REFERENCE, --reference REFERENCE
                        Reference fastq file which table was build
  -o OUTFILE, --outfile OUTFILE
                        Vcf output file
  -cov MIN_COV, --min_COV MIN_COV
                        Minimum Coverage
  -ac MIN_AC, --min_AC MIN_AC
                        Minimum reads supporting alternative allele
  -variant_dist MIN_DIST, --min_DIST MIN_DIST
                        Minimum distance allowed between variants (to avoid
                        clustered errors)
  -dup1 DUPLICATE1_FILTER, --duplicate1_FILTER DUPLICATE1_FILTER
                        Minimum (mean) number of duplicates 1
  -dup2 DUPLICATE2_FILTER, --duplicate2_FILTER DUPLICATE2_FILTER
                        Minimum (mean) number of duplicates 2
  -str {0,1}, --strand {0,1}
                        Strand bias test (Fisher test). 0 for turn it off
  -af MIN_AF, --min_AF MIN_AF
                        Minimum allele frequency allowed
  -end {0,1}, --end_read_filter {0,1}
                        Filtering for variants found at end or begining of the
                        read. 0 for turn it off. Default: Activated
  -som {single,paired}, --somatic_type {single,paired}
                        If analysis is based on two paired samples (paired) or
                        based on only one sample (single). Default: paired
```

Once we get the vcf for all samples, we collapse all PASS variant sites in one zero-based bed file (i.e using bedtools merge), which will be necessary for Step 2.


* Step 2. Re-genotyping variant sites in all samples.

Once we have created the bed file with all variant sites, if possible, we should remove those sites which are known to be germline variants from the analysis. 

At this point, we take the previous `Final.bam` files of all samples (see Workflow) and re-calculate the mpileup file using the variant bed file with the option `-l variant.bed`. Again, we transform this pileup files to tsv files using the `pileup2tsv.py` script. The resulting tsv is then used as input of the R script ` tsv2genotype.r`, which requires:

```
usage: tsv2genotype.r [-h] -i file -o file -cov COV

optional arguments:
  -h, --help            show this help message and exit
  -i file, --var_file file
                        input read count file
  -o file, --output_file file
                        output_file
  -cov COV, --cov COV   Minimum coverage
```

!! Output files shoul be called with next pattern: `Individual.Tissue.tsv`, where Individual is the id used for the individual, and Tissue is the analysed tissue without dots (i.e Brain, Skin_sun_exposed...)


*Step 3. 3D-Matrix: multi-tissue, multi-individual

In this step we collect all tsv per individual obtained in Step 2 to create a multi-tissue per individual matrix. It is performed with the python script ` tsv2matrix.py`, which requires (tsv files should be sorted by genomic coordinates):

```
usage: tsv2matrix.py [-h] --input I --headerlines H --individual ID --output O

optional arguments:
  -h, --help            show this help message and exit
  --input I, -i I       input files
  --headerlines H, -hd H
                        header line number
  --individual ID, -id ID
                        Id of the individual
  --output O, -o O      output file: Matrix per sample with all tissue
                        information
```

Once matrixes per individual are created, they must be merged taking into account the column order (R function `rbind.fill` from `plyr` package)


