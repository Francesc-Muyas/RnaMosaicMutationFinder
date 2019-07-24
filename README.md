# RnaMosaicMutationFinder
Somatic variant calling tools for RNAseq data

RnaMosaicMutationFinder provides a workflow and individual custom scripts to perform somatic variant calling in RNA-seq data. 
As observed in the workflow figure, the workflow is divided in 4 parts:

* Step 1. Single sample somatic calling. 
* Step 2. Re-genotyping variant sites in all samples.
* Step 3. 3D-Matrix: multi-tissue, multi-individual
* Step 4. 3D-Matrix: Filtering
