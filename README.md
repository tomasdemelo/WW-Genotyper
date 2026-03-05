Thank you for downloading these scripts.

If you have any questions, reach out to me at tomas.demelo@ontariotechu.ca.

While these scripts were originally used for genotyping Enterovirus and Norovirus GII reads it can easily be adapted for other viral targets.

Additionally, it utilizes BLAST+ executables, which can be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

1. The first step is obtaining reference sequences of known genotypes from NCBI Virus or similar repositories, and the sequences names must ONLY be the genotype.
    a. Before using this genotyping framework, some manual work is required (options are being explored to automate this). The user must standardize the names of genotypes so that they follow a consistent format.

2. Before proceeding to raw reads, the database must be validated (using the WW_Genotyper-v2_Validation.rmd). This R script requires partial sequences of the desired target. Similarly, these sequence names must also be standardized to match the genotype names in the consensus sequence database. 

4. Once that is done, a consensus sequence of each genotype can be generated using the FASTA-seq_condenser.sh which will act as the backbone for the remainder of the pipeline.

5. After this, ensure the generated consensus sequence is in the same path as the raw reads and raw read pipeline (Raw-read_pipeline.sh) after which the raw read pipeline can be run.
    a. The raw read pipeline also utilizes a fasta file with common host/undesired sequences (cross-contam.fasta) to be removed from the reads.

6. Finally, once the raw reads are processed, the final step can commence (WW_Genotyper-v4). This script will process the reads and using BLAST estimate the genotype based on the database from step 1. In addition to estimating the relative abundance of all genotypes consensus sequences for the observed genotypes are also generated for subsequence phylogenetic assessment -- currently work is being done to perform this within the WW_Genotyper and the code will be updated when that is available.
7. 

