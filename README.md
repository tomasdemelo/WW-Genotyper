# WW-Genotyper

A comprehensive framework for genotyping Enterovirus (EV) and Norovirus (NV) GII from wastewater sequences. Version 11 introduces integrated iVar mutation calling, ESM-2 structural predictions for single nucleotide variants, absolute abundance normalization, and robust longitudinal statistical modeling.

If you have any questions, reach out at tomas.demelo@ontariotechu.ca.

---

## Installation & Dependencies

To test the full functionality of the pipeline on your own machine, ensure the following dependencies are installed:

### 1. R and Required Packages
The core genotyping and plotting scripts are built in R. You will need:
- R (>= 4.1.0)
- The following R packages:
  ```R
  install.packages(c("tidyverse", "furrr", "future", "openxlsx", "ggplot2", "grid", "pheatmap"))
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("ComplexHeatmap", "Biostrings"))
  ```
- **reticulate**: Required to bridge R to Python for the ESM-2 triage.
  ```R
  install.packages("reticulate")
  ```

### 2. Python & ESM-2 Model
The script `score_mutation.py` uses the ESM-2 structural model to score nonsynonymous variants.
- Python 3.9+ is required.
- Install the required Python libraries using pip:
  ```bash
  pip install torch transformers pandas biopython
  ```
*Note: Make sure your `reticulate` in R is configured to use the python environment where these are installed.*

### 3. BLAST+ Executables
The framework utilizes BLAST+ for aligning sequences to reference databases.
- Download from: [NCBI BLAST+ FTP](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- Ensure `blastn` and `makeblastdb` are available in your system's `PATH`.

---

## Directory Structure

This repository provides all the essential pieces:
- `WW_Genotyper-v11.R`: The main analysis script for Mac/Linux environments.
- `WW_Genotyper-v11_Windows.R`: Modified paths/execution tailored for Windows/WSL environments.
- `score_mutation.py`: The ESM-2 triage integration.
- `bash_scripts/`: Contains the preprocessing bash scripts (`FASTA-seq_condenser.sh`, `new-new_EVminimap2.sh`, `new-new_NVminimap2.sh`) required to process raw FASTQs into `_combined.fasta` sequences.
- `references/`: Essential databases, including reference sequences for both EV and NV.
- `example_data/`: Contains a test dataset (`EV/Ajax/` and `EV/Barrie/`) to quickly verify the pipeline.
- `RT-qPCR_data.xlsx`: A template for formatting your longitudinal viral load tracking.

---

## Running a Test Drive

We have included an `example_data` directory to test the pipeline right out of the box and verify all your dependencies are correctly set up.

To test the functionality locally:

1. Clone the repository.
2. Open `WW_Genotyper-v11.R` in RStudio or run it directly.
3. Set your `workspace_root` variable to the directory where you cloned this repository.
4. Set `tsv_dirs <- c("example_data/EV/Ajax", "example_data/EV/Barrie")`.
5. Ensure `RT-qPCR_data.xlsx` is correctly located in the root (or update its path in the script).
6. Run the script!

The output should automatically generate a PDF (`Output/WW_Genotyper_Summary_Report_v11.pdf`) showcasing the relative/absolute abundance, regression models, Tajima's D, pN/pS tracking, and more!
