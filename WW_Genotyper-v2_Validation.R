# install.packages(c("BiocManager","dplyr", "openxlsx", "ggplot2"))
# BiocManager::install("Biostrings")

library(tidyr)
library(dplyr)
library(openxlsx)
library(Biostrings)
library(ggplot2)
library(patchwork)

output_dir <- "validation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Colour palette
cArray <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70",
  "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown", "black", "forestgreen", "red2", "orange", "cornflowerblue",
  "magenta", "darkolivegreen4", "indianred1", "tan4", "darkblue",
  "mediumorchid1", "firebrick4", "yellowgreen", "lightsalmon", "tan3",
  "tan1", "darkgray", "wheat4", "#DDAD4B", "chartreuse",
  "seagreen1", "moccasin", "mediumvioletred", "seagreen", "cadetblue1",
  "darkolivegreen1", "tan2", "tomato3", "#7CE3D8", "gainsboro", "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70",
  "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown", "black", "forestgreen", "red2", "orange", "cornflowerblue",
  "magenta", "darkolivegreen4", "indianred1", "tan4", "darkblue",
  "mediumorchid1", "firebrick4", "yellowgreen", "lightsalmon", "tan3",
  "tan1", "darkgray", "wheat4", "#DDAD4B", "chartreuse",
  "seagreen1", "moccasin", "mediumvioletred", "seagreen", "cadetblue1",
  "darkolivegreen1", "tan2", "tomato3", "#7CE3D8", "gainsboro", "orange", "yellow", "red", "maroon"
)

# Paths to BLAST
blast_bin <- "/Users/tomasdemelo/Desktop/BLAST/ncbi-blast-2.16.0+/bin"
makeblastdb <- file.path(blast_bin, "makeblastdb")
blastn <- file.path(blast_bin, "blastn")

remove_names <- c(
  "EV-sp.", "EV-A", "EV-B", "EV-C", "EV-D",
  "RV-A", "RV-B", "RV-C", "Coxsackievirus"
)

genotype_mapping <- c(
  # EV-A
  "CV-A2" = "EV-A", "CV-A3" = "EV-A", "CV-A4" = "EV-A", "CV-A5" = "EV-A", "CV-A6" = "EV-A", "CV-A7" = "EV-A", "CV-A8" = "EV-A",
  "CV-A10" = "EV-A", "CV-A12" = "EV-A", "CV-A14" = "EV-A", "CV-A16" = "EV-A",
  "EV-A71" = "EV-A", "EV-A76" = "EV-A", "EV-A89" = "EV-A", "EV-A90" = "EV-A", "EV-A91" = "EV-A", "EV-A92" = "EV-A",
  "EV-A114" = "EV-A", "EV-A119" = "EV-A", "EV-A120" = "EV-A", "EV-A121" = "EV-A",
  "SV19" = "EV-A", "SV43" = "EV-A", "SV46" = "EV-A", "BabEV-A13" = "EV-A",

  # EV-B
  "CV-B1" = "EV-B", "CV-B2" = "EV-B", "CV-B3" = "EV-B", "CV-B4" = "EV-B", "CV-B5" = "EV-B", "CV-B6" = "EV-B", "CV-A9" = "EV-B",
  "Echovirus-E1" = "EV-B", "Echovirus-E2" = "EV-B", "Echovirus-E3" = "EV-B", "Echovirus-E4" = "EV-B", "Echovirus-E5" = "EV-B",
  "Echovirus-E6" = "EV-B", "Echovirus-E7" = "EV-B", "Echovirus-E9" = "EV-B",
  "Echovirus-E11" = "EV-B", "Echovirus-E12" = "EV-B", "Echovirus-E13" = "EV-B", "Echovirus-E14" = "EV-B", "Echovirus-E15" = "EV-B",
  "Echovirus-E16" = "EV-B", "Echovirus-E17" = "EV-B", "Echovirus-E18" = "EV-B", "Echovirus-E19" = "EV-B", "Echovirus-E20" = "EV-B", "Echovirus-E21" = "EV-B",
  "Echovirus-E24" = "EV-B", "Echovirus-E25" = "EV-B", "Echovirus-E26" = "EV-B",
  "Echovirus-E27" = "EV-B", "Echovirus-E29" = "EV-B", "Echovirus-E30" = "EV-B", "Echovirus-E31" = "EV-B", "Echovirus-E32" = "EV-B", "Echovirus-E33" = "EV-B",
  "EV-B69" = "EV-B", "EV-B73" = "EV-B", "EV-B74" = "EV-B", "EV-B75" = "EV-B", "EV-B77" = "EV-B", "EV-B78" = "EV-B",
  "EV-B79" = "EV-B", "EV-B80" = "EV-B", "EV-B81" = "EV-B", "EV-B82" = "EV-B", "EV-B83" = "EV-B", "EV-B84" = "EV-B",
  "EV-B85" = "EV-B", "EV-B86" = "EV-B", "EV-B87" = "EV-B", "EV-B88" = "EV-B", "EV-B93" = "EV-B", "EV-B97" = "EV-B",
  "EV-B98" = "EV-B", "EV-B100" = "EV-B", "EV-B101" = "EV-B", "EV-B106" = "EV-B", "EV-B107" = "EV-B", "EV-B110" = "EV-B",
  "EV-B111" = "EV-B", "EV-B112" = "EV-B", "EV-B113" = "EV-B", "SA5" = "EV-B",

  # EV-C
  "CV-A1" = "EV-C", "CV-A11" = "EV-C", "CV-A13" = "EV-C", "CV-A17" = "EV-C", "CV-A19" = "EV-C", "CV-A20" = "EV-C", "CV-A21" = "EV-C",
  "CV-A22" = "EV-C", "CV-A24" = "EV-C",
  "EV-C95" = "EV-C", "EV-C96" = "EV-C", "EV-C99" = "EV-C", "EV-C102" = "EV-C",
  "EV-C104" = "EV-C", "EV-C105" = "EV-C", "EV-C109" = "EV-C", "EV-C113" = "EV-C", "EV-C119" = "EV-C",
  "EV-C116" = "EV-C", "EV-C117" = "EV-C", "EV-C118" = "EV-C", "Poliovirus-1" = "EV-C", "Poliovirus-2" = "EV-C", "Poliovirus-3" = "EV-C",

  # EV-D
  "EV-D68" = "EV-D", "EV-D70" = "EV-D", "EV-D94" = "EV-D", "EV-D111" = "EV-D", "EV-D120" = "EV-D",

  # RV-A
  "RV-A1" = "RV-A", "RV-A1B" = "RV-A", "RV-A2" = "RV-A", "RV-A7" = "RV-A", "RV-A8" = "RV-A", "RV-A9" = "RV-A", "RV-A10" = "RV-A",
  "RV-A11" = "RV-A", "RV-A12" = "RV-A", "RV-A13" = "RV-A", "RV-A15" = "RV-A", "RV-A16" = "RV-A", "RV-A18" = "RV-A", "RV-A19" = "RV-A",
  "RV-A20" = "RV-A", "RV-A21" = "RV-A", "RV-A22" = "RV-A", "RV-A23" = "RV-A", "RV-A24" = "RV-A", "RV-A25" = "RV-A", "RV-A28" = "RV-A",
  "RV-A29" = "RV-A", "RV-A30" = "RV-A", "RV-A31" = "RV-A", "RV-A32" = "RV-A", "RV-A33" = "RV-A", "RV-A34" = "RV-A", "RV-A36" = "RV-A",
  "RV-A38" = "RV-A", "RV-A39" = "RV-A", "RV-A40" = "RV-A", "RV-A41" = "RV-A", "RV-A43" = "RV-A", "RV-A45" = "RV-A", "RV-A46" = "RV-A",
  "RV-A47" = "RV-A", "RV-A49" = "RV-A", "RV-A50" = "RV-A", "RV-A51" = "RV-A", "RV-A53" = "RV-A", "RV-A54" = "RV-A", "RV-A55" = "RV-A",
  "RV-A56" = "RV-A", "RV-A57" = "RV-A", "RV-A58" = "RV-A", "RV-A59" = "RV-A", "RV-A89" = "RV-A", "RV-A73" = "RV-A", "RV-A75" = "RV-A", "RV-A82" = "RV-A",
  "RV-A106" = "RV-A", "RV-A107" = "RV-A", "RV-A108" = "RV-A", "RV-A94" = "RV-A", "RV-A80" = "RV-A", "RV-A77" = "RV-A",

  # RV-B
  "RV-B3" = "RV-B", "RV-B4" = "RV-B", "RV-B5" = "RV-B", "RV-B6" = "RV-B", "RV-B7" = "RV-B", "RV-B8" = "RV-B",
  "RV-B9" = "RV-B", "RV-B14" = "RV-B", "RV-B15" = "RV-B", "RV-B16" = "RV-B", "RV-B17" = "RV-B", "RV-B20" = "RV-B",
  "RV-B26" = "RV-B", "RV-B27" = "RV-B", "RV-B28" = "RV-B", "RV-B29" = "RV-B", "RV-B30" = "RV-B", "RV-B52" = "RV-B", "RV-B42" = "RV-B",

  # RV-C
  "RV-C1" = "RV-C", "RV-C2" = "RV-C", "RV-C3" = "RV-C", "RV-C4" = "RV-C", "RV-C5" = "RV-C", "RV-C6" = "RV-C",
  "RV-C7" = "RV-C", "RV-C8" = "RV-C", "RV-C9" = "RV-C", "RV-C10" = "RV-C", "RV-C11" = "RV-C", "RV-C12" = "RV-C",
  "RV-C13" = "RV-C", "RV-C14" = "RV-C", "RV-C15" = "RV-C", "RV-C16" = "RV-C", "RV-C17" = "RV-C", "RV-C18" = "RV-C",
  "RV-C19" = "RV-C", "RV-C20" = "RV-C", "RV-C21" = "RV-C", "RV-C22" = "RV-C", "RV-C36" = "RV-C", "RV-C40" = "RV-C", "RV-C41" = "RV-C",
  "RV-C45" = "RV-C", "RV-C46" = "RV-C", "RV-C26" = "RV-C"
)

run_validation_pipeline <- function(virus_type, reads_path, db_path, db_name) {
  cat("\n============================================\n")
  cat("Running Validation Pipeline for", virus_type, "\n")
  cat("============================================\n")

  reads <- readDNAStringSet(reads_path)
  db_seqs <- readDNAStringSet(db_path)

  # 1. Filter out reads that don't map to the reference sequences
  reads <- reads[names(reads) %in% names(db_seqs)]

  # 2. Filter out generic names
  WW_reads_filtered <- reads[!names(reads) %in% remove_names]

  seq_names <- names(WW_reads_filtered)
  uniqueID <- sprintf("ID%09d", seq_along(seq_names))
  names(WW_reads_filtered) <- paste0(uniqueID, "_", seq_names)

  temp_dir <- tempdir()
  clean_fasta <- file.path(temp_dir, paste0("clean_", virus_type, ".fasta"))
  writeXStringSet(WW_reads_filtered, clean_fasta)

  noisy_fastq <- file.path(temp_dir, paste0("noisy_reads_", virus_type, ".fastq"))
  noisy_fasta <- file.path(temp_dir, paste0("noisy_reads_", virus_type, ".fasta"))

  cat("Simulating noisy reads with Badread (this may take a while)...\n")
  cmd_badread <- sprintf("/opt/homebrew/Caskroom/miniforge/base/envs/ww-genotyper-sim/bin/badread simulate --reference \"%s\" --quantity 100x --error_model nanopore2023 --qscore_model nanopore2023 --identity 99,100,1 > \"%s\"", clean_fasta, noisy_fastq)
  system(cmd_badread, ignore.stdout = FALSE, ignore.stderr = FALSE)

  if (!file.exists(noisy_fastq) || file.info(noisy_fastq)$size == 0) {
    stop(paste("Error: badread failed to generate noisy reads for", virus_type))
  }

  cat("Converting noisy FASTQ to FASTA and restoring headers...\n")
  awk_script <- "awk 'NR%4==1 {match($0, /ID[0-9]+_[^,]+/); if (RSTART>0) {name = substr($0, RSTART, RLENGTH)} else {name = \"\"}} NR%4==2 {if (name != \"\") {print \">\" name; print $0}}'"
  cmd_convert <- sprintf("%s '%s' > '%s'", awk_script, noisy_fastq, noisy_fasta)
  system(cmd_convert)

  WW_reads <- noisy_fasta

  if (!file.exists(paste0(db_name, ".nin"))) {
    cat("Creating BLAST database from reference file...\n")
    system(paste(shQuote(makeblastdb), "-in", db_path, "-dbtype nucl -out", db_name))
  } else {
    cat("BLAST database already exists. Skipping creation.\n")
  }

  cat("Running BLAST for", WW_reads, "...\n")
  out_file <- sub("\\.fasta$", "_blast_results.txt", WW_reads)
  system(paste("blastn -query", WW_reads, "-db", db_name, "-outfmt 6 -out ", out_file, "-task dc-megablast -num_threads 8 -evalue 1e-3 -word_size 11"))

  if (file.info(out_file)$size == 0) {
    warning("No hits found for ", WW_reads)
    results <- data.frame(query = character(), subject = character(), pident = numeric(), length = numeric(), mismatch = numeric(), gapopen = numeric(), qstart = numeric(), qend = numeric(), sstart = numeric(), send = numeric(), evalue = numeric(), bitscore = numeric(), stringsAsFactors = FALSE)
  } else {
    results <- read.table(out_file,
      header = FALSE, sep = "\t",
      col.names = c("query", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    )
  }

  results$source_file <- WW_reads
  uniqueID_blast <- sprintf("ID%09d", 1:nrow(results))
  results$uniqueID <- uniqueID_blast

  if (nrow(results) > 0) {
    top_hits <- results %>%
      group_by(query) %>%
      top_n(1, bitscore) %>%
      ungroup()
  } else {
    top_hits <- data.frame()
  }

  if (nrow(top_hits) > 0) {
    top_hits$query <- substring(top_hits$query, 13)
    top_hits <- top_hits %>% filter(!query %in% remove_names)
    top_hits <- top_hits %>% filter(!is.na(subject), !is.na(query))

    # Species fallback logic for hits < 80% identity
    for (i in seq_len(nrow(top_hits))) {
      if (!is.na(top_hits$pident[i]) && top_hits$pident[i] < 80) {
        if (top_hits$subject[i] %in% names(genotype_mapping)) {
          top_hits$subject[i] <- unname(genotype_mapping[top_hits$subject[i]])
        }
        if (top_hits$query[i] %in% names(genotype_mapping)) {
          top_hits$query[i] <- unname(genotype_mapping[top_hits$query[i]])
        }
      }
    }

    top_hits$is_right <- (top_hits$query == top_hits$subject)
    rightcalls <- top_hits %>%
      filter(is_right) %>%
      select(Actual = query, call = subject, PercentIdentity = pident, ID = uniqueID)
    wrongcalls <- top_hits %>%
      filter(!is_right) %>%
      select(Actual = query, call = subject, PercentIdentity = pident, ID = uniqueID)
  } else {
    rightcalls <- data.frame(Actual = character(), call = character(), PercentIdentity = numeric(), ID = character())
    wrongcalls <- data.frame(Actual = character(), call = character(), PercentIdentity = numeric(), ID = character())
  }

  total_calls <- nrow(rightcalls) + nrow(wrongcalls)
  percentright <- if (total_calls > 0) (nrow(rightcalls) / total_calls) * 100 else 0

  right_df <- if (nrow(rightcalls) > 0) rightcalls %>% mutate(Type = "Right", PIdent = round(as.numeric(PercentIdentity), 0)) else data.frame(Actual = character(), Type = character(), PIdent = numeric(), stringsAsFactors = FALSE)
  wrong_df <- if (nrow(wrongcalls) > 0) wrongcalls %>% mutate(Type = "Wrong", PIdent = round(as.numeric(PercentIdentity), 0)) else data.frame(Actual = character(), Type = character(), PIdent = numeric(), stringsAsFactors = FALSE)

  stacked_data <- rbind(right_df, wrong_df)

  if (nrow(wrongcalls) > 0) {
    relative_abundance <- wrongcalls %>%
      group_by(Actual, call) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Actual) %>%
      mutate(prop = (Count / sum(Count)) * 100)
  } else {
    relative_abundance <- data.frame(Actual = character(), call = character(), Count = numeric(), prop = numeric())
  }

  if (nrow(stacked_data) > 0) {
    stacked_percent <- stacked_data %>%
      group_by(Actual, Type) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Actual) %>%
      mutate(Percent = (Count / sum(Count)) * 100)
    stacked_data_PIdent <- stacked_data %>%
      group_by(Type, PIdent) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Type) %>%
      mutate(Percent = (Count / total_calls) * 100)
  } else {
    stacked_percent <- data.frame(Actual = character(), Type = character(), Count = numeric(), Percent = numeric())
    stacked_data_PIdent <- data.frame(Type = character(), PIdent = numeric(), Count = numeric(), Percent = numeric())
  }

  percentright_rounded <- round(percentright, 1)
  pdf_path <- file.path(output_dir, paste0(virus_type, "_BLAST-Genotyping_Validations.pdf"))
  pdf(pdf_path, width = 15, height = 10)

  if (nrow(stacked_data_PIdent) > 0) {
    p_histo <- ggplot(stacked_data_PIdent, aes(x = PIdent, y = Percent, fill = Type)) +
      geom_col(color = "black") +
      scale_fill_manual(values = c("Right" = "forestgreen", "Wrong" = "firebrick"), labels = c("Right" = paste0("Right (", percentright_rounded, "%)"), "Wrong" = paste0("Wrong (", round(100 - percentright_rounded, 1), "%)"))) +
      labs(title = paste("BLAST Hit Percent Identity by Outcome (Noisy Reads) -", virus_type), x = "Percent Identity (%)", y = "Percentage of Total Calls (%)", fill = "Outcome") +
      theme_minimal() +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_histo)
  }

  if (nrow(stacked_percent) > 0) {
    p_stacked <- ggplot(stacked_percent, aes(x = Actual, y = Percent, fill = Type)) +
      geom_col(color = "black") +
      scale_fill_manual(values = c("Right" = "forestgreen", "Wrong" = "firebrick")) +
      labs(title = paste("Accuracy of Genotype Calling for Each Variant (Noisy Reads) -", virus_type), x = "Designated Genotype", y = "Percentage of Calls (%)", fill = "Outcome") +
      theme_minimal() +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_stacked)
  }

  if (nrow(relative_abundance) > 0) {
    p_abund <- ggplot(relative_abundance, aes(x = Actual, y = prop, fill = call)) +
      geom_col(color = "black") +
      scale_fill_manual(values = cArray) +
      labs(title = paste("Expected vs. Detected Relative Abundance (Noisy Reads) -", virus_type), x = "Engineered Abundance Profile", y = "Relative Abundance (%)", fill = "Called Virus Type") +
      theme_minimal() +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6))
    print(p_abund)
  }
  dev.off()

  cat("Generated plots for", virus_type, "at", pdf_path, "\n")
  return(list(
    top_hits = top_hits,
    stacked_percent = stacked_percent,
    relative_abundance = relative_abundance,
    percentright = percentright
  ))
}

# Run for EV
ev_hits <- run_validation_pipeline(
  virus_type = "EV",
  reads_path = "Validation/EV-VP1-seq_for_validation.fasta",
  db_path = "EV/Barrie/all_EV_types_combined.fasta",
  db_name = "EV_Reference_Sequences"
)

# Run for NV
nv_hits <- run_validation_pipeline(
  virus_type = "NV",
  reads_path = "Validation/NV-ORF1-ORF2-seq_for_validation.fasta",
  db_path = "NV/Barrie/v3-all_NV_types_combined.fasta",
  db_name = "NV-v3_Reference_Sequences"
)

# Write Excel
excel_file <- file.path(output_dir, "blast_top_hits_by_sample.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "EV Validation")
writeData(wb, sheet = "EV Validation", ev_hits$top_hits)
addWorksheet(wb, "NV Validation")
writeData(wb, sheet = "NV Validation", nv_hits$top_hits)
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat("Excel workbook written to", excel_file, "\n")




output_dir <- "validation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Setup & Parameters
Sys.setenv(PATH = paste0("/opt/homebrew/Caskroom/miniforge/base/envs/ww-genotyper-sim/bin:", Sys.getenv("PATH")))
badread_bin <- "/opt/homebrew/Caskroom/miniforge/base/envs/ww-genotyper-sim/bin/badread"
depth <- 10000
intervals <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)
part2_dir <- "validation/seqs"
dir.create(part2_dir, showWarnings = FALSE, recursive = TRUE)

viruses <- list(
  "EV" = list(
    fasta = "Complex_FASTAs/Reference_CV-A16_VP1_NC_001612.1_nucleotide.fasta",
    gff = "Complex_FASTAs/NC_001612.1.gff",
    cds_start = 751,
    cds_end = 7332,
    g_len = 7413
  ),
  "NV" = list(
    fasta = "Complex_FASTAs/Reference_Norovirus_GII_NC_044932.1_nucleotide.fasta",
    gff = "Complex_FASTAs/NC_044932.1.gff",
    cds_start = 5,
    cds_end = 5095,
    g_len = 7525
  )
)

# Synthetic Generation Function
generate_synthetic_reads <- function(virus_id, v_params, var_frequency) {
  wt_fasta <- v_params$fasta
  out_var_fasta <- file.path(part2_dir, sprintf("%s_variant_%.2f.fasta", virus_id, var_frequency))
  out_mixed_fastq <- file.path(part2_dir, sprintf("%s_mixed_%.2f.fastq", virus_id, var_frequency))

  ref_seqs <- readDNAStringSet(wt_fasta)
  seq_name <- names(ref_seqs)[1]
  seq_char <- as.character(ref_seqs[[1]])
  seq_vec <- unlist(strsplit(seq_char, ""))

  cds_seq <- seq_vec[v_params$cds_start:v_params$cds_end]
  codon_table <- GENETIC_CODE
  bases <- c("A", "C", "G", "T")

  syn_target <- 5
  nonsyn_target <- 10
  syn_muts <- 0
  nonsyn_muts <- 0
  attempts <- 0

  set.seed(1867) # For reproducibility

  while ((syn_muts < syn_target || nonsyn_muts < nonsyn_target) && attempts < 10000) {
    attempts <- attempts + 1
    codon_idx <- sample(1:(length(cds_seq) / 3), 1)
    pos <- (codon_idx - 1) * 3 + 1

    codon <- paste(cds_seq[pos:(pos + 2)], collapse = "")
    if (!(codon %in% names(codon_table))) next
    orig_aa <- codon_table[[codon]]
    if (orig_aa == "*" || is.na(orig_aa)) next

    mut_pos <- sample(0:2, 1)
    orig_base <- substr(codon, mut_pos + 1, mut_pos + 1)
    alt_bases <- setdiff(bases, orig_base)
    mut_base <- sample(alt_bases, 1)

    new_codon <- codon
    substr(new_codon, mut_pos + 1, mut_pos + 1) <- mut_base
    if (!(new_codon %in% names(codon_table))) next
    new_aa <- codon_table[[new_codon]]
    if (new_aa == "*" || is.na(new_aa)) next

    is_syn <- (orig_aa == new_aa)
    if (is_syn && syn_muts < syn_target) {
      cds_seq[pos + mut_pos] <- mut_base
      syn_muts <- syn_muts + 1
    } else if (!is_syn && nonsyn_muts < nonsyn_target) {
      cds_seq[pos + mut_pos] <- mut_base
      nonsyn_muts <- nonsyn_muts + 1
    }
  }

  mutated_seq_char <- seq_char
  substr(mutated_seq_char, v_params$cds_start, v_params$cds_end) <- paste(cds_seq, collapse = "")
  var_dnas <- DNAStringSet(mutated_seq_char)
  names(var_dnas) <- paste0(seq_name, "_variant")
  writeXStringSet(var_dnas, out_var_fasta)

  # Simulate reads using badread
  wt_abundance <- 1.0 - var_frequency
  wt_reads <- as.integer(depth * wt_abundance)
  var_reads <- as.integer(depth * var_frequency)

  wt_out <- file.path(part2_dir, "wt_sim_tmp.fastq")
  var_out <- file.path(part2_dir, "var_sim_tmp.fastq")

  cmd_wt <- sprintf("%s simulate --reference '%s' --quantity %dx --error_model nanopore2023 --qscore_model nanopore2023 --identity 99,100,1 --chimeras 1 > '%s'", badread_bin, wt_fasta, wt_reads, wt_out)
  cmd_var <- sprintf("%s simulate --reference '%s' --quantity %dx --error_model nanopore2023 --qscore_model nanopore2023 --identity 99,100,1 --chimeras 1 > '%s'", badread_bin, out_var_fasta, var_reads, var_out)

  system(cmd_wt)
  system(cmd_var)

  cmd_combine <- sprintf("cat '%s' '%s' > '%s'", wt_out, var_out, out_mixed_fastq)
  system(cmd_combine)

  unlink(wt_out)
  unlink(var_out)

  return(out_mixed_fastq)
}

# Variant Calling Function
call_variants <- function(virus_id, v_params, var_frequency, mixed_fastq) {
  sorted_bam <- file.path(part2_dir, sprintf("%s_sorted_%.2f.bam", virus_id, var_frequency))
  tsv_prefix <- file.path(part2_dir, sprintf("%s_variants_%.2f", virus_id, var_frequency))

  cmd_map <- sprintf("minimap2 -ax map-ont '%s' '%s' | samtools sort -o '%s'", v_params$fasta, mixed_fastq, sorted_bam)
  system(cmd_map, ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(sprintf("samtools index '%s'", sorted_bam))

  cmd_ivar <- sprintf(
    "samtools mpileup -aa -A -d 0 -B -Q 0 --reference '%s' '%s' | ivar variants -p '%s' -q 0 -t 0.05 -r '%s' -g '%s'",
    v_params$fasta, sorted_bam, tsv_prefix, v_params$fasta, v_params$gff
  )
  system(cmd_ivar, ignore.stdout = TRUE, ignore.stderr = TRUE)

  return(paste0(tsv_prefix, ".tsv"))
}

# Orchestration Loop
run_validation <- function() {
  for (virus_id in names(viruses)) {
    v_params <- viruses[[virus_id]]
    for (val in intervals) {
      cat(sprintf("Processing %s at %.2f frequency...\n", virus_id, val))
      mixed_fastq <- generate_synthetic_reads(virus_id, v_params, val)
      call_variants(virus_id, v_params, val, mixed_fastq)
    }
  }
}

# Run the simulation and mapping loop
run_validation()

# Metrics Calculation & Plotting
results <- data.frame()

for (v in names(viruses)) {
  v_params <- viruses[[v]]
  for (val in intervals) {
    file <- file.path(part2_dir, sprintf("%s_variants_%.2f.tsv", v, val))
    if (!file.exists(file)) next

    df <- read.delim(file, stringsAsFactors = FALSE)
    df <- df %>% filter(PASS == TRUE)

    translated_vars <- df %>% filter(REF_AA != "NA" & ALT_AA != "NA")
    recovery_count <- nrow(translated_vars)
    recovery_rate <- (recovery_count / 15) * 100

    pnps_df <- translated_vars %>%
      mutate(
        var_type = case_when(
          REF_AA == ALT_AA ~ "synonymous",
          REF_AA != ALT_AA ~ "nonsynonymous",
          TRUE ~ "other"
        )
      ) %>%
      filter(var_type %in% c("synonymous", "nonsynonymous")) %>%
      summarise(
        n_S = sum(var_type == "synonymous"),
        n_NS = sum(var_type == "nonsynonymous"),
        pN_pS = ifelse(n_S > 0, n_NS / n_S, NA_real_)
      )

    pi_df <- df %>%
      mutate(site_pi = 2 * ALT_FREQ * (1 - ALT_FREQ)) %>%
      summarise(
        mean_pi = sum(site_pi, na.rm = TRUE) / v_params$g_len
      )

    results <- rbind(results, data.frame(
      Virus = v,
      Frequency = as.numeric(val) * 100,
      Recovery_Rate = recovery_rate,
      pN_pS = pnps_df$pN_pS,
      Pi = pi_df$mean_pi
    ))
  }
}

if (nrow(results) > 0) {
  p1 <- ggplot(results, aes(x = Frequency, y = Recovery_Rate, color = Virus, group = Virus)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
    theme_bw() +
    labs(
      title = "Mutation Recovery Rate vs Variant Frequency",
      x = "Engineered Variant Frequency (%)",
      y = "Recovery Rate (%)"
    ) +
    scale_y_continuous(limits = c(0, 105))

  p2 <- ggplot(results, aes(x = Frequency, y = pN_pS, color = Virus, group = Virus)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "red") +
    theme_bw() +
    labs(
      title = "Calculated pN/pS vs Variant Frequency",
      subtitle = "Dashed line = Ground Truth Raw Ratio (2.0)",
      x = "Engineered Variant Frequency (%)",
      y = "Calculated pN/pS Ratio"
    )

  p3 <- ggplot(results, aes(x = Frequency, y = Pi, color = Virus, group = Virus)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    theme_bw() +
    labs(
      title = "Calculated Nucleotide Diversity (\u03c0) vs Variant Frequency",
      x = "Engineered Variant Frequency (%)",
      y = "Calculated \u03c0 Diversity"
    )

  combined_plot <- (p1 / p2 / p3) + plot_annotation(tag_levels = 'A')
  ggsave(file.path(part2_dir, "WW_Genotyper_Validation_Combined.tiff"), plot = combined_plot, width = 10, height = 15, dpi = 300, compression = "lzw")

  pdf(file.path(part2_dir, "WW_Genotyper_Validation_Plots.pdf"), width = 10, height = 15)
  print(combined_plot)
  dev.off()

  cat(sprintf("Plots saved to %s\n", file.path(part2_dir, "WW_Genotyper_Validation_Combined.tiff")))
} else {
  cat("No TSV files found to plot. Uncomment `run_validation()` to generate the data.\n")
}

# Write Summary Excel covering ALL modules
summary_file <- file.path("validation", "validation_summary_metrics.xlsx")
wb_sum <- createWorkbook()

# Genotyper Module - EV Summary
addWorksheet(wb_sum, "EV_Overall")
writeData(wb_sum, "EV_Overall", data.frame(Total_Reads = sum(ev_hits$stacked_percent$Count), Overall_Accuracy = ev_hits$percentright))
addWorksheet(wb_sum, "EV_Accuracy_Per_Genotype")
writeData(wb_sum, "EV_Accuracy_Per_Genotype", ev_hits$stacked_percent)
addWorksheet(wb_sum, "EV_Miscalls")
writeData(wb_sum, "EV_Miscalls", ev_hits$relative_abundance)

# Genotyper Module - NV Summary
addWorksheet(wb_sum, "NV_Overall")
writeData(wb_sum, "NV_Overall", data.frame(Total_Reads = sum(nv_hits$stacked_percent$Count), Overall_Accuracy = nv_hits$percentright))
addWorksheet(wb_sum, "NV_Accuracy_Per_Genotype")
writeData(wb_sum, "NV_Accuracy_Per_Genotype", nv_hits$stacked_percent)
addWorksheet(wb_sum, "NV_Miscalls")
writeData(wb_sum, "NV_Miscalls", nv_hits$relative_abundance)

# Mutation Validation Module
if (nrow(results) > 0) {
  addWorksheet(wb_sum, "Mutation_Validation")
  writeData(wb_sum, "Mutation_Validation", results)
}

saveWorkbook(wb_sum, summary_file, overwrite = TRUE)
cat("\nSummary Excel workbook covering all modules written to", summary_file, "\n")

