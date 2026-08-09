# Define global configuration options

# Pipeline Execution Toggles
# Set these to FALSE to skip specific, time-consuming parts of the pipeline
RUN_PHASE1_GENOTYPER <- TRUE
RUN_PHASE2_IVAR <- TRUE
RUN_PHASE3_COMBINED <- TRUE
BUILD_PDF_REPORT <- TRUE
MIN_ALT_FREQ_POLYMORPHIC <- 0.05 # Minimum variant frequency for a polymorphic site
MIN_INDEL_FREQ <- 0.50 # Minimum variant frequency for indels (Nanopore artifact filter)
MIN_DEPTH <- 400 # Minimum read depth cutoff (can be relaxed to e.g., 50 for low-coverage runs)
WINDOW_SIZE <- 50 # Sliding window size for diversity/Tajima's D calculation (bp)
STEP_SIZE <- 10 # Step size for sliding windows (bp)
STRAND_BIAS_ALPHA <- 0.05 # FDR-adjusted p-value significance threshold for Fisher's Exact strand bias
STRAND_BIAS_FILT_FREQ <- 0.20 # Cutoff below which strand-biased variants are filtered out
RUN_PARALLEL <- FALSE # Set to TRUE to process sample suites in parallel (faster, uses more RAM). Set to FALSE to process sequentially.
RUN_COLABFOLD_LOCAL <- TRUE # Set to TRUE to fold structures locally via ColabFold. Set to FALSE to skip folding and use pre-folded models instead.
MIN_PLDDT_STRUCTURAL <- 70 # Minimum mean pLDDT to trust structural analysis (RMSD/ipTM flagging)

# Paths
workspace_root <- "/mnt/c/Users/Gru/Desktop/Tomas/Windows_Folding/Windows_Transfer"
ivar_workspace <- file.path(workspace_root, "iVar")
blast_bin <- file.path(workspace_root, "BLAST", "ncbi-blast-2.16.0+", "bin")
CHIMERAX_BIN <- "/mnt/c/Program Files/ChimeraX 1.12/bin/ChimeraX.exe"
COLABFOLD_BIN <- "/home/wastewater/colabfold-env/bin/colabfold_batch --msa-mode mmseqs2_uniref_env --num-models 1 --num-recycle 3 --max-seq 128 --max-extra-seq 512"

# Setup Python and plotting tools to run without a display
Sys.setenv(KERAS_HOME = file.path(workspace_root, ".keras"))
Sys.setenv(MPLCONFIGDIR = file.path(workspace_root, ".matplotlib"))
Sys.setenv(PYTHONPATH = file.path(workspace_root, "tmp"))

# Manage GPU memory for ColabFold to prevent crashes
# Stop JAX from taking all GPU memory at once
# This prevents crashes caused by hidden memory limits
Sys.setenv(XLA_PYTHON_CLIENT_PREALLOCATE = "false") # Allocate GPU memory on-demand instead of upfront
Sys.setenv(XLA_PYTHON_CLIENT_MEM_FRACTION = "0.80") # Cap JAX to 80% of VRAM (leaves headroom for driver)
Sys.setenv(TF_FORCE_GPU_ALLOW_GROWTH = "true") # TensorFlow memory growth (ColabFold imports TF)
Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0") # Disable oneDNN custom ops (eliminates numerical warnings)

# Terminal Output Auto-Logging
terminal_log_file <- file.path(workspace_root, "Output", paste0("WW_Genotyper_Run_Log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
if (!dir.exists(dirname(terminal_log_file))) dir.create(dirname(terminal_log_file), recursive = TRUE, showWarnings = FALSE)
sink(terminal_log_file, split = TRUE)

# Global Error/Warning Logging Setup
error_log_file <- file.path(workspace_root, "pipeline_errors.log")
if (file.exists(error_log_file)) file.remove(error_log_file)

log_warning_error <- function(...) {
  msg <- paste(..., sep = " ")
  cat(msg)
  msg_clean <- gsub("\n$", "", msg)
  cat(msg_clean, "\n", file = error_log_file, append = TRUE)
}

# Local caching directories
cds_cache_dir <- file.path(workspace_root, "ncbi_cds_cache")
alignments_cache_dir <- file.path(workspace_root, "ncbi_alignments_cache")
if (!dir.exists(cds_cache_dir)) {
  dir.create(cds_cache_dir, showWarnings = FALSE, recursive = TRUE)
}
if (!dir.exists(alignments_cache_dir)) {
  dir.create(alignments_cache_dir, showWarnings = FALSE, recursive = TRUE)
}

# Packages
suppressPackageStartupMessages({
  required_pkgs <- c(
    "Rsamtools", "dplyr", "openxlsx", "Biostrings", "ggplot2", "tidyr",
    "vegan", "tibble", "RVAideMemoire", "tidyverse", "ComplexHeatmap", "circlize",
    "jsonlite", "scales", "ape", "png", "grid", "rentrez", "future", "future.apply", "furrr",
    "prophet", "reticulate"
  )
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing missing package:", pkg, "\n")
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }

  library(Rsamtools)
  library(dplyr)
  library(openxlsx)
  library(Biostrings)
  library(prophet)
  library(reticulate)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(vegan)
  library(tibble)
  library(RVAideMemoire)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(jsonlite)
  library(scales)
  library(ape)
  library(png)
  library(grid)
  library(rentrez)
  library(future)
  library(future.apply)
  library(furrr)
})

# Initialize Python Bridge for ESM-2 Triage
if (file.exists(file.path(workspace_root, ".venv", "bin", "python"))) {
  use_python(file.path(workspace_root, ".venv", "bin", "python"), required = TRUE)
} else if (file.exists("/home/wastewater/miniconda3/bin/python")) {
  use_python("/home/wastewater/miniconda3/bin/python", required = TRUE)
} else {
  use_python("/opt/homebrew/Caskroom/miniforge/base/bin/python3", required = TRUE)
}
tryCatch(
  {
    source_python(file.path(workspace_root, "score_mutation.py"))
    ESM_READY <- TRUE
  },
  error = function(e) {
    log_warning_error("Warning: Could not load score_mutation.py or Python environment. ESM-2 triage will be skipped.\n")
    log_warning_error("  -> Detailed Error:", e$message, "\n")
    ESM_READY <- FALSE
  }
)
# Setup Parallel Backend
num_workers <- max(1, future::availableCores() - 2)
if (Sys.info()[["sysname"]] %in% c("Windows")) {
  plan(multisession, workers = num_workers)
} else if (Sys.info()[["sysname"]] == "Darwin") {
  # macOS background workers hang during parallel graphics rendering (tiff). Run sequentially.
  plan(sequential)
} else {
  plan(multisession, workers = num_workers)
}

# Wait between NCBI requests to avoid getting blocked
# Increase the wait time if using multiple background workers
NCBI_DELAY <- max(0.4, (num_workers / 3.0) + 0.1)

# Receptor Mapping for Co-folding
receptor_map <- list(
  "CV-A16" = "Q14108",
  "CV-A6" = "Q96MU8",
  "CV-B3" = "P78310",
  "CV-B4" = "P78310",
  "CV-B5" = "P78310",
  "ECHOVIRUS-E11" = "P55899",
  "ECHOVIRUS-E24" = "P55899",
  "ECHOVIRUS-E3" = "P55899",
  "ECHOVIRUS-E7" = "P08174",
  "ECHOVIRUS-E9" = "P08174",
  "CV-A9" = "P06756",
  "CV-A1" = "P05362",
  "CV-A19" = "P05362"
)

# Load receptor sequences immediately
receptor_fasta_path <- file.path(workspace_root, "receptors.fasta")
receptor_seqs <- list()
if (file.exists(receptor_fasta_path)) {
  rec_lines <- readLines(receptor_fasta_path, warn = FALSE)
  current_id <- NULL
  current_seq <- c()
  for (line in rec_lines) {
    if (grepl("^>", line)) {
      if (!is.null(current_id)) receptor_seqs[[current_id]] <- paste(current_seq, collapse = "")
      parts <- strsplit(line, "\\\\|")[[1]]
      if (length(parts) >= 2) {
        current_id <- parts[2]
      } else {
        current_id <- sub("^>", "", line)
      }
      current_seq <- c()
    } else if (nchar(trimws(line)) > 0) {
      current_seq <- c(current_seq, trimws(line))
    }
  }
  if (!is.null(current_id)) receptor_seqs[[current_id]] <- paste(current_seq, collapse = "")
}

fetch_receptor_sequence <- function(region) {
  for (gt in names(receptor_map)) {
    if (grepl(gt, region, ignore.case = TRUE)) {
      rec_id <- receptor_map[[gt]]
      if (!is.null(receptor_seqs[[rec_id]])) {
        return(receptor_seqs[[rec_id]])
      }
    }
  }
  return(NULL)
}

# Helper Functions

# Download EV VP1 sequences from NCBI automatically
build_ev_vp1_database <- function(ev_genotypes, output_fasta, output_csv) {
  if (file.exists(output_fasta) && file.exists(output_csv)) {
    return(TRUE)
  }

  cat("\n=== Building Dynamic EV VP1 Reference Database ===\n")
  all_vp1_seqs <- list()
  mapping_df <- data.frame(REGION = character(), Accession = character(), stringsAsFactors = FALSE)

  for (gt in ev_genotypes) {
    cat("Searching NCBI for VP1 sequence for", gt, "...\n")

    # Dynamically generate the organism search name (e.g. CV-A16 -> Coxsackievirus A16)
    org_name <- gt
    if (grepl("^CV-A", gt)) {
      org_name <- gsub("CV-A", "Coxsackievirus A", gt)
    } else if (grepl("^CV-B", gt)) {
      org_name <- gsub("CV-B", "Coxsackievirus B", gt)
    } else if (grepl("^ECHOVIRUS-E", gt)) {
      org_name <- gsub("ECHOVIRUS-E", "Echovirus ", gt)
    } else if (grepl("^EV-", gt)) {
      org_name <- gsub("EV-", "Enterovirus ", gt)
    }

    search_term <- paste0('"', org_name, '"[Organism] AND VP1[Gene] AND 800:1000[Sequence Length]')

    search_res <- tryCatch(
      {
        res <- entrez_search(db = "nuccore", term = search_term, retmax = 1)
        Sys.sleep(NCBI_DELAY)
        res
      },
      error = function(e) {
        log_warning_error("  -> Error searching:", e$message, "\n")
        return(NULL)
      }
    )

    if (is.null(search_res) || length(search_res$ids) == 0) {
      cat("  -> No isolated VP1 sequences found for", gt, ".\n")
      next
    }

    seq_id <- search_res$ids[1]
    cat("  -> Found accession:", seq_id, "- Fetching FASTA...\n")

    fasta_data <- tryCatch(
      {
        res <- entrez_fetch(db = "nuccore", id = seq_id, rettype = "fasta")
        Sys.sleep(NCBI_DELAY)
        res
      },
      error = function(e) {
        log_warning_error("  -> Error fetching FASTA:", e$message, "\n")
        return(NULL)
      }
    )

    if (is.null(fasta_data)) next

    seq_lines <- strsplit(fasta_data, "\n")[[1]]
    seq_lines <- seq_lines[seq_lines != ""]
    raw_header <- seq_lines[1]
    sequence <- paste(seq_lines[-1], collapse = "")

    accession <- gsub("^>([^ ]+).*", "\\1", raw_header)

    all_vp1_seqs[[gt]] <- sequence
    mapping_df <- rbind(mapping_df, data.frame(REGION = gt, Accession = accession, stringsAsFactors = FALSE))
  }

  if (length(all_vp1_seqs) > 0) {
    if (!dir.exists(dirname(output_fasta))) {
      dir.create(dirname(output_fasta), recursive = TRUE, showWarnings = FALSE)
    }
    combined_set <- Biostrings::DNAStringSet(unlist(all_vp1_seqs))
    Biostrings::writeXStringSet(combined_set, filepath = output_fasta)
    write.csv(mapping_df, file = output_csv, row.names = FALSE)
    cat("\nSuccess! Saved", length(combined_set), "VP1 sequences to", output_fasta, "\n")
    return(TRUE)
  } else {
    cat("\nFailed to find any VP1 sequences.\n")
    return(FALSE)
  }
}

# Fetch PDB from ColabFold
run_colabfold <- function(sequence, output_dir, sample_name, colabfold_bin) {
  out_pdb <- file.path(output_dir, paste0(sample_name, ".pdb"))

  # Check if pre-folded on Windows (or other external machine)
  base_sample_name <- gsub("_complex_with_receptor$|_complex$|_VP1_only$", "", sample_name)
  windows_folded_path <- file.path(workspace_root, "Windows_Folding", "Protein_Folded", paste0(base_sample_name, ".pdb"))

  if (!file.exists(windows_folded_path)) {
    barcode <- sub("_.*", "", base_sample_name)
    acc_match <- regexpr("[A-Z0-9]+_?[0-9]+\\.[0-9]+$", base_sample_name)
    if (acc_match != -1) {
      accession <- regmatches(base_sample_name, acc_match)
      pattern <- paste0("^", barcode, "_.*", accession, ".*\\.pdb$")

      # Check both Protein_Folded and Protein_Folded_Raw
      for (dir_name in c("Protein_Folded", "Protein_Folded_Raw")) {
        folded_dir <- file.path(workspace_root, "Windows_Folding", dir_name)
        if (dir.exists(folded_dir)) {
          # Use recursive = TRUE for Protein_Folded_Raw since ColabFold creates subfolders
          matches <- list.files(folded_dir, pattern = pattern, full.names = TRUE, recursive = (dir_name == "Protein_Folded_Raw"))
          if (length(matches) > 0) {
            # Find the best match (rank_001 or rank_1 if present)
            best_match <- matches[1]
            for (m in matches) {
              if (grepl("rank_001|rank_1", m)) {
                best_match <- m
                break
              }
            }
            windows_folded_path <- best_match
            break
          }
        }
      }
    }
  }

  cat("  [DEBUG] Looking for pre-folded model at:", windows_folded_path, "\n")
  if (file.exists(windows_folded_path)) {
    file.copy(windows_folded_path, out_pdb, overwrite = TRUE)
    cat("  [ColabFold] Pre-folded model FOUND! Copied to:", out_pdb, "\n")
    # For JSON, we look for the corresponding score file if it came from Protein_Folded_Raw
    windows_json_path <- sub("\\.pdb$", ".json", windows_folded_path)
    if (grepl("Protein_Folded_Raw", windows_folded_path)) {
      windows_json_path <- sub("_unrelaxed_.*\\.pdb$", "_scores_rank_001_alphafold2_ptm_model_1_seed_000.json", windows_folded_path)
      # Fuzzy match for json if exact replacement fails
      if (!file.exists(windows_json_path)) {
        json_pattern <- paste0("^", barcode, "_.*", accession, ".*\\.json$")
        json_matches <- list.files(dirname(windows_folded_path), pattern = json_pattern, full.names = TRUE)
        if (length(json_matches) > 0) windows_json_path <- json_matches[1]
      }
    }
    if (file.exists(windows_json_path)) {
      file.copy(windows_json_path, file.path(output_dir, paste0(sample_name, ".json")), overwrite = TRUE)
    }
  }

  if (file.exists(out_pdb) && file.info(out_pdb)$size > 100) {
    return(TRUE)
  }

  # If we have specified not to run ColabFold locally, skip it
  if (!RUN_COLABFOLD_LOCAL) {
    cat("  [ColabFold] Pre-folded model not found for", sample_name, ". Skipping local folding.\n")
    return(FALSE)
  }

  # Otherwise, execute local ColabFold
  cat("  [ColabFold] Folding sequence length", nchar(sequence), "aa for", sample_name, "...\n")

  # Create a temp directory for colabfold
  cf_tmp <- file.path(output_dir, paste0("cf_", sample_name))
  if (dir.exists(cf_tmp)) {
    unlink(cf_tmp, recursive = TRUE)
  }
  dir.create(cf_tmp, showWarnings = FALSE, recursive = TRUE)

  # Write fasta (stripping trailing/internal stop codons)
  fa_file <- file.path(cf_tmp, paste0(sample_name, ".fasta"))
  clean_seq <- gsub("\\*", "", as.character(sequence))
  writeLines(c(paste0(">", sample_name), clean_seq), fa_file)

  # Run colabfold headlessly
  # Fix settings so ColabFold's Python doesn't crash
  old_preload <- Sys.getenv("LD_PRELOAD")
  if (old_preload != "") Sys.unsetenv("LD_PRELOAD")

  old_path <- Sys.getenv("PATH")
  colabfold_dir <- dirname(strsplit(colabfold_bin, " ")[[1]][1])
  Sys.setenv(PATH = paste0(colabfold_dir, ":", old_path))

  # Build the base environment for GPU limits
  gpu_env_base <- paste(
    "XLA_PYTHON_CLIENT_PREALLOCATE=false",
    "XLA_PYTHON_CLIENT_MEM_FRACTION=0.80",
    "TF_FORCE_GPU_ALLOW_GROWTH=true",
    "TF_ENABLE_ONEDNN_OPTS=0"
  )

  # Retry logic: 3-tier fallback strategy for WSL2
  max_attempts <- 3
  cf_out <- NULL
  for (attempt in seq_len(max_attempts)) {
    if (attempt > 1) {
      cat(sprintf("  [ColabFold] Attempt %d/%d — clearing GPU cache and retrying...\n", attempt, max_attempts))
      # Clear any leftover ColabFold temp files from failed attempt
      unlink(list.files(cf_tmp, pattern = "\\.(a3m|hhr|atab)$", full.names = TRUE))
      # Brief pause to let GPU driver release memory
      Sys.sleep(10)
    }

    cf_log_file <- file.path(cf_tmp, paste0("colabfold_run_attempt", attempt, ".log"))

    # Configure command based on attempt number
    if (attempt == 1) {
      cmd <- paste(gpu_env_base, colabfold_bin, shQuote(fa_file), shQuote(cf_tmp), "2>&1 | tee", shQuote(cf_log_file))
    } else if (attempt == 2) {
      # Attempt 2: Try alternative JAX allocator for WSL2
      cmd <- paste(gpu_env_base, "XLA_PYTHON_CLIENT_ALLOCATOR=platform", colabfold_bin, shQuote(fa_file), shQuote(cf_tmp), "2>&1 | tee", shQuote(cf_log_file))
    } else {
      # Attempt 3: Guaranteed CPU Fallback
      cat("  [ColabFold] WARNING: GPU inference failed twice. Falling back to CPU. This will be slow but will complete.\n")
      cmd <- paste("CUDA_VISIBLE_DEVICES=\"\" JAX_PLATFORMS=cpu", colabfold_bin, shQuote(fa_file), shQuote(cf_tmp), "2>&1 | tee", shQuote(cf_log_file))
    }

    cat("  [DEBUG] Executing ColabFold command...\n")
    system(cmd)

    if (file.exists(cf_log_file)) {
      cf_out <- readLines(cf_log_file, warn = FALSE)
    } else {
      cf_out <- ""
    }
    # Check if a PDB was produced (success) or if we need to retry
    pdb_check <- list.files(cf_tmp, pattern = "\\.pdb$", full.names = TRUE)
    if (length(pdb_check) > 0) break
    # Check for segfault in output
    if (!any(grepl("Segmentation fault|core dumped|SIGSEGV|CUDA_ERROR", cf_out, ignore.case = TRUE))) break
  }

  Sys.setenv(PATH = old_path)
  if (old_preload != "") Sys.setenv(LD_PRELOAD = old_preload)

  # Find the rank 1 unrelaxed model or fallback to any generated pdb files
  pdb_files <- list.files(cf_tmp, pattern = "rank_001.*\\.pdb$", full.names = TRUE)
  if (length(pdb_files) == 0) {
    pdb_files <- list.files(cf_tmp, pattern = "rank_1.*\\.pdb$", full.names = TRUE)
  }
  if (length(pdb_files) == 0) {
    pdb_files <- list.files(cf_tmp, pattern = "_unrelaxed_.*\\.pdb$", full.names = TRUE)
  }
  if (length(pdb_files) == 0) {
    pdb_files <- list.files(cf_tmp, pattern = "\\.pdb$", full.names = TRUE)
  }

  if (length(pdb_files) > 0) {
    file.copy(pdb_files[1], out_pdb, overwrite = TRUE)

    # Try to find the corresponding JSON scores file to save ipTM
    json_files <- list.files(cf_tmp, pattern = ".*_scores_.*\\\\.json$", full.names = TRUE)
    if (length(json_files) > 0) {
      file.copy(json_files[1], file.path(output_dir, paste0(sample_name, ".json")), overwrite = TRUE)
    }

    unlink(cf_tmp, recursive = TRUE)
    return(TRUE)
  }

  return(FALSE)
}

# Run ChimeraX headlessly
run_chimerax_alignment <- function(ref_pdb, sample_pdb, output_cxs, output_log, chimerax_bin, mutated_residues = NULL) {
  if (file.exists(output_cxs) && file.exists(output_log)) {
    return(TRUE)
  }

  # Turn off the Chromium sandbox so it doesn't crash from permission errors
  Sys.setenv(QTWEBENGINE_DISABLE_SANDBOX = "1")
  Sys.setenv(ELECTRON_DISABLE_SANDBOX = "1")

  # Create a temporary python script to run in ChimeraX
  py_file <- paste0(tempfile(), ".py")

  # Python script wrapper to ensure ChimeraX exits on error
  # Paths are enclosed in double-quotes to handle any spaces in the paths
  py_cmds <- c(
    "from chimerax.core.commands import run",
    "try:",
    "    run(session, 'windowsize 1600 1200')",
    "    run(session, 'wait 30')",
    paste0("    run(session, 'open \"", ref_pdb, "\"')"),
    paste0("    run(session, 'open \"", sample_pdb, "\"')"),
    "    res = run(session, 'matchmaker #2 to #1')",
    "    run(session, 'color bfactor')",
    "    run(session, 'view')",
    "    run(session, 'wait 60')"
  )

  # Highlight mutations if any exist
  if (!is.null(mutated_residues) && length(mutated_residues) > 0 && any(!is.na(mutated_residues))) {
    res_list <- paste(mutated_residues, collapse = ",")
    py_cmds <- c(
      py_cmds,
      paste0("    run(session, 'color #2:", res_list, " hotpink')"),
      paste0("    run(session, 'show #2:", res_list, " sidechain')")
    )
  }

  py_cmds <- c(
    py_cmds,
    paste0("    run(session, 'save \"", output_cxs, "\"')"),
    "    rmsd_val = res[0]['final RMSD']",
    paste0("    with open(\"", output_log, "\", \"w\") as f:"),
    "        f.write(f\"RMSD between pruned atom pairs is {rmsd_val:.3f} angstroms\\n\")",
    "except Exception as e:",
    "    print(f'Error occurred during ChimeraX alignment: {e}')",
    "finally:",
    "    run(session, 'quit')"
  )
  writeLines(py_cmds, py_file)

  cmd <- paste(shQuote(chimerax_bin), "--nogui --script", shQuote(py_file), "> /dev/null 2>&1")
  system(cmd)

  # Clean up temporary script
  if (file.exists(py_file)) {
    unlink(py_file)
  }

  if (file.exists(output_cxs)) {
    return(TRUE)
  }
  return(FALSE)
}

# Figure out the season (Winter, Spring, etc.) from the date
getSeason <- function(dates) {
  m <- as.numeric(format(dates, "%m"))
  dplyr::case_when(
    m %in% 3:5 ~ "Spring",
    m %in% 6:8 ~ "Summer",
    m %in% 9:11 ~ "Fall",
    TRUE ~ "Winter"
  )
}

# Reads an iVar output TSV, calculates total depth, and extracts the sample ID from the filename.
load_tsv_file <- function(file) {
  sample_name <- basename(file) |> stringr::str_remove("_variants\\.tsv$")
  readr::read_tsv(file, show_col_types = FALSE) |>
    dplyr::mutate(
      dplyr::across(dplyr::contains("POS_AA"), as.character),
      dplyr::across(c(
        POS, ALT_DP, REF_DP, REF_RV, REF_QUAL,
        ALT_RV, ALT_QUAL, ALT_FREQ, TOTAL_DP, PVAL
      ), as.numeric),
      PASS = as.logical(PASS),
      sample_id = sample_name
    ) |>
    dplyr::filter(PASS == TRUE)
}

# Match positions between our custom reference and the official NCBI reference
# by looking at how they align and handling missing pieces
build_coord_map <- function(aln_custom_str, aln_ncbi_str) {
  chars_c <- strsplit(aln_custom_str, "")[[1]]
  chars_n <- strsplit(aln_ncbi_str, "")[[1]]
  is_gap_c <- chars_c == "-"
  is_gap_n <- chars_n == "-"
  cum_n <- cumsum(!is_gap_n)
  cum_n[is_gap_n] <- NA_integer_
  mapping <- cum_n[!is_gap_c]
  return(as.integer(mapping))
}

# Get the NCBI ID from the FASTA header
extract_accession <- function(header) {
  header <- sub("^>", "", header)
  first_word <- strsplit(header, "\\s+")[[1]][1]
  if (grepl("\\|", first_word)) {
    parts <- strsplit(first_word, "\\|")[[1]]
    clean_parts <- parts[
      parts != "" &
        !grepl("^[0-9]+$", parts) &
        !tolower(parts) %in% c("gb", "emb", "dbj", "ref", "gi")
    ]
    if (length(clean_parts) > 0) {
      return(clean_parts[1])
    }
  }
  return(first_word)
}

# Calculate Tajima's D using mutation frequencies
# Calculate a score comparing actual diversity to expected diversity
compute_tajima_d <- function(freqs, n_seq) {
  freqs_seg <- freqs[freqs > 0 & freqs < 1] # Filter out fixed or absent variants
  S <- length(freqs_seg) # Count number of segregating sites
  if (S < 2 || n_seq < 3) {
    return(NA_real_)
  }

  n <- n_seq
  a1 <- sum(1 / seq_len(n - 1)) # Calculate Tajima's a1 coefficient
  a2 <- sum(1 / seq_len(n - 1)^2) # Calculate Tajima's a2 coefficient

  theta_w <- S / a1 # Watterson's estimator (expected variance)
  pi_val <- sum(2 * freqs_seg * (1 - freqs_seg)) # Nucleotide diversity (observed variance)

  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - 1 / a1
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)

  denom <- sqrt(e1 * S + e2 * S * (S - 1)) # Standard deviation of (pi - theta)
  if (denom == 0) {
    return(NA_real_)
  }
  return((pi_val - theta_w) / denom)
}

# Run statistical test to compare forward and reverse reads for mutations
# to detect sequencing strand bias artifacts at a single position.
fisher_strand_bias <- function(ref_rv, ref_dp, alt_rv, alt_dp) {
  fwd_ref <- ref_dp - ref_rv
  fwd_alt <- alt_dp - alt_rv

  if (is.na(fwd_ref) || is.na(ref_rv) || is.na(fwd_alt) || is.na(alt_rv) ||
    fwd_ref < 0 || ref_rv < 0 || fwd_alt < 0 || alt_rv < 0) {
    return(NA_real_)
  }

  # Contingency table:
  #            Forward   Reverse
  # Reference  fwd_ref   ref_rv
  # Alternate  fwd_alt   alt_rv
  mat <- matrix(c(fwd_ref, ref_rv, fwd_alt, alt_rv), nrow = 2, byrow = TRUE)

  # Perform two-tailed Fisher's exact test
  res <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
  return(res)
}

# Run the strand bias test quickly on many columns of data
# Speed up testing by saving previous math results
strand_bias_pvals <- function(ref_rv, ref_dp, alt_rv, alt_dp) {
  if (length(ref_rv) == 0) {
    return(numeric(0))
  }
  df <- data.frame(ref_rv = ref_rv, ref_dp = ref_dp, alt_rv = alt_rv, alt_dp = alt_dp)
  unique_df <- unique(df)
  unique_df$pval <- mapply(fisher_strand_bias, unique_df$ref_rv, unique_df$ref_dp, unique_df$alt_rv, unique_df$alt_dp, USE.NAMES = FALSE)

  key <- paste(df$ref_rv, df$ref_dp, df$alt_rv, df$alt_dp, sep = "_")
  unique_key <- paste(unique_df$ref_rv, unique_df$ref_dp, unique_df$alt_rv, unique_df$alt_dp, sep = "_")
  pvals <- unique_df$pval[match(key, unique_key)]
  return(pvals)
}

# Figure out if a mutation is growing, shrinking, or gone
# based on its variant allele frequency changes over chronologically ordered dates.
classify_mutation_dynamics <- function(freqs_ordered) {
  freqs_ordered <- freqs_ordered[!is.na(freqs_ordered)]
  n <- length(freqs_ordered)
  if (n == 0) {
    return(NA_character_)
  }
  present <- sum(freqs_ordered > 0.05)

  if (present / n >= 0.8) {
    return("persistent")
  }
  if (any(freqs_ordered > 0.95)) {
    return("fixed")
  }
  if (n >= 3 && all(diff(tail(freqs_ordered, 3)) > 0) &&
    freqs_ordered[1] < 0.05) {
    return("emerging")
  }
  if (any(freqs_ordered > 0.10) &&
    tail(freqs_ordered, 1) < 0.05) {
    return("lost")
  }
  if (present <= 2) {
    return("transient")
  }
  return("other")
}

# Download protein boundaries and products from NCBI
# Save files locally to avoid asking NCBI for the same data twice
get_ncbi_cds <- function(accession, cache_dir = NULL) {
  if (!is.null(cache_dir)) {
    cache_path <- file.path(cache_dir, paste0(accession, "_cds.csv"))
    if (file.exists(cache_path)) {
      cat("Loading cached CDS for", accession, "from local file...\n")
      # Return cached dataframe, ensuring correct column types
      df <- read.csv(cache_path, stringsAsFactors = FALSE)
      if (nrow(df) > 0) {
        df$start <- as.integer(df$start)
        df$end <- as.integer(df$end)
        df$product <- as.character(df$product)
      }
      return(df)
    }
  }

  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=", accession, "&rettype=ft&retmode=text")
  Sys.sleep(NCBI_DELAY)
  lines <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
  if (is.null(lines) || length(lines) == 0) {
    return(NULL)
  }

  cds_df <- data.frame(start = integer(), end = integer(), product = character(), stringsAsFactors = FALSE)

  current_start <- NA
  current_end <- NA
  current_feature <- NA

  for (line in lines) {
    if (grepl("^>Feature", line)) next

    parts <- strsplit(line, "\t")[[1]]
    val1 <- if (length(parts) >= 1) gsub("[<>]", "", parts[1]) else ""
    val2 <- if (length(parts) >= 2) gsub("[<>]", "", parts[2]) else ""
    if (length(parts) >= 3 && grepl("^[0-9]+$", val1) && grepl("^[0-9]+$", val2)) {
      current_start <- as.integer(val1)
      current_end <- as.integer(val2)
      current_feature <- parts[3]
    } else if (length(parts) >= 5 && !is.na(current_feature) && current_feature %in% c("CDS", "mat_peptide") && parts[4] == "product") {
      product <- parts[5]
      s <- min(current_start, current_end)
      e <- max(current_start, current_end)
      cds_df <- rbind(cds_df, data.frame(start = s, end = e, product = product, stringsAsFactors = FALSE))
      current_feature <- NA
    }
  }

  if (!is.null(cache_dir) && nrow(cds_df) > 0) {
    write.csv(cds_df, file = cache_path, row.names = FALSE)
    cat("Saved CDS cache for", accession, "to", cache_path, "\n")
  }

  return(cds_df)
}

# Dynamically fetch VP4, VP2, and VP3 sequences for EV complex prediction
fetch_vp_complex <- function(region) {
  # Local cache to prevent redundant fetching
  cache_dir <- file.path(workspace_root, "ncbi_vp_complex_cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0(region, "_vp_complex.json"))

  if (file.exists(cache_file)) {
    cached <- jsonlite::fromJSON(cache_file)
    if (length(cached$VP4) > 0 && length(cached$VP2) > 0 && length(cached$VP3) > 0) {
      return(list(VP4 = cached$VP4, VP2 = cached$VP2, VP3 = cached$VP3))
    } else {
      # If cache is {} or NULL, remove it to force re-fetch
      file.remove(cache_file)
    }
  }

  surrogate_map <- list(
    "CV-A1" = "PQ119794.1",
    "CV-A6" = "PZ182137.1",
    "CV-A16" = "PQ786279.1",
    "CV-A19" = "PP756357.1",
    "CV-B3" = "PZ117059.1",
    "CV-B4" = "PQ553589.1",
    "CV-B5" = "PZ117055.1",
    "ECHOVIRUS-E3" = "PV137773.2",
    "ECHOVIRUS-E9" = "PX441827.1",
    "ECHOVIRUS-E11" = "PV667667.1",
    "ECHOVIRUS-E25" = "OQ842428.1",
    "ECHOVIRUS-E30" = "PX441828.1"
  )

  # Try the primary accession first
  acc <- search_ncbi_accession(region)

  fetch_from_acc <- function(accession) {
    if (is.null(accession)) {
      return(list(VP4 = NULL, VP2 = NULL, VP3 = NULL))
    }
    cds_df <- get_ncbi_cds(accession)
    if (is.null(cds_df)) {
      return(list(VP4 = NULL, VP2 = NULL, VP3 = NULL))
    }
    fasta_data <- tryCatch(
      {
        res <- rentrez::entrez_fetch(db = "nuccore", id = accession, rettype = "fasta")
        Sys.sleep(NCBI_DELAY)
        res
      },
      error = function(e) NULL
    )
    if (is.null(fasta_data)) {
      return(list(VP4 = NULL, VP2 = NULL, VP3 = NULL))
    }
    seq_lines <- strsplit(fasta_data, "\n")[[1]]
    sequence <- paste(seq_lines[-1], collapse = "")
    vp4 <- NULL
    vp2 <- NULL
    vp3 <- NULL
    for (i in seq_len(nrow(cds_df))) {
      s <- cds_df$start[i]
      e <- cds_df$end[i]
      prod <- cds_df$product[i]
      if (grepl("VP4", prod, ignore.case = TRUE) && is.null(vp4)) {
        vp4 <- suppressWarnings(tryCatch(as.character(Biostrings::translate(Biostrings::DNAString(substr(sequence, s, e)))), error = function(e) NULL))
      }
      if (grepl("VP2", prod, ignore.case = TRUE) && is.null(vp2)) {
        vp2 <- suppressWarnings(tryCatch(as.character(Biostrings::translate(Biostrings::DNAString(substr(sequence, s, e)))), error = function(e) NULL))
      }
      if (grepl("VP3", prod, ignore.case = TRUE) && is.null(vp3)) {
        vp3 <- suppressWarnings(tryCatch(as.character(Biostrings::translate(Biostrings::DNAString(substr(sequence, s, e)))), error = function(e) NULL))
      }
    }
    return(list(VP4 = vp4, VP2 = vp2, VP3 = vp3))
  }

  res <- fetch_from_acc(acc)

  # If primary accession failed (likely partial sequence), try the surrogate complete genome
  if ((is.null(res$VP4) || is.null(res$VP2) || is.null(res$VP3)) && region %in% names(surrogate_map)) {
    cat(sprintf("    [Info] Primary accession lacks VP4/VP2/VP3. Falling back to surrogate complete genome for %s: %s\n", region, surrogate_map[[region]]))
    res <- fetch_from_acc(surrogate_map[[region]])
  }

  # If the backup genome is also missing specific protein labels
  # (which happens often), use saved species-level protein references instead
  # This makes sure ColabFold has a valid starting structure
  if (is.null(res$VP4) || is.null(res$VP2) || is.null(res$VP3)) {
    cat(sprintf("    [Info] Surrogate also lacked granular annotation. Using reference proxy for %s.\n", region))
    if (region %in% c("CV-A1", "CV-A6", "CV-A16", "CV-A19")) {
      # Enterovirus Species A proxies (Based on CV-A16 Q9QF31)
      res$VP4 <- "GAQVSTQKTGAHENSNSGPNIIHYTNINYYKDAASNSANRQDFTQDPGKFTEPIKDMLIKSMQGLNSP"
      res$VP2 <- "SPSAEACGYSDRVAQLTIGNSTITTQEAANIVIAYGEWPEYCPDTDATAVDKPTRPDVSVNRFFTLDTKSWAKDSKGWYWKFPDVLTEVGVFGQNAQFHYLYRSGFCVHVQCNASKFHQGALLVAVLPEYVLGTIAGGTGNENSHPPYATTQPGQVGAVLTHPYVLDAGIPLSQLTVCPHQWINLRTNNCATIIVPYMNTVPFDSALNHCNFGLLVVPVVPLDFNAGATSEIPITVTIAPMCAEFAGLRQAVKQ"
      res$VP3 <- "GIPTELKPGTNQFLTTDDGVSAPILPGFHPTPPIHIPGEVHNLLEICRVETILEVNNLKTNETTPMQRLCFPVSVQSKTGELCAAFRADPGRDGPWQSTILGQLCRYYTQWSGSLEVTFMFAGSFMATGKMLIAYTPPGGNVPADRITAMLGTHVIWDFGLQSSVTLVVPWISNTHYRAHARAGYFDYYTTGIITIWYQTNYVVPIGAPTTAYIVALAAAQDNFTMKLCKDTEDIEQTANIQ"
    } else {
      # Enterovirus Species B proxies (Based on CV-B3 P03313)
      res$VP4 <- "GAQVSTQKTGAHETRLNASGNSIIHYTNINYYKDAASNSANRQDFTQDPGKFTEPVKDIMIKSLPALN"
      res$VP2 <- "SPTVEECGYSDRARSITLGNSTITTQECANVVVGYGVWPDYLKDSEATAEDQPTQPDVATCRFYTLDSVQWQKTSPGWWWKLPDALSNLGLFGQNMQYHYLGRTGYTVHVQCNASKFHQGCLLVVCVPEAEMGCATLDNTPSSAELLGGDTAKEFADKPVASGSNKLVQRVVYNAGMGVGVGNLTIFPHQWINLRTNNSATIVMPYTNSVPMDNMFRHNNVTLMVIPFVPLDYCPGSTTYVPITVTIAPMCAEYNGLRLAGHQ"
      res$VP3 <- "GLPTMNTPGSCQFLTSDDFQSPSAMPQYDVTPEMRIPGEVKNLMEIAEVDSVVPVQNVGEKVNSMEAYQIPVRSNEGSGTQVFGFPLQPGYSSVFSRTLLGEILNYYTHWSGSIKLTFMFCGSAMATGKFLLAYSPPGAGAPTKRVDAMLGTHVIWDVGLQSSCVLCIPWISQTHYRFVASDEYTAGGFITCWYQTNIVVPADAQSSCYIMCFVSACNDFSVRLLKDTPFISQQNFFQ"
    }
  }

  jsonlite::write_json(res, cache_file, null = "null")
  return(res)
}

# Determines the Amino Acid substitution caused by a single nucleotide mutation.
# Extracts the local codon from the reference sequence and translates the mutated codon.
get_ncbi_aa_change <- function(seq_str, cds_start, cds_end, ncbi_pos, ncbi_alt) {
  if (is.na(ncbi_pos) || ncbi_pos < cds_start || ncbi_pos > cds_end) {
    return(list(pos_aa = NA_integer_, ref_aa = NA_character_, alt_aa = NA_character_))
  }
  cds_pos <- ncbi_pos - cds_start + 1 # Relative nucleotide position in CDS
  pos_aa <- ceiling(cds_pos / 3) # Amino acid position in protein
  codon_start <- cds_start + (pos_aa - 1) * 3 # Genomic start of codon

  if (codon_start + 2 > length(seq_str)) {
    return(list(pos_aa = pos_aa, ref_aa = NA_character_, alt_aa = NA_character_))
  }
  ref_codon <- seq_str[codon_start:(codon_start + 2)] # Extract 3bp REF codon
  ref_codon_str <- paste(ref_codon, collapse = "")

  if (nchar(ref_codon_str) == 3 && !grepl("NA", ref_codon_str)) {
    ref_aa <- suppressWarnings(tryCatch(
      {
        as.character(Biostrings::translate(Biostrings::DNAString(ref_codon_str)))
      },
      error = function(e) NA_character_
    ))
  } else {
    ref_aa <- NA_character_
  }

  if (is.na(ncbi_alt) || nchar(ncbi_alt) > 1 || grepl("[^ACGTUacgtuNn]", ncbi_alt)) {
    alt_aa <- NA_character_
  } else {
    alt_codon <- ref_codon
    codon_idx <- cds_pos - (pos_aa - 1) * 3 # Find mutated pos within codon (1, 2, or 3)
    alt_codon[codon_idx] <- ncbi_alt # Substitute ALT nucleotide
    alt_codon_str <- paste(alt_codon, collapse = "")

    if (nchar(alt_codon_str) == 3 && !grepl("NA", alt_codon_str)) {
      alt_aa <- suppressWarnings(tryCatch(
        {
          as.character(Biostrings::translate(Biostrings::DNAString(alt_codon_str)))
        },
        error = function(e) NA_character_
      ))
    } else {
      alt_aa <- NA_character_
    }
  }
  return(list(pos_aa = pos_aa, ref_aa = ref_aa, alt_aa = alt_aa))
}

# Build a search query to find the right reference genome on NCBI
# accession number based on a generic viral genotype name (e.g., "CV-A16").
search_ncbi_accession <- function(region_name) {
  term_base <- ""
  if (grepl("^(GI|GII)-", region_name, ignore.case = TRUE)) {
    genotype <- gsub("-", ".", region_name, ignore.case = TRUE)
    term_base <- paste0('("Norovirus ', genotype, '" OR "Norovirus ', region_name, '")')
  } else if (grepl("^CV-B", region_name, ignore.case = TRUE)) {
    type_num <- gsub("^CV-B", "", region_name, ignore.case = TRUE)
    term_base <- paste0('"Coxsackievirus B', type_num, '"')
  } else if (grepl("^CV-A", region_name, ignore.case = TRUE)) {
    type_num <- gsub("^CV-A", "", region_name, ignore.case = TRUE)
    term_base <- paste0('"Coxsackievirus A', type_num, '"')
  } else if (grepl("^Echovirus-", region_name, ignore.case = TRUE)) {
    type_num <- gsub("^Echovirus-[Ee]?", "", region_name, ignore.case = TRUE)
    term_base <- paste0('("Echovirus ', type_num, '" OR "Echovirus E', type_num, '")')
  } else if (grepl("^EV-", region_name, ignore.case = TRUE)) {
    type_name <- gsub("^EV-", "Enterovirus ", region_name, ignore.case = TRUE)
    term_base <- paste0('"', type_name, '"')
  } else if (grepl("^RV-", region_name, ignore.case = TRUE)) {
    type_name <- gsub("^RV-", "Rhinovirus ", region_name, ignore.case = TRUE)
    term_base <- paste0('"', type_name, '"')
  } else {
    term_base <- paste0('"', region_name, '"')
  }

  queries <- c(
    paste0(term_base, ' AND "complete genome" AND refseq[Filter]'),
    paste0(term_base, ' AND "complete genome"'),
    paste0(term_base, ' AND "complete"')
  )

  for (query in queries) {
    encoded_query <- URLencode(query)
    search_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=",
      encoded_query, "&retmode=json&retmax=1"
    )
    res <- tryCatch(fromJSON(search_url), error = function(e) NULL)
    if (!is.null(res)) {
      ids <- res$esearchresult$idlist
      if (length(ids) > 0 && ids[1] != "") {
        id <- ids[1]
        summary_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=", id, "&retmode=json")
        sum_res <- tryCatch(fromJSON(summary_url), error = function(e) NULL)
        if (!is.null(sum_res)) {
          accession <- sum_res$result[[id]]$accessionversion
          if (!is.null(accession) && accession != "") {
            return(accession)
          }
        }
        return(id)
      }
    }
    Sys.sleep(NCBI_DELAY) # respect rate limit
  }
  return(NULL)
}

# Find the variant files in the local folder first, then check the main folder
# if they aren't in the local folder
locate_tsv_files <- function(subfolder, suite_name, ivar_workspace) {
  local_paths <- c(
    file.path(subfolder, "iVar"),
    subfolder
  )
  for (p in local_paths) {
    if (file.exists(p)) {
      tsvs <- list.files(path = p, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
      tsvs <- tsvs[grepl("barcode0?([1-9]|[1-8][0-9]|9[0-6])(_|\\.)", basename(tsvs))]
      if (length(tsvs) > 0) {
        cat("Found TSV files locally at:", p, "\n")
        return(list(dir = p, files = tsvs))
      }
    }
  }

  if (!is.null(ivar_workspace) && file.exists(ivar_workspace)) {
    parts <- strsplit(subfolder, "/")[[1]]
    virus <- parts[length(parts) - 1]
    run_name <- parts[length(parts)]

    remote_paths <- c(
      # Standard suite names
      file.path(ivar_workspace, "with iVar", suite_name, "Combined Run", "iVar"),
      file.path(ivar_workspace, "with iVar", suite_name, "iVar"),
      file.path(ivar_workspace, "with iVar", suite_name),
      file.path(ivar_workspace, suite_name, "Combined Run", "iVar"),
      file.path(ivar_workspace, suite_name, "iVar"),
      file.path(ivar_workspace, suite_name),
      # Structural nested paths (e.g., EV/Ajax)
      file.path(ivar_workspace, virus, run_name, "Combined Run", "iVar"),
      file.path(ivar_workspace, virus, run_name, "iVar"),
      file.path(ivar_workspace, virus, run_name)
    )
    for (p in remote_paths) {
      if (file.exists(p)) {
        tsvs <- list.files(path = p, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
        tsvs <- tsvs[grepl("barcode0?([1-9]|[1-8][0-9]|9[0-6])(_|\\.)", basename(tsvs))]
        if (length(tsvs) > 0) {
          cat("Found TSV files in iVar workspace at:", p, "\n")
          return(list(dir = p, files = tsvs))
        }
      }
    }
  }
  return(NULL)
}

# Setup and Initialization
cArray <- c(
  "dodgerblue2", "#E31A1C", "darkorange", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70",
  "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "forestgreen",
  "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown", "black", "forestgreen", "red2", "orange", "cornflowerblue",
  "magenta", "darkolivegreen4", "indianred1", "tan4", "darkblue",
  "mediumorchid1", "firebrick4", "yellowgreen", "lightsalmon", "tan3",
  "tan1", "darkgray", "wheat4", "#DDAD4B", "chartreuse",
  "seagreen1", "moccasin", "mediumvioletred", "seagreen", "cadetblue1",
  "darkolivegreen1", "tan2", "tomato3", "#7CE3D8", "gainsboro", "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1"
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

default_fallback_mappings <- c(
  "GII-1" = "NC_039474.1", "GII-2" = "NC_039476.1", "GII-3" = "U02030.1",
  "GII-4" = "NC_039477.1", "GII-6" = "AJ277614.1", "GII-7" = "AJ277615.1",
  "GII-8" = "AF195848.1", "GII-10" = "AF504671.1", "GII-11" = "AB074893.1",
  "GII-12" = "NC_029646.1", "GII-13" = "AY113106.1", "GII-14" = "AY134748.1",
  "GII-15" = "AY113105.1", "GII-16" = "AY502010.1", "GII-17" = "NC_039475.1",
  "GII-20" = "AJ626814.1", "GII-21" = "AY675554.1", "GII-22" = "AB083780.1"
)

# Genotyper Pipeline
# PHASE 1: Scan sequences to identify Enterovirus/Norovirus strains
# Make plots showing mutation amounts over time
run_genotyper_pipeline <- function(subfolder, output_dir, blast_bin, cArray, genotype_mapping, qpcr_dfs = NULL, data_only = FALSE, plot_only = FALSE) {
  writeLines(paste("DEBUG: top of run_genotyper_pipeline. qpcr_dfs is null?", is.null(qpcr_dfs)), file.path(output_dir, "debug_top.txt"))
  cat("\n==================================================\n")
  cat("Running Genotyper pipeline on:", subfolder, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("==================================================\n")

  parts <- strsplit(subfolder, "/")[[1]]
  virus <- parts[length(parts) - 1]
  city <- parts[length(parts)]
  p_pcoa <- NULL

  orig_dir <- getwd()
  setwd(subfolder)
  on.exit(setwd(orig_dir), add = TRUE)

  xlsx_files <- list.files(pattern = "\\.xlsx$")
  xlsx_files <- xlsx_files[!grepl("Analysis_Final|Frequencies|blast_top_hits|genotype_pc|abund_mat|bc_scores", xlsx_files)]
  if (length(xlsx_files) == 0) {
    stop("No metadata Excel file found in ", subfolder)
  }
  metadata_file <- xlsx_files[1]
  cat("Using metadata file:", metadata_file, "\n")

  df_meta <- read.xlsx(metadata_file, sheet = "Sheet1")
  fasta_file_mapping <- setNames(df_meta$date, df_meta$file)

  if (!plot_only) {

  makeblastdb <- if (blast_bin == "") "makeblastdb" else file.path(blast_bin, "makeblastdb")
  blastn <- if (blast_bin == "") "blastn" else file.path(blast_bin, "blastn")

  if (virus == "EV") {
    types <- "all_EV_types_combined.fasta"
    db_name <- "EV_Reference_Sequences"
  } else {
    types <- "v3-all_NV_types_combined.fasta"
    db_name <- "NV-v3_Reference_Sequences"
  }

  if (!file.exists(types)) {
    src_fasta <- NULL
    if (file.exists(file.path("..", types))) {
      src_fasta <- file.path("..", types)
    } else if (file.exists(file.path("../..", types))) {
      src_fasta <- file.path("../..", types)
    } else if (file.exists(file.path(ivar_workspace, types))) {
      src_fasta <- file.path(ivar_workspace, types)
    }

    if (!is.null(src_fasta)) {
      file.copy(src_fasta, types, overwrite = TRUE)
    } else {
      stop("Reference FASTA file ", types, " not found!")
    }
  }

  if (!file.exists(paste0(db_name, ".nin"))) {
    cat("Creating BLAST database from reference file...\n")
    real_tmp <- "/tmp"
    if (!dir.exists(real_tmp)) real_tmp <- tempdir()

    unique_suffix <- paste0("_", city, "_", sample(1e6:9e6, 1))
    tmp_ref_fasta <- file.path(real_tmp, paste0(basename(types), unique_suffix))
    file.copy(types, tmp_ref_fasta, overwrite = TRUE)

    tmp_db_prefix <- file.path(real_tmp, paste0(db_name, unique_suffix))

    cmd <- paste(shQuote(makeblastdb), "-in", shQuote(tmp_ref_fasta), "-dbtype nucl -out", shQuote(tmp_db_prefix))
    system(cmd)

    db_exts <- c(".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto")
    for (ext in db_exts) {
      src_db <- paste0(tmp_db_prefix, ext)
      dest_db <- paste0(db_name, ext)
      if (file.exists(src_db)) {
        file.copy(src_db, dest_db, overwrite = TRUE)
        file.remove(src_db)
      }
    }
    file.remove(tmp_ref_fasta)
  } else {
    cat("BLAST database already exists. Skipping creation.\n")
  }

  WW_reads <- list.files(pattern = "^barcode[0-96]+\\_combined.fasta$")
  if (length(WW_reads) == 0) {
    log_warning_error("Warning: No barcode fasta files found in", subfolder, "\n")
    return(NULL)
  }

  wb <- createWorkbook()
  total_queries_list <- list()

  for (query in WW_reads) {
    cat("Checking if", query, "results exist...")
    out_file <- sub("\\.fasta$", "_blast_results.txt", query)

    if (file.exists(out_file)) {
      cat(" Found (cached). Loading results...\n")
    } else {
      cat(" Not found. Running BLAST...\n")
      system(paste(
        shQuote(blastn), "-query", shQuote(query), "-db", shQuote(db_name),
        "-outfmt 6 -out", shQuote(out_file),
        "-task dc-megablast -num_threads 8 -evalue 1e-3 -word_size 11"
      ))
      cat("BLAST query completed.\n")
    }

    query_seqs <- readDNAStringSet(query)
    total_queries <- length(query_seqs)
    total_queries_list[[query]] <- total_queries

    if (file.info(out_file)$size == 0) {
      warning("No hits found for ", query)
      results <- data.frame()
    } else {
      results <- read.table(out_file,
        header = FALSE, sep = "\t",
        col.names = c(
          "query", "subject", "pident", "length",
          "mismatch", "gapopen", "qstart", "qend",
          "sstart", "send", "evalue", "bitscore"
        )
      )
    }

    if (nrow(results) > 0) {
      results$source_file <- query
      attr(results, "total_queries") <- total_queries

      top_hits <- results %>%
        group_by(query) %>%
        top_n(1, bitscore) %>%
        ungroup()

      low_id_idx <- which(top_hits$pident < 80)
      if (length(low_id_idx) > 0) {
        mapped_subjects <- genotype_mapping[top_hits$subject[low_id_idx]]
        top_hits$subject[low_id_idx] <- ifelse(!is.na(mapped_subjects), mapped_subjects, top_hits$subject[low_id_idx])
      }
    } else {
      top_hits <- data.frame()
    }

    sheet_name <- gsub("\\.fasta$", "", query)
    sheet_name <- substr(sheet_name, 1, 31)

    addWorksheet(wb, sheet_name)
    if (nrow(top_hits) > 0) {
      writeData(wb, sheet = sheet_name, top_hits)
    }
    cat("Finished processing", query, "\n")
  }

  excel_file <- file.path(output_dir, "blast_top_hits_by_sample.xlsx")
  saveWorkbook(wb, excel_file, overwrite = TRUE)

  sheets <- getSheetNames(excel_file)
  summary_list <- list()

  for (sheet in sheets) {
    data <- read.xlsx(excel_file, sheet = sheet)
    fasta_file <- paste0(sheet, ".fasta")
    sample_date <- fasta_file_mapping[[fasta_file]]
    if (is.null(sample_date)) sample_date <- sheet

    total_queries <- total_queries_list[[fasta_file]]
    if (is.null(total_queries)) {
      total_queries <- 0
    }

    if (!is.null(data) && nrow(data) > 0) {
      summary <- data %>%
        group_by(subject) %>%
        summarise(count = n(), avg_pident = mean(pident), .groups = "drop")
      assigned_count <- sum(summary$count)
    } else {
      summary <- tibble(subject = character(), count = integer(), avg_pident = numeric())
      assigned_count <- 0
    }

    unassigned_count <- total_queries - assigned_count
    if (unassigned_count > 0) {
      summary <- bind_rows(
        summary,
        tibble(subject = "Could not assign", count = unassigned_count, avg_pident = NA)
      )
    }

    summary <- summary %>%
      mutate(sample = sample_date, total_hits = sum(count))
    summary_list[[sheet]] <- summary
  }

  summary_df <- bind_rows(summary_list)

  genotype_percentage_results <- summary_df %>%
    group_by(sample) %>%
    mutate(Percentage = (count / sum(count)) * 100) %>%
    select(sample, genotype = subject, Percentage, total_hits, avg_pident)

  desired_order <- fasta_file_mapping %>%
    unname() %>%
    unique()
  genotype_percentage_results$sample <- factor(genotype_percentage_results$sample, levels = desired_order)

  # Create an unfiltered version for regression
  genotype_percentage_results_unfiltered <- genotype_percentage_results[genotype_percentage_results$total_hits >= 50, ]

  # Export unfiltered Excel
  wb_unfiltered <- createWorkbook()
  addWorksheet(wb_unfiltered, virus)
  writeData(wb_unfiltered, sheet = virus, genotype_percentage_results_unfiltered)
  saveWorkbook(wb_unfiltered, file.path(output_dir, "genotype_pc_results_unfiltered.xlsx"), overwrite = TRUE)

  # Pivot unfiltered matrix and save
  abund_mat_unfiltered <- genotype_percentage_results_unfiltered %>%
    select(sample, genotype, Percentage) %>%
    pivot_wider(
      names_from = genotype,
      values_fn = sum,
      values_from = Percentage
    )
  if (nrow(abund_mat_unfiltered) > 0) {
    abund_mat_unfiltered <- abund_mat_unfiltered %>%
      column_to_rownames("sample") %>%
      as.matrix()
    abund_mat_unfiltered[is.na(abund_mat_unfiltered)] <- 0
    write.csv(abund_mat_unfiltered, file = file.path(output_dir, paste0(city, "_", virus, "_abund_mat_unfiltered.csv")), row.names = TRUE)
  }

  genotype_percentage_results <- genotype_percentage_results %>%
    mutate(genotype = if_else(Percentage < 5 & genotype != "Could not assign", "Other", genotype))

  genotype_percentage_results <- genotype_percentage_results[genotype_percentage_results$total_hits >= 50, ]

  wb_pc <- createWorkbook()
  addWorksheet(wb_pc, virus)
  writeData(wb_pc, sheet = virus, genotype_percentage_results)
  saveWorkbook(wb_pc, file.path(output_dir, "genotype_pc_results.xlsx"), overwrite = TRUE)

  } else {
    # Load data from previous data_only phase
    pc_file <- file.path(output_dir, "genotype_pc_results.xlsx")
    if (!file.exists(pc_file)) {
      cat("Data file not found for", subfolder, "- skipping plot phase.\n")
      return(NULL)
    }
    genotype_percentage_results <- openxlsx::read.xlsx(pc_file, sheet = 1)
    # Restore factor levels for samples
    desired_order <- fasta_file_mapping %>% unname() %>% unique()
    genotype_percentage_results$sample <- factor(genotype_percentage_results$sample, levels = desired_order)
  } # End of !plot_only block

  if (data_only) return(NULL)

  p_abund <- ggplot(genotype_percentage_results, aes(x = sample, y = Percentage, fill = genotype)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    scale_fill_manual(values = cArray) +
    labs(x = "Samples", y = "Relative Abundance", fill = "Genotype") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(output_dir, paste0(virus, "_", city, "_relative_abundance.tiff")),
    plot = p_abund, device = "tiff", dpi = 600, width = 10, height = 6, compression = "lzw"
  )

  abund_mat <- genotype_percentage_results %>%
    select(sample, genotype, Percentage) %>%
    pivot_wider(
      names_from = genotype,
      values_fn = sum,
      values_from = Percentage
    )

  if (nrow(abund_mat) > 0) {
    abund_mat <- abund_mat %>%
      column_to_rownames("sample") %>%
      as.matrix()
    abund_mat[is.na(abund_mat)] <- 0

    meta <- data.frame(
      sample_id = rownames(abund_mat),
      sample_date = as.Date(rownames(abund_mat), format = "%d/%m/%Y")
    )
    meta$season <- getSeason(meta$sample_date)
    meta$season <- as.factor(meta$season)

    write.csv(abund_mat, file = file.path(output_dir, paste0(city, "_", virus, "_abund_mat.csv")), row.names = TRUE)

    if (nrow(abund_mat) >= 2) {
      bc_dist <- vegdist(abund_mat, method = "bray")
      bc_pcoa <- cmdscale(bc_dist, eig = TRUE, k = 2)

      bc_scores <- as.data.frame(bc_pcoa$points)
      colnames(bc_scores) <- c("PCoA1", "PCoA2")
      bc_scores$sample_id <- rownames(bc_scores)
      bc_scores$sample_date <- as.Date(bc_scores$sample_id, format = "%d/%m/%Y")
      bc_scores$sample_id <- factor(bc_scores$sample_id, levels = bc_scores$sample_id)

      n_samples <- length(unique(bc_scores$sample_id))
      rainbow_colors <- rainbow(n_samples, start = 0, end = 0.9)
      names(rainbow_colors) <- levels(bc_scores$sample_id)

      pcoa1_var <- round(100 * bc_pcoa$eig[1] / sum(bc_pcoa$eig), 1)
      pcoa2_var <- round(100 * bc_pcoa$eig[2] / sum(bc_pcoa$eig), 1)

      bc_scores <- merge(bc_scores, meta[, c("sample_id", "season")], by = "sample_id")

      write.csv(bc_scores, file = file.path(output_dir, paste0(city, "_", virus, "_bc_scores.csv")), row.names = FALSE)

      season_shapes <- c("Spring" = 21, "Summer" = 22, "Fall" = 23, "Winter" = 24)
      level_season_lookup <- bc_scores %>%
        distinct(sample_id, season) %>%
        mutate(sample_id = factor(sample_id, levels = levels(bc_scores$sample_id))) %>%
        arrange(sample_id)

      # Identify top 5% outliers for season
      bc_scores$outlier_label <- NA
      valid_idx <- !is.na(bc_scores$season)
      if (sum(valid_idx) > 0) {
        group_cents <- bc_scores[valid_idx, ] %>% 
          dplyr::group_by(season) %>% 
          dplyr::summarize(cx = mean(PCoA1), cy = mean(PCoA2), .groups = "drop")
        
        p_dist_data <- bc_scores[valid_idx, ] %>% 
          dplyr::left_join(group_cents, by = "season") %>%
          dplyr::mutate(dist_to_cent = sqrt((PCoA1 - cx)^2 + (PCoA2 - cy)^2))
          
        dist_thresh <- quantile(p_dist_data$dist_to_cent, 0.95, na.rm = TRUE)
        bc_scores$outlier_label[valid_idx] <- ifelse(p_dist_data$dist_to_cent > dist_thresh, 
                                                     as.character(p_dist_data$sample_id), NA)
      }

      p_pcoa <- ggplot(bc_scores, aes(x = PCoA1, y = PCoA2, fill = season, shape = season, color = season)) +
        geom_label_repel(aes(label = outlier_label), max.overlaps = 50, box.padding = 0.5, fill = alpha("white", 0.7), color = "black", size = 3, show.legend = FALSE, na.rm = TRUE) +
        geom_point(size = 3, alpha = 0.7, stroke = 0.8) +
        scale_shape_manual(values = season_shapes, name = "Season", na.value = 25) +
        scale_color_manual(
          values = c(
            "Spring" = "#00BA38",
            "Summer" = "#F8766D",
            "Fall" = "#FF7F00",
            "Winter" = "#619CFF"
          ),
          name = "Season"
        ) +
        scale_fill_manual(
          values = c(
            "Spring" = "#00BA38",
            "Summer" = "#F8766D",
            "Fall" = "#FF7F00",
            "Winter" = "#619CFF"
          ),
          name = "Season"
        ) +
        labs(
          x = paste0("PCoA1 (", pcoa1_var, "%)"),
          y = paste0("PCoA2 (", pcoa2_var, "%)")
        ) +
        theme_bw()

      if (length(unique(bc_scores$season)) > 1 && any(table(bc_scores$season) >= 3)) {
        p_pcoa <- p_pcoa + stat_ellipse(aes(group = season),
          type = "t",
          linetype = 2,
          linewidth = 0.6,
          level = 0.95
        ) +
          guides(
            fill = guide_legend(override.aes = list(shape = 21, color = "black")),
            color = guide_legend(override.aes = list(fill = "gray80", linetype = "dashed"))
          )
      }

      ggsave(file.path(output_dir, paste0(virus, "_", city, "_pcoa_plot.tiff")),
        plot = p_pcoa, device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw"
      )

      if (length(unique(meta$season)) >= 2 && nrow(meta) >= 4) {
        tryCatch(
          {
            perm_result <- adonis2(bc_dist ~ season, data = meta, permutations = 1000, method = "bray")
            sink(file.path(output_dir, "permanova_results.txt"))
            print(perm_result)
            sink()
          },
          error = function(e) {
            cat("PERMANOVA failed for", subfolder, ":", e$message, "\n")
          }
        )
      }
    } else {
      cat("Too few samples to perform PCoA/PERMANOVA for", subfolder, "\n")
    }
  }

  p_abs_abund <- NULL
  p_qc <- NULL
  if (!is.null(qpcr_dfs) && city %in% names(qpcr_dfs)) {
    tryCatch(
      {
        cat("Integrating RT-qPCR data for absolute abundance and QC...\n")
        qpcr_df <- qpcr_dfs[[city]]

        if ("Date" %in% colnames(qpcr_df)) {
          date_col_idx <- which(colnames(qpcr_df) == "Date")
          if (date_col_idx < ncol(qpcr_df)) {
            if (inherits(qpcr_df[[date_col_idx + 1]], "POSIXct") || inherits(qpcr_df[[date_col_idx + 1]], "Date")) {
              cat("  -> Detected shifted columns in RT-qPCR data. Auto-correcting...\n")
              for (i in date_col_idx:(ncol(qpcr_df) - 1)) {
                qpcr_df[[i]] <- qpcr_df[[i + 1]]
              }
              qpcr_df[[ncol(qpcr_df)]] <- NA
            }
          }
        }

        qpcr_df$Date <- as.Date(qpcr_df$Date)

        # Set up full traceback on error
        old_opt <- options(error = function() {
          sink(file.path(output_dir, "traceback.txt"))
          print(sys.calls())
          sink()
        })
        on.exit(options(old_opt), add = TRUE)
        
        abund_df <- genotype_percentage_results %>%
          select(sample, genotype, Percentage) %>%
          rename(Date = sample) %>%
          mutate(
            Date = as.Date(Date, format = "%d/%m/%Y"),
            Rel_Abund = Percentage / 100
          )

        # Add an 'Other' category so the total always equals 100%
        other_df <- abund_df %>%
          dplyr::group_by(Date) %>%
          dplyr::summarize(
            Rel_Abund = max(0, 1 - sum(Rel_Abund, na.rm = TRUE)),
            genotype = "Other",
            .groups = "drop"
          ) %>%
          dplyr::filter(Rel_Abund > 0)
        
        abund_df <- dplyr::bind_rows(abund_df, other_df)
        
        # Complete the grid to avoid geom_area rendering gaps for missing genotypes
        abund_df <- abund_df %>%
          dplyr::ungroup() %>%
          tidyr::complete(Date, genotype, fill = list(Rel_Abund = 0, Percentage = 0))

        merged_df <- abund_df %>% dplyr::left_join(qpcr_df, by = "Date")

        qpcr_cols_available <- colnames(merged_df)
        qpcr_col <- NA

        if (grepl("NV|NoV|Norovirus", virus, ignore.case = TRUE)) {
          possible_cols <- qpcr_cols_available[grepl("^Norovirus.*\\(gc/mL\\)$", qpcr_cols_available, ignore.case = TRUE)]
        } else if (grepl("EV|Enterovirus", virus, ignore.case = TRUE)) {
          possible_cols <- qpcr_cols_available[grepl("^EV.*\\(gc/mL\\)$|^Enterovirus.*\\(gc/mL\\)$", qpcr_cols_available, ignore.case = TRUE)]
        } else {
          possible_cols <- character(0)
        }

        raw_cols <- possible_cols[!grepl("7-day|average|SD", possible_cols, ignore.case = TRUE) & grepl("gc/mL", possible_cols, ignore.case = TRUE)]

        if (length(raw_cols) > 0) {
          qpcr_col <- raw_cols[1]
        } else if (length(possible_cols) > 0) {
          qpcr_col <- possible_cols[1]
        }

        debug_file <- file.path(output_dir, "debug_abs_abund.txt")
        sink(debug_file)
        cat("DEBUG [", city, "]: qpcr_col =", qpcr_col, "\n")
        cat("DEBUG [", city, "]: colnames =", paste(colnames(merged_df), collapse = ", "), "\n")
        cat("DEBUG [", city, "]: nrow before filter =", nrow(merged_df), "\n")
        sink()

        if (qpcr_col %in% colnames(merged_df)) {
          merged_df <- merged_df %>% dplyr::filter(!is.na(.data[[qpcr_col]]))

          merged_df$Absolute_Abundance <- merged_df$Rel_Abund * merged_df[[qpcr_col]]

          total_vl_df <- merged_df %>% dplyr::distinct(Date, .keep_all = TRUE)

          virus_display <- if (grepl("NV|NoV|Norovirus", virus, ignore.case = TRUE)) "Norovirus GII" else if (grepl("EV|Enterovirus", virus, ignore.case = TRUE)) "Enterovirus" else virus

          p_abs_abund <- ggplot(merged_df, aes(x = Date)) +
            geom_area(aes(y = Absolute_Abundance, fill = genotype), alpha = 0.8) +
            geom_line(data = total_vl_df, aes(y = .data[[qpcr_col]]), color = "black", linewidth = 1, alpha = 0.7) +
            geom_point(data = total_vl_df, aes(y = .data[[qpcr_col]]), color = "black", size = 2) +
            theme_minimal() +
            scale_fill_manual(values = cArray) +
            labs(
              title = paste(virus_display, "-", city),
              x = "Date", y = "Viral Load (gc/mL)"
            ) +
            theme(legend.position = "bottom")

          qpcr_out_dir <- file.path(dirname(output_dir), "RT-qPCR")
          if (!dir.exists(qpcr_out_dir)) dir.create(qpcr_out_dir, recursive = TRUE, showWarnings = FALSE)

          ggsave(file.path(qpcr_out_dir, paste0(virus, "_", city, "_absolute_abundance.tiff")),
            plot = p_abs_abund, device = "tiff", dpi = 600, width = 10, height = 6, compression = "lzw"
          )

          calc_shannon <- function(x) {
            x <- x[x > 0]
            -sum(x * log(x))
          }
          diversity_mat <- abund_mat / 100

          div_res <- data.frame(
            Date = as.Date(rownames(diversity_mat), format = "%d/%m/%Y"),
            Shannon_Diversity = apply(diversity_mat, 1, calc_shannon)
          ) %>%
            dplyr::left_join(qpcr_df, by = "Date") %>%
            dplyr::filter(!is.na(.data[[qpcr_col]]))

          log_x <- log10(div_res[[qpcr_col]])
          y <- div_res$Shannon_Diversity
          valid_idx <- is.finite(log_x) & is.finite(y)

          if (sum(valid_idx) > 2) {
            model <- lm(y[valid_idx] ~ log_x[valid_idx])
            r2 <- round(summary(model)$r.squared, 3)
            pval <- summary(model)$coefficients[2, 4]
            pval_str <- if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)
            subtitle_comb <- paste0("R² = ", r2, " | p-value = ", pval_str, " | n = ", sum(valid_idx))
          } else {
            subtitle_comb <- "Insufficient points for regression"
          }

          # Save shannon diversity df for combined analysis
          div_res$log10_viral_load <- log_x
          write.csv(div_res, file.path(qpcr_out_dir, paste0(virus, "_", city, "_shannon_div.csv")), row.names = FALSE)


          p_qc <- ggplot(div_res, aes(x = .data[[qpcr_col]], y = Shannon_Diversity)) +
            geom_point(alpha = 0.7, color = "blue", size = 3) +
            geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
            scale_x_log10() +
            theme_bw() +
            labs(
              title = paste(virus_display, "-", city),
              x = paste(qpcr_col, "[Log10 Scale]"), y = "Shannon Diversity Index"
            ) +
            annotate("label", x = Inf, y = Inf, label = subtitle_comb, hjust = 1, vjust = 1, fill = "white", color = "black")

          ggsave(file.path(qpcr_out_dir, paste0(virus, "_", city, "_qc_scatter.tiff")),
            plot = p_qc, device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw"
          )
            merged_df$Rel_Abund <- merged_df$Rel_Abund * 100
        }
      },
      error = function(e) {
        cat("Failed to integrate RT-qPCR data:", e$message, "\n")
        sink(file.path(output_dir, "qPCR_error.txt"))
        cat("Error message:\n", e$message, "\n")
        print(sys.calls())
        cat("\n--- abund_df before error ---\n")
        try(print(str(abund_df)), silent = TRUE)
        sink()
      }
    )
  }

  p_lol <- NULL

  ret_plots <- list(abundance = p_abund, abs_abund = p_abs_abund, qc = p_qc)
  if (!is.null(p_pcoa)) ret_plots$pcoa <- p_pcoa
  return(ret_plots)
}

# iVar Pipeline
# PHASE 2: Analyze variant files from iVar
# Calculate diversity metrics and plot mutation graphs
run_ivar_pipeline <- function(subfolder, output_dir, ivar_workspace, virus, suite_name) {
  parts <- strsplit(subfolder, "/")[[1]]
  city <- parts[length(parts)]
  virus_display <- ifelse(virus == "NV", "Norovirus GII", "Enterovirus")

  cat("\n==================================================\n")
  cat("Running iVar pipeline on suite:", suite_name, "\n")
  cat("==================================================\n")

  # Restart Python in the background worker so it doesn't crash
  if (exists("ESM_READY") && ESM_READY) {
    if (file.exists(file.path(workspace_root, ".venv", "bin", "python"))) {
      reticulate::use_python(file.path(workspace_root, ".venv", "bin", "python"), required = TRUE)
    } else if (file.exists("/home/wastewater/miniconda3/bin/python")) {
      reticulate::use_python("/home/wastewater/miniconda3/bin/python", required = TRUE)
    } else {
      reticulate::use_python("/opt/homebrew/Caskroom/miniforge/base/bin/python3", required = TRUE)
    }
    tryCatch(
      {
        reticulate::source_python(file.path(workspace_root, "score_mutation.py"))
      },
      error = function(e) {
        cat("Worker Python Init Warning:", e$message, "\n")
      }
    )
  }

  # Locate TSV files
  tsv_info <- locate_tsv_files(subfolder, suite_name, ivar_workspace)
  if (is.null(tsv_info)) {
    cat("No iVar TSV files found for suite:", suite_name, ". Skipping iVar pipeline.\n")
    return(NULL)
  }

  ivar_plots <- list(heatmaps = list(), lollipops = list())
  p_depth <- NULL
  genome_plot <- NULL
  traj_plot <- NULL
  div_plot <- NULL
  pnps_plot <- NULL
  tajima_plot <- NULL
  strand_plot <- NULL
  plot_df <- NULL
  mutation_dynamics <- NULL
  reg_plot <- NULL
  reg_plot_faceted <- NULL

  # Find metadata file
  meta_search_dirs <- c(tsv_info$dir, dirname(tsv_info$dir), subfolder)
  metadata_file <- NULL
  for (d in meta_search_dirs) {
    if (!file.exists(d)) next
    files <- list.files(path = d, pattern = "\\.xlsx$", full.names = TRUE)
    files <- files[!grepl("Analysis_Final|Frequencies|blast_top_hits|genotype_pc|abund_mat|bc_scores", basename(files))]
    if (length(files) > 0) {
      metadata_file <- files[1]
      break
    }
  }

  if (is.null(metadata_file)) {
    log_warning_error("Warning: No metadata Excel file found for iVar analysis. Skipping.\n")
    return(NULL)
  }
  cat("Using iVar metadata file:", metadata_file, "\n")

  df_meta <- read.xlsx(metadata_file, sheet = "Sheet1")
  df_meta$file <- sub("\\_combined.fasta$", "", df_meta$file)
  df_meta$date_obj <- as.Date(df_meta$date, format = "%d/%m/%Y")

  barcode_date <- setNames(as.character(df_meta$date), trimws(as.character(df_meta$file)))
  chronological_dates <- df_meta %>%
    arrange(date_obj) %>%
    pull(date) %>%
    unique()

  # Load TSV files
  cat("Loading", length(tsv_info$files), "TSV files in parallel...\n")
  snv_df <- furrr::future_map_dfr(tsv_info$files, load_tsv_file, .options = furrr::furrr_options(seed = TRUE))
  cat("Loaded", nrow(snv_df), "PASS-filtered rows across", length(tsv_info$files), "samples.\n")

  if (nrow(snv_df) == 0) {
    cat("Warning: All loaded TSVs were empty or filtered out. Skipping iVar pipeline.\n")
    return(NULL)
  }

  # Filter iVar analyses by genotypes shown in the relative abundance plot
  genotype_pc_file <- file.path(dirname(output_dir), "genotype_pc_results.xlsx")
  if (file.exists(genotype_pc_file)) {
    genotype_pc <- read.xlsx(genotype_pc_file, sheet = 1)
    valid_genotypes <- unique(genotype_pc$genotype[
      !genotype_pc$genotype %in% c("Other", "Could not assign")
    ])
    if (length(valid_genotypes) > 0) {
      norm_valid <- toupper(gsub("\\.", "-", valid_genotypes))
      snv_df <- snv_df %>%
        filter(toupper(gsub("\\.", "-", REGION)) %in% norm_valid)
      cat("Filtered iVar variants to", length(valid_genotypes), "valid genotypes present in abundance plot.\n")
    }
  }

  # Fix an iVar bug: Recalculate true depth and frequency without gaps
  snv_df <- snv_df %>%
    group_by(sample_id, POS) %>%
    mutate(TOTAL_DP = REF_DP + sum(ALT_DP, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(ALT_FREQ = ALT_DP / TOTAL_DP)

  # Filter indel artifacts
  snv_df <- snv_df %>%
    filter(!(grepl("[+-]", ALT) & ALT_FREQ < MIN_INDEL_FREQ))

  # Add sample dates mapping prior to filtering
  snv_df$date_str <- barcode_date[as.character(snv_df$sample_id)]
  snv_df$sample <- factor(snv_df$date_str, levels = chronological_dates)

  # Generate read depth summary plot BEFORE filtering by depth
  plot_df <- snv_df %>% filter(!is.na(sample))
  if (nrow(plot_df) > 0) {
    cat("Generating read depth summary plot...\n")
    tryCatch(
      {
        tiff(file.path(output_dir, "iVar_Read_Depth_Summary.tiff"), units = "in", width = 12, height = 8, res = 300)
        p_depth <- ggplot(plot_df, aes(x = sample, y = TOTAL_DP, fill = REGION)) +
          geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.alpha = 0.5) +
          geom_hline(yintercept = MIN_DEPTH, linetype = "dashed", color = "#E31A1C", linewidth = 0.8) +
          annotate("label", x = 1.2, y = MIN_DEPTH, label = paste("Min Depth Cutoff (", MIN_DEPTH, "x)"), color = "white", fill = "#E31A1C", fontface = "bold", size = 3, hjust = 0, vjust = -0.5) +
          scale_y_log10(labels = scales::comma) +
          scale_fill_manual(values = rep(cArray, length.out = length(unique(plot_df$REGION)))) +
          theme_bw() +
          labs(
            title = paste(virus_display, "-", city),
            subtitle = "Distribution of total read depth (log10 scale) at variant positions across samples prior to filtering.",
            x = "Sample Date",
            y = "Read Depth (log10 scale)",
            fill = "Genotype"
          ) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", size = 14),
            legend.position = "bottom"
          )
        print(p_depth)
        dev.off()
      },
      error = function(e) {
        log_warning_error("Warning: Failed to generate read depth summary plot:", e$message, "\n")
        if (!is.null(dev.list())) dev.off()
      }
    )
  }

  snv_df <- snv_df %>% filter(TOTAL_DP >= MIN_DEPTH)
  cat("Filtered to", nrow(snv_df), "rows with read depth >=", MIN_DEPTH, "x\n")

  snv_df <- snv_df %>%
    filter(!is.na(sample)) %>%
    mutate(
      mutation = paste0(REF, POS, ALT),
      date_obj = as.Date(as.character(sample), format = "%d/%m/%Y"),
      site_pi = 2 * ALT_FREQ * (1 - ALT_FREQ),
      var_type = case_when(
        is.na(REF_AA) | is.na(ALT_AA) ~ "noncoding",
        REF_AA == ALT_AA ~ "synonymous",
        ALT_AA == "*" ~ "stop_gain",
        REF_AA == "*" ~ "stop_loss",
        TRUE ~ "nonsynonymous"
      )
    )

  # Check for strand bias and correct for multiple testing
  snv_df <- snv_df %>%
    group_by(sample_id) %>%
    mutate(
      strand_bias_pval = strand_bias_pvals(REF_RV, REF_DP, ALT_RV, ALT_DP),
      # Apply FDR multiple testing correction within each sample
      strand_bias_adj_pval = p.adjust(strand_bias_pval, method = "BH"),
      strand_bias_flag = !is.na(strand_bias_adj_pval) & strand_bias_adj_pval < STRAND_BIAS_ALPHA
    ) %>%
    ungroup()

  cat(sprintf(
    "Strand bias: %d / %d variants flagged after FDR correction (%.1f%%)\n",
    sum(snv_df$strand_bias_flag, na.rm = TRUE), nrow(snv_df),
    100 * mean(snv_df$strand_bias_flag, na.rm = TRUE)
  ))

  # Remove rare mutations that only appear on one strand
  orig_count <- nrow(snv_df)
  snv_df <- snv_df %>%
    filter(!(strand_bias_flag & ALT_FREQ < STRAND_BIAS_FILT_FREQ))
  cat(sprintf(
    "Filtered out %d low-frequency strand-biased variants (ALT_FREQ < %.2f). Remaining: %d rows.\n",
    orig_count - nrow(snv_df), STRAND_BIAS_FILT_FREQ, nrow(snv_df)
  ))

  # Reference sequences
  if (virus == "EV") {
    types <- "all_EV_types_combined_VP1.fasta"
  } else {
    types <- "v3-all_NV_types_combined.fasta"
  }

  # Locate reference sequence FASTA
  custom_fasta_path <- NULL
  fasta_search_paths <- c(
    file.path(subfolder, types),
    file.path(tsv_info$dir, types),
    file.path(dirname(tsv_info$dir), types),
    file.path(dirname(dirname(tsv_info$dir)), types),
    file.path(ivar_workspace, "with iVar", types),
    file.path(ivar_workspace, types)
  )
  for (p in fasta_search_paths) {
    if (file.exists(p)) {
      custom_fasta_path <- p
      break
    }
  }

  if (is.null(custom_fasta_path)) {
    log_warning_error("Warning: Reference FASTA file not found for iVar analysis. Skipping.\n")
    return(NULL)
  }
  cat("Using custom FASTA reference:", custom_fasta_path, "\n")

  custom_seqs <- readDNAStringSet(custom_fasta_path)
  names(custom_seqs) <- trimws(names(custom_seqs))

  genome_lengths <- setNames(as.integer(width(custom_seqs)), names(custom_seqs))
  partial_seqs <- names(genome_lengths)[genome_lengths < 3000]
  if (length(partial_seqs) > 0) {
    warning("Potentially partial sequences (< 3000 bp): ", paste(partial_seqs, collapse = ", "))
  }

  # Accession mappings
  mapping_cache_file <- file.path(output_dir, "ncbi_accession_mappings.csv")

  if (virus == "EV") {
    vp1_csv <- file.path(ivar_workspace, "EV_VP1_mappings.csv")
    if (file.exists(vp1_csv)) {
      file.copy(vp1_csv, mapping_cache_file, overwrite = TRUE)
    }
  }

  ncbi_mappings <- c()
  if (file.exists(mapping_cache_file)) {
    ncbi_mappings_df <- read.csv(mapping_cache_file, stringsAsFactors = FALSE)
    if (nrow(ncbi_mappings_df) > 0) {
      ncbi_mappings_df <- ncbi_mappings_df[!grepl("^[0-9]+$", ncbi_mappings_df$Accession), ]
      ncbi_mappings <- setNames(ncbi_mappings_df$Accession, ncbi_mappings_df$REGION)
    }
  } else {
    orig_cache <- file.path(subfolder, "ncbi_accession_mappings.csv")
    if (file.exists(orig_cache)) {
      ncbi_mappings_df <- read.csv(orig_cache, stringsAsFactors = FALSE)
      if (nrow(ncbi_mappings_df) > 0) {
        ncbi_mappings_df <- ncbi_mappings_df[!grepl("^[0-9]+$", ncbi_mappings_df$Accession), ]
        ncbi_mappings <- setNames(ncbi_mappings_df$Accession, ncbi_mappings_df$REGION)
      }
      if (length(ncbi_mappings) > 0) {
        df_to_save <- data.frame(
          REGION = names(ncbi_mappings),
          Accession = as.character(ncbi_mappings),
          stringsAsFactors = FALSE
        )
        write.csv(df_to_save, file = mapping_cache_file, row.names = FALSE)
      }
    }
  }

  cache_file <- file.path(output_dir, "ncbi_refseqs_cached.fasta")
  orig_refseq_cache <- file.path(subfolder, "ncbi_refseqs_cached.fasta")
  if (file.exists(orig_refseq_cache) && !file.exists(cache_file)) {
    file.copy(orig_refseq_cache, cache_file)
  }

  if (file.exists(cache_file)) {
    cached_seqs <- readDNAStringSet(cache_file)
  } else {
    cached_seqs <- DNAStringSet()
  }

  coord_mappings <- list()
  ncbi_seqs_loaded <- list()
  ncbi_cds_loaded <- list()

  get_ncbi_accession_for_region_local <- function(region_name) {
    match_idx <- match(toupper(gsub("\\.", "-", region_name)), toupper(gsub("\\.", "-", names(ncbi_mappings))))
    if (!is.na(match_idx)) {
      return(ncbi_mappings[match_idx])
    }

    cat("Searching NCBI for reference sequence for", region_name, "...\n")
    accession <- search_ncbi_accession(region_name)
    if (is.null(accession)) {
      match_fallback <- match(toupper(region_name), toupper(names(default_fallback_mappings)))
      if (!is.na(match_fallback)) {
        accession <- default_fallback_mappings[match_fallback]
      }
    }

    if (!is.null(accession)) {
      ncbi_mappings[region_name] <<- accession
      df_to_save <- data.frame(
        REGION = names(ncbi_mappings),
        Accession = as.character(ncbi_mappings),
        stringsAsFactors = FALSE
      )
      write.csv(df_to_save, file = mapping_cache_file, row.names = FALSE)
      return(accession)
    }
    return(NULL)
  }

  get_ncbi_sequence_local <- function(accession) {
    if (accession %in% names(cached_seqs)) {
      seq <- cached_seqs[[accession]]
      if (length(seq) > 0 && nchar(seq) <= 20000) {
        return(seq)
      }
    }
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=", accession, "&rettype=fasta")
    cat("Downloading NCBI reference sequence:", accession, "\n")
    tryCatch(
      {
        new_seq <- readDNAStringSet(url)
        if (length(new_seq) > 0 && nchar(new_seq[[1]]) <= 20000) {
          raw_header <- names(new_seq)[[1]]
          extracted_acc <- extract_accession(raw_header)
          names(new_seq) <- extracted_acc
          cached_seqs <<- c(cached_seqs[names(cached_seqs) != accession], new_seq)
          writeXStringSet(cached_seqs, file = cache_file)
          return(new_seq[[1]])
        } else {
          log_warning_error("Warning: Sequence for accession ", accession, " is empty or too large.")
          return(NULL)
        }
      },
      error = function(e) {
        warning("Could not download sequence for ", accession, ": ", e$message)
        return(NULL)
      }
    )
  }

  # Align and map coordinates
  for (region_name in unique(snv_df$REGION)) {
    accession <- get_ncbi_accession_for_region_local(region_name)
    if (is.null(accession)) next

    # Check local cache first to avoid pairwiseAlignment overhead
    alignment_cache_file <- file.path(alignments_cache_dir, paste0(region_name, "_vs_", accession, "_alignment_map.csv"))
    if (file.exists(alignment_cache_file)) {
      cat("Loading cached alignment mapping for", region_name, "vs", accession, "from local file...\n")
      map_df <- read.csv(alignment_cache_file, stringsAsFactors = FALSE)
      coord_mappings[[region_name]] <- as.integer(map_df$NCBI_POS)

      # Still need to load NCBI seq and CDS for downstream protein translation
      seq_ncbi <- get_ncbi_sequence_local(accession)
      if (!is.null(seq_ncbi)) {
        ncbi_seqs_loaded[[region_name]] <- seq_ncbi
        ncbi_cds_loaded[[region_name]] <- get_ncbi_cds(accession, cache_dir = cds_cache_dir)
      }
      next
    }

    match_idx <- match(toupper(gsub("\\.", "-", region_name)), toupper(gsub("\\.", "-", names(custom_seqs))))
    if (is.na(match_idx)) {
      warning("Genotype ", region_name, " not found in custom FASTA.")
      next
    }

    seq_ncbi <- get_ncbi_sequence_local(accession)
    if (is.null(seq_ncbi)) {
      match_fallback <- match(toupper(region_name), toupper(names(default_fallback_mappings)))
      if (!is.na(match_fallback)) {
        fallback_acc <- default_fallback_mappings[match_fallback]
        cat("Online search returned an invalid reference. Falling back to default for", region_name, ":", fallback_acc, "\n")
        seq_ncbi <- get_ncbi_sequence_local(fallback_acc)
        if (!is.null(seq_ncbi)) {
          accession <- fallback_acc
          ncbi_mappings[region_name] <- fallback_acc
          df_to_save <- data.frame(
            REGION = names(ncbi_mappings),
            Accession = as.character(ncbi_mappings),
            stringsAsFactors = FALSE
          )
          write.csv(df_to_save, file = mapping_cache_file, row.names = FALSE)
        }
      }
    }

    if (is.null(seq_ncbi)) next

    seq_custom <- custom_seqs[[match_idx]]
    ncbi_seqs_loaded[[region_name]] <- seq_ncbi
    ncbi_cds_loaded[[region_name]] <- get_ncbi_cds(accession, cache_dir = cds_cache_dir)

    cat("Aligning custom", region_name, "to NCBI", accession, "...\n")
    aln <- pairwiseAlignment(seq_custom, seq_ncbi, type = "global")

    mapping_vector <- build_coord_map(
      as.character(Biostrings::pattern(aln)),
      as.character(IRanges::subject(aln))
    )
    coord_mappings[[region_name]] <- mapping_vector

    # Save to local cache
    map_df <- data.frame(
      POS = seq_along(mapping_vector),
      NCBI_POS = as.integer(mapping_vector)
    )
    write.csv(map_df, file = alignment_cache_file, row.names = FALSE)
    cat("Saved alignment mapping cache for", region_name, "vs", accession, "to", alignment_cache_file, "\n")
  }

  # Build coordinate lookup table
  if (length(coord_mappings) > 0) {
    coord_df <- imap_dfr(coord_mappings, function(mapping, region) {
      tibble(REGION = region, POS = seq_along(mapping), NCBI_POS = as.integer(mapping))
    })

    snv_df <- snv_df %>%
      left_join(coord_df, by = c("REGION", "POS")) %>%
      mutate(NCBI_Accession = ncbi_mappings[REGION]) %>%
      group_by(REGION) %>%
      mutate(NCBI_REF = {
        seq_obj <- ncbi_seqs_loaded[[REGION[1]]]
        if (is.null(seq_obj)) {
          NA_character_
        } else {
          seq_str <- strsplit(as.character(seq_obj), "")[[1]]
          ifelse(!is.na(NCBI_POS) & NCBI_POS >= 1 & NCBI_POS <= length(seq_str),
            seq_str[NCBI_POS], NA_character_
          )
        }
      }) %>%
      ungroup() %>%
      mutate(
        NCBI_ALT = ALT,
        NCBI_Mutation = if_else(
          !is.na(NCBI_POS) & !is.na(NCBI_REF),
          paste0(NCBI_Accession, ":", NCBI_REF, NCBI_POS, NCBI_ALT),
          NA_character_
        )
      )
  } else {
    snv_df$NCBI_POS <- NA_integer_
    snv_df$NCBI_Accession <- NA_character_
    snv_df$NCBI_REF <- NA_character_
    snv_df$NCBI_ALT <- NA_character_
    snv_df$NCBI_Mutation <- NA_character_
  }

  # NCBI Amino Acid annotations
  snv_df <- snv_df %>%
    group_by(REGION) %>%
    mutate(
      NCBI_GENE = map_chr(NCBI_POS, function(pos) {
        cds <- ncbi_cds_loaded[[REGION[1]]]
        if (is.null(cds) || is.na(pos)) {
          return(NA_character_)
        }
        idx <- which(cds$start <= pos & cds$end >= pos)
        if (length(idx) > 0) {
          return(cds$product[idx[1]])
        }
        return(NA_character_)
      }),
      ncbi_aa_res = purrr::pmap(list(NCBI_POS, NCBI_ALT, NCBI_GENE), function(pos, alt, gene) {
        if (is.na(gene)) {
          return(list(pos_aa = NA_integer_, ref_aa = NA_character_, alt_aa = NA_character_))
        }
        cds <- ncbi_cds_loaded[[REGION[1]]]
        cds_row <- cds[cds$product == gene, ][1, ]
        seq_obj <- ncbi_seqs_loaded[[REGION[1]]]
        seq_str <- strsplit(as.character(seq_obj), "")[[1]]
        get_ncbi_aa_change(seq_str, cds_row$start, cds_row$end, pos, alt)
      }),
      NCBI_POS_AA = map_int(ncbi_aa_res, "pos_aa"),
      NCBI_REF_AA = map_chr(ncbi_aa_res, "ref_aa"),
      NCBI_ALT_AA = map_chr(ncbi_aa_res, "alt_aa"),
      NCBI_AA_Change = ifelse(
        !is.na(NCBI_POS_AA) & !is.na(NCBI_REF_AA) & !is.na(NCBI_ALT_AA),
        paste0(NCBI_GENE, ":", NCBI_REF_AA, NCBI_POS_AA, NCBI_ALT_AA),
        NA_character_
      ),
      PLOT_POS = as.integer(POS - min(POS, na.rm = TRUE) + 1),
      mutation = paste0(REF, PLOT_POS, ALT)
    ) %>%
    select(-ncbi_aa_res) %>%
    ungroup()

  total_samples <- length(unique(snv_df$sample_id))

  genome_div <- snv_df %>%
    group_by(REGION) %>%
    group_modify(~ {
      match_idx <- match(tolower(gsub("\\.", "-", .y$REGION)), tolower(gsub("\\.", "-", names(genome_lengths))))
      gl <- if (!is.na(match_idx)) genome_lengths[match_idx] else NA_integer_
      if (is.na(gl) || gl == 0) gl <- 7500L
      max_pos <- max(.x$PLOT_POS, na.rm = TRUE)
      if (max_pos <= WINDOW_SIZE) {
        return(tibble())
      }
      windows <- seq(1, max_pos - WINDOW_SIZE, by = STEP_SIZE)
      map_df(windows, function(w_start) {
        w_end <- w_start + WINDOW_SIZE
        window_data <- .x %>% filter(PLOT_POS >= w_start & PLOT_POS < w_end)
        tibble(
          Window_Mid = w_start + (WINDOW_SIZE / 2),
          Mean_Pi    = sum(window_data$site_pi, na.rm = TRUE) / (WINDOW_SIZE * total_samples),
          Tajima_D   = compute_tajima_d(window_data$ALT_FREQ, total_samples)
        )
      })
    }) %>%
    ungroup()

  # Diversity Plot
  if (nrow(genome_div) > 0) {
    tiff(file.path(output_dir, "Diversity_Across_Alignment.tiff"), units = "in", width = 12, height = 10, res = 300)
    genome_plot <- ggplot(genome_div, aes(x = Window_Mid, y = Mean_Pi, color = REGION)) +
      geom_area(aes(fill = REGION), alpha = 0.3) +
      geom_line(linewidth = 0.8) +
      facet_grid(REGION ~ ., scales = "fixed") +
      scale_color_manual(values = rep(cArray, length.out = length(unique(genome_div$REGION)))) +
      scale_fill_manual(values = rep(cArray, length.out = length(unique(genome_div$REGION)))) +
      theme_minimal() +
      labs(title = paste(virus_display, "-", city), x = "VP1 Position", y = "Nucleotide Diversity (π)") +
      theme(legend.position = "none", strip.text.y = element_text(angle = 0, size = 8))
    suppressWarnings(print(genome_plot))
    dev.off()
  }

  polymorphic_sites <- snv_df %>%
    group_by(REGION, PLOT_POS) %>%
    filter(any(ALT_FREQ > 0.05 & ALT_FREQ < 0.95)) %>%
    ungroup()

  # Mutation Dynamics
  if (nrow(polymorphic_sites) > 0) {
    mutation_dynamics <- polymorphic_sites %>%
      arrange(REGION, mutation, date_obj) %>%
      group_by(REGION, mutation) %>%
      summarise(
        dyn_class = classify_mutation_dynamics(ALT_FREQ),
        .groups = "drop"
      )

    polymorphic_sites <- polymorphic_sites %>%
      left_join(mutation_dynamics, by = c("REGION", "mutation"))

    dyn_colors <- c(
      emerging    = "#FF0000", persistent = "#0000FF", fixed   = "#000000",
      lost        = "#FF8C00", transient  = "#AAAAAA", other   = "#00CC00"
    )

    tiff(file.path(output_dir, "Mutation_Trajectories.tiff"), units = "in", width = 14, height = 12, res = 300)
    traj_plot <- ggplot(polymorphic_sites, aes(x = date_obj, y = ALT_FREQ, group = mutation, color = dyn_class)) +
      geom_line(alpha = 0.7, linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_wrap(~REGION, scales = "free") +
      scale_color_manual(values = dyn_colors, name = "Dynamics", na.value = "grey80") +
      theme_minimal() +
      labs(title = paste(virus_display, "-", city), x = "Date", y = "ALT Allele Frequency")
    suppressWarnings(print(traj_plot))
    dev.off()

    # Prophet Forecasting for Emerging Variants
    if (nrow(polymorphic_sites) > 0) {
      cat("  [Prophet] Forecasting trajectories for emerging mutations...\n")

      emerging_muts <- polymorphic_sites %>% filter(dyn_class == "emerging")
      prophet_plots <- list()

      for (mut in unique(emerging_muts$mutation)) {
        mut_data <- emerging_muts %>%
          filter(mutation == mut) %>%
          select(ds = date_obj, y = ALT_FREQ) %>%
          arrange(ds)

        if (nrow(mut_data) >= 2) {
          m <- suppressMessages(prophet(mut_data, yearly.seasonality = FALSE, weekly.seasonality = FALSE, daily.seasonality = FALSE))
          future <- make_future_dataframe(m, periods = 14)
          forecast <- predict(m, future)

          forecast_plot_df <- forecast %>%
            select(ds, yhat, yhat_lower, yhat_upper) %>%
            mutate(ds = as.Date(ds)) %>%
            left_join(mut_data, by = "ds")

          p_for <- ggplot(forecast_plot_df, aes(x = ds)) +
            geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper), fill = "red", alpha = 0.2) +
            geom_line(aes(y = yhat), color = "red", linetype = "dashed", linewidth = 1) +
            geom_point(aes(y = y), color = "black", size = 2) +
            theme_minimal() +
            labs(
              title = paste(virus_display, "-", city),
              subtitle = "Shaded region represents 80% confidence interval",
              x = "Date", y = "Predicted ALT Frequency"
            )

          prophet_plots[[mut]] <- p_for
        }
      }

      if (length(prophet_plots) > 0) {
        ivar_plots$forecasts <- prophet_plots
        pdf_path_prophet <- file.path(output_dir, "Prophet_Emerging_Forecasts.pdf")
        pdf(pdf_path_prophet, width = 10, height = 8)
        for (p in prophet_plots) suppressWarnings(print(p))
        dev.off()
        cat(sprintf("  => [Prophet] Forecasts generated and saved to: %s\n", pdf_path_prophet))
      } else {
        cat("  => [Prophet] No emerging mutations with sufficient data for forecasting.\n")
      }
    }
  }

  # Mean Nucleotide Diversity Over Time
  diversity_df <- snv_df %>%
    group_by(REGION, date_obj) %>%
    summarise(
      mean_pi = {
        match_idx <- match(tolower(gsub("\\.", "-", REGION[1])), tolower(gsub("\\.", "-", names(genome_lengths))))
        g_len <- if (!is.na(match_idx)) genome_lengths[match_idx] else NA_integer_
        if (is.na(g_len) || g_len == 0) g_len <- 7500L
        sum(site_pi, na.rm = TRUE) / (g_len * length(unique(sample_id)))
      },
      .groups = "drop"
    )

  if (nrow(diversity_df) > 0) {
    tiff(file.path(output_dir, "Mean_Nucleotide_Diversity.tiff"), units = "in", width = 10, height = 8, res = 300)
    div_plot <- ggplot(diversity_df, aes(x = date_obj, y = mean_pi, color = REGION)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_color_manual(values = rep(cArray, length.out = length(unique(diversity_df$REGION)))) +
      theme_bw() +
      labs(title = paste(virus_display, "-", city), x = "Date", y = "Mean Diversity (π)")
    suppressWarnings(print(div_plot))
    dev.off()
  }

  # Heatmaps
  col_fun <- colorRamp2(c(0, 0.5, 1), c("purple", "white", "green"))
  genotypes_to_plot <- unique(snv_df$REGION)

  for (gt in genotypes_to_plot) {
    gt_subset <- snv_df %>%
      filter(REGION == gt) %>%
      group_by(date_str, mutation, PLOT_POS) %>%
      summarise(ALT_FREQ = mean(ALT_FREQ), .groups = "drop")

    mutation_order <- gt_subset %>%
      select(mutation, PLOT_POS) %>%
      distinct() %>%
      arrange(PLOT_POS) %>%
      pull(mutation)

    gt_mat <- gt_subset %>%
      select(-PLOT_POS) %>%
      pivot_wider(names_from = mutation, values_from = ALT_FREQ, values_fill = NA_real_) %>%
      column_to_rownames("date_str") %>%
      as.matrix()

    existing_dates <- chronological_dates[chronological_dates %in% rownames(gt_mat)]
    gt_mat <- gt_mat[existing_dates, mutation_order, drop = FALSE]

    if (ncol(gt_mat) > 0 && nrow(gt_mat) > 0) {
      clean_gt <- gsub("[^[:alnum:]]", "_", gt)
      acc <- if (gt %in% names(ncbi_mappings)) ncbi_mappings[gt] else "NA"
      mean_dp <- round(mean(snv_df$TOTAL_DP[snv_df$REGION == gt], na.rm = TRUE))
      title_str <- paste0(gt, " (NCBI RefSeq: ", acc, ", n = ", mean_dp, ")")

      tiff(file.path(output_dir, paste(clean_gt, "Heatmap.tiff", sep = "_")), units = "in", width = 12, height = 10, res = 300)
      h <- Heatmap(gt_mat,
        name = "ALT Freq",
        col = col_fun,
        na_col = "grey",
        column_title = paste(virus_display, "-", city),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 9),
        border = TRUE
      )
      draw(h)
      dev.off()
      ivar_plots$heatmaps[[clean_gt]] <- h
    }
  }

  # pN/pS Selection Pressure
  pnps_plot_df <- snv_df %>%
    filter(var_type %in% c("synonymous", "nonsynonymous")) %>%
    group_by(REGION, date_obj) %>%
    summarise(
      n_S = sum(var_type == "synonymous"),
      n_NS = sum(var_type == "nonsynonymous"),
      pN_pS = ifelse(n_S > 0, n_NS / n_S, NA_real_),
      mean_depth = mean(TOTAL_DP, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(pN_pS))

  if (nrow(pnps_plot_df) > 0) {
    tiff(file.path(output_dir, "pN_pS_Over_Time.tiff"), units = "in", width = 12, height = 8, res = 300)
    pnps_plot <- ggplot(pnps_plot_df, aes(x = date_obj, y = pN_pS, color = REGION)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = rep(cArray, length.out = length(unique(pnps_plot_df$REGION)))) +
      theme_bw() +
      labs(
        title = paste(virus_display, "-", city),
        subtitle = "Dashed line = neutral expectation (pN/pS = 1). Above = positive selection, below = purifying.",
        x = "Date", y = "pN/pS Ratio"
      )
    suppressWarnings(print(pnps_plot))
    dev.off()
  }

  # Relative Abundance vs. pN/pS Linear Regression
  # Load the unfiltered Phase 1 abundance matrix to merge with pN/pS
  phase1_out_dir <- dirname(output_dir)
  abund_mat_file <- file.path(phase1_out_dir, paste0(city, "_", virus, "_abund_mat_unfiltered.csv"))

  if (file.exists(abund_mat_file) && nrow(pnps_plot_df) > 0) {
    cat("Merging Phase 1 abundance data with pN/pS for regression analysis...\n")
    abund_df <- read.csv(abund_mat_file, row.names = 1, check.names = FALSE)

    if (nrow(abund_df) > 0) {
      # Convert abundance matrix to long format and normalize REGION
      abund_long <- abund_df %>%
        dplyr::mutate(date_obj = as.Date(rownames(abund_df), format = "%d/%m/%Y")) %>%
        tidyr::pivot_longer(
          cols = -date_obj,
          names_to = "REGION",
          values_to = "Relative_Abundance"
        ) %>%
        dplyr::mutate(REGION = toupper(gsub("\\.", "-", REGION)))

      pnps_plot_df_norm <- pnps_plot_df %>%
        dplyr::mutate(REGION = toupper(gsub("\\.", "-", REGION)))

      # Merge datasets
      regression_df <- dplyr::inner_join(pnps_plot_df_norm, abund_long, by = c("date_obj", "REGION")) %>%
        dplyr::filter(!is.na(Relative_Abundance), !is.na(pN_pS), Relative_Abundance >= 1)

      if (nrow(regression_df) > 3) { # Need enough points for a regression
        # Save dataframe for combined analysis
        write.csv(regression_df, file.path(output_dir, paste0(city, "_", virus, "_regression_df.csv")), row.names = FALSE)

        # 1. Combined Plot (Regression controlling for Depth)
        fit_comb <- lm(pN_pS ~ Relative_Abundance + log10(mean_depth), data = regression_df)
        fit_red <- lm(pN_pS ~ log10(mean_depth), data = regression_df)
        r2_full <- summary(fit_comb)$r.squared
        r2_red <- summary(fit_red)$r.squared
        r2_comb <- round((r2_full - r2_red) / (1 - r2_red), 3)

        # Safely extract p-value for the abundance term
        coefs <- summary(fit_comb)$coefficients
        if ("Relative_Abundance" %in% rownames(coefs)) {
          pval_comb <- coefs["Relative_Abundance", "Pr(>|t|)"]
        } else {
          pval_comb <- NA_real_
        }
        pval_comb_str <- if (is.na(pval_comb)) "NA" else if (pval_comb < 0.001) "< 0.001" else sprintf("%.3f", pval_comb)
        subtitle_comb <- paste0("R² = ", r2_comb, " | p-value = ", pval_comb_str, " | n = ", nrow(regression_df))

        # Save a log of all model terms to Output for manuscript writing
        write.csv(as.data.frame(coefs), file.path(output_dir, "Abundance_vs_pN_pS_Regression_Terms.csv"))
        write.csv(data.frame(n = nrow(regression_df), R2 = r2_comb, pval = pval_comb), file.path(output_dir, "Abundance_vs_pN_pS_Summary.csv"), row.names = FALSE)
        write.csv(regression_df, file.path(output_dir, paste0(city, "_", virus, "_regression_df.csv")), row.names = FALSE)


        tiff(file.path(output_dir, "Abundance_vs_pN_pS_Regression.tiff"), units = "in", width = 10, height = 8, res = 300)
        reg_plot <- ggplot(regression_df, aes(x = Relative_Abundance, y = pN_pS)) +
          geom_point(aes(color = REGION), size = 3, alpha = 0.7) +
          geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
          scale_color_manual(values = rep(cArray, length.out = length(unique(regression_df$REGION)))) +
          theme_bw() +
          labs(
            title = paste(virus_display, "-", city),
            x = "Relative Genotype Abundance (%)",
            y = "pN/pS Ratio"
          ) +
          annotate("label", x = Inf, y = Inf, label = subtitle_comb, hjust = 1, vjust = 1, fill = "white", color = "black")
        suppressWarnings(print(reg_plot))
        dev.off()

        # 2. Faceted Grid Plot
        # Calculate stats for each genotype facet
        stats_df <- regression_df %>%
          dplyr::group_by(REGION) %>%
          dplyr::filter(dplyr::n() > 2, var(Relative_Abundance) > 0) %>%
          dplyr::summarize(
            n_obs = dplyr::n(),
            r2 = {
              r_full <- summary(lm(pN_pS ~ Relative_Abundance + log10(mean_depth)))$r.squared
              r_red <- summary(lm(pN_pS ~ log10(mean_depth)))$r.squared
              round((r_full - r_red) / (1 - r_red), 3)
            },
            pval = {
              coefs <- summary(lm(pN_pS ~ Relative_Abundance + log10(mean_depth)))$coefficients
              if ("Relative_Abundance" %in% rownames(coefs)) coefs["Relative_Abundance", 4] else NA_real_
            },
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            pval_str = ifelse(!is.na(pval) & pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
            label = paste0("R² = ", r2, "\np = ", pval_str, "\nn = ", n_obs)
          )

        tiff(file.path(output_dir, "Abundance_vs_pN_pS_Regression_Faceted.tiff"), units = "in", width = 12, height = 10, res = 300)
        reg_plot_faceted <- ggplot(regression_df, aes(x = Relative_Abundance, y = pN_pS)) +
          geom_point(aes(color = REGION), size = 3, alpha = 0.7) +
          geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
          facet_wrap(~REGION, scales = "free") +
          geom_label(data = stats_df, aes(x = Inf, y = Inf, label = label), hjust = 1, vjust = 1, size = 3.5, inherit.aes = FALSE, fill = "white", color = "black") +
          scale_color_manual(values = rep(cArray, length.out = length(unique(regression_df$REGION)))) +
          theme_bw() +
          theme(legend.position = "right") +
          labs(
            title = paste(virus_display, "-", city),
            x = "Relative Genotype Abundance (%)",
            y = "pN/pS Ratio",
            color = "Genotype"
          )
        suppressWarnings(print(reg_plot_faceted))
        dev.off()
      }
    }
  }

  # Tajima's D Plot
  tajima_plot_df <- genome_div %>% filter(!is.na(Tajima_D))
  if (nrow(tajima_plot_df) > 0) {
    tiff(file.path(output_dir, "Tajima_D_Across_Alignment.tiff"), units = "in", width = 12, height = 10, res = 300)
    tajima_plot <- ggplot(tajima_plot_df, aes(x = Window_Mid, y = Tajima_D)) +
      geom_ribbon(aes(ymin = pmin(Tajima_D, 0), ymax = pmax(Tajima_D, 0), fill = REGION), alpha = 0.3) +
      geom_line(aes(color = REGION), linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      facet_grid(REGION ~ ., scales = "free_y") +
      scale_color_manual(values = rep(cArray, length.out = length(unique(tajima_plot_df$REGION)))) +
      scale_fill_manual(values = rep(cArray, length.out = length(unique(tajima_plot_df$REGION)))) +
      theme_minimal() +
      labs(
        title = paste(virus_display, "-", city),
        subtitle = "D < 0: purifying selection / sweep. D > 0: balancing selection / structure.",
        x = "VP1 Position", y = "Tajima's D"
      ) +
      theme(legend.position = "none", strip.text.y = element_text(angle = 0, size = 8))
    suppressWarnings(print(tajima_plot))
    dev.off()
  }

  # Strand Bias Summary Bar Chart
  strand_plot_df <- snv_df %>%
    group_by(REGION, sample_id) %>%
    summarise(pct_flagged = 100 * sum(strand_bias_flag, na.rm = TRUE) / n(), .groups = "drop")

  if (nrow(strand_plot_df) > 0) {
    tiff(file.path(output_dir, "Strand_Bias_Summary.tiff"), units = "in", width = 12, height = 8, res = 300)
    strand_plot <- ggplot(strand_plot_df, aes(x = reorder(sample_id, pct_flagged), y = pct_flagged, fill = REGION)) +
      geom_col(width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = rep(cArray, length.out = length(unique(strand_plot_df$REGION)))) +
      theme_bw() +
      labs(
        title = paste(virus_display, "-", city),
        subtitle = paste0("Flagged = Fisher's exact strand-bias FDR-adjusted p < ", STRAND_BIAS_ALPHA),
        x = "Sample (Barcode)", y = "% Variants Flagged"
      )
    suppressWarnings(print(strand_plot))
    dev.off()
  }

  # Excel Output Workbook
  export_wb <- createWorkbook()
  for (gt in genotypes_to_plot) {
    gt_data <- snv_df %>%
      filter(REGION == gt) %>%
      group_by(date_str, mutation, PLOT_POS) %>%
      summarise(ALT_FREQ = mean(ALT_FREQ), .groups = "drop") %>%
      arrange(match(date_str, chronological_dates), PLOT_POS)

    if (nrow(gt_data) > 0) {
      gt_data <- gt_data %>%
        select(-PLOT_POS) %>%
        pivot_wider(names_from = mutation, values_from = ALT_FREQ, values_fill = NA_real_) %>%
        dplyr::rename(sample = date_str)

      clean_name <- gsub("[^[:alnum:]]", "_", substr(gt, 1, 31))
      addWorksheet(export_wb, sheetName = clean_name)
      writeData(export_wb, sheet = clean_name, x = gt_data)
    }
  }

  translation_map_df <- snv_df %>%
    select(REGION, POS, REF, ALT, mutation, NCBI_Accession, NCBI_POS, NCBI_REF, NCBI_ALT, NCBI_Mutation, NCBI_GENE, NCBI_POS_AA, NCBI_REF_AA, NCBI_ALT_AA, NCBI_AA_Change) %>%
    distinct() %>%
    filter(!is.na(NCBI_Mutation)) %>%
    arrange(REGION, POS)

  addWorksheet(export_wb, sheetName = "NCBI_Translation_Map")
  writeData(export_wb, sheet = "NCBI_Translation_Map", x = translation_map_df)

  pnps_df <- snv_df %>%
    filter(var_type %in% c("synonymous", "nonsynonymous")) %>%
    group_by(REGION, date_obj) %>%
    summarise(
      n_S = sum(var_type == "synonymous"),
      n_NS = sum(var_type == "nonsynonymous"),
      pN_pS = ifelse(n_S > 0, n_NS / n_S, NA_real_),
      .groups = "drop"
    )

  addWorksheet(export_wb, sheetName = "pN_pS")
  writeData(export_wb, sheet = "pN_pS", x = pnps_df)

  strand_summary <- snv_df %>%
    group_by(REGION, sample_id) %>%
    summarise(
      n_total = n(),
      n_flagged = sum(strand_bias_flag, na.rm = TRUE),
      pct_flagged = round(100 * n_flagged / n_total, 1),
      .groups = "drop"
    )

  addWorksheet(export_wb, sheetName = "Strand_Bias_Summary")
  writeData(export_wb, sheet = "Strand_Bias_Summary", x = strand_summary)

  addWorksheet(export_wb, sheetName = "Tajima_D")
  writeData(export_wb, sheet = "Tajima_D", x = genome_div)

  # Read Depth Summary Export in Excel
  if (!is.null(plot_df) && nrow(plot_df) > 0) {
    depth_summary_table <- plot_df %>%
      group_by(sample, REGION) %>%
      summarise(
        Mean_Depth = round(mean(TOTAL_DP, na.rm = TRUE), 1),
        Median_Depth = median(TOTAL_DP, na.rm = TRUE),
        Min_Depth = min(TOTAL_DP, na.rm = TRUE),
        Max_Depth = max(TOTAL_DP, na.rm = TRUE),
        Total_Variants = n(),
        .groups = "drop"
      )

    addWorksheet(export_wb, sheetName = "Read_Depth_Summary")
    writeData(export_wb, sheet = "Read_Depth_Summary", x = depth_summary_table)

    depth_tiff <- file.path(output_dir, "iVar_Read_Depth_Summary.tiff")
    if (file.exists(depth_tiff)) {
      insertImage(export_wb, sheet = "Read_Depth_Summary", file = depth_tiff, width = 10, height = 6.67, startRow = 2, startCol = 8, units = "in")
    }
  }

  if (!is.null(mutation_dynamics) && nrow(mutation_dynamics) > 0) {
    emerging_df <- snv_df %>%
      filter(var_type == "nonsynonymous") %>%
      left_join(mutation_dynamics, by = c("REGION", "mutation")) %>%
      filter(dyn_class == "emerging") %>%
      select(REGION, mutation, POS, REF_AA, ALT_AA, POS_AA, GFF_FEATURE, NCBI_Accession, NCBI_Mutation, NCBI_AA_Change, date_obj, ALT_FREQ) %>%
      arrange(REGION, mutation, date_obj)

    addWorksheet(export_wb, sheetName = "Emerging_Mutations")
    writeData(export_wb, sheet = "Emerging_Mutations", x = emerging_df)
  }

  convergent_df <- snv_df %>%
    filter(!is.na(POS_AA), var_type == "nonsynonymous") %>%
    mutate(cds_short = gsub("cds-", "", GFF_FEATURE)) %>%
    group_by(cds_short, POS_AA, REF_AA, ALT_AA) %>%
    summarise(
      n_genotypes = n_distinct(REGION),
      genotypes = paste(sort(unique(REGION)), collapse = ", "),
      NCBI_Accessions = paste(sort(unique(NCBI_Accession)), collapse = " / "),
      mean_freq = round(mean(ALT_FREQ, na.rm = TRUE), 4),
      n_timepoints = n_distinct(sample_id),
      .groups = "drop"
    ) %>%
    filter(n_genotypes >= 2) %>%
    arrange(desc(n_genotypes), desc(mean_freq))

  addWorksheet(export_wb, sheetName = "Convergent_Evolution")
  writeData(export_wb, sheet = "Convergent_Evolution", x = convergent_df)

  # Consensus protein sequence generation
  cat("Generating consensus protein sequences...\n")
  protein_seqs <- list()
  for (region in unique(snv_df$REGION)) {
    ncbi_seq_obj <- ncbi_seqs_loaded[[region]]
    cds_df <- ncbi_cds_loaded[[region]]
    if (is.null(ncbi_seq_obj) || is.null(cds_df) || nrow(cds_df) == 0) next

    seq_chars <- strsplit(as.character(ncbi_seq_obj), "")[[1]]
    ref_seq_str <- paste(seq_chars, collapse = "")

    for (i in seq_len(nrow(cds_df))) {
      s <- cds_df$start[i]
      e <- cds_df$end[i]
      gene_name <- cds_df$product[i]
      if (e <= nchar(ref_seq_str)) {
        cds_seq <- substr(ref_seq_str, s, e)
        prot <- suppressWarnings(tryCatch(
          {
            as.character(Biostrings::translate(Biostrings::DNAString(cds_seq)))
          },
          error = function(e) NA_character_
        ))
        if (!is.na(prot)) {
          # Cut the protein sequence to match the part we sequenced
          map_pos <- coord_mappings[[region]]
          if (!is.null(map_pos)) {
            valid_pos <- map_pos[!is.na(map_pos)]
            if (length(valid_pos) > 0) {
              start_nt <- min(valid_pos)
              end_nt <- max(valid_pos)

              # Calculate amino acid boundaries relative to the start of this CDS
              # If start_nt is before CDS start, we crop from AA 1.
              start_aa <- max(1, floor((start_nt - s) / 3) + 1)
              end_aa <- min(nchar(prot), ceiling((end_nt - s + 1) / 3))

              if (start_aa <= end_aa && start_aa <= nchar(prot)) {
                prot <- substr(prot, start_aa, end_aa)
              }

              # Clean up EV polyprotein name if necessary
              if (virus == "EV" && grepl("polyprotein", gene_name, ignore.case = TRUE)) {
                gene_name <- "VP1"
              }
            }
          }
          acc <- ncbi_mappings[region]
          header <- paste0("Reference_", region, "_", gsub("[^a-zA-Z0-9]", "_", gene_name), "_", acc)
          protein_seqs[[header]] <- prot
        }
      }
    }
  }

  for (samp in unique(snv_df$sample_id)) {
    samp_df <- snv_df %>% filter(sample_id == samp)
    for (region in unique(samp_df$REGION)) {
      ncbi_seq_obj <- ncbi_seqs_loaded[[region]]
      cds_df <- ncbi_cds_loaded[[region]]
      if (is.null(ncbi_seq_obj) || is.null(cds_df) || nrow(cds_df) == 0) next

      ncbi_seq_chars <- strsplit(as.character(ncbi_seq_obj), "")[[1]]
      custom_seq_obj <- custom_seqs[[region]]
      map_pos <- coord_mappings[[region]]

      if (!is.null(custom_seq_obj) && !is.null(map_pos)) {
        custom_seq_chars <- strsplit(as.character(custom_seq_obj), "")[[1]]
        maj_vars <- samp_df %>% filter(REGION == region, ALT_FREQ > 0.5)
        for (i in seq_len(nrow(maj_vars))) {
          row <- maj_vars[i, ]
          if (!is.na(row$POS) && row$POS >= 1 && row$POS <= length(custom_seq_chars)) {
            alt_base <- as.character(row$ALT)
            if (nchar(alt_base) == 1 && grepl("^[ACGTacgtNn]$", alt_base)) {
              custom_seq_chars[row$POS] <- alt_base
            }
          }
        }
        for (j in seq_along(map_pos)) {
          n_pos <- map_pos[j]
          if (!is.na(n_pos) && n_pos >= 1 && n_pos <= length(ncbi_seq_chars)) {
            ncbi_seq_chars[n_pos] <- custom_seq_chars[j]
          }
        }
      } else {
        maj_vars <- samp_df %>% filter(REGION == region, ALT_FREQ > 0.5)
        for (i in seq_len(nrow(maj_vars))) {
          row <- maj_vars[i, ]
          if (!is.na(row$NCBI_POS) && row$NCBI_POS >= 1 && row$NCBI_POS <= length(ncbi_seq_chars)) {
            alt_base <- as.character(row$NCBI_ALT)
            if (nchar(alt_base) == 1 && grepl("^[ACGTacgtNn]$", alt_base)) {
              ncbi_seq_chars[row$NCBI_POS] <- alt_base
            }
          }
        }
      }

      cons_seq_str <- paste(ncbi_seq_chars, collapse = "")
      for (i in seq_len(nrow(cds_df))) {
        s <- cds_df$start[i]
        e <- cds_df$end[i]
        gene_name <- cds_df$product[i]
        if (e <= nchar(cons_seq_str)) {
          cds_seq <- substr(cons_seq_str, s, e)
          prot <- suppressWarnings(tryCatch(
            {
              as.character(Biostrings::translate(Biostrings::DNAString(cds_seq)))
            },
            error = function(e) NA_character_
          ))
          if (!is.na(prot)) {
            # Cut the protein sequence to match the part we sequenced
            map_pos <- coord_mappings[[region]]
            if (!is.null(map_pos)) {
              valid_pos <- map_pos[!is.na(map_pos)]
              if (length(valid_pos) > 0) {
                start_nt <- min(valid_pos)
                end_nt <- max(valid_pos)

                start_aa <- max(1, floor((start_nt - s) / 3) + 1)
                end_aa <- min(nchar(prot), ceiling((end_nt - s + 1) / 3))

                if (start_aa <= end_aa && start_aa <= nchar(prot)) {
                  prot <- substr(prot, start_aa, end_aa)
                }

                if (virus == "EV" && grepl("polyprotein", gene_name, ignore.case = TRUE)) {
                  gene_name <- "VP1"
                }
              }
            }
            acc <- ncbi_mappings[region]
            header <- paste0(samp, "_", region, "_", gsub("[^a-zA-Z0-9]", "_", gene_name), "_", acc)
            protein_seqs[[header]] <- prot
          }
        }
      }
    }
  }

  if (length(protein_seqs) > 0) {
    prot_set <- AAStringSet(unlist(protein_seqs))
    writeXStringSet(prot_set, filepath = file.path(output_dir, "All_Samples_Consensus_Proteins.fasta"))
    # 3D Structure Pipeline
    struct_dir <- file.path(output_dir, "Structures")
    dir.create(struct_dir, showWarnings = FALSE, recursive = TRUE)

    cat("Generating 3D structures and computing RMSD via ColabFold & ChimeraX...\n")
    rmsd_results <- list()

    for (region in unique(snv_df$REGION)) {
      ref_headers <- names(protein_seqs)[grepl(paste0("^Reference_", region, "_"), names(protein_seqs))]
      if (length(ref_headers) == 0) next

      for (ref_h in ref_headers) {
        gene_suffix <- sub(paste0("^Reference_", region, "_"), "", ref_h)

        # Only perform structural alignment for VP1 / viral_protein_1
        # For EV, the whole sequence is VP1, so skip checking the gene name
        if (virus == "NV" && !grepl("VP1|viral_protein_1", gene_suffix, ignore.case = TRUE)) {
          next
        }

        ref_seq <- protein_seqs[[ref_h]]
        ref_seq_clean <- gsub("\\*", "", as.character(ref_seq))
        if (nchar(ref_seq_clean) > 2000) {
          cat("  -> Reference sequence (", nchar(ref_seq_clean), "aa) is too long for local CPU folding (>2000 aa). Skipping.\n")
          next
        }
        ref_pdb <- file.path(struct_dir, paste0(ref_h, ".pdb"))

        complex_ref_seq <- ref_seq_clean
        suffix <- "_VP1_only"
        if (virus == "NV") {
          complex_ref_seq <- paste(ref_seq_clean, ref_seq_clean, sep = ":")
          suffix <- "_complex"
        } else {
          vp_res <- fetch_vp_complex(region)
          if (!is.null(vp_res$VP2) && !is.null(vp_res$VP3)) {
            rec_seq <- fetch_receptor_sequence(region)
            if (!is.null(rec_seq)) {
              complex_ref_seq <- paste(vp_res$VP2, vp_res$VP3, ref_seq_clean, rec_seq, sep = ":")
              suffix <- "_complex_with_receptor"
            } else {
              complex_ref_seq <- paste(vp_res$VP2, vp_res$VP3, ref_seq_clean, sep = ":")
              suffix <- "_complex"
            }
          }
        }
        fasta_out_dir <- file.path(workspace_root, "Complex_FASTAs")
        dir.create(fasta_out_dir, showWarnings = FALSE)
        ref_h_suffix <- paste0(ref_h, suffix)
        writeLines(c(paste0(">", ref_h_suffix), complex_ref_seq), file.path(fasta_out_dir, paste0(ref_h_suffix, ".fasta")))

        samp_headers <- names(protein_seqs)[grepl(paste0("_", region, "_", gene_suffix, "$"), names(protein_seqs)) & !grepl("^Reference_", names(protein_seqs))]

        # ALWAYS generate sample FASTAs first!
        for (samp_h in samp_headers) {
          samp_seq <- protein_seqs[[samp_h]]
          samp_seq_clean <- gsub("\\*", "", as.character(samp_seq))
          if (nchar(samp_seq_clean) > 2000) next

          complex_samp_seq <- samp_seq_clean
          suffix_samp <- "_VP1_only"
          if (virus == "NV") {
            complex_samp_seq <- paste(samp_seq_clean, samp_seq_clean, sep = ":")
            suffix_samp <- "_complex"
          } else {
            vp_res <- fetch_vp_complex(region)
            if (!is.null(vp_res$VP2) && !is.null(vp_res$VP3)) {
              rec_seq <- fetch_receptor_sequence(region)
              if (!is.null(rec_seq)) {
                complex_samp_seq <- paste(vp_res$VP2, vp_res$VP3, samp_seq_clean, rec_seq, sep = ":")
                suffix_samp <- "_complex_with_receptor"
              } else {
                complex_samp_seq <- paste(vp_res$VP2, vp_res$VP3, samp_seq_clean, sep = ":")
                suffix_samp <- "_complex"
              }
            }
          }

          samp_h_suffix <- paste0(samp_h, suffix_samp)
          writeLines(c(paste0(">", samp_h_suffix), complex_samp_seq), file.path(fasta_out_dir, paste0(samp_h_suffix, ".fasta")))
        }

        # NOW check if we can perform structural alignment (requires the PDB models)
        if (!run_colabfold(complex_ref_seq, struct_dir, ref_h_suffix, COLABFOLD_BIN)) next

        for (samp_h in samp_headers) {
          samp_seq <- protein_seqs[[samp_h]]
          samp_seq_clean <- gsub("\\*", "", as.character(samp_seq))
          if (nchar(samp_seq_clean) > 2000) next

          highest_risk_score <- NA
          # ESM-2 Triage Step
          samp_id_raw <- sub("_.*", "", samp_h)
          samp_muts <- snv_df %>%
            filter(sample_id == samp_id_raw, REGION == region, var_type == "nonsynonymous", !is.na(NCBI_POS_AA))

          skip_colabfold <- FALSE
          if (RUN_COLABFOLD_LOCAL && ESM_READY && nrow(samp_muts) > 0) {
            highest_risk_score <- 0.0
            for (i in seq_len(nrow(samp_muts))) {
              pos <- as.integer(samp_muts$NCBI_POS_AA[i])
              mut_aa <- samp_muts$NCBI_ALT_AA[i]
              risk_score <- calculate_mutation_risk(ref_seq_clean, pos, mut_aa)
              if (risk_score < highest_risk_score) highest_risk_score <- risk_score
            }

            if (highest_risk_score >= -2.0) {
              cat(sprintf("  [ESM-2] Skipping ColabFold for %s. Mutations tolerated (Score: %.2f)\n", samp_h, highest_risk_score))
              skip_colabfold <- TRUE
            } else {
              cat(sprintf("  [ESM-2] High risk mutation detected in %s (Score: %.2f). Proceeding to ColabFold...\n", samp_h, highest_risk_score))
            }
          }

          if (skip_colabfold) next
          # End Triage Step

          samp_pdb <- file.path(struct_dir, paste0(samp_h, ".pdb"))

          suffix_samp <- "_VP1_only"
          if (virus == "NV") {
            suffix_samp <- "_complex"
          } else {
            vp_res <- fetch_vp_complex(region)
            if (!is.null(vp_res$VP2) && !is.null(vp_res$VP3)) {
              rec_seq <- fetch_receptor_sequence(region)
              if (!is.null(rec_seq)) {
                suffix_samp <- "_complex_with_receptor"
              } else {
                suffix_samp <- "_complex"
              }
            }
          }
          samp_h_suffix <- paste0(samp_h, suffix_samp)

          if (!run_colabfold(complex_samp_seq, struct_dir, samp_h_suffix, COLABFOLD_BIN)) next

          out_cxs <- file.path(struct_dir, paste0(samp_h, "_alignment.cxs"))
          out_log <- file.path(struct_dir, paste0(samp_h, "_rmsd.log"))

          ref_chars <- strsplit(ref_seq_clean, "")[[1]]
          samp_chars <- strsplit(samp_seq_clean, "")[[1]]
          min_len <- min(length(ref_chars), length(samp_chars))
          mut_pos <- which(ref_chars[1:min_len] != samp_chars[1:min_len])

          if (run_chimerax_alignment(ref_pdb, samp_pdb, out_cxs, out_log, CHIMERAX_BIN, mutated_residues = mut_pos)) {
            log_lines <- tryCatch(readLines(out_log, warn = FALSE), error = function(e) NULL)
            rmsd_line <- grep("RMSD between", log_lines, value = TRUE)

            ref_iptm <- NA
            ref_plddt <- NA
            samp_iptm <- NA
            samp_plddt <- NA

            ref_json <- file.path(struct_dir, paste0(ref_h_suffix, ".json"))
            if (file.exists(ref_json)) {
              jdata <- tryCatch(jsonlite::fromJSON(ref_json), error = function(e) list())
              if (!is.null(jdata$iptm)) ref_iptm <- as.numeric(jdata$iptm)
              if (!is.null(jdata$plddts)) ref_plddt <- mean(unlist(jdata$plddts), na.rm = TRUE)
            }

            samp_json <- file.path(struct_dir, paste0(samp_h_suffix, ".json"))
            if (file.exists(samp_json)) {
              jdata <- tryCatch(jsonlite::fromJSON(samp_json), error = function(e) list())
              if (!is.null(jdata$iptm)) samp_iptm <- as.numeric(jdata$iptm)
              if (!is.null(jdata$plddts)) {
                samp_plddt <- mean(unlist(jdata$plddts), na.rm = TRUE)
                cat(sprintf(
                  "  [pLDDT] %s: mean pLDDT = %.1f (%s)\n", samp_h, samp_plddt,
                  ifelse(samp_plddt >= MIN_PLDDT_STRUCTURAL, "PASS", "LOW CONFIDENCE")
                ))
              }
            }

            if (length(rmsd_line) > 0) {
              rmsd_val <- as.numeric(sub(".*is\\s+([0-9.]+)\\s+angstroms.*", "\\1", rmsd_line[1]))
              rmsd_results[[samp_h]] <- data.frame(
                Sample_Protein = samp_h,
                Region = region,
                Gene = gene_suffix,
                RMSD = rmsd_val,
                Ref_ipTM = ref_iptm,
                Sample_ipTM = samp_iptm,
                Ref_pLDDT = ref_plddt,
                Sample_pLDDT = samp_plddt,
                ESM2_Score = highest_risk_score,
                CXS_Path = out_cxs,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }

    if (length(rmsd_results) > 0) {
      rmsd_df <- dplyr::bind_rows(rmsd_results)
      rmsd_df <- rmsd_df[order(-rmsd_df$RMSD), ]

      cat("\n=======================================================\n")
      non_zero <- rmsd_df[rmsd_df$RMSD > 0, ]
      if (nrow(non_zero) > 0) {
        cat("Samples with Structural Deviations (Highest RMSD First):\n")
        print(head(non_zero[, c("Sample_Protein", "RMSD", "Ref_ipTM", "Sample_ipTM", "Sample_pLDDT", "CXS_Path")], 20))
      } else {
        cat("All aligned samples had 0.000 RMSD (no deviations from reference).\n")
      }
      cat("=======================================================\n")

      csv_path <- file.path(struct_dir, "Structural_Deviation_Batch.csv")
      write.csv(rmsd_df, csv_path, row.names = FALSE)
      cat("\n  => Saved local results to", csv_path, "\n")


      addWorksheet(export_wb, sheetName = "Structural_Deviation")
      writeData(export_wb, sheet = "Structural_Deviation", x = rmsd_df[, c("Sample_Protein", "Region", "Gene", "RMSD", "Sample_pLDDT", "Ref_pLDDT")])
      ivar_plots$structures <- rmsd_df

      # Generate ipTM binding affinity plot if data is available
      iptm_df <- rmsd_df[!is.na(rmsd_df$Ref_ipTM) & !is.na(rmsd_df$Sample_ipTM), ]
      if (nrow(iptm_df) > 0) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          text_geom <- ggrepel::geom_text_repel(aes(label = Sample_Protein), size = 3.5, max.overlaps = 15)
        } else {
          text_geom <- geom_text(aes(label = Sample_Protein), size = 3.5, vjust = -1)
        }

        p_iptm <- ggplot(iptm_df, aes(x = Ref_ipTM, y = Sample_ipTM, color = Region)) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
          geom_point(size = 4, alpha = 0.8) +
          text_geom +
          labs(
            x = "Reference ipTM (Baseline Binding Confidence)",
            y = "Variant ipTM (Mutated Binding Confidence)",
            title = paste(virus_display, "-", city),
            subtitle = "Points ABOVE the dashed line indicate variants with potentially INCREASED receptor binding affinity.",
            color = "Genotype"
          ) +
          theme_bw(base_size = 14)

        ivar_plots$iptm <- p_iptm

        p_rmsd_iptm <- ggplot(iptm_df, aes(x = RMSD, y = Sample_ipTM, color = Region)) +
          geom_point(size = 4, alpha = 0.8) +
          text_geom +
          labs(
            x = "RMSD (Structural Deviation in Angstroms)",
            y = "Variant ipTM (Mutated Binding Confidence)",
            title = paste(virus_display, "-", city),
            subtitle = "Comparing magnitude of structural change (RMSD) to predicted receptor binding confidence.",
            color = "Genotype"
          ) +
          theme_bw(base_size = 14)

        ivar_plots$rmsd_iptm <- p_rmsd_iptm
      }

      # Generate pLDDT Summary Bar Plot
      plddt_df <- rmsd_df[!is.na(rmsd_df$Sample_pLDDT), ]
      if (nrow(plddt_df) > 0) {
        p_plddt_summary <- ggplot(plddt_df, aes(x = reorder(Sample_Protein, -Sample_pLDDT), y = Sample_pLDDT, fill = Region)) +
          geom_col(color = "black", alpha = 0.8) +
          geom_hline(yintercept = MIN_PLDDT_STRUCTURAL, linetype = "dashed", color = "red", linewidth = 1) +
          labs(
            x = "Sample (Barcode)",
            y = "Mean pLDDT (Model Confidence)",
            title = paste(virus_display, "-", city),
            subtitle = "Dashed red line indicates minimum confidence threshold.",
            fill = "Genotype"
          ) +
          theme_bw(base_size = 14) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ivar_plots$plddt_summary <- p_plddt_summary
      }

      # Generate ipTM Summary Bar Plot
      iptm_df_bar <- rmsd_df[!is.na(rmsd_df$Sample_ipTM), ]
      if (nrow(iptm_df_bar) > 0) {
        p_iptm_summary <- ggplot(iptm_df_bar, aes(x = reorder(Sample_Protein, -Sample_ipTM), y = Sample_ipTM, fill = Region)) +
          geom_col(color = "black", alpha = 0.8) +
          labs(
            x = "Sample (Barcode)",
            y = "Predicted Binding Affinity (ipTM)",
            title = paste(virus_display, "-", city),
            fill = "Genotype"
          ) +
          theme_bw(base_size = 14) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ivar_plots$iptm_summary <- p_iptm_summary
      }
    }
  }

  # Save DNA consensus sequences and construct phylogenetic trees
  dna_seqs <- list()
  for (region in unique(snv_df$REGION)) {
    ncbi_seq_obj <- ncbi_seqs_loaded[[region]]
    if (is.null(ncbi_seq_obj)) next

    # Save reference sequence
    acc <- ncbi_mappings[region]
    ref_header <- paste0("Reference_", region, "_", acc)
    dna_seqs[[ref_header]] <- as.character(ncbi_seq_obj)
  }

  for (samp in unique(snv_df$sample_id)) {
    samp_df <- snv_df %>% filter(sample_id == samp)
    for (region in unique(samp_df$REGION)) {
      ncbi_seq_obj <- ncbi_seqs_loaded[[region]]
      if (is.null(ncbi_seq_obj)) next

      ncbi_seq_chars <- strsplit(as.character(ncbi_seq_obj), "")[[1]]
      custom_seq_obj <- custom_seqs[[region]]
      map_pos <- coord_mappings[[region]]

      if (!is.null(custom_seq_obj) && !is.null(map_pos)) {
        custom_seq_chars <- strsplit(as.character(custom_seq_obj), "")[[1]]
        maj_vars <- samp_df %>% filter(REGION == region, ALT_FREQ > 0.5)
        for (i in seq_len(nrow(maj_vars))) {
          row <- maj_vars[i, ]
          if (!is.na(row$POS) && row$POS >= 1 && row$POS <= length(custom_seq_chars)) {
            alt_base <- as.character(row$ALT)
            if (nchar(alt_base) == 1 && grepl("^[ACGTacgtNn]$", alt_base)) {
              custom_seq_chars[row$POS] <- alt_base
            }
          }
        }
        for (j in seq_along(map_pos)) {
          n_pos <- map_pos[j]
          if (!is.na(n_pos) && n_pos >= 1 && n_pos <= length(ncbi_seq_chars)) {
            ncbi_seq_chars[n_pos] <- custom_seq_chars[j]
          }
        }
      } else {
        maj_vars <- samp_df %>% filter(REGION == region, ALT_FREQ > 0.5)
        for (i in seq_len(nrow(maj_vars))) {
          row <- maj_vars[i, ]
          if (!is.na(row$NCBI_POS) && row$NCBI_POS >= 1 && row$NCBI_POS <= length(ncbi_seq_chars)) {
            alt_base <- as.character(row$NCBI_ALT)
            if (nchar(alt_base) == 1 && grepl("^[ACGTacgtNn]$", alt_base)) {
              ncbi_seq_chars[row$NCBI_POS] <- alt_base
            }
          }
        }
      }
      cons_seq_str <- paste(ncbi_seq_chars, collapse = "")
      acc <- ncbi_mappings[region]
      header <- paste0(samp, "_", region, "_", acc)
      dna_seqs[[header]] <- cons_seq_str
    }
  }

  if (length(dna_seqs) > 0) {
    dna_set <- DNAStringSet(unlist(dna_seqs))
    writeXStringSet(dna_set, filepath = file.path(output_dir, "All_Samples_Consensus_DNA.fasta"))
  }

  # Lollipop plots
  lollipop_base <- snv_df %>%
    filter(
      ALT_FREQ > 0.05, !is.na(GFF_FEATURE), !is.na(POS_AA), POS_AA != "NA",
      var_type %in% c("synonymous", "nonsynonymous", "stop_gain", "stop_loss")
    ) %>%
    mutate(
      cds_short = gsub("cds-", "", GFF_FEATURE),
      POS_AA_num = as.integer(POS_AA)
    )

  if (nrow(lollipop_base) > 0) {
    region_labels <- snv_df %>%
      group_by(REGION) %>%
      summarise(mean_dp_genotype = round(mean(TOTAL_DP, na.rm = TRUE)), .groups = "drop") %>%
      mutate(REGION_LABEL = paste0(REGION, " (n=", mean_dp_genotype, "x)"))

    lollipop_base <- lollipop_base %>% left_join(region_labels, by = "REGION")

    lollipop_data <- lollipop_base %>%
      group_by(REGION_LABEL, cds_short, POS_AA_num, var_type) %>%
      summarise(mean_freq = mean(ALT_FREQ, na.rm = TRUE), locus_mean_dp = mean(TOTAL_DP, na.rm = TRUE), .groups = "drop")

    gff_accessions <- unique(lollipop_data$cds_short)
    gene_name_lookup <- setNames(gff_accessions, gff_accessions)

    for (acc in gff_accessions) {
      tryCatch(
        {
          url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=", acc, "&rettype=gp&retmode=text")
          lines <- readLines(url, warn = FALSE)
          prod_line <- grep("product=", lines, value = TRUE)[1]
          if (!is.na(prod_line)) {
            gene_name <- sub('.*product="([^"]+)".*', "\\1", prod_line)
            gene_name_lookup[acc] <- gene_name
          }
        },
        error = function(e) {}
      )
    }

    var_type_shapes <- c(nonsynonymous = 16, synonymous = 1, stop_gain = 17, stop_loss = 15)

    for (gene_acc in gff_accessions) {
      gene_data <- lollipop_data %>% filter(cds_short == gene_acc)
      if (nrow(gene_data) == 0) next
      gene_name <- gene_name_lookup[gene_acc]
      clean_name <- gsub("[^[:alnum:]_]", "_", gene_name)
      clean_acc <- gsub("[^[:alnum:]_]", "_", gene_acc)

      display_name <- paste0(gene_name, " [", gene_acc, "]")
      file_name <- file.path(output_dir, paste0("Lollipop_", clean_name, "_", clean_acc, ".tiff"))

      max_aa <- max(gene_data$POS_AA_num, na.rm = TRUE)

      tiff(file_name, units = "in", width = 18, height = 8, res = 300)
      p_lol <- ggplot(gene_data, aes(x = POS_AA_num, y = mean_freq, color = REGION_LABEL)) +
        geom_segment(aes(x = POS_AA_num, xend = POS_AA_num, y = 0, yend = mean_freq), color = "grey60", linewidth = 0.5) +
        geom_point(aes(size = locus_mean_dp, shape = var_type), alpha = 0.8) +
        scale_x_continuous(breaks = seq(0, max_aa, by = 500), minor_breaks = seq(0, max_aa, by = 100)) +
        scale_color_manual(values = rep(cArray, length.out = length(unique(gene_data$REGION_LABEL)))) +
        scale_shape_manual(values = var_type_shapes, name = "Variant Type") +
        theme_bw() +
        labs(
          title = paste(virus_display, "-", city),
          x = "Amino Acid Position", y = "Mean ALT Frequency", color = "Genotype", size = "Locus Read Depth"
        ) +
        theme(legend.position = "bottom", panel.grid.minor.x = element_line(color = "grey85", linewidth = 0.3))
      suppressWarnings(print(p_lol))
      dev.off()
      ivar_plots$lollipops[[clean_acc]] <- p_lol
    }
  }

  saveWorkbook(export_wb, file = file.path(output_dir, "Genotype_Diversity_Analysis_Final_v8.xlsx"), overwrite = TRUE)
  cat("Completed iVar pipeline on suite:", suite_name, "\n")

  if (!is.null(p_depth)) ivar_plots$read_depth <- p_depth
  if (!is.null(genome_plot)) ivar_plots$diversity_across_alignment <- genome_plot
  if (!is.null(traj_plot)) ivar_plots$mutation_trajectories <- traj_plot
  if (!is.null(div_plot)) ivar_plots$mean_diversity_over_time <- div_plot
  if (!is.null(pnps_plot)) ivar_plots$pn_ps_over_time <- pnps_plot
  if (!is.null(reg_plot)) ivar_plots$abundance_vs_pnps <- reg_plot
  if (!is.null(reg_plot_faceted)) ivar_plots$abundance_vs_pnps_faceted <- reg_plot_faceted
  if (!is.null(tajima_plot)) ivar_plots$tajima_d <- tajima_plot
  if (!is.null(strand_plot)) ivar_plots$strand_bias <- strand_plot

  return(ivar_plots)
}

# Combined Analysis
# PHASE 3: Combine data from all locations
# Run Beta-Diversity tests and make PCoA plots to compare differences over time and location
run_combined_analysis <- function(virus, output_dirs, combined_out_dir, cArray) {
  cat("\n==================================================\n")
  cat("Running COMBINED analysis for virus:", virus, "\n")
  cat("Output directory:", combined_out_dir, "\n")
  cat("==================================================\n")

  dir.create(combined_out_dir, showWarnings = FALSE, recursive = TRUE)

  virus_full_name <- if (virus == "EV") "Enterovirus" else "Norovirus"
  fit_all <- NULL

  abund_files <- c()
  bc_files <- c()
  cities <- c()

  for (out_dir in output_dirs) {
    base_dir_name <- basename(out_dir)
    city <- sub("^(EV|Noro|NV)_", "", base_dir_name)

    abund_f <- file.path(out_dir, paste0(city, "_", virus, "_abund_mat.csv"))
    bc_f <- file.path(out_dir, paste0(city, "_", virus, "_bc_scores.csv"))

    if (file.exists(abund_f) && file.exists(bc_f)) {
      abund_files <- c(abund_files, abund_f)
      bc_files <- c(bc_files, bc_f)
      cities <- c(cities, city)
    }
  }

  if (length(abund_files) < 2) {
    cat("Skipping combined analysis for", virus, ": fewer than 2 cities with abundance matrix files found.\n")
    return(NULL)
  }

  abund_list <- lapply(seq_along(abund_files), function(i) {
    df <- read.csv(abund_files[i])
    colnames(df)[1] <- "date"
    df$city <- cities[i]
    df
  })

  combined <- bind_rows(abund_list)
  combined[is.na(combined)] <- 0
  combined$date <- as.Date(combined$date, format = "%d/%m/%Y")
  combined <- combined[order(combined$date), ]
  combined$date_str <- format(combined$date, "%d/%m/%Y")
  combined$sample_id <- paste(combined$city, combined$date_str, sep = " - ")

  exclude_cols <- c("date", "city", "date_str", "sample_id", "Other", "Could.not.assign")
  genotype_cols <- setdiff(colnames(combined), exclude_cols)

  meta_list <- lapply(seq_along(bc_files), function(i) {
    df <- read.csv(bc_files[i])
    df$city <- cities[i]
    df$sample_date <- as.Date(df$sample_id, format = "%d/%m/%Y")
    df$sample_id <- paste(cities[i], format(df$sample_date, "%d/%m/%Y"), sep = " - ")
    df$collection_year <- format(df$sample_date, "%Y")
    df$month_year <- format(df$sample_date, "%Y-%m")

    keep_cols <- c("sample_id", "season", "city", "collection_year", "month_year")
    df[, intersect(colnames(df), keep_cols)]
  })

  meta_combined <- do.call(rbind, meta_list)
  meta_combined <- meta_combined[match(combined$sample_id, meta_combined$sample_id), ]

  abund_matrix <- combined[, genotype_cols]
  bc_dist <- vegdist(abund_matrix, method = "bray")

  bc_pcoa <- cmdscale(bc_dist, eig = TRUE, k = 2)
  pcoa_scores <- as.data.frame(bc_pcoa$points)
  colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
  pcoa_scores$sample_id <- combined$sample_id

  pcoa_data <- dplyr::inner_join(pcoa_scores, meta_combined, by = "sample_id")

  pcoa_data$temp_date <- as.Date(gsub("^.* - ", "", pcoa_data$sample_id), format = "%d/%m/%Y")
  pcoa_data <- pcoa_data[order(pcoa_data$city, pcoa_data$temp_date), ]
  pcoa_data$sample_id <- factor(pcoa_data$sample_id, levels = unique(pcoa_data$sample_id))
  pcoa_data$temp_date <- NULL

  pcoa1_var <- round(100 * bc_pcoa$eig[1] / sum(bc_pcoa$eig), 1)
  pcoa2_var <- round(100 * bc_pcoa$eig[2] / sum(bc_pcoa$eig), 1)

  cat("\n=== Combined Statistical Testing Results ===\n")
  candidate_vars <- c("city", "season", "collection_year", "month_year")
  valid_vars <- c()

  for (v in candidate_vars) {
    if (v %in% colnames(pcoa_data)) {
      non_na_vals <- na.omit(pcoa_data[[v]])
      unique_count <- length(unique(non_na_vals))
      if (unique_count >= 2) {
        valid_vars <- c(valid_vars, v)
      } else {
        cat("Skipping variable", v, "for PERMANOVA/betadisper: it has only", unique_count, "unique level(s).\n")
      }
    }
  }

  permanova_results_list <- list()
  combined_vars <- setdiff(valid_vars, "month_year")

  if (length(combined_vars) > 0) {
    cat("\n--- Combined PERMANOVA (adonis2) Results (excluding month_year) ---\n")
    set.seed(42)
    formula_str <- as.formula(paste("bc_dist ~", paste(combined_vars, collapse = " + ")))
    fit_all <- adonis2(formula_str, data = pcoa_data, permutations = 1000, by = "margin")
    print(fit_all)

    cat("\n--- Independent PERMANOVA tests ---\n")
    for (v in valid_vars) {
      cat("\nPERMANOVA for variable:", v, "\n")
      ind_formula <- as.formula(paste("bc_dist ~", v))
      res <- adonis2(ind_formula, data = pcoa_data, permutations = 1000)
      print(res)

      df_res <- as.data.frame(res)
      df_res$Variable <- rownames(df_res)
      df_res <- df_res[, c(ncol(df_res), 1:(ncol(df_res) - 1))]
      permanova_results_list[[v]] <- df_res
    }
  } else {
    cat("\nNo metadata variables had enough levels (>= 2) to perform PERMANOVA.\n")
  }

  cat("\n--- Beta-Dispersion (betadisper) Results ---\n")
  valid_categorical <- intersect(candidate_vars, valid_vars)

  for (v in valid_categorical) {
    cat("\nDispersion by", v, ":\n")
    non_na_idx <- !is.na(pcoa_data[[v]])

    tryCatch(
      {
        if (sum(non_na_idx) == nrow(pcoa_data)) {
          group_factor <- as.factor(pcoa_data[[v]])
          bd_res <- betadisper(bc_dist, group_factor)
        } else {
          sub_dist <- as.dist(as.matrix(bc_dist)[non_na_idx, non_na_idx])
          sub_group <- factor(pcoa_data[[v]][non_na_idx])
          bd_res <- betadisper(sub_dist, sub_group)
        }
        print(anova(bd_res))
        print(permutest(bd_res))
      },
      error = function(e) {
        cat("Failed to run betadisper for variable", v, ":", e$message, "\n")
      }
    )
  }

  get_perm_stats_local <- function(var_name, results_list) {
    if (var_name %in% names(results_list)) {
      df <- results_list[[var_name]]
      if ("R2" %in% colnames(df) && "Pr(>F)" %in% colnames(df)) {
        r2_val <- round(df$R2[1], 3)
        p_val <- df[["Pr(>F)"]][1]
        p_str <- if (is.na(p_val)) "p = NA" else if (p_val < 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
        return(paste0("PERMANOVA: R² = ", r2_val, ", ", p_str))
      }
    }
    return(NULL)
  }

  temp_plots <- list()
  plot_filenames <- list()

  plot_vars <- intersect(candidate_vars, colnames(pcoa_data))

  for (v in plot_vars) {
    clean_var_name <- tolower(v)

    plot_subtitle <- dplyr::case_when(
      v == "city" ~ "Points show individual samples; ellipses group samples by City",
      v == "season" ~ "Samples grouped/colored by sampling season",
      v == "collection_year" ~ "Samples grouped/colored by sample collection year",
      v == "month_year" ~ "Samples grouped/colored by sample month and year",
      TRUE ~ paste0("Samples analyzed by ", gsub("_", " ", v))
    )

    if (v == "city") {
      season_shapes <- c("Spring" = 21, "Summer" = 22, "Fall" = 23, "Winter" = 24)


      # Identify top 5% outliers for this grouping variable
      pcoa_data$outlier_label <- NA
      valid_idx <- !is.na(pcoa_data$city)
      if (sum(valid_idx) > 0) {
        group_cents <- pcoa_data[valid_idx, ] %>% 
          dplyr::group_by(city) %>% 
          dplyr::summarize(cx = mean(PCoA1), cy = mean(PCoA2), .groups = "drop")
        
        p_dist_data <- pcoa_data[valid_idx, ] %>% 
          dplyr::left_join(group_cents, by = "city") %>%
          dplyr::mutate(dist_to_cent = sqrt((PCoA1 - cx)^2 + (PCoA2 - cy)^2))
          
        dist_thresh <- quantile(p_dist_data$dist_to_cent, 0.95, na.rm = TRUE)
        pcoa_data$outlier_label[valid_idx] <- ifelse(p_dist_data$dist_to_cent > dist_thresh, 
                                                     as.character(p_dist_data$sample_id), NA)
      }

      p_temp <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, fill = city, shape = season)) +
        geom_label_repel(aes(label = outlier_label), max.overlaps = 50, box.padding = 0.5, fill = alpha("white", 0.7), color = "black", size = 3, show.legend = FALSE, na.rm = TRUE) +
        geom_point(size = 5, alpha = 0.8, stroke = 0.8, color = "black") +
        scale_fill_manual(values = cArray, name = "City") +
        scale_shape_manual(values = season_shapes, name = "Season") +
        labs(
          x = paste0("PCoA1 (", pcoa1_var, "%)"),
          y = paste0("PCoA2 (", pcoa2_var, "%)"),
          title = paste(virus_full_name, "-", paste(sort(unique(combined$city)), collapse = " and ")),
          subtitle = plot_subtitle
        ) +
        theme_bw(base_size = 14) +
        theme(legend.position = "right")

      if ("city" %in% colnames(pcoa_data) && length(unique(pcoa_data$city)) >= 2) {
        p_temp <- p_temp +
          stat_ellipse(aes(x = PCoA1, y = PCoA2, group = city, color = city),
            inherit.aes = FALSE, type = "t", linetype = 2, linewidth = 0.6, level = 0.95
          ) +
          scale_color_manual(values = cArray, name = "City")
      }

      plot_filenames[[v]] <- paste0(virus_full_name, "_pcoa_city_samples")
    } else {

      legend_name <- if (v == "season") "Season" else gsub("_", " ", v)

      # Identify top 5% outliers for this grouping variable
      pcoa_data$outlier_label <- NA
      valid_idx <- !is.na(pcoa_data[[v]])
      if (sum(valid_idx) > 0) {
        group_cents <- pcoa_data[valid_idx, ] %>% 
          dplyr::group_by(.data[[v]]) %>% 
          dplyr::summarize(cx = mean(PCoA1), cy = mean(PCoA2), .groups = "drop")
        
        p_dist_data <- pcoa_data[valid_idx, ] %>% 
          dplyr::left_join(group_cents, by = v) %>%
          dplyr::mutate(dist_to_cent = sqrt((PCoA1 - cx)^2 + (PCoA2 - cy)^2))
          
        dist_thresh <- quantile(p_dist_data$dist_to_cent, 0.95, na.rm = TRUE)
        pcoa_data$outlier_label[valid_idx] <- ifelse(p_dist_data$dist_to_cent > dist_thresh, 
                                                     as.character(p_dist_data$sample_id), NA)
      }

      p_temp <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, fill = .data[[v]], shape = season)) +
        geom_label_repel(aes(label = outlier_label), max.overlaps = 50, box.padding = 0.5, fill = alpha("white", 0.7), color = "black", size = 3, show.legend = FALSE, na.rm = TRUE) +
        geom_point(size = 5, alpha = 0.8, stroke = 0.8, color = "black") +
        scale_shape_manual(values = c("Spring" = 21, "Summer" = 22, "Fall" = 23, "Winter" = 24), name = "Season") +
        scale_fill_manual(values = cArray, name = legend_name) +
        scale_color_manual(values = cArray, name = legend_name)

      if (v == "season") {
        p_temp <- p_temp + guides(
          fill = guide_legend(override.aes = list(color = "black")),
          color = "none"
        )
      } else {
        p_temp <- p_temp + guides(
          fill = guide_legend(override.aes = list(shape = 21, color = "black")),
          shape = guide_legend(override.aes = list(fill = "gray70"))
        )
      }

      if (length(unique(pcoa_data[[v]])) >= 2) {
        p_temp <- p_temp +
          stat_ellipse(aes(x = PCoA1, y = PCoA2, group = .data[[v]], color = .data[[v]]),
            inherit.aes = FALSE, type = "t", linetype = 2, linewidth = 0.6, level = 0.95
          )
      }

      p_temp <- p_temp +
        labs(
          x = paste0("PCoA1 (", pcoa1_var, "%)"),
          y = paste0("PCoA2 (", pcoa2_var, "%)"),
          fill = legend_name,
          color = legend_name,
          title = paste(virus_full_name, "-", paste(sort(unique(combined$city)), collapse = " and ")),
          subtitle = plot_subtitle
        ) +
        theme_bw(base_size = 14) +
        theme(legend.position = "right")

      plot_filenames[[v]] <- paste0(virus_full_name, "_pcoa_", clean_var_name)
    }

    temp_plots[[v]] <- p_temp
  }

  plot_lim <- max(1.0, max(abs(c(pcoa_data$PCoA1, pcoa_data$PCoA2)), na.rm = TRUE) * 1.05)

  final_plots <- list()
  for (v in names(temp_plots)) {
    p_final <- temp_plots[[v]] + coord_fixed(xlim = c(-plot_lim, plot_lim), ylim = c(-plot_lim, plot_lim))

    stat_text <- get_perm_stats_local(v, permanova_results_list)
    if (!is.null(stat_text)) {
      p_final <- p_final + annotate("label",
        x = plot_lim * 0.95, y = plot_lim * 0.95,
        label = stat_text, hjust = 1, vjust = 1,
        fill = "white", color = "black", alpha = 0.8,
        fontface = "italic", size = 4.5
      )
    }
    final_plots[[v]] <- p_final
  }

  for (v in names(final_plots)) {
    base_fn <- plot_filenames[[v]]
    tiff_fn <- file.path(combined_out_dir, paste0(base_fn, ".tiff"))
    ggsave(tiff_fn, plot = final_plots[[v]], device = "tiff", dpi = 600, width = 16, height = 10, compression = "lzw")
  }

  wb <- createWorkbook()
  addWorksheet(wb, "PCoA_Metadata_Scores")
  writeData(wb, "PCoA_Metadata_Scores", pcoa_data)

  if (!is.null(fit_all)) {
    addWorksheet(wb, "PERMANOVA_Combined")
    writeData(wb, "PERMANOVA_Combined", as.data.frame(fit_all), rowNames = TRUE)
  }

  for (v in names(permanova_results_list)) {
    sheet_name <- paste0("PERMANOVA_", substr(v, 1, 10))
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, permanova_results_list[[v]])

    tiff_path <- file.path(combined_out_dir, paste0(plot_filenames[[v]], ".tiff"))
    if (file.exists(tiff_path)) {
      insertImage(wb, sheet = sheet_name, file = tiff_path, width = 10, height = 7.5, startRow = 2, startCol = 9, units = "in")
    }
  }

  cat("\nGenerating superimposed regression plots and minor genotype dynamics...\n")
  
  # Set consistent colors for Ajax and Barrie
  unique_cities_plot <- sort(unique(combined$city))
  default_city_colors <- c("blue", "red", "green", "purple", "orange", "cyan", "magenta")
  city_cols <- setNames(default_city_colors[1:length(unique_cities_plot)], unique_cities_plot)
  
  # 1. pN/pS combined regression
  pnps_list <- list()
  shannon_list <- list()
  
  for (i in seq_along(output_dirs)) {
    odir <- output_dirs[i]
    c_city <- cities[i]
    
    pnps_file <- file.path(odir, "iVar", paste0(c_city, "_", virus, "_regression_df.csv"))
    if (file.exists(pnps_file)) {
      df <- read.csv(pnps_file)
      df$city <- c_city
      pnps_list[[c_city]] <- df
    }
    
    shannon_file <- file.path(dirname(odir), "RT-qPCR", paste0(virus, "_", c_city, "_shannon_div.csv"))
    if (file.exists(shannon_file)) {
      df <- read.csv(shannon_file)
      df$city <- c_city
      shannon_list[[c_city]] <- df
    }
  }

  if (length(pnps_list) > 0) {
    comb_pnps <- bind_rows(pnps_list)
    stats_pnps <- comb_pnps %>%
      dplyr::group_by(city) %>%
      dplyr::filter(dplyr::n() > 3) %>%
      dplyr::summarize(
        r2 = {
          r_full = summary(lm(pN_pS ~ Relative_Abundance + log10(mean_depth)))$r.squared
          r_red = summary(lm(pN_pS ~ log10(mean_depth)))$r.squared
          round((r_full - r_red) / (1 - r_red), 3)
        },
        pval = {
           coefs <- summary(lm(pN_pS ~ Relative_Abundance + log10(mean_depth)))$coefficients
           if ("Relative_Abundance" %in% rownames(coefs)) coefs["Relative_Abundance", 4] else NA_real_
        },
        n_obs = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        pval_str = ifelse(!is.na(pval) & pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
        label = paste0(city, ": R² = ", r2, ", p = ", pval_str, ", n = ", n_obs)
      )
    
    subtitle_str <- paste(stats_pnps$label, collapse = "\n")
    
    p_comb_pnps <- ggplot(comb_pnps, aes(x = Relative_Abundance, y = pN_pS, color = city)) +
      geom_point(alpha = 0.7, size = 3) +
      geom_smooth(method = "lm", aes(fill = city), linetype = "dashed", se = TRUE, alpha = 0.15) +
      scale_color_manual(values = city_cols) +
      scale_fill_manual(values = city_cols) +
      theme_bw() +
      labs(
        title = paste(virus_full_name, "-", paste(sort(unique(combined$city)), collapse = " and ")),
        x = "Relative Genotype Abundance (%)",
        y = "pN/pS Ratio",
        color = "City", fill = "City"
      ) +
      annotate("label", x = Inf, y = Inf, label = subtitle_str, hjust = 1, vjust = 1, fill = "white", color = "black")
    ggsave(file.path(combined_out_dir, paste0(virus, "_Combined_pN_pS_Regression.tiff")), plot = p_comb_pnps, width = 10, height = 8, dpi = 600, compression="lzw")
    final_plots$comb_pnps <- p_comb_pnps
  }

  if (length(shannon_list) > 0) {
    comb_shannon <- bind_rows(shannon_list)
    stats_shannon <- comb_shannon %>%
      dplyr::group_by(city) %>%
      dplyr::filter(sum(is.finite(log10_viral_load) & is.finite(Shannon_Diversity)) > 2) %>%
      dplyr::summarize(
        r2 = round(summary(lm(Shannon_Diversity ~ log10_viral_load))$r.squared, 3),
        pval = summary(lm(Shannon_Diversity ~ log10_viral_load))$coefficients[2, 4],
        n_obs = sum(is.finite(log10_viral_load) & is.finite(Shannon_Diversity)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        pval_str = ifelse(!is.na(pval) & pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
        label = paste0(city, ": R² = ", r2, ", p = ", pval_str, ", n = ", n_obs)
      )
    
    subtitle_str_s <- paste(stats_shannon$label, collapse = "\n")
    
    p_comb_shannon <- ggplot(comb_shannon, aes(x = log10_viral_load, y = Shannon_Diversity, color = city)) +
      geom_point(alpha = 0.7, size = 3) +
      geom_smooth(method = "lm", aes(fill = city), linetype = "dashed", se = TRUE, alpha = 0.15) +
      scale_color_manual(values = city_cols) +
      scale_fill_manual(values = city_cols) +
      theme_bw() +
      labs(
        title = paste(virus_full_name, "-", paste(sort(unique(combined$city)), collapse = " and ")),
        subtitle = subtitle_str_s,
        x = "Viral Load [Log10 Scale]",
        y = "Shannon Diversity Index",
        color = "City", fill = "City"
      )
    ggsave(file.path(combined_out_dir, paste0(virus, "_Combined_Shannon_Regression.tiff")), plot = p_comb_shannon, width = 10, height = 8, dpi = 600, compression="lzw")
    final_plots$comb_shannon <- p_comb_shannon
  }

  # Minor Genotype Plot (GII.15 and CV-A1)
  if (nrow(combined) > 0) {
    minor_gens <- c()
    if (virus == "NV" && "GII-15" %in% colnames(combined)) minor_gens <- c(minor_gens, "GII-15")
    if (virus == "EV" && "CV-A1" %in% colnames(combined)) minor_gens <- c(minor_gens, "CV-A1")
    
    if (length(minor_gens) > 0) {
      minor_abund <- combined %>%
        dplyr::select(date, city, dplyr::all_of(minor_gens)) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(minor_gens), names_to = "Genotype", values_to = "Abundance") %>%
        dplyr::filter(Abundance > 0)
        
      if (length(shannon_list) > 0) {
        vl_df <- bind_rows(shannon_list) %>% 
          dplyr::mutate(date = as.Date(Date)) %>%
          dplyr::select(date, city, log10_viral_load)
          
        minor_plot_df <- dplyr::left_join(minor_abund, vl_df, by = c("date", "city")) %>% dplyr::filter(!is.na(log10_viral_load))
        
        if (nrow(minor_plot_df) > 0) {
          max_abund <- max(minor_plot_df$Abundance, na.rm = TRUE)
          max_vl <- max(minor_plot_df$log10_viral_load, na.rm = TRUE)
          scale_fac <- max_abund / max_vl
          
          p_minor <- ggplot(minor_plot_df, aes(x = date)) +
            geom_col(aes(y = Abundance, fill = Genotype), position = "stack", alpha = 0.7) +
            geom_line(aes(y = log10_viral_load * scale_fac), color = "black", linewidth = 1) +
            geom_point(aes(y = log10_viral_load * scale_fac), color = "black", size = 2) +
            scale_y_continuous(
              name = "Relative Abundance (%)",
              sec.axis = sec_axis(~ . / scale_fac, name = "Viral Load (Log10)")
            ) +
            facet_wrap(~city, ncol = 1) +
            theme_bw() +
            labs(
              title = paste(virus_full_name, "-", paste(sort(unique(combined$city)), collapse = " and ")),
              x = "Date"
            ) +
            theme(legend.position = "bottom") +
            scale_fill_manual(values = cArray)
          
          ggsave(file.path(combined_out_dir, paste0(virus, "_Minor_Genotype_Dynamics.tiff")), plot = p_minor, width = 10, height = 8, dpi = 600, compression="lzw")
          final_plots$minor_genotypes <- p_minor
        }
      }
    }
  }

  excel_out <- file.path(combined_out_dir, paste0(virus_full_name, "_PCoA_Metadata_Results.xlsx"))
  saveWorkbook(wb, excel_out, overwrite = TRUE)
  cat("Excel exported successfully to '", excel_out, "'\n", sep = "")
  return(final_plots)
}

# Main Execution Loop

if (Sys.info()[["sysname"]] %in% c("Windows")) {
  plan(multisession)
} else if (Sys.info()[["sysname"]] == "Darwin") {
  plan(sequential)
} else {
  plan(multicore)
}

subfolders <- c(
  file.path(workspace_root, "EV", "Ajax"),
  file.path(workspace_root, "EV", "Barrie"),
  file.path(workspace_root, "NV", "Ajax"),
  file.path(workspace_root, "NV", "Barrie")
)

output_directories <- list()
global_reports <- list()
meta_names <- list()

# First loop: Setup directories
for (subf in subfolders) {
  xlsx_files <- list.files(path = subf, pattern = "\\.xlsx$")
  xlsx_files <- xlsx_files[!grepl("Analysis_Final|Frequencies|blast_top_hits|genotype_pc|abund_mat|bc_scores", xlsx_files)]
  if (length(xlsx_files) > 0) {
    meta_name <- sub("\\.xlsx$", "", xlsx_files[1])
    meta_name_nv <- sub("^Noro_", "NV_", meta_name)
    out_dir <- file.path(workspace_root, "Output", meta_name_nv)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    output_directories[[subf]] <- out_dir
    meta_names[[subf]] <- meta_name_nv
  }
}

# PHASE 1: Genotyper (Parallel)
if (RUN_PHASE1_GENOTYPER) {
  main_qpcr_dfs <- NULL
  qpcr_file_main <- file.path(workspace_root, "RT-qPCR_data.xlsx")
  if (file.exists(qpcr_file_main)) {
    main_qpcr_dfs <- list()
    available_sheets <- openxlsx::getSheetNames(qpcr_file_main)
    for (sh in available_sheets) {
      main_qpcr_dfs[[sh]] <- openxlsx::read.xlsx(qpcr_file_main, sheet = sh, detectDates = TRUE, check.names = FALSE)
    }
  }

  cat("Beginning Phase 1: Sequential Genotyper processing of", length(subfolders), "suites...\n")
  res_p1 <- lapply(subfolders, function(subf) {
    out_dir <- output_directories[[subf]]
    if (is.null(out_dir)) {
      return(NULL)
    }
    tryCatch(
      {
        return(run_genotyper_pipeline(subf, out_dir, blast_bin, cArray, genotype_mapping, qpcr_dfs = main_qpcr_dfs))
      },
      error = function(e) {
        log_warning_error(subf, "Phase 1 Genotyper", e$message)
        cat("Error running Genotyper on", subf, ":", e$message, "\n")
        return(NULL)
      }
    )
  })

  for (i in seq_along(subfolders)) {
    subf <- subfolders[i]
    if (!is.null(meta_names[[subf]]) && !is.null(res_p1[[i]])) {
      if (is.null(global_reports[[meta_names[[subf]]]])) {
        global_reports[[meta_names[[subf]]]] <- list()
      }
      global_reports[[meta_names[[subf]]]]$genotyper <- res_p1[[i]]
    }
  }
}

# DYNAMIC DATABASE BUILDING
# Scan the exact genotypes identified by Phase 1 before running Phase 2
cat("Scanning Genotyper outputs to dynamically identify valid EV genotypes...\n")
dynamic_ev_genotypes <- c()
for (subf in subfolders) {
  if (grepl("/EV/", subf)) {
    out_dir <- output_directories[[subf]]
    if (!is.null(out_dir)) {
      pc_file <- file.path(out_dir, "genotype_pc_results.xlsx")
      if (file.exists(pc_file)) {
        tryCatch(
          {
            df <- openxlsx::read.xlsx(pc_file, sheet = 1)
            if ("genotype" %in% names(df)) {
              valid <- unique(df$genotype[!df$genotype %in% c("Other", "Could not assign")])
              dynamic_ev_genotypes <- c(dynamic_ev_genotypes, valid)
            }
          },
          error = function(e) NULL
        )
      }
    }
  }
}

ev_genotypes_pre <- toupper(gsub("\\.", "-", unique(dynamic_ev_genotypes)))
if (length(ev_genotypes_pre) > 0) {
  cat("Identified", length(ev_genotypes_pre), "valid EV genotypes from Phase 1.\n")
  cat("Pre-building dynamic EV VP1 database sequentially to avoid NCBI rate limits...\n")
  vp1_fasta_pre <- file.path(ivar_workspace, "all_EV_types_combined_VP1.fasta")
  vp1_csv_pre <- file.path(ivar_workspace, "EV_VP1_mappings.csv")
  build_ev_vp1_database(ev_genotypes_pre, vp1_fasta_pre, vp1_csv_pre)
} else {
  cat("No valid EV genotypes identified from Phase 1 outputs.\n")
}

# PHASE 2: iVar
if (RUN_PHASE2_IVAR) {
  if (RUN_COLABFOLD_LOCAL) {
    cat("Beginning Phase 2: Sequential iVar processing of", length(subfolders), "suites (ColabFold enabled - preventing memory exhaustion)...\n")
    res_p2 <- lapply(subfolders, function(subf) {
      out_dir <- output_directories[[subf]]
      if (is.null(out_dir)) {
        return(NULL)
      }
      parts <- strsplit(subf, "/")[[1]]
      virus <- parts[length(parts) - 1]
      meta_name_nv <- meta_names[[subf]]
      tryCatch(
        {
          ivar_out_dir <- file.path(out_dir, "iVar")
          dir.create(ivar_out_dir, showWarnings = FALSE, recursive = TRUE)
          return(run_ivar_pipeline(subf, ivar_out_dir, ivar_workspace, virus, meta_name_nv))
        },
        error = function(e) {
          log_warning_error("Error running iVar pipeline on", subf, ":", e$message, "\n")
          return(NULL)
        }
      )
    })
  } else {
    cat("Beginning Phase 2: Parallel iVar processing of", length(subfolders), "suites...\n")
    res_p2 <- future.apply::future_lapply(subfolders, function(subf) {
      out_dir <- output_directories[[subf]]
      if (is.null(out_dir)) {
        return(NULL)
      }
      parts <- strsplit(subf, "/")[[1]]
      virus <- parts[length(parts) - 1]
      meta_name_nv <- meta_names[[subf]]
      tryCatch(
        {
          ivar_out_dir <- file.path(out_dir, "iVar")
          dir.create(ivar_out_dir, showWarnings = FALSE, recursive = TRUE)
          return(run_ivar_pipeline(subf, ivar_out_dir, ivar_workspace, virus, meta_name_nv))
        },
        error = function(e) {
          log_warning_error("Error running iVar pipeline on", subf, ":", e$message, "\n")
          return(NULL)
        }
      )
    }, future.seed = TRUE)
  }

  for (i in seq_along(subfolders)) {
    subf <- subfolders[i]
    if (!is.null(meta_names[[subf]]) && !is.null(res_p2[[i]])) {
      if (is.null(global_reports[[meta_names[[subf]]]])) {
        global_reports[[meta_names[[subf]]]] <- list()
      }
      global_reports[[meta_names[[subf]]]]$ivar <- res_p2[[i]]
    }
  }
}

# Fix the folder list so the combined analysis still works
# It expects a named list with names corresponding to meta_name_nv, not subf
clean_output_directories <- list()
for (subf in names(output_directories)) {
  clean_output_directories[[meta_names[[subf]]]] <- output_directories[[subf]]
}
output_directories <- clean_output_directories

# Execute Combined Analysis
if (RUN_PHASE3_COMBINED) {
  cat("\nRunning combined analysis across cities...\n")

  # Combined EV Analysis (using EV_Ajax and EV_Barrie output directories)
  ev_dirs <- output_directories[names(output_directories) %in% c("EV_Ajax", "EV_Barrie")]
  if (length(ev_dirs) >= 2) {
    tryCatch(
      {
        ev_plots <- run_combined_analysis("EV", unlist(ev_dirs), file.path(workspace_root, "Output", "EV_Combined"), cArray)
        global_reports$combined_EV <- ev_plots
      },
      error = function(e) {
        log_warning_error("Error in combined EV analysis:", e$message, "\n")
      }
    )
  }

  # Combined NV Analysis (using NV_Ajax and NV_Barrie output directories)
  nv_dirs <- output_directories[names(output_directories) %in% c("NV_Ajax", "NV_Barrie")]
  if (length(nv_dirs) >= 2) {
    tryCatch(
      {
        nv_plots <- run_combined_analysis("NV", unlist(nv_dirs), file.path(workspace_root, "Output", "NV_Combined"), cArray)
        global_reports$combined_NV <- nv_plots
      },
      error = function(e) {
        log_warning_error("Error in combined NV analysis:", e$message, "\n")
      }
    )
  }
}

# PDF Report Helpers
draw_title_page <- function(title, subtitle) {
  # Draws on the automatically initialized first page of the PDF device (no grid.newpage)
  grid::grid.rect(gp = grid::gpar(fill = "#F3F4F6", col = NA))
  grid::grid.text(title, y = 0.6, gp = grid::gpar(fontsize = 24, fontface = "bold", col = "#1F2937"))
  grid::grid.text(subtitle, y = 0.5, gp = grid::gpar(fontsize = 14, fontface = "italic", col = "#4B5563"))
  grid::grid.text(paste("Generated on:", Sys.time()), y = 0.2, gp = grid::gpar(fontsize = 10, col = "#9CA3AF"))
}

draw_section_page <- function(title) {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "#1F2937", col = NA))
  grid::grid.text(title, gp = grid::gpar(fontsize = 20, fontface = "bold", col = "white"))
}

# Compile Final PDF Report
if (BUILD_PDF_REPORT) {
  cat("\nCompiling final PDF Summary Report...\n")
  pdf_path <- file.path(workspace_root, "Output", "WW_Genotyper_Summary_Report_v11.pdf")
  if (!dir.exists(dirname(pdf_path))) dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    {
      pdf(pdf_path, width = 11, height = 8.5)

      # Title Page
      draw_title_page(
        title = "Wastewater Genotyper Summary Report (v11)",
        subtitle = "Comprehensive Analysis of Enterovirus & Norovirus Wastewater Sequences"
      )

      # Loop over each subfolder results
      for (suite in names(output_directories)) {
        # Section Page
        draw_section_page(paste("Run Suite:", suite))

        suite_plots <- global_reports[[suite]]
        if (is.null(suite_plots)) next

        # 1. Genotyper plots
        if (!is.null(suite_plots$genotyper)) {
          if (!is.null(suite_plots$genotyper$abundance)) {
            suppressWarnings(print(suite_plots$genotyper$abundance))
          }
          if (!is.null(suite_plots$genotyper$abs_abund)) {
            suppressWarnings(print(suite_plots$genotyper$abs_abund))
          }
          if (!is.null(suite_plots$genotyper$qc)) {
            suppressWarnings(print(suite_plots$genotyper$qc))
          }
          if (!is.null(suite_plots$genotyper$pcoa)) {
            suppressWarnings(print(suite_plots$genotyper$pcoa))
          }
        }

        # 2. iVar plots
        if (!is.null(suite_plots$ivar)) {
          ivar_p <- suite_plots$ivar

          if (!is.null(ivar_p$read_depth)) {
            suppressWarnings(print(ivar_p$read_depth))
          }
          if (!is.null(ivar_p$diversity_across_alignment)) {
            suppressWarnings(print(ivar_p$diversity_across_alignment))
          }
          if (!is.null(ivar_p$mutation_trajectories)) {
            suppressWarnings(print(ivar_p$mutation_trajectories))
          }
          if (!is.null(ivar_p$mean_diversity_over_time)) {
            suppressWarnings(print(ivar_p$mean_diversity_over_time))
          }
          if (!is.null(ivar_p$pn_ps_over_time)) {
            suppressWarnings(print(ivar_p$pn_ps_over_time))
          }
          if (!is.null(ivar_p$abundance_vs_pnps)) {
            suppressWarnings(print(ivar_p$abundance_vs_pnps))
          }
          if (!is.null(ivar_p$abundance_vs_pnps_faceted)) {
            suppressWarnings(print(ivar_p$abundance_vs_pnps_faceted))
          }
          if (!is.null(ivar_p$tajima_d)) {
            suppressWarnings(print(ivar_p$tajima_d))
          }
          if (!is.null(ivar_p$strand_bias)) {
            suppressWarnings(print(ivar_p$strand_bias))
          }

          if (!is.null(ivar_p$forecasts)) {
            for (p in ivar_p$forecasts) {
              suppressWarnings(print(p))
            }
          }

          # Heatmaps
          if (!is.null(ivar_p$heatmaps) && length(ivar_p$heatmaps) > 0) {
            for (h_name in names(ivar_p$heatmaps)) {
              ComplexHeatmap::draw(ivar_p$heatmaps[[h_name]])
            }
          }

          # Lollipops
          if (!is.null(ivar_p$lollipops) && length(ivar_p$lollipops) > 0) {
            for (l_name in names(ivar_p$lollipops)) {
              suppressWarnings(print(ivar_p$lollipops[[l_name]]))
            }
          }

          # 3D Structure Alignments
          if (!is.null(ivar_p$structures) && nrow(ivar_p$structures) > 0) {
            # Print the summary bar plots if they were generated
            if (!is.null(ivar_p$plddt_summary)) {
              suppressWarnings(print(ivar_p$plddt_summary))
            }
            if (!is.null(ivar_p$iptm_summary)) {
              suppressWarnings(print(ivar_p$iptm_summary))
            }

            # Print the ipTM scatter plot if it was generated
            if (!is.null(ivar_p$iptm)) {
              suppressWarnings(print(ivar_p$iptm))
            }

            # Print the RMSD vs ipTM scatter plot if it was generated
            if (!is.null(ivar_p$rmsd_iptm)) {
              suppressWarnings(print(ivar_p$rmsd_iptm))
            }

            for (i in seq_len(nrow(ivar_p$structures))) {
              row <- ivar_p$structures[i, ]
              if (file.exists(row$CXS_Path)) {
                draw_section_page(paste("3D Structural Alignment:", row$Sample_Protein))
                grid::grid.newpage()
                grid::grid.text(paste("RMSD (Deviation from Reference):", row$RMSD, "Å\n\n(See .cxs file for 3D alignment in ChimeraX)"),
                  y = 0.5, gp = grid::gpar(fontsize = 16, fontface = "bold", col = "red")
                )
              }
            }
          }
        }
      }

      # 3. Combined Analysis
      if (!is.null(global_reports$combined_EV) && length(global_reports$combined_EV) > 0) {
        draw_section_page("Combined Enterovirus Analysis")
        for (p_name in names(global_reports$combined_EV)) {
          suppressWarnings(print(global_reports$combined_EV[[p_name]]))
        }
      }

      if (!is.null(global_reports$combined_NV) && length(global_reports$combined_NV) > 0) {
        draw_section_page("Combined Norovirus Analysis")
        for (p_name in names(global_reports$combined_NV)) {
          suppressWarnings(print(global_reports$combined_NV[[p_name]]))
        }
      }

      # 4. Summary of Flags for Manual Investigation
      flagged_samples <- data.frame(Sample = character(), Reason = character(), stringsAsFactors = FALSE)

      for (suite in names(output_directories)) {
        suite_plots <- global_reports[[suite]]
        if (is.null(suite_plots$ivar)) next

        ivar_p <- suite_plots$ivar

        # Check Structures for high RMSD or ipTM shift
        if (!is.null(ivar_p$structures) && nrow(ivar_p$structures) > 0) {
          for (i in seq_len(nrow(ivar_p$structures))) {
            row <- ivar_p$structures[i, ]

            # Check pLDDT quality gate first
            if (!is.na(row$Sample_pLDDT) && row$Sample_pLDDT < MIN_PLDDT_STRUCTURAL) {
              flagged_samples <- rbind(flagged_samples, data.frame(Sample = row$Sample_Protein, Reason = paste0("Low Model Confidence (pLDDT: ", round(row$Sample_pLDDT, 1), ") - structural flags suppressed"), stringsAsFactors = FALSE))
            } else {
              # Only check structure scores if the model is good enough
              if (!is.na(row$RMSD) && row$RMSD > 1.0) {
                flagged_samples <- rbind(flagged_samples, data.frame(Sample = row$Sample_Protein, Reason = paste0("High Structural Deviation (RMSD: ", row$RMSD, " Å)"), stringsAsFactors = FALSE))
              }
              if (!is.na(row$Ref_ipTM) && !is.na(row$Sample_ipTM)) {
                iptm_diff <- row$Sample_ipTM - row$Ref_ipTM
                if (iptm_diff <= -0.05) {
                  flagged_samples <- rbind(flagged_samples, data.frame(Sample = row$Sample_Protein, Reason = paste0("Decreased Binding Affinity (ipTM Drop: ", round(iptm_diff, 3), ")"), stringsAsFactors = FALSE))
                } else if (iptm_diff >= 0.05) {
                  flagged_samples <- rbind(flagged_samples, data.frame(Sample = row$Sample_Protein, Reason = paste0("Increased Binding Affinity (ipTM Rise: ", round(iptm_diff, 3), ")"), stringsAsFactors = FALSE))
                }
              }
            }

            # The ESM2 score only uses the sequence, not the 3D structure score
            if ("ESM2_Score" %in% names(row) && !is.na(row$ESM2_Score) && row$ESM2_Score <= -5.0) {
              flagged_samples <- rbind(flagged_samples, data.frame(Sample = row$Sample_Protein, Reason = paste0("Extreme Risk Mutation (ESM-2 Score: ", round(row$ESM2_Score, 2), ")"), stringsAsFactors = FALSE))
            }
          }
        }

        # Check Forecasts for emerging status
        if (!is.null(ivar_p$forecasts) && length(ivar_p$forecasts) > 0) {
          for (mut_name in names(ivar_p$forecasts)) {
            flagged_samples <- rbind(flagged_samples, data.frame(Sample = suite, Reason = paste0("Emerging Variant Forecasted: ", mut_name), stringsAsFactors = FALSE))
          }
        }
      }

      if (nrow(flagged_samples) > 0) {
        draw_section_page("Samples Flagged for Manual Investigation")
        grid::grid.newpage()

        flag_text <- "The following samples/variants met criteria for manual review:\n(Criteria: RMSD > 1.0, |ipTM shift| > 0.05, ESM-2 Score <= -5.0, or Emerging Status)\n\n"
        max_lines <- 35
        for (i in seq_len(min(nrow(flagged_samples), max_lines))) {
          line <- paste0("• [", flagged_samples$Sample[i], "] ", flagged_samples$Reason[i])
          flag_text <- paste0(flag_text, line, "\n")
        }
        if (nrow(flagged_samples) > max_lines) {
          flag_text <- paste0(flag_text, "\n... and ", nrow(flagged_samples) - max_lines, " more (see raw data).")
        }

        grid::grid.text(flag_text, x = 0.05, y = 0.95, just = c("left", "top"), gp = grid::gpar(fontsize = 10, col = "darkred"))
      }

      # Generate Structural Pass Summary PDF
      struct_pass_list <- list()
      for (suite in names(global_reports)) {
        if (!is.null(global_reports[[suite]]$ivar$structures) && nrow(global_reports[[suite]]$ivar$structures) > 0) {
          tmp <- global_reports[[suite]]$ivar$structures
          tmp$Suite <- suite
          struct_pass_list[[suite]] <- tmp
        }
      }

      struct_pdf_path <- file.path(workspace_root, "Output", "WW_Genotyper_Structural_Pass_Summary.pdf")
      if (length(struct_pass_list) > 0) {
        cat("\nGenerating Structural Pass Summary PDF...\n")
        struct_combined <- dplyr::bind_rows(struct_pass_list)
        if (!"Sample_ipTM" %in% colnames(struct_combined)) struct_combined$Sample_ipTM <- NA
        if (!"Sample_pLDDT" %in% colnames(struct_combined)) struct_combined$Sample_pLDDT <- NA

        struct_summary <- data.frame(
          Suite = struct_combined$Suite,
          Sample = struct_combined$Sample_Protein,
          Model_Confidence = ifelse(!is.na(struct_combined$Sample_pLDDT), sprintf("%.1f", struct_combined$Sample_pLDDT), "N/A"),
          RMSD = sprintf("%.3f", struct_combined$RMSD),
          Binding_Affinity = ifelse(!is.na(struct_combined$Sample_ipTM), sprintf("%.3f", struct_combined$Sample_ipTM), "N/A"),
          stringsAsFactors = FALSE
        )
        struct_summary <- struct_summary[order(-struct_combined$RMSD), ]

        pdf(struct_pdf_path, width = 8.5, height = 11)
        grid::grid.newpage()
        grid::grid.text("Samples Successfully Processed for Structural Analysis", y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))

        table_text <- paste(capture.output(print(struct_summary, row.names = FALSE, right = FALSE)), collapse = "\n")
        grid::grid.text(table_text, x = 0.05, y = 0.90, just = c("left", "top"), gp = grid::gpar(fontsize = 10, fontfamily = "mono"))
        dev.off()
      } else {
        if (file.exists(struct_pdf_path)) file.remove(struct_pdf_path)
      }

      dev.off()
      cat("PDF Summary Report generated successfully at:", pdf_path, "\n")
    },
    error = function(e) {
      log_warning_error("Warning: Failed to generate PDF Summary Report:", e$message, "\n")
      if (!is.null(dev.list())) dev.off()
    }
  )
}
if (file.exists(error_log_file)) {
  cat("\n==================================================\n")
  cat("               PIPELINE WARNINGS & ERRORS         \n")
  cat("==================================================\n")
  log_lines <- readLines(error_log_file)
  if (length(log_lines) > 0) {
    for (line in unique(log_lines)) {
      cat(line, "\n")
    }
  } else {
    cat("No warnings or errors detected.\n")
  }
  cat("==================================================\n")
  file.remove(error_log_file)
}

cat("\n==================================================\n")
cat("               PIPELINE OUTPUT SUMMARY            \n")
cat("==================================================\n")

if (length(output_directories) > 0) {
  for (suite_name in names(output_directories)) {
    out_dir <- output_directories[[suite_name]]
    cat(sprintf("\n[%s]\n", suite_name))

    # Phase 1: Genotyper (outputs are saved directly in out_dir)
    if (dir.exists(out_dir)) {
      # See if Phase 1 output files exist yet
      phase1_files <- list.files(out_dir, pattern = "genotype_pc_results\\.xlsx|abund_mat\\.csv|bc_scores\\.csv|blast_top_hits", full.names = TRUE)
      if (length(phase1_files) > 0) {
        cat("  -> Phase 1 (Genotyping & Mapping)\n")
        cat("       Location: ", out_dir, "\n")
        cat("       Contents: FASTAs, mapping stats, primer trimming, initial typing CSVs/Excels\n")
      } else {
        cat("  -> Phase 1 (Genotyping & Mapping) : [Not Generated/Missing]\n")
      }
    } else {
      cat("  -> Phase 1 (Genotyping & Mapping) : [Not Generated/Missing]\n")
    }

    # Phase 2: iVar & Population Genetics
    ivar_dir <- file.path(out_dir, "iVar")
    if (dir.exists(ivar_dir)) {
      ivar_xlsx <- file.path(ivar_dir, paste0(suite_name, ".xlsx"))
      if (file.exists(ivar_xlsx)) {
        cat("  -> Phase 2 (Variant & PopGen Excel)\n")
        cat("       Location: ", ivar_xlsx, "\n")
        cat("       Contents: iVar SNVs/Indels, Coverage, Tajima's D, Pi, Strand Bias, ESM-2 Risk Scores\n")
      } else {
        cat("  -> Phase 2 (iVar & PopGen) Outputs:\n")
        cat("       Location: ", ivar_dir, "\n")
      }

      # Structural Analysis
      struct_dir <- file.path(ivar_dir, "Structures")
      if (dir.exists(struct_dir) && length(list.files(struct_dir)) > 0) {
        cat("  -> Phase 2 (3D Struct. / ColabFold / ChimeraX)\n")
        cat("       Location: ", struct_dir, "\n")
        cat("       Contents: .pdb folds, .json pLDDT/ipTM scores, .cxs alignments, Structural_Deviation_Batch.csv\n")
      }

      # Prophet
      prophet_pdf <- file.path(ivar_dir, "Prophet_Emerging_Forecasts.pdf")
      if (file.exists(prophet_pdf)) {
        cat("  -> Phase 2 (Emerging Variant Forecasts)\n")
        cat("       Location: ", prophet_pdf, "\n")
        cat("       Contents: ARIMA/Prophet timeseries forecasting plots for variant frequencies\n")
      }
    } else {
      cat("  -> Phase 2 (iVar & PopGen)        : [Not Generated/Missing]\n")
    }
  }
}

if (RUN_PHASE3_COMBINED) {
  cat("\n[Phase 3: Global Combined Analysis]\n")

  pdf_path <- file.path(workspace_root, "Output", "WW_Genotyper_Summary_Report_v11.pdf")
  if (file.exists(pdf_path)) {
    cat("  -> Global PDF Summary Report\n")
    cat("       Location: ", pdf_path, "\n")
    cat("       Contents: Lollipop plots, diversity timeseries, heatmaps, Structural Confidence Bar Plots, and flagged samples\n")
  } else {
    cat("  -> Global PDF Summary Report      : [Not Generated/Missing]\n")
  }

  struct_pdf_path <- file.path(workspace_root, "Output", "WW_Genotyper_Structural_Pass_Summary.pdf")
  if (file.exists(struct_pdf_path)) {
    cat("  -> Global Structural Confidence Summary\n")
    cat("       Location: ", struct_pdf_path, "\n")
    cat("       Contents: Table of all structurally analyzed samples, their mean pLDDT, RMSD, and binding affinities\n")
  } else {
    cat("  -> Global Structural Confidence Summary: [None Passed/Missing]\n")
  }

  ev_excel <- file.path(workspace_root, "Output", "EV_Combined", "EV_PCoA_Metadata_Results.xlsx")
  if (file.exists(ev_excel)) {
    cat("  -> Enterovirus PCoA & Metadata\n")
    cat("       Location: ", ev_excel, "\n")
    cat("       Contents: PCoA coordinates, metadata correlations, and combined genetic diversity statistics\n")
  }

  nv_excel <- file.path(workspace_root, "Output", "NV_Combined", "NV_PCoA_Metadata_Results.xlsx")
  if (file.exists(nv_excel)) {
    cat("  -> Norovirus PCoA & Metadata\n")
    cat("       Location: ", nv_excel, "\n")
    cat("       Contents: PCoA coordinates, metadata correlations, and combined genetic diversity statistics\n")
  }

  if (file.exists(error_log_file) && file.info(error_log_file)$size > 0) {
    cat("  -> Error/Warning Log\n")
    cat("       Location: ", error_log_file, "\n")
    cat("       Contents: Detailed list of warnings and errors generated during this run\n")
  }
}

# Close the terminal logging connection
sink()
