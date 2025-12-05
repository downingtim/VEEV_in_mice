#!/usr/bin/env Rscript

# VCF Quality Filtering Script
# Filters noisy SNPs from VCF files based on multiple quality metrics

library(tidyverse)

config <- list(
  min_qual = 30,           # Minimum QUAL score
  min_depth = 5,          # Minimum read depth (DP)
  max_depth = 1500,         # Maximum read depth (remove extremely high coverage)
  min_mq = 40,             # Minimum mapping quality
  max_mq0f = 0.1,          # Maximum fraction of MQ0 reads
  min_af = 0.05,           # Minimum allele frequency for haploid calls
  min_fq = -50,            # Minimum FQ (more negative = more confident)
  vdb_threshold = 0.01     # Minimum VDB (Variant Distance Bias)
)

# Function to parse a single VCF file
parse_vcf <- function(vcf_file) {
  
  if (grepl("\\.gz$", vcf_file)) {
    lines <- readLines(gzfile(vcf_file))
  } else {
    lines <- readLines(vcf_file)
  }
  
  header_end <- max(which(grepl("^##", lines)))
  col_line <- lines[header_end + 1]
  data_lines <- lines[(header_end + 2):length(lines)]
  
  if (length(data_lines) == 0) {
    message(paste("No variants in:", basename(vcf_file)))
    return(NULL)
  }
  
  # Parse column names - remove leading # and handle spaces
  col_names <- str_split(str_remove(col_line, "^#"), "\t")[[1]]
  col_names <- str_trim(col_names)
  
  # Standard VCF columns (first 9 are always the same)
  std_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  
  # Parse data - use fill=TRUE to handle variable column counts
  vcf_data <- read.table(text = paste(data_lines, collapse = "\n"),
                         sep = "\t",
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote = "",
                         fill = TRUE,
                         col.names = c(std_cols, paste0("SAMPLE_", seq_len(length(col_names) - 9))))
  
  # Keep only the columns we need (first 10 columns at most)
  n_cols <- min(ncol(vcf_data), 10)
  vcf_data <- vcf_data[, 1:n_cols]
  
  # Convert to tibble
  vcf_data <- as_tibble(vcf_data)
  
  # Extract INFO field components
  vcf_data <- vcf_data %>%
    mutate(
      QUAL = as.numeric(as.character(QUAL)),
      # Extract DP (depth)
      DP = as.numeric(str_extract(INFO, "DP=([0-9]+)") %>% str_remove("DP=")),
      # Extract MQ (mapping quality)
      MQ = as.numeric(str_extract(INFO, "MQ=([0-9]+)") %>% str_remove("MQ=")),
      # Extract MQ0F (fraction of MQ0 reads)
      MQ0F = as.numeric(str_extract(INFO, "MQ0F=([0-9.]+)") %>% str_remove("MQ0F=")),
      # Extract AF1 (allele frequency)
      AF1 = as.numeric(str_extract(INFO, "AF1=([0-9.]+)") %>% str_remove("AF1=")),
      # Extract FQ (quality)
      FQ = as.numeric(str_extract(INFO, "FQ=(-?[0-9.]+)") %>% str_remove("FQ=")),
      # Extract VDB (variant distance bias)
      VDB = as.numeric(str_extract(INFO, "VDB=([0-9.]+)") %>% str_remove("VDB=")),
      # Extract DP4 (strand counts)
      DP4 = str_extract(INFO, "DP4=([0-9,]+)") %>% str_remove("DP4="),
      # Check if INDEL
      IS_INDEL = str_detect(INFO, "INDEL")    )
  
  # Calculate strand bias if DP4 is available
  vcf_data <- vcf_data %>%
    rowwise() %>%
    mutate(
      strand_bias = if_else(
        !is.na(DP4) && DP4 != "",
        {
          dp4_vals <- as.numeric(str_split(DP4, ",")[[1]])
          
          if (length(dp4_vals) == 4) {
            ref_fwd <- dp4_vals[1]
            ref_rev <- dp4_vals[2]
            alt_fwd <- dp4_vals[3]
            alt_rev <- dp4_vals[4]
            
            # Calculate strand bias ratio
            ref_ratio <- if_else(ref_fwd + ref_rev > 0, 
                                 min(ref_fwd, ref_rev) / (ref_fwd + ref_rev),
                                 NA_real_)
            alt_ratio <- if_else(alt_fwd + alt_rev > 0,
                                 min(alt_fwd, alt_rev) / (alt_fwd + alt_rev),
                                 NA_real_)
            
            min(ref_ratio, alt_ratio, na.rm = TRUE)
          } else {
            NA_real_
          }
        },
        NA_real_
      )
    ) %>%
    ungroup()
  
  return(vcf_data)
}

# Function to apply quality filters
filter_vcf <- function(vcf_data, sample_name, cfg = config) {
  
  if (is.null(vcf_data)) return(NULL)
  
  original_count <- nrow(vcf_data)
  
  filtered <- vcf_data %>%
    mutate(
      pass_qual = QUAL >= cfg$min_qual,
      pass_depth = DP >= cfg$min_depth & DP <= cfg$max_depth,
      pass_mq = is.na(MQ) | MQ >= cfg$min_mq,
      pass_mq0f = is.na(MQ0F) | MQ0F <= cfg$max_mq0f,
      pass_af = is.na(AF1) | AF1 >= cfg$min_af,
      pass_fq = is.na(FQ) | FQ <= cfg$min_fq,
      pass_vdb = is.na(VDB) | VDB >= cfg$vdb_threshold,
      pass_strand = is.na(strand_bias) | strand_bias >= 0.2,  # At least 20% on minority strand
      
      # Overall pass
      PASS = pass_qual & pass_depth & pass_mq & pass_mq0f & 
             pass_af & pass_fq & pass_vdb & pass_strand
    )
  
  passed <- filtered %>% filter(PASS)
  
  # Summary statistics
  summary_stats <- tibble(
    sample = sample_name,
    original_variants = original_count,
    passed_variants = nrow(passed),
    pct_passed = round(100 * nrow(passed) / original_count, 2),
    failed_qual = sum(!filtered$pass_qual),
    failed_depth = sum(!filtered$pass_depth),
    failed_mq = sum(!filtered$pass_mq),
    failed_mq0f = sum(!filtered$pass_mq0f),
    failed_af = sum(!filtered$pass_af),
    failed_vdb = sum(!filtered$pass_vdb),
    failed_strand = sum(!filtered$pass_strand)
  )
  
  return(list(
    filtered_data = passed,
    all_data = filtered,
    summary = summary_stats
  ))
}

# Main execution
main <- function(vcf_dir = "BCF_VCF_FILES", output_dir = "FILTERED_VCF") {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all VCF files
  vcf_files <- list.files(vcf_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE)
  
  if (length(vcf_files) == 0) {
    stop("No VCF files found in directory: ", vcf_dir)
  }
  
  message(paste("Processing", length(vcf_files), "VCF files..."))
  
  # Process all files
  all_summaries <- list()
  
  for (vcf_file in vcf_files) {
    sample_name <- basename(vcf_file) %>% str_remove("\\.vcf\\.gz$")
    message(paste("\nProcessing:", sample_name))
    
    # Parse and filter
    vcf_data <- parse_vcf(vcf_file)
    results <- filter_vcf(vcf_data, sample_name)
    
    if (!is.null(results)) {
      # Save filtered VCF
      output_file <- file.path(output_dir, paste0(sample_name, ".filtered.tsv"))
      write_tsv(results$filtered_data, output_file)
      
      # Store summary
      all_summaries[[sample_name]] <- results$summary
      
      message(paste("  Original:", results$summary$original_variants,
                    "| Passed:", results$summary$passed_variants,
                    paste0("(", results$summary$pct_passed, "%)")))
    }
  }
  
  # Combine and save summaries
  summary_df <- bind_rows(all_summaries)
  write_csv(summary_df, file.path(output_dir, "filtering_summary.csv"))
  
  # Create visual summary
  summary_df <- summary_df %>%
    separate(sample, into = c("status", "dose", "timepoint", "replicate"), 
             sep = "_", remove = FALSE, fill = "right")
  
  # Plot filtering results
  p1 <- ggplot(summary_df, aes(x = sample, y = passed_variants, fill = status)) +
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Variants Passing Quality Filters",
         x = "Sample", y = "Number of Variants") +
    scale_fill_brewer(palette = "Set2")
  
  ggsave(file.path(output_dir, "variants_passed.png"), p1, width = 14, height = 6)
  
  # Plot pass rate
  p2 <- ggplot(summary_df, aes(x = status, y = pct_passed, fill = timepoint)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Variant Pass Rate by Sample Type",
         x = "Sample Status", y = "% Variants Passed")
  
  ggsave(file.path(output_dir, "pass_rate.png"), p2, width = 8, height = 6)
  
  message("\n=== Filtering Complete ===")
  message(paste("Results saved to:", output_dir))
  message(paste("Summary saved to:", file.path(output_dir, "filtering_summary.csv")))
  
  return(summary_df)
}

# Run the analysis
if (!interactive()) {
  # If running as script
  args <- commandArgs(trailingOnly = TRUE)
  vcf_dir <- ifelse(length(args) > 0, args[1], "BCF_VCF_FILES")
  output_dir <- ifelse(length(args) > 1, args[2], "FILTERED_VCF")
  
  results <- main(vcf_dir, output_dir)
} else {
  # If running interactively
  message("Run main(vcf_dir, output_dir) to process VCF files")
  message("Example: results <- main('BCF_VCF_FILES', 'FILTERED_VCF')")
}