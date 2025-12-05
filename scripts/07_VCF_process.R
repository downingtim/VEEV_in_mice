#!/usr/bin/env Rscript
# annotate_variants_with_AA.R
# Requires: tidyverse, data.table, purrr, stringr
# Optional (for AA effects): Biostrings

library(tidyverse)
library(data.table)
library(purrr)
library(stringr)

# --- IMPORTANT: SET PATHS TO YOUR REFERENCE FASTAs ---
# Update these paths to point to your actual FASTA files
ref_fasta_std  <- "../68U201_PROKKA/68U201.fasta"
ref_fasta_7141 <- "../68U201_7141_PROKKA/68U201_7141.fasta"

have_fasta <- !is.null(ref_fasta_std) && file.exists(ref_fasta_std) &&
              !is.null(ref_fasta_7141) && file.exists(ref_fasta_7141)

if (have_fasta) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings required for AA annotation. Install with: BiocManager::install('Biostrings')")
  }
  # Load Biostrings with library() to ensure all methods are loaded
  suppressPackageStartupMessages(library(Biostrings))
  
  fasta_std  <- readDNAStringSet(ref_fasta_std)
  fasta_7141 <- readDNAStringSet(ref_fasta_7141)
  message("Loaded FASTA files for AA prediction")
  
  # Print sequence names for debugging
  message("Sequences in std FASTA: ", paste(names(fasta_std), collapse=", "))
  message("Sequences in 7141 FASTA: ", paste(names(fasta_7141), collapse=", "))

  # helper to get seq by id: try matching several id styles
  get_ref_seq <- function(fasta_obj, seqid) {
    # try direct name match
    if (seqid %in% names(fasta_obj)) return(fasta_obj[[seqid]])
    
    # try removing _v suffix (common in VCF chromosome names)
    seqid_no_v <- str_replace(seqid, "_v$", "")
    if (seqid_no_v %in% names(fasta_obj)) return(fasta_obj[[seqid_no_v]])
    
    # try matching substring (for partial matches)
    idx <- which(str_detect(names(fasta_obj), fixed(seqid)))
    if (length(idx) >= 1) return(fasta_obj[[idx[1]]])
    
    # try matching without _v suffix
    idx <- which(str_detect(names(fasta_obj), fixed(seqid_no_v)))
    if (length(idx) >= 1) return(fasta_obj[[idx[1]]])
    
    return(NULL)
  }
} else {
  message("WARNING: FASTA files not found. AA prediction will be skipped (all NA).")
  message("Set ref_fasta_std and ref_fasta_7141 paths to enable AA prediction.")
}

# ---------------- helpers for reading variant tables safely ----------------

read_vcf_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)

  # cut off FASTA section if present
  fasta_start <- grep("^##FASTA", lines)
  if (length(fasta_start) > 0) lines <- lines[seq_len(fasta_start - 1)]

  # remove # header lines (keep header if present without #)
  lines <- lines[!startsWith(lines, "#")]

  if (length(lines) == 0) {
    message("No variant rows in: ", path)
    return(tibble())
  }

  df <- tryCatch({
    fread(text = lines, sep = "\t", header = TRUE, data.table = FALSE, fill = TRUE, showProgress = FALSE)
  }, error = function(e) {
    message("Could not parse: ", path)
    return(tibble())
  })

  as_tibble(df)
}

load_gff <- function(gff_path) {
  raw <- readLines(gff_path, warn = FALSE)
  raw <- raw[!grepl("^#", raw)]
  if (length(raw) == 0) return(tibble())

  gff <- fread(text = raw, sep = "\t", header = FALSE, data.table = FALSE, fill = TRUE, showProgress = FALSE)
  # ensure colnames
  colnames(gff)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  gff <- as_tibble(gff)

  gff %>%
    filter(type %in% c("gene", "CDS")) %>%
    mutate(
      locus_tag = str_extract(attributes, "locus_tag=[^;]+") %>% str_replace("^locus_tag=", "") %>% replace_na(NA_character_),
      gene_id   = str_extract(attributes, "ID=[^;]+")        %>% str_replace("^ID=", "")        %>% replace_na(NA_character_),
      product   = str_extract(attributes, "product=[^;]+")   %>% str_replace("^product=", "")   %>% replace_na(NA_character_)
    )
}

# load GFFs (your paths)
gff_std   <- load_gff("../68U201_PROKKA/68U201.gff")
gff_7141  <- load_gff("../68U201_7141_PROKKA/68U201_7141.gff")

# ---------------- sample name cleaning ----------------
clean_sample_name <- function(filepath) {
  nm <- basename(filepath) %>%
    # robustly strip the long suffix and variations
    str_remove("_v\\.fa(\\.norm)?(\\.bcf)?(\\.filtered)?(\\.tsv)?$") %>%
    str_remove("\\.fa(\\.norm)?(\\.bcf)?(\\.filtered)?(\\.tsv)?$") %>%
    str_replace_all("_RNA_?", "_") %>%
    str_remove("\\.tsv$") %>%
    str_remove("\\.$")
  nm
}

# ---------------- annotation helpers ----------------

# returns tibble with gene/locus/product and CDS start/end/strand/phase for the hit (first CDS hit)
find_cds_hit <- function(chrom, pos, reftag = c("std","7141")) {
  reftag <- reftag[1]
  gff_use <- if (str_detect(chrom, "7141") || reftag == "7141") gff_7141 else gff_std
  
  # Filter for CDS features that overlap the position
  hit <- gff_use %>% 
    filter(type == "CDS", start <= pos, end >= pos)
  
  if (nrow(hit) == 0) {
    return(tibble(Gene = NA_character_, Locus = NA_character_, Product = NA_character_,
                  CDS_start = NA_integer_, CDS_end = NA_integer_, Strand = NA_character_, Phase = NA_integer_))
  }
  hit1 <- hit[1,]
  tibble(
    Gene = hit1$gene_id,
    Locus = hit1$locus_tag,
    Product = hit1$product,
    CDS_start = as.integer(hit1$start),
    CDS_end = as.integer(hit1$end),
    Strand = hit1$strand,
    Phase = as.integer(ifelse(is.na(hit1$phase) | hit1$phase == "." | hit1$phase == "", 0, hit1$phase))
  )
}

# compute codon / AA change if FASTA provided
compute_aa_change <- function(chrom, pos, ref_nt, alt_nt, cds_start, cds_end, strand, phase, reftag = "std") {
  
  na_result <- tibble(Consequence = NA_character_, AA_pos = NA_integer_, 
                      Ref_AA = NA_character_, Alt_AA = NA_character_, 
                      AA_change = NA_character_)
  
  if (!have_fasta) {
    return(na_result)
  }

  # pick FASTA set
  fasta_obj <- if (str_detect(chrom, "7141") || reftag == "7141") fasta_7141 else fasta_std

  # try to get the matching seq
  seq_obj <- get_ref_seq(fasta_obj, chrom)
  if (is.null(seq_obj)) {
    # Try first sequence as fallback
    if (length(fasta_obj) > 0) {
      seq_obj <- fasta_obj[[1]]
    } else {
      return(na_result)
    }
  }
  
  # Validate that we got a sequence
  if (is.null(seq_obj)) {
    return(na_result)
  }

  # Validate inputs
  if (is.na(cds_start) || is.na(cds_end) || is.na(strand)) {
    return(na_result)
  }
  
  # Validate coordinates are reasonable
  if (cds_start < 1 || cds_end < cds_start) {
    return(na_result)
  }

  # extract CDS sequence as DNAString with error handling
  cds_seq <- tryCatch({
    subseq(seq_obj, start = cds_start, end = cds_end)
  }, error = function(e) {
    message("Error extracting sequence at ", chrom, ":", pos, " (CDS ", cds_start, "-", cds_end, ") - ", e$message)
    NULL
  })
  
  # Check if extraction failed
  if (is.null(cds_seq)) {
    return(na_result)
  }
  
  # Handle reverse strand
  if (strand == "-") {
    cds_seq <- tryCatch({
      reverseComplement(cds_seq)
    }, error = function(e) {
      message("Error reverse complementing at ", chrom, ":", pos, " - ", e$message)
      NULL
    })
    if (is.null(cds_seq)) return(na_result)
  }

  # Get sequence length - use direct slot access or nchar as fallback
  seq_length <- tryCatch({
    # First try: direct nchar conversion (most reliable)
    as.integer(nchar(as.character(cds_seq)))
  }, error = function(e) {
    message("Error getting sequence length at ", chrom, ":", pos, " - ", e$message)
    NULL
  })
  
  if (is.null(seq_length) || seq_length == 0) {
    return(na_result)
  }
  
  # compute position within CDS (1-based)
  if (strand == "+") {
    cds_nt_pos <- pos - cds_start + 1
  } else {
    cds_nt_pos <- cds_end - pos + 1
  }
  
  if (cds_nt_pos < 1 || cds_nt_pos > seq_length) {
    return(na_result)
  }

  # Handle phase (reading frame offset)
  frame_offset <- ifelse(is.na(phase), 0L, as.integer(phase) %% 3L)
  
  # Adjust CDS position for phase
  adjusted_pos <- cds_nt_pos - frame_offset
  if (adjusted_pos < 1) {
    return(na_result)
  }
  
  # compute codon index (0-based) and codon start within CDS (1-based)
  codon0_index <- (adjusted_pos - 1) %/% 3
  codon_start_pos <- codon0_index * 3 + 1 + frame_offset

  # ensure codon within bounds
  if (codon_start_pos + 2 > seq_length || codon_start_pos < 1) {
    return(na_result)
  }

  codon_seq <- tryCatch({
    as.character(Biostrings::subseq(cds_seq, start = codon_start_pos, width = 3))
  }, error = function(e) {
    message("Error extracting codon at ", chrom, ":", pos)
    return(NULL)
  })
  
  if (is.null(codon_seq) || nchar(codon_seq) != 3) {
    return(na_result)
  }
  
  # Position within the codon (1-3)
  pos_in_codon <- ((adjusted_pos - 1) %% 3) + 1
  
  # Build alternate codon
  codon_chars <- strsplit(codon_seq, "")[[1]]
  
  # For reverse strand, we need to complement the alt base since cds_seq is already rev-complemented
  if (strand == "-") {
    alt_base <- tryCatch({
      as.character(complement(DNAString(alt_nt)))
    }, error = function(e) {
      return(alt_nt)  # fallback to original if complement fails
    })
  } else {
    alt_base <- alt_nt
  }
  
  codon_chars[pos_in_codon] <- alt_base
  alt_codon <- paste0(codon_chars, collapse = "")

  # translate codons
  ref_aa <- tryCatch({
    as.character(translate(DNAString(codon_seq)))
  }, error = function(e) {
    "?"
  })
  
  alt_aa <- tryCatch({
    as.character(translate(DNAString(alt_codon)))
  }, error = function(e) {
    "?"
  })
  
  aa_pos <- codon0_index + 1  # 1-based amino acid position

  # Determine consequence
  consequence <- if (ref_aa == alt_aa) {
    "synonymous"
  } else if (alt_aa == "*" && ref_aa != "*") {
    "stop_gain"
  } else if (ref_aa == "*" && alt_aa != "*") {
    "stop_loss"
  } else {
    "nonsynonymous"
  }

  tibble(
    Consequence = consequence,
    AA_pos = aa_pos,
    Ref_AA = ref_aa,
    Alt_AA = alt_aa,
    AA_change = paste0(ref_aa, aa_pos, alt_aa)
  )
}

# ---------------- parse single variant file ----------------

parse_variant_file <- function(path) {
  df <- read_vcf_tsv(path)
  if (nrow(df) == 0) return(tibble())

  sample <- clean_sample_name(path)
  ref_type <- ifelse(str_detect(basename(path), "_7141_"), "7141", "std")

  # normalise column names
  names(df) <- names(df) %>% str_replace_all("^#", "")
  if (!"CHROM" %in% names(df) && "chrom" %in% names(df)) names(df)[names(df)=="chrom"] <- "CHROM"
  if (!"POS" %in% names(df) && "pos" %in% names(df)) names(df)[names(df)=="pos"] <- "POS"

  if (!("CHROM" %in% names(df) && "POS" %in% names(df))) {
    message("Missing CHROM or POS in: ", path)
    return(tibble())
  }

  # Prepare basic table
  tbl <- df %>%
    mutate(
      Sample = sample,
      RefType = ref_type,
      CHROM = as.character(CHROM),
      POS = as.integer(POS)
    ) %>%
    select(any_of(c("CHROM","POS","ID","REF","ALT","Sample","RefType")), everything())

  # For each row, find CDS info and compute AA change
  message("Processing ", nrow(tbl), " variants from ", basename(path))
  
  results <- pmap_dfr(list(tbl$CHROM, tbl$POS, tbl$REF, tbl$ALT), function(chrom, pos, ref_nt, alt_nt) {
    # Find CDS
    cds_info <- find_cds_hit(chrom, pos, reftag = ref_type)
    
    # Compute AA change
    if (!is.na(cds_info$CDS_start)) {
      aa_info <- compute_aa_change(chrom, pos, ref_nt, alt_nt, 
                                   cds_info$CDS_start, cds_info$CDS_end, 
                                   cds_info$Strand, cds_info$Phase, 
                                   reftag = ref_type)
    } else {
      aa_info <- tibble(Consequence = NA_character_, AA_pos = NA_integer_, 
                       Ref_AA = NA_character_, Alt_AA = NA_character_, 
                       AA_change = NA_character_)
    }
    
    bind_cols(cds_info, aa_info)
  })
  
  out <- bind_cols(tbl, results)

  # produce final output
  final <- out %>%
    transmute(
      CHROM = CHROM,
      POS = POS,
      ID = if("ID" %in% names(out)) out$ID else NA_character_,
      REF = if("REF" %in% names(out)) out$REF else NA_character_,
      ALT = if("ALT" %in% names(out)) out$ALT else NA_character_,
      Sample = Sample,
      Ref = RefType,
      Type = if("type" %in% names(out)) out$type else NA_character_,
      Gene = Gene,
      Locus = Locus,
      Product = Product,
      In_CDS = !is.na(CDS_start),
      CDS_nt_pos = ifelse(!is.na(CDS_start), 
                          if_else(Strand == "+", POS - CDS_start + 1, CDS_end - POS + 1),
                          NA_integer_),
      Consequence = Consequence,
      AA_pos = AA_pos,
      Ref_AA = Ref_AA,
      Alt_AA = Alt_AA,
      AA_change = AA_change
    )

  final
}

# ------------------- run on files -------------------

input_pattern <- "v\\.fa\\.norm\\.bcf\\.filtered\\.tsv$"
files <- list.files(path = ".", pattern = input_pattern, full.names = TRUE)
if (length(files) == 0) stop("No input files found matching pattern: ", input_pattern)
message("Found ", length(files), " input files — processing...")

all_variants <- map_df(files, parse_variant_file)

out_fn <- "all_variants_minimal_annotated.tsv"
write_tsv(all_variants, out_fn)
message("Wrote: ", out_fn, " with ", nrow(all_variants), " variants")