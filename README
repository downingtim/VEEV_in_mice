# VEEV RNA-Seq Analysis in Mouse Cell Lines

## Overview

This repository contains a complete bioinformatics pipeline for analyzing RNA-seq data from mouse cell lines infected with Venezuelan Equine Encephalitis Virus (VEEV). The analysis includes quality control, read filtering, viral/host discrimination, variant calling, and differential gene expression analysis.

**Experiment Design:** The dataset consists of 9 samples across three treatment groups:
- **Control** (3 samples): Mock-infected cells
- **WT** (3 samples): Cells infected with wild-type VEEV
- **4X** (3 samples): Cells infected with a 4X-mutant VEEV variant

## Repository Structure

```
VEEV_in_mice/
├── FIGURES/                 # Output figures and visualizations
├── TABLES/                  # Processed count and expression tables
├── VCFS/                    # Variant call files in TSV format
├── REFERENCES/              # Reference genome sequences
├── scripts/                 # Analysis scripts (detailed below)
├── LICENSE
└── README.md
```

## Reference Genomes

The `REFERENCES/` folder contains multiple VEEV reference sequences:
- **68U201.fa / 68U201.gff**: Original 68U201 VEEV strain reference and annotation
- **68U201_v.fa**: Variant version of 68U201 strain
- **68U201_7141.fa / 68U201_7141.gff**: 68U201 strain with 7141 variant
- **68U201_7141_v.fa**: Combined variant of 68U201 and 7141 strain
- **.fai files**: FASTA index files for each reference sequence

## Analysis Pipeline

### Phase 1: Quality Control and Read Preprocessing

**00_Snakemake.qc** (Snakemake workflow)
- Automated workflow orchestrating quality control steps
- Reads configuration from `config.yaml`
- Performs FastQC on raw and trimmed reads
- Outputs: FastQC reports and MultiQC aggregated quality reports

**01_runqc.illumina.sh** (Bash script)
- Quality control for individual Illumina sequencing runs
- Runs `fastp` for adapter detection, poly-G/poly-X trimming, and quality filtering
- Filters reads with mean quality < 20 and length < 30 bp
- Applies additional quality trimming using `fastq_quality_trimmer` (quality threshold: 30, minimum length: 60 bp)
- Outputs: HTML and JSON reports for each sample pair

**02_kraken.illumina.sh** (Bash script)
- Taxonomic classification using Kraken2
- Classifies reads against comprehensive reference database
- Processes paired-end reads for each sample
- Generates Kraken reports with species abundance
- Outputs: Classification files and abundance reports in `KRAKEN_FILES/` and `KRAKEN_REPORT/`

**03_kraken.illumina.extract.sh** (Bash script)
- Extracts reads belonging to specific taxa from Kraken2 output
- Separates VEEV reads (taxon ID: 11036) and mouse reads (taxon ID: 862507)
- Uses `extract_kraken_reads.py` with `--include-children` and `--include-parents` flags
- Outputs: Classified FASTQ files in `KRAKEN_VALID_FILES/` directory

### Phase 2: Read Mapping and Variant Calling

**05_mapping.sh** (Bash script)
- Maps quality-filtered VEEV reads to reference genomes
- Processes samples against two FASTA references: `68U201_v.fa` and `68U201_7141_v.fa`
- Implements comprehensive SAM/BAM processing pipeline:
  - **Mapping:** Minimap2 with short-read mode (`-ax sr`)
  - **SAM to BAM conversion:** samtools view with compression
  - **Quality control:** Sorts by query name, fixes mate information, sorts by coordinate
  - **Duplicate marking:** `samtools markdup` with removal of secondary and supplementary alignments
  - **Indexing and statistics:** Generates BAM index, flagstat, and coverage reports
- **Variant calling:** Dual approach using both FreeBayes and BCFtools
  - FreeBayes: Sensitive variant calling (allele frequency threshold: 0.1%, minimum 1 alternate read)
  - BCFtools: `mpileup` → `call` with haploid ploidy
  - Variant normalization: `bcftools norm` with `-m-both` to split complex variants
- **Consensus sequences:** Generates two consensus variants (LA and LR haplotypes)
- Outputs: BAM files, flagstat reports, coverage files, and normalized VCF files

### Phase 3: SNP Processing and Differential Expression

**06_get_snps.R** (R script)
- Extracts SNP information from normalized BCFtools VCF files
- Filters variants for quality and functional impact
- Processes variants for all samples and reference combinations
- Outputs: SNP annotations and summaries

**07_VCF_process.R** (R script)
- Converts processed SNP data into tab-separated format
- Formats variant information for downstream analysis
- Creates clean, queryable variant tables
- Outputs: TSV files in `VCFS/` directory

**code.R** (R script - Primary Analysis)
- **Differential expression analysis** using two complementary methods:
  - **Limma-voom:** Analyzes host (mouse) gene expression
    - Generates DEG lists for three comparisons: Control vs WT, Control vs 4X, WT vs 4X
    - Outputs: `limma.DE.c_4x.csv`, `limma.DE.wt_4x.csv`, `limma.DE.wt_c.csv`
  - **Sleuth:** Transcript-level analysis with bootstrap resampling for uncertainty quantification
    - Outputs: `sleuth.DE.genes.csv`
- Reads host read counts from `count_host.txt` (generated by 04_viz.levels.R)
- Handles read count normalization and filtering
- Outputs: Log-fold change estimates, p-values, and statistical summaries

**04_viz.levels.R** (R script)
- Calculates VEEV/host read ratio from Kraken2 classification results
- Processes raw FASTQ read counts from `KRAKEN_VALID_FILES/`
- Groups samples by treatment (Control, WT, 4X) and replicate number
- **Statistical analysis:** ANOVA with replicate as blocking factor, post-hoc emmeans comparisons
- **Visualization:** 
  - Box plots with individual points showing viral percentage per group
  - Generates both standard and y-axis limited versions for clarity
- Outputs:
  - `count_host.txt`: Sample-level viral and host read counts (input for code.R)
  - `VEEV_mouse_summary.txt`: Group-level summary statistics
  - `VEEV_mouse_boxplot.pdf` and `VEEV_mouse_boxplot.ylim.pdf`: Comparative plots

### Phase 4: Visualization and Reporting

**08_depth.py** (Python script)
- Analyzes sequencing depth across reference sequences
- Generates depth-of-coverage visualizations
- Helps identify regions with low/high coverage

**09_viz.py.R** (R script)
- Creates publication-quality figures from differential expression results
- Generates volcano plots, MA plots, and heatmaps
- Visualizes log fold-change correlations between comparisons

**10_viz.hist.R** (R script)
- Generates histograms and distribution plots
- Visualizes depth-of-coverage distributions
- Summarizes variant burden across samples

## Output Files

### TABLES/
- **count_host.txt**: Read counts (viral vs mouse) per sample
- **samples.txt**: Sample metadata and group assignments
- **limma.DE.c_4x.csv**: Limma differential expression (Control vs 4X)
- **limma.DE.wt_4x.csv**: Limma differential expression (WT vs 4X)
- **limma.DE.wt_c.csv**: Limma differential expression (WT vs Control)
- **sleuth.DE.genes.csv**: Sleuth transcript-level DEG analysis

### VCFS/
- **{Sample}.{Reference}.norm.bcf.filtered.tsv**: Normalized variant calls in TSV format
  - Files for each sample-reference combination
  - Two references per sample: 68U201_v and 68U201_7141_v

### FIGURES/
- **VEEV_mouse_boxplot.pdf**: Viral read percentage by group (with outliers)
- **VEEV_mouse_boxplot.ylim.pdf**: Same plot with y-axis limited for clarity
- **volcano_combined.pdf**: Volcano plots for differential expression
- **logFC_correlation_plots.pdf/png**: Correlation matrices of log fold-changes
- **plot_sample_heatmap.pdf**: Hierarchical clustering heatmap of samples
- **rdaf_*.pdf**: Read depth and variant frequency visualizations

## Usage

### Running the Complete Pipeline

1. **Prepare inputs:**
   - Raw Illumina FASTQ files in `ILLUMINA_FASTQS2/` directory
   - Reference sequences in `REFERENCES/` directory

2. **Execute QC and preprocessing:**
   ```bash
   snakemake -s 00_Snakemake.qc --cores 8
   for sample in $(ls ILLUMINA_FASTQS2 | sort | uniq); do
     bash 01_runqc.illumina.sh $sample
   done
   ```

3. **Run taxonomic classification:**
   ```bash
   bash 02_kraken.illumina.sh
   bash 03_kraken.illumina.extract.sh
   ```

4. **Map reads and call variants:**
   ```bash
   for sample in Control_RNA_1_S3 Control_RNA_2_S6 Control_RNA_3_S9 VEEV_WT_RNA_1_S1 VEEV_WT_RNA_2_S4 VEEV_WT_RNA_3_S7 VEEV_4X_RNA_1_S2 VEEV_4X_RNA_2_S5 VEEV_4X_RNA_3_S8; do
     bash 05_mapping.sh $sample
   done
   ```

5. **Process variants:**
   ```bash
   Rscript 06_get_snps.R
   Rscript 07_VCF_process.R
   ```

6. **Analyze host gene expression:**
   ```bash
   Rscript 04_viz.levels.R
   Rscript code.R
   ```

7. **Generate visualizations:**
   ```bash
   Rscript 08_depth.py
   Rscript 09_viz.py.R
   Rscript 10_viz.hist.R
   ```

## Key Parameters

- **Quality filtering:** Phred score ≥ 30, minimum length 60 bp
- **Variant calling:** Minimum allele frequency 0.1% (FreeBayes), 1% (BCFtools)
- **Mapping:** Minimap2 short-read mode optimized for Illumina sequencing
- **Consensus:** Both LA (latest-allele) and LR (latest-read) haplotypes generated

## Dependencies

- **QC:** FastQC, fastp, FastQ quality trimmer, MultiQC
- **Taxonomic:** Kraken2, KrakenTools
- **Mapping/Variant:** Minimap2, SAMtools, FreeBayes, BCFtools
- **Analysis:** R (limma, sleuth, tidyverse, ggplot2, emmeans), Python

## References

- VEEV reference: 68U201 strain and variants
- Mouse reference: Integrated from Kraken2 mouse database
- Variant calling methodology: BCFtools and FreeBayes best practices
- Statistical methods: Limma-voom for bulk RNA-seq, Sleuth for transcript-level analysis

## License

See LICENSE file for details.