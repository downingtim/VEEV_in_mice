#!/bin/bash

sample=$1

# All FASTA files to process
FASTA_FILES=( "68U201_v.fa"
	      "68U201_7141_v.fa")

# Create necessary directories
mkdir -p SAM_FILES BAM_FILES BAM_FINAL COVERAGE FB_VCF_FILES BCF_VCF_FILES CONSENSUS ERROR_FILES

# Function to check if command succeeded
check_command() {
    if [ $? -ne 0 ]; then
        echo "ERROR: $1 failed for sample $sample with fasta $fasta"
        return 1
    fi
    return 0
}

# Function to check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Expected file $1 not found for sample $sample with fasta $fasta"
        return 1
    fi
    echo "SUCCESS: File $1 created"
    return 0
}

# Function to determine report folder based on fasta filename
get_report_folder() {
    local fasta_file="$1"
    local base_name=$(basename "$fasta_file" | sed 's/\.[^.]*$//')  # Remove file extension
    # FFASTA_FILES=( "68U201.fa", "68U201-7141.fa")
    
    case "$base_name" in
	"68U201")
            echo "REPORT_68U201"
            ;;
	"68U201_7141")
            echo "REPORT_68U201_7141"
            ;;
        *)
            echo "REPORT_Unknown"  # Default fallback
            ;;
    esac
}

echo "Processing sample: $sample against all FASTA references..."

# Process the sample against each FASTA file
for fasta in "${FASTA_FILES[@]}"; do
    echo "DEBUG: Processing FASTA file: $fasta"
    
    # Check if FASTA file exists
    if [ ! -f "$fasta" ]; then
        echo "WARNING: FASTA file '$fasta' not found - skipping"
        continue
    fi
    
    # Get the appropriate report folder
    report_folder=$(get_report_folder "$fasta")
    echo "DEBUG: Report folder determined as: $report_folder"
    
    # Check if report folder exists, create if it doesn't
    if [ ! -d "$report_folder" ]; then
        echo "Creating report folder: $report_folder"
        mkdir -p "$report_folder"
    fi
    
    echo "=== Processing $sample with $fasta ==="
    
    # Mapping Illumina reads with validation
    echo "Mapping Illumina reads..."
     minimap2 -ax sr "$fasta" "../KRAKEN_VALID_FILES/${sample}_1.VEEV.fq.gz" "../KRAKEN_VALID_FILES/${sample}_2.VEEV.fq.gz" > "SAM_FILES/${sample}.${fasta}.illumina.sam" 2> "ERROR_FILES/${sample}.illumina.sam.errors.txt"
    if ! check_command "Illumina mapping"; then continue; fi
    if ! check_file "SAM_FILES/${sample}.${fasta}.illumina.sam"; then continue; fi
    
    ls -lt "SAM_FILES/${sample}.${fasta}.illumina.sam"
    
    samtools view -bS "SAM_FILES/${sample}.${fasta}.illumina.sam" > "BAM_FILES/${sample}.${fasta}.illumina.bam"
    if ! check_file "BAM_FILES/${sample}.${fasta}.illumina.bam"; then continue; fi
    rm -rf "SAM_FILES/${sample}.${fasta}.illumina.sam"
    
    # Sort by queryname first
    echo "Sorting Illumina BAM by queryname..."
    samtools sort -n -o "BAM_FILES/${sample}.${fasta}.illumina.queryname.bam" "BAM_FILES/${sample}.${fasta}.illumina.bam"
    if ! check_file "BAM_FILES/${sample}.${fasta}.illumina.queryname.bam"; then continue; fi
    
    # Fix mate information
    echo "Fixing mate information for Illumina reads..."
    samtools fixmate -m "BAM_FILES/${sample}.${fasta}.illumina.queryname.bam" "BAM_FILES/${sample}.${fasta}.illumina.fixmate.bam"
    if ! check_file "BAM_FILES/${sample}.${fasta}.illumina.fixmate.bam"; then continue; fi

    # Sort by coordinate for markdup
    echo "Sorting Illumina BAM by coordinate..."
    samtools sort -o "BAM_FILES/${sample}.${fasta}.illumina.coordsort.bam" "BAM_FILES/${sample}.${fasta}.illumina.fixmate.bam"
    if ! check_file "BAM_FILES/${sample}.${fasta}.illumina.coordsort.bam"; then continue; fi

    # Mark duplicates
    echo "Marking duplicates in Illumina reads..."
    samtools markdup -r -s "BAM_FILES/${sample}.${fasta}.illumina.coordsort.bam" "BAM_FINAL/${sample}.${fasta}.merged.bam"
    if ! check_file "BAM_FINAL/${sample}.${fasta}.merged.bam"; then continue; fi
    rm -rf "BAM_FILES/${sample}.${fasta}.illumina.coordsort.bam" "BAM_FILES/${sample}.${fasta}.illumina.queryname.bam"
    
    samtools index "BAM_FINAL/${sample}.${fasta}.merged.bam"
    samtools flagstat "BAM_FINAL/${sample}.${fasta}.merged.bam" > "BAM_FINAL/${sample}.${fasta}.flagstat"
    samtools coverage "BAM_FINAL/${sample}.${fasta}.merged.bam" > "COVERAGE/${sample}.${fasta}.coverage.txt"
    
    # Variant calling with freebayes
    echo "Calling variants with freebayes..."
    freebayes -f "$fasta" -F 0.01 -p 1 --min-alternate-count 1 --min-alternate-fraction 0.001 "BAM_FINAL/${sample}.${fasta}.merged.bam" > "FB_VCF_FILES/${sample}.${fasta}.fb.vcf" 2> "ERROR_FILES/${sample}.${fasta}.fb.errors.txt"
    
    # Compress and index
    bgzip -f "FB_VCF_FILES/${sample}.${fasta}.fb.vcf"
    tabix -f -p vcf "FB_VCF_FILES/${sample}.${fasta}.fb.vcf.gz"
    
    # Normalize variants
    bcftools norm -c w -f "$fasta" -m-both -Oz -o "FB_VCF_FILES/${sample}.${fasta}.norm.fb.vcf.gz" "FB_VCF_FILES/${sample}.${fasta}.fb.vcf.gz"
    tabix -f -p vcf "FB_VCF_FILES/${sample}.${fasta}.norm.fb.vcf.gz"
    
    # BCFtools variant calling
    echo "Calling variants with BCFtools..."
    bcftools mpileup -A -Ob -f "$fasta" "BAM_FINAL/${sample}.${fasta}.merged.bam" | bcftools call -cvO v --ploidy 1 -o "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf"
    bgzip -f "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf"
    tabix -f -p vcf "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf.gz"
    
    # Normalize BCF variants
    bcftools norm -c w -f "$fasta" -m-both -Oz -o "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz" "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf.gz"
    tabix -f -p vcf "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz"
    
    # Generate consensus
    echo "Generating consensus sequences..."
    bcftools consensus -H LA -f "$fasta" "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz" -o "CONSENSUS/${sample}.${fasta}.LA.fasta"
    bcftools consensus -H LR -f "$fasta" "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz" -o "CONSENSUS/${sample}.${fasta}.LR.fasta"
    
    # Copy results to report folder
    cp "BAM_FINAL/${sample}.${fasta}.flagstat" "$report_folder/"
    cp "COVERAGE/${sample}.${fasta}.coverage.txt" "$report_folder/"
    cp "CONSENSUS/${sample}.${fasta}.LA.fasta" "$report_folder/"
    cp "CONSENSUS/${sample}.${fasta}.LR.fasta" "$report_folder/"
    echo "=== Completed $fasta for sample $sample ==="
    
done
echo "Completed processing $sample against all FASTA references"
