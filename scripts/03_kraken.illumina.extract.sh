#!/bin/bash

fastp_files_dir="ILLUMINA_FASTQS2"

base_command="python /mnt/lustre/RDS-live/downing/KrakenTools/extract_kraken_reads.py -k "
output_dir_base="KRAKEN_FILES"
output_dir_base2="KRAKEN_REPORT"
sample_names=$(ls $fastp_files_dir | sort | uniq)

for sample_name in $sample_names; do
    output_dir="${output_dir_base}/${sample_name}"
    output_dir2="${output_dir_base2}/${sample_name}"
    out_file=${sample_name}

    taxon=11036  # VEEV
#    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name} -o KRAKEN_VALID_FILES/${out_file} -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents & "
#    echo "${command}"
#    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name} -o KRAKEN_VALID_FILES/${out_file} -t $taxon --fastq-output -r $output_dir_base2/${sample_name}  --include-children --include-parents "
#    echo "${command}"

    taxon=862507 # mouse
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name} -o KRAKEN_VALID_FILES/${out_file}_1_mouse.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents & "
    echo "${command}"
#    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name} -o KRAKEN_VALID_FILES/${out_file}_2_mouse.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name}  --include-children --include-parents "
#    echo "${command}"

done
