#!/bin/bash

module load tools
file1=$1
file2=$1
trim1=$1
trim2=$1
fastp1=$1
fastp2=$1
fastphtml1=$1
fastphtml2=$1
stats1=$1
stats2=$1
file1+="_R1_001.fastq.gz" # S9_S9_R1_001
file2+="_R2_001.fastq.gz"
trim1+="_1_fastp_trim.fastq"
trim2+="_2_fastp_trim.fastq"
fastp1+="_1_fastp.fastq"
fastp2+="_2_fastp.fastq"
fastp1c=$1
fastp2c=$1
fastp1c+=".1.cut.fastq"
fastp2c+=".2.cut.fastq"
fastphtml1+="_1_fastp.html"
fastphtml2+="_2_fastp.html"
stats1+="_1.trim.txt"
stats2+="_2.trim.txt"
echo $file1 $file2 ;

if test -f $file1; then
#    fastqc $file1 &
#    fastqc $file2 &
    echo "fastqc $file1 "

    fastp --overrepresentation_analysis --detect_adapter_for_pe --trim_poly_g \
	  --cut_mean_quality 20 --trim_poly_x --length_required 30 --thread 4 \
	  --cut_window_size 4 --cut_front --cut_tail  --html $fastphtml1 --json ${fastphtml1%.html}.json \
	  -i $file1 -I $file2 -o $fastp1 -O $fastp2
    
#    echo "cutadapt -a "A{10}" -a "T{10}" -A "A{10}" -A "T{10}" -o $fastp1c -p $fastp2c $fastp1 $fastp2"
#    cutadapt -a "A{10}" -a "T{10}" -A "A{10}" -A "T{10}" -o $fastp1c -p $fastp2c $fastp1 $fastp2
#    rm -rf $fastp1 $fastp2
    
    fastq_quality_trimmer -Q 33 -t 30 -l 60 -i $fastp1 -o $trim1 -v  $stats1 &
    echo "~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 60 -i $fastp1 -o $trim1 -v  $stats1"
    echo "~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 60 -i $fastp2 -o $trim2 -v  $stats2"
    fastq_quality_trimmer -Q 33 -t 30 -l 60 -i $fastp2 -o $trim2 -v  $stats2
    rm -rf $fastp1c $fastp2c 
    
#    fastqc $trim1 &
    echo "fastqc $trim1"
#    fastqc $trim2
fi

echo $file1 $file2 $trim1 $trim2 $fastp1 $fastp2 $fastphtml1 $fastphtml2
