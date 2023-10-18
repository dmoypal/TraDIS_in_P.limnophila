#!/bin/bash
# Author: Ildefonso Cases
echo "WARNING: Make sure you have the environment with all the dependencies installed." 
#Concatenate all the FASTQ files in the current directory into a single file.
echo "Concatenating all the .fastq..."
cat *.fastq > all_data.fastq
echo -e "\033[32mDONE\033[0m" && sleep 8 
echo -e "\n\n"


#Use cutadapt to trim specific adapters from the sequences in the FASTQ file.
echo "Removing adapters with cutadapt..."
cutadapt -g XNNNNNGTTCGAAATGAGATGTGTATAAGAGACAG -e 0 -O 34 -o all_data.only_trans.trimmer.fastq --untrimmed-output all_data.only_trans.untrimmed.fastq all_data.fastq --cores=8

#Trim another adapter with cutadapt and remove read sequences that are shorter than 15 bases.
cutadapt --cores=10 -m 15 -a ATGGAATTCTCGGGTGCCAAGG -o all_data.only_trans.trimmer.no_adapt.fastq all_data.only_trans.trimmer.fastq
echo -e "\n\n\n"
echo -e "\033[32mDONE\033[0m"
echo -e "\n\n"


#Align the trimmed read sequences to the reference genome with Bowtie2
echo "Aligning sequences with bowtie2..."
bowtie2 -p 10 -x genomebowtieindex/GCA_000092105.1_ASM9210v1_genomic -U all_data.only_trans.trimmer.no_adapt.fastq -S trimmed_t2_b2_all.sam
echo -e "\n\n\n"
echo -e "\033[32mDONE\033[0m"
echo -e "\n\n"


#Convert the SAM file to BAM format using samtools
samtools view -b trimmed_t2_b2_all.sam -o trimmed_t2_b2_all.bam

#Sort the read alignments in the BAM file by their reference positions using samtools
samtools sort trimmed_t2_b2_all.bam  -o  trimmed_t2_b2_all.sorted.bam

#Create an index for the sorted BAM file
samtools index trimmed_t2_b2_all.sorted.bam

#Remove duplicates with picard
picard MarkDuplicates I=trimmed_t2_b2_all.sorted.bam O=trimmed2_t2_b2_all.sorted.dedup.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=trimmed_t2_b2_all.picard.txt

#Generate an index for the duplicate-free BAM file using samtools
samtools index trimmed2_t2_b2_all.sorted.dedup.bam  

