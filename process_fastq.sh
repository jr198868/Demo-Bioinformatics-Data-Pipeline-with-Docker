#!/bin/bash

# Preprocess the FastQ files with Trimmomatic
trimmomatic PE -threads 4 read1.fastq read2.fastq \
    output_forward_paired.fastq output_forward_unpaired.fastq \
    output_reverse_paired.fastq output_reverse_unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Perform taxonomic classification with Kraken 2
kraken2 --db standard --paired \
    output_forward_paired.fastq output_reverse_paired.fastq \
    --output kraken2_output.txt
