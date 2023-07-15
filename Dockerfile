FROM ubuntu:latest

# Install dependencies
RUN apt-get update && apt-get install -y openjdk-11-jre-headless wget unzip pigz perl libfindbin-libs-perl bwa samtools default-jdk bedtools seqtk mash 

# Download and install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    ln -s /FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.11.9.zip

# Download Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    rm Trimmomatic-0.39.zip

# Download Picard Tools
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.0/picard.jar -O /picard.jar

# Install bedtools

# Create the working directory
RUN mkdir /app

# Set the working directory
WORKDIR /app
# Set the output directory as a volume
VOLUME /output

# Copy the gzipped paired-end FASTQ files to the container
COPY /input_data/read1.fastq /app/read1.fastq
COPY /input_data/read2.fastq /app/read2.fastq
COPY /input_data/TruSeq3-PE.fa /app/TruSeq3-PE.fa
COPY /input_data/GCF_000155815.1_ASM15581v1_genomic.fna /app/reference.fa
COPY /input_data/refseq.genomes.k21s1000.msh /app/refseq.genomes.k21s1000.msh
COPY /generate_report.py /app/generate_report.py

# Run FastQC on the paired-end FASTQ files
RUN fastqc /app/read1.fastq /app/read2.fastq

# Run Trimmomatic on the paired-end FASTQ files
RUN java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 \
     /app/read1.fastq /app/read2.fastq /app/trim_read1_paired.fastq \
     /app/trim_read1_unpaired.fastq /app/trim_read2_paired.fastq \
     /app/trim_read2_unpaired.fastq \
     ILLUMINACLIP:/app/TruSeq3-PE.fa:2:30:20 LEADING:3 TRAILING:3 \
     SLIDINGWINDOW:4:15 MINLEN:36


#####################################################################################################################


# Align the trimmed paired-end FASTQ files using BWA
RUN bwa index /app/reference.fa
RUN bwa mem -t 4 /app/reference.fa /app/trim_read1_paired.fastq /app/trim_read2_paired.fastq > /app/trim_read1_read2_alig.sam

# Convert SAM to BAM
RUN samtools view -bS /app/trim_read1_read2_alig.sam > /app/trim_read1_read2_alig.bam


# Convert read1.fastq and read2.fastq to fasta
RUN seqtk mergepe /app/trim_read1_paired.fastq /app/trim_read2_paired.fastq > /app/trim_read1_read2_paired.fasta

# Generate MASH sketches for the paired-end fasta files
RUN mash sketch -p 4 -o /app/isolate.msh /app/trim_read1_read2_paired.fasta

# create mash sketch for the reference genome
# mash sketch --refseq GCF_000155815.1_ASM15581v1_genomic.fna.gz

# Compare isolate's MASH sketch to the reference database 
Run mash dist -p 4 isolate.msh refseq.genomes.k21s1000.msh > /app/mash_isolate_output.txt

# sort the mash output file based on 3th column, less number = higher similar distance
Run sort -k3,3 -n /app/mash_isolate_output.txt > /app/mash_isolate_output_sorted.txt

# get the RefSeq accession numbers of the sequenced organism.
Run head -1 /app/mash_isolate_output_sorted.txt > RefSeq_accession_numbers.txt

# Calculate the insert sizes
RUN bedtools bamtobed -i /app/trim_read1_read2_alig.bam | awk '{print $3 - $2}' > /app/insert_sizes.txt

# Calculate the mean and standard deviation
RUN awk '{ sum += $1; sumsq += ($1)^2 } END { print "Mean:", sum/NR, "Standard Deviation:", sqrt(sumsq/NR - (sum/NR)^2) }' /app/insert_sizes.txt > /app/insert_size_metrics.txt


# Get the number of molecules sequenced using samtools flagstat
RUN samtools flagstat /app/trim_read1_read2_alig.bam > /app/flagstat.txt

# Copy the BAM file to the output directory
RUN mkdir /output
RUN cp /app/trim_read1_read2_alig.bam /output/trim_read1_read2_alig.bam

# run generate_report.py
RUN python3 generate_report.py

# Copy the entire /app directory from the container to the host machine
VOLUME /app

# EXPOSE port 80 for HTTP traffic
EXPOSE 80

# Set the entrypoint to an interactive shell
CMD ["/bin/bash"]
