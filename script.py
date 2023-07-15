import gzip
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Input files
read1_file = 'input/read1.fastq.gz'
read2_file = 'input/read2.fastq.gz'

# Output files
report_file = 'output/report.csv'
plot_file = 'output/quality_scores_plot.png'

# Constants
PHRED_OFFSET = 33

# Function to calculate quality scores
def calculate_quality_scores(qualities):
    return [-10 * np.log10(10 ** ((ord(q) - PHRED_OFFSET) / 10)) for q in qualities]

# Function to parse FASTQ file and extract information
def parse_fastq(file):
    with gzip.open(file, 'rt') as f:
        lines = f.readlines()
        sequence = lines[1].strip()
        qualities = lines[3].strip()
        return sequence, qualities

# Function to perform basic adapter trimming
def trim_sequence(sequence, adapter):
    adapter_len = len(adapter)
    if sequence.startswith(adapter):
        sequence = sequence[adapter_len:]
    if sequence.endswith(adapter):
        sequence = sequence[:-adapter_len]
    return sequence

# Parse read 1
with gzip.open(read1_file, 'rt') as f:
    lines = f.readlines()
    read1_sequence = lines[1].strip()
    read1_qualities = lines[3].strip()
    read1_sequence = trim_sequence(read1_sequence, 'ADAPTER_SEQ')
    read1_qualities = read1_qualities[:len(read1_sequence)]
    read1_quality_scores = calculate_quality_scores(read1_qualities)

# Parse read 2
with gzip.open(read2_file, 'rt') as f:
    lines = f.readlines()
    read2_sequence = lines[1].strip()
    read2_qualities = lines[3].strip()
    read2_sequence = trim_sequence(read2_sequence, 'ADAPTER_SEQ')
    read2_qualities = read2_qualities[:len(read2_sequence)]
    read2_quality_scores = calculate_quality_scores(read2_qualities)

# Ensure sequences and qualities have the same length
min_length = min(len(read1_sequence), len(read2_sequence))
read1_sequence = read1_sequence[:min_length]
read2_sequence = read2_sequence[:min_length]
read1_qualities = read1_qualities[:min_length]
read2_qualities = read2_qualities[:min_length]
read1_quality_scores = read1_quality_scores[:min_length]
read2_quality_scores = read2_quality_scores[:min_length]

# Calculate statistics
molecules_sequenced = len(read1_sequence)
read1_length = len(read1_sequence)
read2_length = len(read2_sequence)
mean_insert_size = 0  # Calculate mean insert size
std_dev_insert_size = 0  # Calculate standard deviation of insert size
library_type = 'Nextera'  # Replace with appropriate library type
refseq_accession = 'ABC123'  # Replace with RefSeq or Genbank accession number

# Generate report
report_data = {
    'Number of molecules sequenced': [molecules_sequenced],
    'Sequencing configuration': ['Read 1', 'Index 1', 'Index 2', 'Read 2'],
    'Read lengths': [read1_length, 0, 0, read2_length],
    'Library type': [library_type],
    'Mean insert size': [mean_insert_size],
    'Standard deviation of insert size': [std_dev_insert_size],
    'RefSeq or Genbank accession numbers': [refseq_accession]
}

df_report = pd.DataFrame(report_data)
df_report.to_csv(report_file, index=False)

# Generate plot
positions = range(1, read1_length + 1)
plt.plot(positions, read1_quality_scores, label='Instrument-reported Phred scores')
plt.plot(positions, read1_quality_scores, label='Alignment inferred quality scores')
plt.xlabel('Position (Base Pairs)')
plt.ylabel('Quality Score')
plt.legend()
plt.savefig(plot_file)

# Display the plot
plt.show()
