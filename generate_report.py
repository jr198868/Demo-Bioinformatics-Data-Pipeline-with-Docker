import subprocess

# Paths to the input files
insert_size_file = '/app/insert_size_metrics.txt'
isolate_output = '/app/mash_isolate_output_sorted.txt'
flagstat = '/app/flagstat.txt'
fastq_file_read1 = '/app/read1.fastq'
fastq_file_read2 = '/app/read2.fastq'
report_file = '/output/report.csv'

# parse the samtools flagstat results
molecules_sequenced_data, molecules_sequenced = [], ''
with open(flagstat) as f:
    for line in f:
        molecules_sequenced_data.append(line.strip())
        
for i in molecules_sequenced_data:     
    if '(QC-passed reads + QC-failed reads)' in i:
        molecules_sequenced = i

for i in molecules_sequenced_data:  
    if 'read1' in i:
        print(i)
        molecules_sequenced =  molecules_sequenced + '; '+ 'read1: ' + i.split('+')[0]
        
for i in molecules_sequenced_data:     
    if 'read2' in i:
        molecules_sequenced = molecules_sequenced + '; '+ 'read2:' + i.split('+')[0]

# parse fastq data to get sequencing configuration
def get_fastq_index_length(fastq_file):
    with open(fastq_file, "r") as f:
        first_line = f.readline().strip()
        index_length = len(first_line.split(":")[-1])
    return index_length


def get_fastq_read_length(fastq_file):
    with open(fastq_file, "r") as f:
        for line in f:
            if line.startswith("@"):
                read_length = len(f.readline().strip())
                return read_length

read1_length = get_fastq_read_length(fastq_file_read1)
index1_length = get_fastq_index_length(fastq_file_read1)
read2_length = get_fastq_read_length(fastq_file_read2)
index2_length = get_fastq_index_length(fastq_file_read2)


# parse the Mean library insert size 
insertsize = ''
with open(insert_size_file) as f:
    for line in f:
        insertsize = line.strip().split(',')
        
# parse the NCBI RefSeq assembly of sequenced organism 
RefSeq_genome = ''
mash_distance = ''
with open (isolate_output) as f:
    for line in f:
        RefSeq_genome = '_'.join(line.strip().split()[1].split('_')[0:3])
        mash_distance = line.strip().split()[2]
        break
        



# Perform necessary computations and generate the report
def generate_report(read1_length, index1_length, read2_length, index2_length, molecules_sequenced, RefSeq_genome, mash_distance, insertsize):
    data = [
        [
            "Number of molecules sequenced",
            "Sequencing configuration",
            "Library type",
            "Mean library insert size",
            "Standard deviation of library insert size",
            "RefSeq or Genbank accession numbers"
        ],
        [
            molecules_sequenced,
            "Read1: {}bp; Index1: {}bp; Read2: {}bp; Index2: {}bp".format(read1_length, index1_length, read2_length, index2_length),
            "TruSeq",
            insertsize[0].split('Standard Deviation:')[0].split(':')[1],
            insertsize[0].split('Standard Deviation:')[1],
            "RefSeq genome: {}; MASH distance: {}".format(RefSeq_genome, mash_distance)
        ]
    ]

    # Save the data as a CSV file
    with open(report_file, 'w') as file:
        for row in data:
            file.write(','.join(str(cell) for cell in row) + '\n')

# Generate the report
generate_report(read1_length, index1_length, read2_length, index2_length, molecules_sequenced, RefSeq_genome, mash_distance, insertsize)
