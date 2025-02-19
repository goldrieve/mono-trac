#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <reads_file>"
    exit 1
fi

# Define the paths to the input files and output directory
READS_FILE=$1
REFERENCE_FILE="/Users/goldriev/mono-trac/primer_design/targets_sequence.fasta"
OUTPUT_DIR="/Users/goldriev/mono-trac/test_reads"
SAM_FILE="${OUTPUT_DIR}/mapped_reads.sam"
BAM_FILE="${OUTPUT_DIR}/mapped_reads.bam"
SORTED_BAM_FILE="${OUTPUT_DIR}/mapped_reads_sorted.bam"
MAPPED_READS_FILE="${OUTPUT_DIR}/mapped_reads.fastq"

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Step 1: Map the reads using minimap2
minimap2 -ax map-ont ${REFERENCE_FILE} ${READS_FILE} > ${SAM_FILE}

# Step 2: Convert the SAM file to a BAM file and sort it
samtools view -bS ${SAM_FILE} | samtools sort -o ${SORTED_BAM_FILE}

# Step 3: Filter for mapped reads
samtools view -b -F 4 ${SORTED_BAM_FILE} > ${BAM_FILE}

# Step 4: Extract the mapped reads
samtools fastq ${BAM_FILE} > ${MAPPED_READS_FILE}

# Clean up intermediate files
rm ${SAM_FILE} ${SORTED_BAM_FILE}

echo "Mapped reads have been extracted to ${MAPPED_READS_FILE}"
