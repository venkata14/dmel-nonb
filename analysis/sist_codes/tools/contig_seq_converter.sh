#!/bin/bash

# Purpose of script
# Take in a .fasta file and contig and outputs a continous string of the sequences with out the header
# -f: path to fasta file
# -c: name of contig

while getopts c:f: flag
do
    case "${flag}" in
        c) CONTIG=${OPTARG};;
        f) FASTA_FILE=${OPTARG};;
    esac
done

module load samtools

# Indexes the file. Creates a $FASTA_FILE.fai
samtools faidx $FASTA_FILE
# Extract the entire contig from the fasta including the header. Want to remove that header
# The sed command removes the first line
SEQUENCES=$(samtools faidx $FASTA_FILE $CONTIG | sed '1d')

SEQUENCES=${SEQUENCES//$'\n'/} # Remove all new line
SEQUENCES=${SEQUENCES%$'\n'} # Remove trailing new lines

# Output SEQUENCES
echo "$SEQUENCES"


