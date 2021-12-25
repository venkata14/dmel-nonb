#!/bin/bash

FULL_PATH = ""

less all_contig_names.txt | sed -e "s/>//g" | while read line
do
BAM_FILE="split-files/$line.bam"
IPD_FILE="ipd/$line.bam.ipd"
python /$FULL_PATH/kineticsTools/kineticsTools/ipdSummary.py $BAM_FILE --reference dm6.fasta --numWorkers 12 --csv $IPD_FILE --useChemistry "P5-C3"
echo "Extracted IPDs from $BAM_FILE"
done
