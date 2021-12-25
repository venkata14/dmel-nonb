#!/bin/bash

(head -n1 sorted.final.merge.pbalign.sam; grep "PL:PACBIO" sorted.final.merge.pbalign.sam; grep "PN:bax2bam" sorted.final.merge.pbalign.sam) > sorted_alignment_metadata.txt

less all_contig_names.txt | sed -e "s/>//g" | while read line
do
echo "Running on $line"
OUTPUT_SAM="split-files/$line.sam"
OUTPUT_BAM="split-files/$line.bam"
( less sorted_alignment_metadata.txt ; grep "$line" sorted.final2.merge.pbalign.sam) > $OUTPUT_SAM
samtools view -@ 2 -h -bS $OUTPUT_SAM > $OUTPUT_BAM
samtools index -@ 2 $OUTPUT_BAM
pbindex $OUTPUT_BAM
done

echo "done"

