#!/bin/bash

# This file is the sort the final merged file from the alignment steps

samtools sort final.merge.pbalign.bam -o sorted.final.merge.pbalign.bam -@ 16
samtools index sorted.final.merge.pbalign.bam -@ 16
samtools view -h -@ 16 -o sorted.final.merge.pbalign.sam sorted.final.merge.pbalign.bam
