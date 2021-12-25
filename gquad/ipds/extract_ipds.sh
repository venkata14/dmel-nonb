#!/bin/bash

FULL_PATH = ""
PATH_TO_REF = ""

python /$FULL_PATH/kineticsTools/kineticsTools/ipdSummary.py /$PWD/final.merge.pbalign.bam --reference /$PATH_TO_REF/dm6.fasta --numWorkers 12 --csv /$PWD/final_ipds.csv --useChemistry "P5-C3"
