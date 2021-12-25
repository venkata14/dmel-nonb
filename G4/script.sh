#!/usr/bin/bash

# G4Hunter Requires Python 2
# You need to make the results directory and put in a random file that starts with the basename of the dna file
python G4Hunter.py -i ref/dm6.fasta -o results -w 25 -s 1
python G4Hunter.py -i ref/dm6.fasta -o results_higher_str -w 25 -s 1.5
