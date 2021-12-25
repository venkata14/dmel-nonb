#!/bin/bash

# create centromere results
CONTROL="control"
DNA="dna"
RESULTS="sist-results"
CURRENT_DIR="$PWD"

declare -a CENS
CENS=("C2" "C3" "C4" "CX" "CY")

# creating the results dirs
RESULTS_DIR="$CONTROL/$RESULTS"
if [ ! -d $RESULTS_DIR ]
then
    mkdir "$RESULTS_DIR"
fi

# Making the subfolders for individual centromeres
for m in ${CENS[@]}
do 
    SUBFOLDER="$RESULTS_DIR/$m"
    if [ ! -d $SUBFOLDER ]
    then
        mkdir "$SUBFOLDER"
    fi
done

DNA_DIR="$CONTROL/$DNA"

for n in $DNA_DIR/*/*.fasta
do 
    # extract the file info from the file name
    BASENAME=$(basename $n)
    
    # This splits the BASENAME variable with the "-" character and it extracts the 1st of the list that was split (this is defined by -f2)
    CENTROMERE=$(cut -d"-" -f2 <<< $BASENAME)
    CONTIG_PT1=$(cut -d"-" -f3 <<< $BASENAME)
    CONTIG_PT2=$(cut -d"-" -f4 <<< $BASENAME)
    CONTIG="$CONTIG_PT1-$CONTIG_PT2"

    # You need the $CURRENT_DIR in there because this command only outputs based on the relative path from sist_codes dir
	CHECK_FILE="control-C3-$CONTIG-.fasta-$CONTIG-18-COMBINED-.sist.csv"
	CHECK_DIR="control/sist-results/C3/T18"
        
	if [ -f "$CHECK_DIR/$CHECK_FILE" ]; then
		echo "true"
	else
		echo "$CHECK_FILE"
		bash sist_codes/tools/split_file_3.sh -f "$CURRENT_DIR/$n" -c $CONTIG -o "$CURRENT_DIR/$CONTROL/$RESULTS/$CENTROMERE"
	fi
done

