#!/bin/bash

# create centromere results
CEN_DIR="cen"
DNA="dna"
RESULTS="sist-results"
CURRENT_DIR="$PWD"

declare -a CENS
CENS=( "C2" "C3" "C4" "CX" "CY")

# creating the results dirs
RESULTS_DIR="$CEN_DIR/$RESULTS"
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

DNA_DIR="$CEN_DIR/$DNA"

for n in $DNA_DIR/*/*.fasta
do 
    # extract the file info from the file name
    BASENAME=$(basename $n)

    echo "-------------------------------"
    echo "Running on this file: $BASENAME"
    echo "-------------------------------"
    
    # This splits the BASENAME variable with the "-" character and it extracts the 1st of the list that was split (this is defined by -f1)
    CENTROMERE=$(cut -d"-" -f1 <<< $BASENAME)
    CONTIG=$(cut -d"-" -f2 <<< $BASENAME)

    bash sist_codes/tools/split_file.sh -f "$CURRENT_DIR/$n" -c $CONTIG -o "$CURRENT_DIR/$CEN_DIR/$RESULTS/$CENTROMERE"
done

