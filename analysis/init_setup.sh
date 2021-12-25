#!/usr/bin/bash

# This script is to initialize the centromere and control folders
# This will put in the centromere sequences into the right centromere folder
# this requires the dm6.fasta folder in the ref/ directory

# First Part: Make all the folderse/subfolders

# Making array of the main folders and the subfolders
declare -a MAINDIR
declare -a CENS

MAINDIR=("cen" "control")
DNADIR="dna"
CENS=("C2" "C3" "C4" "CX" "CY")

# Make all the folders
for n in ${MAINDIR[@]}
do
    # Making the main cen or control folder
    if [ ! -d $n ]
    then
        mkdir "$n"
    fi

    # Making the dna sub folder where all the dan sequences will go
    DNA_DIR="$n/$DNADIR"
    if [ ! -d $DNA_DIR ]
    then
        mkdir "$DNA_DIR"
    fi

    # Making the subfolders for individual centromeres
    for m in ${CENS[@]}
    do 
        SUBFOLDER="$DNA_DIR/$m"
        if [ ! -d $SUBFOLDER ]
        then
            mkdir "$SUBFOLDER"
        fi
    done
done

# Second Part: Put in centromere contig sequences into the respective subfolders

module load samtools

# Making hash table of all contigs
declare -A CONTIGS

CONTIGS=( ["C2"]="tig00057289" ["C3"]="3R_5" ["C4"]="Contig119" ["CX"]="Contig79" ["CY"]="Y_Contig26" )

REF="ref/dm6.fasta"
CEN_DIR="cen/dna"

# need to first index the ref genome
samtools faidx "$REF"

# "${!CONTIGS[@]}" notation to get all keys
# "${CONTIGS[@]}" notation to get all values
for n in "${!CONTIGS[@]}"
do
    CONTIG="${CONTIGS[$n]}"
    FILE="$CEN_DIR/$n/$n-$CONTIG-.fasta"
    samtools faidx "$REF" "$CONTIG" > "$FILE"
done
