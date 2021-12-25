#!/usr/bin/bash

module load samtools

declare -a CENS
CENS=("C2" "C3" "C4" "CX" "CY")
REF='ref'
CONTROL="control"
DNA="$REF/dm6.fasta"

# Making the 
if [ ! -d "$CONTROL/dna" ]; then
    mkdir "$CONTROL/dna"
fi

for cen in ${CENS[@]}; do
    FILENAME="$REF/$cen-controls-.bed"
    DIR="$CONTROL/dna/$cen"

    if [ ! -d "$DIR" ]; then
        mkdir "$DIR"
    fi
    
    for region in $(less $FILENAME); do
        samtools faidx "$DNA" "$region" > "$DIR/control-$cen-$region-.fasta"
    done
done