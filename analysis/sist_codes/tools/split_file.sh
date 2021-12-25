#!/usr/bin/bash

while getopts c:f:o: flag
do
    case "${flag}" in
        c) CONTIG=${OPTARG};;
        f) FASTA_FILE=${OPTARG};;
        # this output is relative to the sist_codes folder
        o) OUT_DIR=${OPTARG};;
    esac
done

module load samtools

# Getting the sequences of the contig
FULL_PATH = ""
SEQUENCE=$(bash /$FULL_PATH/analysis/sist_codes/tools/contig_seq_converter.sh -c $CONTIG -f $FASTA_FILE)


# declaring WINDOW and SLIDE as integers so you you can do math operations to it
declare -i WINDOW
declare -i SLIDE
declare -i INT_SEQLEN
declare -a TEMPS # declaring this an array. This is the 5 temps in Celsius are delared
declare -i NUMWORKERS
declare -i COUNTER

WINDOW=5000
SLIDE=2500
TEMPS=(18 22 25 30 35)
NUMWORKERS=22
COUNTER=0

# The command counts the number of sequences
SEQLEN=${#SEQUENCE}
INT_SEQLEN=$SEQLEN

FIRST_RANGES=$(seq -s " " 0 $WINDOW $SEQLEN)
SECOND_RANGES=$(seq -s " " $SLIDE $WINDOW $SEQLEN)

# Spliting the input fasta file into two windows
# Split into two for loops so that I can differentiate the two when combining them
# The WINDOW1 and WINDOW2 part
for range in $FIRST_RANGES 
do 
    # declaring variables FIRST and SECOND as integers for future computation
    declare -i FIRST
    declare -i SECOND

    FIRST=$range
    SECOND=FIRST+WINDOW
    FIRST=$FIRST+1 # samtool indexing
 

    if [ $SECOND -gt $INT_SEQLEN ]
    then
        SECOND=$INT_SEQLEN
    fi

    # moving sequences into the sist working directory. That is only how it works
    FILENAME="/$FULL_PATH/analysis/sist_codes/$CONTIG-split-WINDOW1-$FIRST-$SECOND-.fasta"
    samtools faidx $FASTA_FILE "$CONTIG:$FIRST-$SECOND" > "$FILENAME"
done

for range in $SECOND_RANGES
do 
    # declaring variables FIRST and SECOND as integers for future computation
    declare -i FIRST
    declare -i SECOND

    FIRST=$range
    SECOND=FIRST+WINDOW
    FIRST=$FIRST+1 # samtool indexing
 

    if [ $SECOND -gt $INT_SEQLEN ]
    then
        SECOND=$INT_SEQLEN
    fi

    # moving sequences into the sist working directory. That is only how it works
    FILENAME="/$FULL_PATH/analysis/sist_codes/$CONTIG-split-WINDOW2-$FIRST-$SECOND-.fasta"
    samtools faidx $FASTA_FILE "$CONTIG:$FIRST-$SECOND" > "$FILENAME"
done



# SIST programs only seems to work when the files are in the same directory as master.pl so Im going to that directory
cd "/$FULL_PATH/analysis/sist_codes" 
TMP_DIR="tmp_sist_results"

# Make this temporary directory if it doesn't exist
if [ ! -d "$TMP_DIR" ] 
then
    mkdir $TMP_DIR
fi

# Looping through the 5 temperatures
for m in ${TEMPS[@]}
do
    # Make sure the directory where to results goes exists
    if [ ! -d "$TMP_DIR/T$m" ] 
    then
        mkdir $TMP_DIR/T$m
    fi

    # Looping though all the fasta files
    for n in *fasta
    do
        # This if statement limits the amount of cores used defined by NUMWORKERS
        if [ $COUNTER -ge $NUMWORKERS ]
        then
            wait
            COUNTER=0
        fi

        OUT_FILE="$TMP_DIR/T$m/$n.T$m.sist"
        # Bash does not do float arithmatic so I did a work around to convert Celsius to Kelvin
        # Put in backend. Parallelizes the opperation 
        perl master.pl -f "$n" -a A -t $(($m + 273)).15 -o $OUT_FILE &

        COUNTER=$((COUNTER+1)) # Increment the counter by 1
    done
    wait # Waits for the backend commands to finish
    rm one_line* # remove these files after each iteration
done
# remove fasta files when done
rm *fasta

# This next part is to merge the files for the 5 temperatures
for m in ${TEMPS[@]}
do
    INP_DIR="$TMP_DIR/T$m"
    BASENAME=$(basename $FASTA_FILE)
    OUT_COMB_DIR="$OUT_DIR/T$m"
    if [ ! -d "$OUT_COMB_DIR" ]; then
        mkdir "$OUT_COMB_DIR"
    fi
    OUT_COMB_FILE="$OUT_COMB_DIR/$BASENAME-$CONTIG-$m-COMBINED-.sist.csv"

    echo "Outputing to: $OUT_COMB_FILE"
    # merge with this python script
    python tools/converter.py -i "$INP_DIR" -o "$OUT_COMB_FILE"
done
rm -r "$TMP_DIR"
cd -