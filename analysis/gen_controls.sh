#!/usr/bin/bash

declare -i NUMWORKERS
declare -i COUNTER

NUMWORKERS=20
COUNTER=0
GENOME="ref/dm6.fasta"
ALL_CONTIGS="/ref/all-contigs-.txt"
NO_CEN_GENOME='ref/dm6-no_cen-.fasta'
REF="ref"

module load samtools
module load bedtools

# Note: You do not need to make a whole ref genome without the centromeric contigs to do this.
# TO make the $NO_CEN_GENOME.fai you can do: less dm6.fasta.fai | sed "/tig00057289\|3R_5\|Contig119\|Contig79\|Y_Contig26/d" > some-name

# Get all contigs
CONTIGS_NO_CEN=$(grep ">" $GENOME | sed -e "s/>//g" | sed "/tig00057289\|3R_5\|Contig119\|Contig79\|Y_Contig26/d" | tr "\n" " ")

samtools faidx $GENOME
samtools faidx $GENOME $CONTIGS_NO_CEN > $NO_CEN_GENOME
samtools faidx $NO_CEN_GENOME


# Generate random sets of genome intervals for each centromere

C2_SEQ=$(bash ./sist_codes/tools/contig_seq_converter.sh -f "./cen/dna/C2/*.fasta" -c "tig00057289" )
C2_SIZE=${#C2_SEQ}

C3_SEQ=$(bash ./sist_codes/tools/contig_seq_converter.sh -f "./cen/dna/C3/*.fasta" -c "3R_5" )
C3_SIZE=${#C3_SEQ}

C4_SEQ=$(bash ./sist_codes/tools/contig_seq_converter.sh -f "./cen/dna/C4/*.fasta" -c "Contig119" )
C4_SIZE=${#C4_SEQ}

CX_SEQ=$(bash ./sist_codes/tools/contig_seq_converter.sh -f "./cen/dna/CX/*.fasta" -c "Contig79" )
CX_SIZE=${#CX_SEQ}

CY_SEQ=$(bash ./sist_codes/tools/contig_seq_converter.sh -f "./cen/dna/CY/*.fasta" -c "Y_Contig26" )
CY_SIZE=${#CY_SEQ}


declare -A CEN_SIZES
CEN_SIZES=( ['C2']=$C2_SIZE ["C3"]=$C3_SIZE ["C4"]=$C4_SIZE ["CX"]=$CX_SIZE ["CY"]=$CY_SIZE )

for n in "${!CEN_SIZES[@]}"
do
    # This command needs a all the contigs with the chromosome sizes. So sending in the .fai works
    # sed command removes the negative strands
    # awk command puts the lines in a format that samtools faidx can use
    # Had to subtract 1 due to weird indexing by samtools
    bedtools random -g $NO_CEN_GENOME.fai -l $((${CEN_SIZES[$n]}-1)) -n 1000000 -seed 42 | sed -e "/-/d" |  awk '{print $1 ":" $2 "-" $3}' > "$REF/$n-potential_controls-.bed"
done

echo "Generated Random Intervals"

# Make this better. Make Function 
C2_GC=$(awk 'BEGIN{RS=">";FS="\n"}NR>1{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq="";for (i=2;i<=NF;i++) seq=seq""$i; k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; print (sumC+sumG)/k*100}' cen/dna/C2/*.fasta)
C3_GC=$(awk 'BEGIN{RS=">";FS="\n"}NR>1{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq="";for (i=2;i<=NF;i++) seq=seq""$i; k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; print (sumC+sumG)/k*100}' cen/dna/C3/*.fasta)
C4_GC=$(awk 'BEGIN{RS=">";FS="\n"}NR>1{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq="";for (i=2;i<=NF;i++) seq=seq""$i; k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; print (sumC+sumG)/k*100}' cen/dna/C4/*.fasta)
CX_GC=$(awk 'BEGIN{RS=">";FS="\n"}NR>1{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq="";for (i=2;i<=NF;i++) seq=seq""$i; k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; print (sumC+sumG)/k*100}' cen/dna/CX/*.fasta)
CY_GC=$(awk 'BEGIN{RS=">";FS="\n"}NR>1{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq="";for (i=2;i<=NF;i++) seq=seq""$i; k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; print (sumC+sumG)/k*100}' cen/dna/CY/*.fasta)

declare -A CEN_GC
CEN_GC=( ['C2']=$C2_GC ["C3"]=$C3_GC ["C4"]=$C4_GC ["CX"]=$CX_GC ["CY"]=$CY_GC )

for m in "${!CEN_GC[@]}"
do
    FILENAME="$REF/$m-potential_controls-.bed"
    OUT_FILENAME="$FILENAME.filtered"

    # clear out the file for debugging purposes
    > $OUT_FILENAME

    for n in $(less $FILENAME)
    do
	# This is for DEBUGGING purposes. This is to see what contig it is running at the moment
	
	# echo "$FILENAME"
	# echo "$n"

        # This if statement limits the amount of cores used defined by NUMWORKERS
        if [ $COUNTER -ge $NUMWORKERS ]
        then
            wait
            COUNTER=0
        fi

        # awk -v tag to put variables in the statement
        samtools faidx 'ref/dm6.fasta' $n | awk -vGC="${CEN_GC[$m]}" 'BEGIN{RS=">";
                                            FS="\n"}NR>1{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq="";
                                            for (i=2;i<=NF;i++) seq=seq""$i; 
                                            k=length(seq); 
                                            for (i=1;i<=k;i++) {
                                                if (substr(seq,i,1)=="T") sumT+=1; 
                                                else if (substr(seq,i,1)=="A") sumA+=1; 
                                                else if (substr(seq,i,1)=="G") sumG+=1; 
                                                else if (substr(seq,i,1)=="C") sumC+=1; 
                                                else if (substr(seq,i,1)=="N") sumN+=1}; 
                                            if ((sumC+sumG)/k*100<(GC+10.0) && (sumC+sumG)/k*100>(GC-10.0)) print $1"\t"k"\t"GC"\t"(sumC+sumG)/k*100}' >> $OUT_FILENAME &
        

        COUNTER=$((COUNTER+1)) # Increment the counter by 1
    done
    wait
done
