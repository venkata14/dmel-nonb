#!/bin/bash
#SBATCH --job-name=run_palindrome
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=6G
#SBATCH --mail-user=venkata.patchigolla@uconn.edu
#SBATCH -o output_run_palindrome_%j.out
#SBATCH -e output_run_palindrome_%j.err

CEN="cen"
CONTROL="control"
DNA="dna"
RESULTS="dyad"
declare -a REGIONS
REGIONS=( "cen" "control" )
declare -a CENS
CENS=( "C2" "C3" "C4" "CX" "CY" )

# make the dyad folder if not there
if [ ! -d "$CEN/$RESULTS" ]; then
	mkdir "$CEN/$RESULTS"
	mkdir "$CONTROL/$RESULTS"
fi

for reg in ${REGIONS[@]}; do
	for cen in ${CENS[@]}; do
		FOLDER="$reg/$DNA/$cen/*.fasta"
		echo "$FOLDER"
		
		# make output folders
		if [ ! -d "$reg/$RESULTS/$cen" ]; then
			mkdir "$reg/$RESULTS/$cen"
		fi

		for file in $FOLDER; do
			CURRENT_FILE=$(basename $file)
			MOD_file=${file/:/-}
			mv $file $MOD_file
			OUTPUT="$reg/$RESULTS/$cen/$CURRENT_FILE-.dyad.txt"
			palindrome -minpallen 5 -maxpallen 100 -gaplimit 20 -nummismatches 0 -overlap -sequence "./$MOD_file" -outfile "./$OUTPUT"
			mv $MOD_file $file
		done
	done
done


