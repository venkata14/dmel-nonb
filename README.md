# dmel-nonb

This repository has all the code for the paper "Enrichment of non-B-form DNA at *D. melanogaster* centromeres"

Note: put the reference genome in the `ref` folders. The genome can be found at [Heterochromatin-Enriched Assemblies Reveal the Sequence and Organization of the Drosophila melanogaster Y Chromosome](https://academic.oup.com/genetics/article/211/1/333/5931168).

Additional Note: You must fill in the `FULL_PATH` variable in multiple Shell scripts: `split_file.sh`, `extract_ipds_from_split_files.sh`, `extract_ipds.sh`,  

Generating Controls
---

In the `analysis` folder run the `gen_controls.ah` and `fill_controls.sh` to create the controls. You must first rin `init_setup.sh` to initialize the centromere folders before you can generate the controls. 

Preperation for Dyad Symmetry
---

```
sudo apt-get install emboss
```

*palindrome* is found in the EMBOSS suite of programs. Follow the commands in the `dyad-symmetry` folder for this part of the analysis.


Preparation for SIST
---

Install `sist-codes` from: https://academic.oup.com/bioinformatics/article/31/3/421/2365978. Follow the instructions for proper installation and place files into `analysis/sist-codes/` folder.

If on a 64-bit system, 32-bit support needs to be enabled.

```
sudo apt-get update
sudo dpkg --add-architecture i386
sudo apt-get update
sudo apt-get install libc6:i386 libncurses5:i386 libstdc++6:i386
```

If this does not work, try this.

```
sudo apt-get install multiarch-support
```

The analysis can be followed in the `analysis` folder

Set up the analysis with `init_setup.sh`. Run SIST on the centromere and controls with `run_cen_sist.sh` and `run_control_sist.sh`. ANalysis of the data can be found in the `analyze_enrichment.ipynb` file


Preperation for GQuad
---

Install R.

```
sudo apt-get install r-base
```

Install GQuad from https://cran.r-project.org/web/packages/gquad/index.html

Split the reference genome file into the respective contig fils and place them into the `ind-contigs` folder. Run `gquad_r_gquad.py`.

Place the respective files in the respective folders
```
mv aphased_* R_gquad_results/aphased/
mv gquad_* R_gquad_results/gquad/
mv hdna_* R_gquad_results/hdna/
mv slipped_* R_gquad_results/slipped/
mv str_* R_gquad_results/str/
mv tfo_* R_gquad_results/tfo/
mv zdna_* R_gquad_results/zdna/
```

Extraction of IPD Values for *D. melanogaster*
---

1. Initialize the environment and setting up `miniconda/smrttools 7.0.0`

```
sudo apt-get update
sudo apt-get install build-essential wget locales
sudo dpkg-reconfigure locales
sudo apt-get install rsync
sudo locale-gen en_US.UTF-8
sudo apt-get install pbbamtools
```

The tool used for the alignment and extraction was `smrttools 7.0.0`

The download can be found here: https://www.pacb.com/support/software-downloads/.
The download instructions can be found here: https://www.pacb.com/wp-content/uploads/SMRT_Link_Installation_v701.pdf

Download `miniconda` here: https://docs.conda.io/en/latest/miniconda.html


2. Install the datasets. Datasets from (https://bergmanlab.uga.edu/high-coverage-pacbio-shotgun-sequences-aligned-to-the-d-melanogaster-genome/)

Install the datasets
```
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro1_24NOV2013_398.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro2_25NOV2013_399.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro3_26NOV2013_400.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro4_28NOV2013_401.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro5_29NOV2013_402.tgz
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro6_1DEC2013_403.tgz
```

Getting md5sums

```
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro1_24NOV2013_398.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro2_25NOV2013_399.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro3_26NOV2013_400.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro4_28NOV2013_401.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro5_29NOV2013_402.md5
wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro6_1DEC2013_403.md5
```

Verify archives have been downloaded properly
```
md5sum -c Dro1_24NOV2013_398.md5 > Dro1_24NOV2013_398.md5.checkresult
md5sum -c Dro2_25NOV2013_399.md5 > Dro2_25NOV2013_399.md5.checkresult
md5sum -c Dro3_26NOV2013_400.md5 > Dro3_26NOV2013_400.md5.checkresult
md5sum -c Dro4_28NOV2013_401.md5 > Dro4_28NOV2013_401.md5.checkresult
md5sum -c Dro5_29NOV2013_402.md5 > Dro5_29NOV2013_402.md5.checkresult
md5sum -c Dro6_1DEC2013_403.md5  > Dro6_1DEC2013_403.md5.checkresult
```


3. Extract files with `tar -xvzf <file.tgz>`


4. Get the reference *D. melanogaster* genome from (HERE REF)


5. Convert the .h5 files to .bam files while extracting IPDs

```
bax2bam -o <output>.bax2ban.bam <input>.1.bax.h5 <input>.2.bax.h5 <input>.3.bax.h5 --subread —pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,PulseWidth, MergeQV,SubstitutionQV,SubstitutionTag —losslessframes
```

Do this for all the subread folders


6. Align the `.bam` files.

```
samtools faidx dm6.fa 
pbalign <input>.bax2bam.bam dm6.fa <output>.bax2bam.pbalign.bam
```


7. Merge the bam files and index it

```
samtools merge final.merge.pbalign.bam <input>*.bam 
pbindex final.merge.pbalign.bam
```

8. Extract IPD values

Install the required tools from: https://github.com/PacificBiosciences/kineticsTools/

This package does not have the IPD model for the chemistry of the dataset. The correct model can be found in the `smrttools 7.0.0` package.

Either this

```
ipdSummary final.merge.pbalign.bam --reference dm6.fa --useChemistry "P5-C3" --ipdModel /path-to-directory/P5-C3.npz.gz --csv kinetics.csv
```

or use `extract_ipds.sh` in the `gquad/ipds` folder. This requires that you put the `final.merge.pbalign.bam` in this folder as well.

However, we ran into RAM issues when using the methods above. As such, we split the final merged file into individual contigs and extracted IPDs from those split files.

In the `gquad/ipds`, `sort_bam.sh` sorts the bam file and outputs a SAM file. `make_split_file.sh` splits the SAM file into individual SAM files for each contig.  `extract_ipds_from_split_files.sh` extracts IPDs from the split SAM files and places them in the `gquad/ipds/ipd` folder. Analysis for this data is in `analyzeRData.ipynb`.


Preperation for G4Hunter
---


Place `G4Hunter.py` from the github found in Re-evaluation of G-quadruplex propensity with G4Hunter(https://academic.oup.com/nar/article/44/4/1746/1854457) in the `G4` folder. Run `script.sh`.


---

The list of controls used in this paper can be found in `analysis/list_of_controls_used`.