#! /bin/bash

# name of experiment
expName="20171027"

# barcodes used in this experiment
barcodesOfInterest=( barcode05 barcode06 barcode07 barcode08 )

# location of reference genome file
genomeFile=/data/projects/p025/Jenny/genomeVer/PCR_wPM28_32/PCR_wPM28_32.fasta

# location of MinION program
MINIONQC=/home/pmeister/software/MinIONQC.R

########## don't edit below this line ###############
# note: the .bashrc should point to the location of nanopolish with the variable
# NANOPOLISH_DIR

export NUMBC=${#barcodesOfInterest[@]}
