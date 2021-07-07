#! /bin/bash

# name of experiment
expName="20190411"

# barcodes used in this experiment
barcodesOfInterest=( barcode01 barcode02 barcode03 barcode04 barcode05 unclassified )
# barcodes with non C.elegans DNA spike in
bcWithSpikeIn=( barcode01 barcode05 unclassified )


# location of reference genome file
#genomeFile=/data/projects/p025/Jenny/genomeVer/PCR_wPM28_32/PCR_wPM28_32.fasta
genomeFile=/mnt/imaging.data/jsemple/genomeVer/ws265/c_elegans.PRJNA13758.WS265.genomic.fa

# location of reference genomes for spiked-in DNA
#lambdaFile="/home/ubelix/izb/semple/genomeVer/lambda/lambdaPhage.fasta"
#phiXfile="home/ubelix/izb/semple/genomeVer/phiX/phiX.fasta"
lambdaFile=/mnt/imaging.data/jsemple/genomeVer/lambda/lambdaPhage.fasta
phiXfile=/mnt/imaging.data/jsemple/genomeVer/phiX/phiX.fasta

dataDir=/mnt/imaging.data/jsemple/20190411_dSMFv021-025np_N2gw
workDir=/home/jsemple/20190411_dSMFv021-025np_N2gw

# location of MinION program
#MINIONQC=/home/pmeister/software/MinIONQC.R

########## don't edit below this line ###############
# note: the .bashrc should point to the location of nanopolish with the variable
# NANOPOLISH_DIR

export NUMBC=${#barcodesOfInterest[@]}
