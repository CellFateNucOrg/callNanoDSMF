#! /bin/bash
## script to basecall all fast5 in the folder called fast5files and its recursive directories with Guppy
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## uses barcoding
## output will be in e.g: ./fastqFiles/pass/firstfastq.fastq
## next the files are sorted by barcode and output to e.g ./bcFastq/pass/barcodeXX.fastq.gz
## pycoQC is run and results are in ./qc/pycoQC/pycoQC.html

source ./varSettings.sh


##################
# basecall
#################

guppy_basecaller --input_path ${dataDir}/fast5Files --save_path ${workDir}/fastqFiles --records_per_fastq 200000 --recursive  --qscore_filtering --min_qscore 3  --device auto --config ${guppyConfigFile} 
