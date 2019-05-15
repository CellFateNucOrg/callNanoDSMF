#! /bin/bash
## script to basecall all fast5 in the folder called fast5files and its recursive directories with Guppy
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## uses barcoding
## output will be in e.g: ./fastqFiles/pass/firstfastq.fastq
## next the files are sorted by barcode and output to e.g ./bcFastq/pass/barcodeXX.fastq.gz
## pycoQC is run and results are in ./qc/pycoQC/pycoQC.html

source ./varSettings.sh
#fullPath=`readlink -f ${relPath}`
fullPath=$PWD

##################
# call barcodes
#################

mkdir -p ${fullPath}/bcFast5
deepbinner classify --native ${fullPath}/fast5Files > classifications 


##################
# basecall
#################


${GUPPYDIR}/guppy_basecaller --input_path ${fullPath}/fast5Files --save_path ${fullPath}/fastqFiles --flowcell FLO-MIN106 --kit SQK-LSK109 --records_per_fastq 200000 --recursive  --cpu_threads_per_caller 8 --qscore_filtering --min_qscore 3 


##################
# bin by barcode
##################

mkdir -p ${fullPath}/bcFastq/pass
mkdir -p ${fullPath}/bcFastq/fail

cat ${fullPath}/fastqFiles/pass/* > ${fullPath}/fastqFiles/pass/passed.fq
deepbinner bin --classes ${fullPath}/classifications --reads ${fullPath}/fastqFiles/pass/passed.fq --out_dir ${fullPath}/bcFastq/pass

cat ${fullPath}/fastqFiles/fail/* > ${fullPath}/fastqFiles/fail/failed.fq
deepbinner bin --classes ${fullPath}/classifications --reads ${fullPath}/fastqFiles/fail/failed.fq --out_dir ${fullPath}/bcFastq/fail

rm ${fullPath}/fastqFiles/pass/passed.fq
rm ${fullPath}/fastqFiles/fail/failed.fq

##################
# run pycoQC
#################

source activate pycoQC

# create qc output directory
qcDir=${fullPath}/qc
mkdir -p $qcDir/pycoQC

#https://github.com/a-slide/pycoQC

pycoQC -f ${fullPath}/fastqFiles/sequencing_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}.html --min_pass_qual 3

source deactivate

