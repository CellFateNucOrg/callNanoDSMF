#! /bin/bash
## script to basecall all fast5 in the folder called fast5files and its recursive directories with Guppy
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## uses barcoding
## output will be in e.g: ./fastqFiles/pass/firstfastq.fastq
## next the files are sorted by barcode and output to e.g ./bcFastq/pass/barcodeXX.fastq.gz
## pycoQC is run and results are in ./qc/pycoQC/pycoQC.html

source ./varSettings.sh

##################
# call barcodes
#################

source ${CONDA_ACTIVATE} MC-HiC-AA

deepbinner classify --native ${dataDir}/fast5Files > ${workDir}/classifications


##################
# bin by barcode
##################

mkdir -p ${workDir}/dbbcFastq/pass
mkdir -p ${workDir}/dbbcFastq/fail

cat ${workDir}/fastqFiles/pass/* > ${workDir}/fastqFiles/pass/passed.fq
cat ${workDir}/fastqFiles/fail/* > ${workDir}/fastqFiles/fail/failed.fq



deepbinner bin --classes ${workDir}/classifications --reads ${workDir}/fastqFiles/pass/passed.fq --out_dir ${workDir}/dbbcFastq/pass
deepbinner bin --classes ${workDir}/classifications --reads ${workDir}/fastqFiles/fail/failed.fq --out_dir ${workDir}/dbbcFastq/fail

rm ${workDir}/fastqFiles/pass/passed.fq
rm ${workDir}/fastqFiles/fail/failed.fq



##################
# run pycoQC
#################

#source ${HOME}/.bashrc
#source ${CONDA_ACTIVATE} pycoQC
#conda activate pycoQC

# create qc output directory
qcDir=${workDir}/qc
mkdir -p ${qcDir}/pycoQC

#https://github.com/a-slide/pycoQC

pycoQC --summary_file ${workDir}/fastqFiles/sequencing_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}.html

conda deactivate
