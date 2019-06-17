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


#guppy_basecaller --input_path ${dataDir}/fast5Files --save_path ${workDir}/fastqFiles --flowcell FLO-MIN106 --kit SQK-LSK109 --records_per_fastq 200000 --recursive  --qscore_filtering --min_qscore 3  --device auto  


##################
# call barcodes
#################

deepbinner classify --native ${dataDir}/fast5Files > ${workDir}/classifications


##################
# bin by barcode
##################

mkdir -p ${workDir}/dbbcFastq/pass
mkdir -p ${workDir}/dbbcFastq/fail

cat ${workDir}/fastqFiles/pass/* > ${workDir}/fastqFiles/pass/passed.fq
cat ${workDir}/fastqFiles/fail/* > ${workDir}/fastqFiles/fail/failed.fq

deepbinner classify --native ${dataDir}/fast5Files > ${workDir}/classifications


deepbinner bin --classes ${workDir}/classifications --reads ${workDir}/fastqFiles/pass/passed.fq --out_dir ${workDir}/dbbcFastq/pass
deepbinner bin --classes ${workDir}/classifications --reads ${workDir}/fastqFiles/fail/failed.fq --out_dir ${workDir}/dbbcFastq/fail

#rm ${workDir}/fastqFiles/pass/passed.fq
#rm ${workDir}/fastqFiles/fail/failed.fq



###################
## call barcodes
##################

#mkdir -p ${workDir}/bcFastq/pass
#guppy_barcoder -i ${workDir}/fastqFiles/pass -s ${workDir}/bcFastq/pass --barcode_kits EXP-NBD104 --device auto  --recursive


#mkdir -p ${workDir}/bcFastq/fail
#guppy_barcoder -i ${workDir}/fastqFiles/fail -s ${workDir}/bcFastq/fail --barcode_kits EXP-NBD104 --device auto  --recursive

##################
# run pycoQC
#################

#source ${HOME}/.bashrc
#source ${CONDA_ACTIVATE}
#conda activate pycoQC

# create qc output directory
qcDir=${workDir}/qc
mkdir -p ${qcDir}/pycoQC

#https://github.com/a-slide/pycoQC

#pycoQC --summary_file ${workDir}/fastqFiles/sequencing_summary.txt --barcode_file ${workDir}/bcFastq/pass/barcoding_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}_pass.html --min_pass_qual 3


#pycoQC --summary_file ${workDir}/fastqFiles/sequencing_summary.txt --barcode_file ${workDir}/bcFastq/fail/barcoding_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}_fail.html

#conda deactivate

#pycoQC --summary_file ${workDir}/fastqFiles/sequencing_summary.txt --barcode_file ${workDir}/bcFastq/pass/barcoding_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}_pass.html --min_pass_qual 3


pycoQC --summary_file ${workDir}/fastqFiles/sequencing_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}.html
