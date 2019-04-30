#! /bin/bash
## script to basecall all fast5 in the folder one level up called fast5files and its recursive directories with Guppy
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## uses barcoding
## outputs 200,000 reads per fastq file (-q option)
## output will be in e.g: ../fastqFiles/pass/firstfastq.fastq
## next the files are sorted by barcode and output to e.g ../bcFastq/pass/barcodeXX/firstfastq.fastq
## pycoQC is run and results are in ../qc/pycoQC/pycoQC.html

source ./varSettings.sh
fullPath=`readlink -f ${relPath}`
fast5path=$fullPath   #`readlink -f ../`

##################
# call barcodes
#################

mkdir -p ${fast5path}/bcFast5
deepbinner classify --native ${fast5path}/fast5Files > classifications 
#deepbinner realtime --in_dir ${fast5path}/fast5Files --out_dir ${fast5path}/bcFast5 --native


###################
## basecall
##################
#
#
#${GUPPYDIR}/guppy_basecaller --input_path ${fast5path}/bcFast5 --save_path ${fullPath}/bcFastq --flowcell FLO-MIN106 --kit SQK-LSK109 --records_per_fastq 200000 --recursive  --cpu_threads_per_caller 8 --qscore_filtering --min_qscore 3 
##--num_callers
#
#
#
###################
## run pycoQC
##################
#
#source activate pycoQC
#
## create qc output directory
#qcDir=${fullPath}/qc
#mkdir -p $qcDir/pycoQC
#
##https://github.com/a-slide/pycoQC
#
#pycoQC -f ${fullPath}/bcFastq/sequencing_summary.txt -b ${fullPath}/bcFastq/pass/barcoding_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}_pass.html --title "nanoDSMF "${expName}" passed reads" --min_pass_qual 3
#
#pycoQC -f ${fullPath}/fastqFiles/sequencing_summary.txt -b ${fullPath}/bcFastq/fail/barcoding_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}_fail.html --title "nanoDSMF "${expName}" failed reads" --min_pass_qual 3 
#
#source deactivate