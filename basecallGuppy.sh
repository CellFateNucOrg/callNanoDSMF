#! /bin/bash
## script to basecall all fast5 in the folder one level up called fast5files and its recursive directories with Guppy
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## uses barcoding
## outputs 200,000 reads per fastq file (-q option)
## output will be in e.g: ../fastqFiles/pass/firstfastq.fastq
## next the files are sorted by barcode and output to e.g ../bcFastq/pass/barcodeXX/firstfastq.fastq
## pycoQC is run and results are in ../qc/pycoQC/pycoQC.html

source ./varSettings.sh
fullPath=`readlink -f ../`

${GUPPYDIR}/guppy_basecaller --input_path ${fullPath}/fast5Files --save_path ${fullPath}/fastqFiles --flowcell FLO-MIN106 --kit SQK-LSK109 --records_per_fastq 200000 --recursive  --cpu_threads_per_caller 8 --qscore_filtering --min_qscore 3 
#--num_callers

#basecall passed files
mkdir -p ../bcFastq/pass
${GUPPYDIR}/guppy_barcoder -i ${fullPath}/fastqFiles/pass -s ${fullPath}/bcFastq/pass --barcode_kits EXP-NBD104 --worker_threads 8 --recursive -q 200000

#basecall failed files
mkdir -p ../bcFastq/fail
${GUPPYDIR}/guppy_barcoder -i ${fullPath}/fastqFiles/fail -s ${fullPath}/bcFastq/fail --barcode_kits EXP-NBD104 --worker_threads 8 --recursive -q 200000

#rm ../fastqFiles


##################
# run pycoQC
#################

source activate pycoQC

# create qc output directory
qcDir=../qc
mkdir -p $qcDir/pycoQC

#https://github.com/a-slide/pycoQC

pycoQC -f ${fullPath}/fastqFiles/sequencing_summary.txt -b ${fullPath}/bcFastq/pass/barcoding_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}_pass.html --title "nanoDSMF "${expName}" passed reads" --min_pass_qual 3

pycoQC -f ${fullPath}/fastqFiles/sequencing_summary.txt -b ${fullPath}/bcFastq/fail/barcoding_summary.txt -o ${qcDir}/pycoQC/pycoQC_${expName}_fail.html --title "nanoDSMF "${expName}" failed reads" 

source deactivate
