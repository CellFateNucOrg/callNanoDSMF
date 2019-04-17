#! /bin/bash
## script to basecall all fast5 in the folder one level up called fast5files and its recursive directories with Guppy
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## uses barcoding
## outputs 200,000 reads per fastq file (-q option)
## output will be in e.g: ../fastq/pass/firstfastq.fastq

mkdir -p fastqFiles

${GUPPYDIR}/guppy_basecaller --input_path ../fast5Files --save_path ../fastqFiles --flowcell FLO-MIN106 --kit SQK-LSK109 --records_per_fastq 200000 --recursive  --cpu_threads_per_caller 8 --qscore_filtering --min_qscore 3 
#--num_callers

mkdir -p bcFastq

${GUPPYDIR}/guppy_barcoder -i ../fastqFiles -s ../bcFastq --barcode_kits EXP-NBD104 --worker_threads 8 --recursive -q 200000

#rm ../fastqFiles


##################
# run pycoQC
#################

source activate pycoQC

# create qc output directory
qcDir=../qc
mkdir -p $qcDir/pycoQC

#https://github.com/a-slide/pycoQC

pycoQC -f ../fastqFiles/sequencing_summary.txt -b ../bcFastq/barcoding_summary.txt -o $qcDir/pycoQC/pycoQC.html --title "nanoDSMF1" --min_pass_qual 3

source deactivate
