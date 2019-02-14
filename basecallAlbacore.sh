#! /bin/bash
## script to basecall all fast5 in the folder one level up called fast5files and its recursive directories with albacore
## uses barcoding
## outputs 200,000 reads per fastq file (-q option)
## does not use batch subdirectories (-n 0 option)
## output will be in e.g: ../workspace/pass/barcode07/firstfastq.fastq


source activate albacore_env

read_fast5_basecaller.py --flowcell FLO-MIN106 --kit SQK-LSK108 --barcoding --output_format fastq -n 0 -q 200000 --input ../fast5files/ --recursive --save_path ../ --worker_threads 8

##################
# run MinIONQC
#################

module load vital-it
module load R/3.5.1

# create qc output directory
qcDir=../fastqQC
mkdir -p $qcDir

echo "doing MinIONQC..."
# run Minion_qc
# source: https://github.com/roblanf/minion_qc

source ./varSettings.sh
mkdir -p ${qcDir}/MinIONQC

Rscript ${MINIONQC} -i ../sequencing_summary.txt -o $qcDir

mv ../*.png ${qcDir}/MinIONQC/
