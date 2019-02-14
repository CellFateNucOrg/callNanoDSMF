#! /bin/bash
## script to basecall all fast5 in the folder one level up called fast5files and its recursive directories with albacore
## uses barcoding
## outputs 200,000 reads per fastq file (-q option)
## does not use batch subdirectories (-n 0 option)
## output will be in e.g: ../workspace/pass/barcode07/firstfastq.fastq


source activate albacore_env

read_fast5_basecaller.py --flowcell FLO-MIN106 --kit SQK-LSK108 --barcoding --output_format fastq -n 0 -q 200000 --input ../fast5files/ --recursive --save_path ../ --worker_threads 8
