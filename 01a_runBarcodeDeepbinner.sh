#! /bin/bash

## Allocate resources
#SBATCH --gres=gpu:1
#SBATCH --time=2-12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=all
#SBATCH --gres=gpu:1

## job name
#SBATCH --job-name="barcodeBin"


# Call barcodes with deepbinner and runs pycoQC
./barcodeDeepbinner.sh 

