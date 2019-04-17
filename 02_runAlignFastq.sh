#! /bin/bash

## Allocate resources
#SBATCH --time=0-03:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-4

## job name
#SBATCH --job-name="align"

# read in the run specific settings
source ./varSettings.sh

# make destination directories

let i=$SLURM_ARRAY_TASK_ID-1

# organise fastq, do QC plots and align with minimap
./alignFastq.sh $expName ${barcodesOfInterest[$i]} $genomeFile
