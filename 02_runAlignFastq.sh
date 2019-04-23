#! /bin/bash

## Allocate resources
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-6
## you should submit as many jobs as there are barcodes in barcodesOfInterest
## (don't forget to include unclassfied in barcodesOfInterst in the varSettings.sh file)

## job name
#SBATCH --job-name="npAlign"

# read in the run specific settings
source ./varSettings.sh

# make destination directories

let i=$SLURM_ARRAY_TASK_ID-1

# organise fastq, do QC plots and align with minimap
./alignFastq.sh $expName ${barcodesOfInterest[$i]} $genomeFile


# if barcode is one with a spikein, align also to other genomes
for bc in ${bcWithSpikeIn[@]}
do
    if [ "${barcodesOfInterest[$i]}" == "$bc" ]; then
        ./alignFastq_spikeIn.sh $expName ${barcodesOfInterest[$i]} ${lambdaFile}
        ./alignFastq_spikeIn.sh $expName ${barcodesOfInterest[$i]} ${phiXfile}
    	exit
    fi
done
