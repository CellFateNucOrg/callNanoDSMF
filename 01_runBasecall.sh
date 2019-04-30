#! /bin/bash

## Allocate resources
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4

## job name
#SBATCH --job-name="guppy"

./basecallGuppy.sh
