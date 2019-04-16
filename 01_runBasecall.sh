#! /bin/bash

## Allocate resources
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

## job name
#SBATCH --job-name="guppy"

./basecallGuppy.sh
