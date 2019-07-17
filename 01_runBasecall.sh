#! /bin/bash

## Allocate resources
#SBATCH --time=3-00:00:00
#SBATCH --partition=all
#SBATCH --gres=gpu:1
##SBATCH --mem-per-cpu=16G
##SBATCH --cpus-per-task=4

## job name
#SBATCH --job-name="dSMFguppy"

./basecallGuppy.sh
