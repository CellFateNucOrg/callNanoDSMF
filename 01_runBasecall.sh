#! /bin/bash

## Allocate resources
#SBATCH --time=0-10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

## job name
#SBATCH --job-name="albacore"

./basecallAlbacore.sh
