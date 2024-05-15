#!/bin/bash

#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-user=nicolas.klutchnikoff@univ-rennes2.fr
#SBATCH --mail-type=ALL
#SBATCH --output=message/POLY-result.txt
#SBATCH --error=message/POLY-error.txt

module load StdEnv/2023 r/4.3.1
Rscript SIMULATE_poly.R
