#!/bin/bash
#SBATCH --job-name=psmc
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH -p schumer

for FASTA in input/*psmcfa; do \
PREFIX=${FASTA%%.psmcfa}; \
./psmc/utils/splitfa $FASTA > $PREFIX.split.psmcfa; done
