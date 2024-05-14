#!/bin/bash
#SBATCH --job-name=cov
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH -p schumer,owners

ml biology bedtools samtools

for file in `ls raw/*.sorted.dedup.realigned.bam`; do \
basename $file; \
samtools depth $file -a | awk '{sum+=$3} END { print "Average = ",sum/NR}'; \
done
