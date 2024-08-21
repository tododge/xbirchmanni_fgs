#!/bin/bash
#SBATCH --job-name=slim_COAC_short
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000
#SBATCH -p schumer,owners

source activate /home/groups/schumer/lab_member_folders/tris/miniconda/envs/slim_env/
#conda install python=3.6.1
#conda install gsl

#pip install tskit --user
#pip install msprime --user
#pip install pyslim --user

ml gcc

for i in {1..25}; do /home/groups/schumer/shared_bin/SLiM/slim neutral_sim_IZAP_64725gen.slim; done
