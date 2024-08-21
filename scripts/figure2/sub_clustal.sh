#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --time=6:00:00
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000
#SBATCH -p schumer

#/home/groups/schumer/shared_bin/clustalo-1.2.4-Ubuntu-x86_64 -i xbirFGS_xbirWT_xmal_xvar.fa -o clustal_xbirFGS_xbirWT_xmal_xvar.fa --outfmt=fasta -v --force --threads=4

/home/groups/schumer/shared_bin/clustalo-1.2.4-Ubuntu-x86_64 -i combined_longread_FGS_haplotypes_formolly_linearized_renamed.fasta -o clustal_combined_longread_FGS_haplotypes_formolly_linearized_renamed.fasta --outfmt=fasta -v --force --threads=16

snp-sites -c clustal_combined_longread_FGS_haplotypes_formolly_linearized_renamed.fasta -o clustal_combined_longread_FGS_haplotypes_formolly_linearized_renamed_snpsites.fasta

/home/groups/schumer/shared_bin/raxmlHPC-PTHREADS -f a -m GTRGAMMA -p 12345 -x 12345 -s clustal_combined_longread_FGS_haplotypes_formolly_linearized_renamed_snpsites.fasta -N 100 -n clustal_combined_longread_FGS_haplotypes_formolly_linearized_renamed_snpsites.fasta

