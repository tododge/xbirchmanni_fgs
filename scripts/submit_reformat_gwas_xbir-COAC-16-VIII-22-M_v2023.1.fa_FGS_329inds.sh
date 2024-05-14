#!/bin/bash
#SBATCH --job-name=reformat
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH -p schumer,owners

gwas_name=gwas_xbir-COAC-16-VIII-22-M_v2023.1.fa_FGS_329inds

for f in `ls $gwas_name.txt*vcf`; do (perl print_alleles_depth_freq_chi_per_site_GWAS_v3.pl $f 0) > $f.summary; done

for f in $gwas_name.txtindivs*summary; do ( tail -n +2 $f) > $f.noheader; done
first_vcf_in_list=`ls $gwas_name*summary | head -n 1`

rm $first_vcf_in_list.noheader
cat $first_vcf_in_list $gwas_name.txtindivs*summary.noheader > $gwas_name.txtindivs.allindiv.allscaff.vcf.summary

awk -F "\t" '{ if(($9 < 0.01 || $9 == "LR1")) { print } }' $gwas_name.txtindivs.allindiv.allscaff.vcf.summary > $gwas_name.txtindivs.allindiv.allscaff.vcf.summary_LR1lessthan0.01

#sed -i 's/chr-//g' $gwas_name.txtindivs.allindiv.allscaff.vcf.summary_LR1lessthan0.01

sed -i 's/ /\t/g' $gwas_name.txtindivs.allindiv.allscaff.vcf.summary_LR1lessthan0.01
sed -i 's/\t\t/\t/g' $gwas_name.txtindivs.allindiv.allscaff.vcf.summary_LR1lessthan0.01
sed -i 's/\t\t/\t/g' $gwas_name.txtindivs.allindiv.allscaff.vcf.summary_LR1lessthan0.01
