genome=/scratch/groups/schumer/tris/falsegravid/gwas/xbir_WT/ref/xbir-COAC-16-VIII-22-M_v2023.1.fa
gwas_name=gwas_xbir-COAC-16-VIII-22-M_v2023.1.fa_FGS_329inds
header=/scratch/groups/schumer/lyle/nezzy_gwas/scripts/sbatch_header_schumer.txt

awk ' $2 > 1000000 ' $genome.fai | cut -f1 > xiphophorus_chroms.txt

chrom_list=xiphophorus_chroms.txt

for line in `cat $chrom_list`
do
    touch submit_mpileup.$gwas_name.$line.sh
    cat $header | while read n
    do
          echo $n >> submit_mpileup.$gwas_name.$line.sh
    done

  echo "perl run_mpileup_bcftools_GWAS.pl $gwas_name.txt $gwas_name.phenotypes.txt $genome /home/groups/schumer/shared_bin/ $line" >> submit_mpileup.$gwas_name.$line.sh
done

 
