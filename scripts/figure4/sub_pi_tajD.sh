#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH -p schumer

#conda create -n "pixy" python=3.8

pixy --stats pi \
--vcf COAC_2016-2018_23_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz \
--populations pops \
--output_prefix COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf_pixy_pi_8794 \
--window_size 8794

ml biology vcftools

zcat COAC_2016-2018_23_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02_28345429_28354223.g.vcf.gz | awk 'BEGIN{OFS="\t"} /^#/{print;next} {$2=$2-28345429;print}' |  bgzip > COAC_2016-2018_23_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02_1_8794.g.vcf.gz 

tabix -p vcf COAC_2016-2018_23_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02_1_8794.g.vcf.gz 
vcftools --gzvcf COAC_2016-2018_23_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02_1_8794.g.vcf.gz --TajimaD 8794 --out tajimasd.inversion
vcftools --gzvcf COAC_2016-2018_23_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.SNP.PASS.chr-02.g.vcf.gz --TajimaD 8794 --out chr-02.tajD
