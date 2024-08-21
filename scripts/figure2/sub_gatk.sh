#!/bin/sh
#SBATCH --ntasks=1 
#SBATCH --mem=64000
#SBATCH --cpus-per-task=1
#SBATCH -p schumer
#SBATCH --time=6:00:00
#SBATCH --job-name=submit

#ls raw/*realigned*bam > bam_list.txt

#perl submit_gatk-indel-hc_jobs_list.pl bam_list.txt submit_gatk-indel-hc.sh /home/groups/schumer/shared_bin/ xbir-COAC-16-VIII-22-M_v2023.1.fa chr-02.list chr-02

#for file in `ls raw/*.chr-02.g.vcf`; do echo "bgzip -c ${file} > ${file}.gz"; \
#bgzip -c ${file} > ${file}.gz; \
#echo "tabix -p vcf ${file}.gz"; \
#tabix -p vcf ${file}.gz; done

####workflow after merging with bcftools####
ml biology bcftools
#bcftools merge raw/*vcf.gz -Oz -o COAC_2016-2018_27_allhighcovinds.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz

#bcftools view COAC_2016-2018_27_allhighcovinds.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz |
#    sed 's/=nan\([;\t]\)/=NaN\1/g' |
#    bgzip > COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz
#tabix -p vcf COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz

#module load java gatk bwa samtools

#gatk SelectVariants \
#-R xbir-COAC-16-VIII-22-M_v2023.1.fa \
#-V COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz \
#--select-type-to-include SNP \
#-O COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.g.vcf.gz

#gatk VariantFiltration \
#    -V COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.g.vcf.gz \
#    -filter "QD < 2.0" --filter-name "QD2" \
#    -filter "QUAL < 30" --filter-name "QUAL30" \
#    -filter "SOR > 3.0" --filter-name "SOR3" \
#    -filter "FS > 60.0" --filter-name "FS60" \
#    -filter "MQ < 40.0" --filter-name "MQ40" \
#    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#    --exclude-filtered TRUE \
#    -O COAC_2016-2018_27_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.PASS.g.vcf.gz

#ls raw/*ma*gz*realigned*bam raw/TETI*gz*realigned*bam > bam_list.txt

#perl submit_gatk-indel-hc_jobs_list.pl bam_list.txt submit_gatk-indel-hc.sh /home/groups/schumer/shared_bin/ xbir-COAC-16-VIII-22-M_v2023.1.fa chr-02.list chr-02

#perl submit_gatk-indel-hc_jobs_list.pl bam_list.txt submit_gatk-indel-hc.sh /home/groups/schumer/shared_bin/ xbir-COAC-16-VIII-22-M_v2023.1.fa chr.list all_chroms

for file in `ls raw/*ma*gz*.chr-02.g.vcf raw/*TETI*gz*.chr-02.g.vcf`; do echo "bgzip -c ${file} > ${file}.gz"; \
bgzip -c ${file} > ${file}.gz; \
echo "tabix -p vcf ${file}.gz"; \
tabix -p vcf ${file}.gz; done


ml biology bcftools

bcftools merge raw/*vcf.gz -Oz -o COAC_malinche_2016-2018_31_allhighcovinds.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz
bcftools view COAC_malinche_2016-2018_31_allhighcovinds.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz |
    sed 's/=nan\([;\t]\)/=NaN\1/g' |
    bgzip > COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz
tabix -p vcf COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz


module load java gatk bwa samtools

gatk SelectVariants \
-R xbir-COAC-16-VIII-22-M_v2023.1.fa \
-V COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.g.vcf.gz \
--select-type-to-include SNP \
-O COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.g.vcf.gz

gatk VariantFiltration \
    -V COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.g.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.PASS.g.vcf.gz

gatk SelectVariants \
-R xbir-COAC-16-VIII-22-M_v2023.1.fa \
-V COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.PASS.g.vcf.gz \
--exclude-filtered TRUE \
-O COAC_malinche_2016-2018_31_allhighcovinds.fixed.fq.gz_2_xbir-COAC-16-VIII-22-M_v2023.1.fa.sorted.dedup.realigned.bam.chr-02.SNP.PASS.only.g.vcf.gz


