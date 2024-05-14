#!/bin/bash
#SBATCH --job-name=longshot
#SBATCH --time=12:00:00
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH -p schumer

ml biology samtools python/3.9.0

#/home/groups/schumer/shared_bin/minimap2 -ax map-hifi -t 16 xbir-COAC-16-VIII-22-M_v2023.1.fa XDOVE_20230106_S64049_PL100275076-1_C01.ccs.fastq.gz XDOVE_20230106_S64049_PL100275076-1_D01.ccs.fastq.gz | samtools sort -@ 4 -m 4G > xbir-COAC-16-VIII-22-M-M01-WT.hifi.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.bam
#samtools index -@ 4 xbir-COAC-16-VIII-22-M-M01-WT.hifi.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.bam

BAM=xbir-COAC-16-VIII-22-M-M01-WT.hifi.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.bam

#BAM_PREFIX=${BAM%%.bam}
#echo $BAM_PREFIX

#samtools view ${BAM} chr-02 -bo ${BAM}_chr-02.bam
samtools index -@ 4 ${BAM}_chr-02.bam

longshot --bam $BAM --ref xbir-COAC-16-VIII-22-M_v2023.1.fa --out ${BAM_PREFIX}.vcf


#bam2fastq -o xbir-COAC-S238-19-IX-23-M01-FGS.hifi m84046_231104_170339_s3.hifi_reads.bc2026.bam
#/home/groups/schumer/shared_bin/minimap2 -ax map-hifi -t 16 xbir-COAC-16-VIII-22-M_v2023.1.fa xbir-COAC-S238-19-IX-23-M01-FGS.hifi.fastq.gz | samtools sort -@ 4 -m 4G > xbir-COAC-S238-19-IX-23-M01-FGS.hifi.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.bam
#samtools index -@ 4 xbir-COAC-S238-19-IX-23-M01-FGS.hifi.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.bam

BAM=xbir-COAC-S238-19-IX-23-M01-FGS.hifi.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.bam

BAM_PREFIX=${BAM%%.bam}
echo $BAM_PREFIX

samtools view ${BAM} chr-02 -bo ${BAM}_chr-02.bam
samtools index -@ 4 ${BAM}_chr-02.bam

#longshot --bam $BAM --ref xbir-COAC-16-VIII-22-M_v2023.1.fa --out ${BAM_PREFIX}.vcf

