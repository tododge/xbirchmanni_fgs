#!/bin/bash
#SBATCH --job-name=psmc
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH -p schumer

#mkdir results/bootstrap/
#VCF=xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-BEJU-S231-9-VIII-23-M01-WT.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-BEJU-S231-9-VIII-23-M02-FGS.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-COAC-S217-22-VI-23-M01-WT.guppy646.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-COAC-S217-22-VI-23-M03-WT.guppy646.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-COAC-S217-9-VIII-23-M01-WT.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-COAC-S217-9-VIII-23-M02-FGS.bonito.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-COAC-S217-9-VIII-23-M03-FGS.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-IZAP-S234-9-VIII-23-M01-FGS.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
#VCF=xbir-IZAP-S234-9-VIII-23-M02-WT.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf
VCF=xbir-IZAP-S234-9-VIII-23-M03-FGS.dorado.0.3.3.fastq_2_xbir-COAC-16-VIII-22-M_v2023.1.fa_sorted.vcf

VCF_PREFIX=${VCF%%.vcf}
echo ${VCF_PREFIX}

seq 51 100 | xargs -i echo ./psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ./results/bootstrap/$VCF_PREFIX.round-{}.psmc ./input/$VCF_PREFIX.split.psmcfa | sh
