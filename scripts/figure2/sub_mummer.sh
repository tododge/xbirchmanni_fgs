mkdir haplotypes

for file in `ls *.fa *fasta`; \
do ./extract_haplotype_results.sh $file barcoding_dup1.fa; done

mv *utg_utg*.fa xmal*ctg_h*.fa *chr-02*.fa *Assembly-Phased_* ../haplotypes/

cd ./haplotypes/

#nanopore read haplotypes for xbir-BEJU-S231-9-VIII-23-M02-FGS
seqkit grep -p 89fe6f45-6dd4-4e38-87f0-a926a308bacc ~/Nanopore_data/xbir-BEJU-S231-9-VIII-23-M02-FGS.dorado.0.3.3.fastq.gz | seqkit fq2fa | seqkit seq -r -p -t DNA -o xbir-BEJU-S231-9-VIII-23-M02-FGS.dorado.0.3.3.fastq.gz_89fe6f45-6dd4-4e38-87f0-a926a308bacc.fa #FGS
seqkit grep -p 61093bf8-78f3-4155-a09d-a44a3df78a45 ~/Nanopore_data/xbir-BEJU-S231-9-VIII-23-M02-FGS.dorado.0.3.3.fastq.gz | seqkit fq2fa -o xbir-BEJU-S231-9-VIII-23-M02-FGS.dorado.0.3.3.fastq.gz_61093bf8-78f3-4155-a09d-a44a3df78a45.fa #WT

#nanopore read haplotypes for xbir-BEJU-S231-3-IV-23-M01-FGS
seqkit grep -f names ~/Nanopore_data/xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3.fastq.gz -o xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3_spanning_reads.fastq.gz
seqkit grep -p 431af464-ea5d-4b01-8e75-232885ef71f6 xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3_spanning_reads.fastq.gz | seqkit fq2fa -o xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3.fastq.gz_431af464-ea5d-4b01-8e75-232885ef71f6.fa #FGS
seqkit grep -p 18802934-b622-4bb3-9363-5becc4d359f1 xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3_spanning_reads.fastq.gz | seqkit fq2fa | seqkit seq -r -p -t DNA -o xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3.fastq.gz_18802934-b622-4bb3-9363-5becc4d359f1.fa #FGS
seqkit grep -p b1b1cc76-078a-49d3-9833-37a2d493114b xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3_spanning_reads.fastq.gz | seqkit fq2fa | seqkit seq -r -p -t DNA -o xbir-BEJU-S231-3-IV-23-M01-FGS.dorado.0.3.3.fastq.gz_b1b1cc76-078a-49d3-9833-37a2d493114b.fa #WT

#nanopore read haplotypes for xbir-BEJU-S231-9-VIII-23-M01-WT
seqkit grep -f names1 ~/Nanopore_data/xbir-BEJU-S231-9-VIII-23-M01-WT.dorado.0.3.3.fastq.gz -o xbir-BEJU-S231-9-VIII-23-M01-WT.dorado.0.3.3_spanning_reads.fastq.gz
seqkit grep -p 1ca2560c-bfed-49b6-ae19-3cffc7d2a090 xbir-BEJU-S231-9-VIII-23-M01-WT.dorado.0.3.3_spanning_reads.fastq.gz | seqkit fq2fa | seqkit seq -r -p -t DNA -o xbir-BEJU-S231-9-VIII-23-M01-WT.dorado.0.3.3.fastq.gz_1ca2560c-bfed-49b6-ae19-3cffc7d2a090.fa #hap1
seqkit grep -p 709a44bf-d1b7-44a6-841c-c6774f6137e4 xbir-BEJU-S231-9-VIII-23-M01-WT.dorado.0.3.3_spanning_reads.fastq.gz | seqkit fq2fa -o xbir-BEJU-S231-9-VIII-23-M01-WT.dorado.0.3.3.fastq.gz_709a44bf-d1b7-44a6-841c-c6774f6137e4.fa #hap2

ln -s ~/pangenome/barcoding/barcoding.fa  .
ml biology viz gnuplot
MUMMERPATH=/home/groups/schumer/shared_bin/mummer-4.0.0rc1/

for ASSEM in `ls x???-????-*.fa`; \
do SEQ_X=barcoding.fa; \
PREFIX=`echo ${SEQ_X}\_2_${ASSEM}`;
echo ${ASSEM}; \
${MUMMERPATH}/nucmer ${SEQ_X} ${ASSEM} --prefix=${PREFIX}; \
${MUMMERPATH}/delta-filter -m ${PREFIX}.delta > ${PREFIX}.delta.m; \
${MUMMERPATH}/show-coords ${PREFIX}.delta.m -T > ${PREFIX}.delta.m.coords; \
done

#remove duplicated haps
rm *utg000057l.fa.delta* #xbir fgs ref
rm *utg001220l.fa.delta* #xbir wt ref
rm xmal-CHIC-XI-20-M_v2023.2_chr-02.fa.delta* #xmal ref
rm *_UR.199.2.fa.delta* #xbir BEJU3
Rscript process_seq.R
