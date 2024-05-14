#!/bin/bash
# usage: ./batch-trim-galore.sh chunk_size *read_1.fastq.gz
# for chunk: analysis on login takes 10 min each pair

ml biology trim_galore py-cutadapt/1.18_py27

CHUNK=$1
FQ="${@:2}"
COUNTER=0
for i in $FQ; do

  if [ $COUNTER -eq 0 ]; then
  echo -e "#!/bin/bash\n#SBATCH -p schumer,hns,normal,owners\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=16\n#SBATCH -t 3:00:00\n#SBATCH --mem 30000" > TEMPBATCH.sbatch; fi
  BASE=$( basename $i R1.fastq.gz )
  echo "ml biology trim_galore py-cutadapt/1.18_py27"
  echo "trim_galore --path_to_cutadapt /share/software/user/open/py-cutadapt/1.18_py27/bin/cutadapt --phred33 --quality 30 -e 0.1 --stringency 1 --length 32 --paired --retain_unpaired ${BASE}R1.fastq.gz ${BASE}R2.fastq.gz" >> TEMPBATCH.sbatch

  let COUNTER=COUNTER+1
  if [ $COUNTER -eq $CHUNK ]; then
    sbatch TEMPBATCH.sbatch
  COUNTER=0; fi
done
if [ $COUNTER -ne 0 ]; then
  sbatch TEMPBATCH.sbatch; fi 
