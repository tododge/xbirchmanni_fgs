#!/bin/bash

# Submit a kallisto run for each pair of reads passed in bulk
# usage: bash ./batch-kallisto.sh transcriptome.idx 2 *.R1_val_1.fq.gz
# usage: bash ./batch-kallisto.sh genome.idx chunk_size reads

ml biology kallisto/0.46.1

IDX_FILE=$1
CHUNK=$2
FQ="${@:3}"
COUNTER=0
for i in $FQ; do

  if [ $COUNTER -eq 0 ]; then
  echo -e "#!/bin/bash\n#SBATCH -p schumer,owners\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=16\n#SBATCH -t 3:00:00\n#SBATCH --mem 48000" > TEMPBATCH.sbatch; fi
  DIR=$( dirname $i)
  BASE=$( basename $i .R1_val_1.fq.gz )
  echo "$BASE"
    echo "kallisto quant -i $IDX_FILE -o ${BASE}_kallisto ${DIR}/${BASE}.R1_val_1.fq.gz ${DIR}/${BASE}.R2_val_2.fq.gz" >> TEMPBATCH.sbatch

  let COUNTER=COUNTER+1
  if [ $COUNTER -eq $CHUNK ]; then
    sbatch TEMPBATCH.sbatch
  COUNTER=0; fi
done
if [ $COUNTER -ne 0 ]; then
  sbatch TEMPBATCH.sbatch; fi 
