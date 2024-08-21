#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: ./extract_hap_results.sh <ref.fa> <query>"
    exit 1
fi
ml biology ncbi-blast+ R
#requires seqkit

FASTA=$1
QUERY=$2

if [[ -n $FASTA ]]; then
  PREFIX=$(echo "$FASTA" | cut -d'_' -f1)
  if [[ -z $PREFIX ]]; then
    PREFIX=$(echo "$FASTA" | cut -d'.' -f1)
  fi
else
  echo "Usage: $0 <FASTA>"
  exit 1
fi

PREFIX2=`basename "$(echo $FASTA | sed -E 's/\.fa$|\.fasta$//')"`

makeblastdb -in ${FASTA} -dbtype nucl
blastn -db ${FASTA} -query ${QUERY} -outfmt 6 | awk '$4 > 1000' > blastresults_kitlga_${PREFIX}_blastout
echo "blast results for "${PREFIX} " > 650bp:"
cat blastresults_kitlga_${PREFIX}_blastout

rm ${FASTA}.n*

Rscript get_coords_blast.R blastresults_kitlga_${PREFIX}_blastout

IFS=$'\n'
for line in $(cat blastresults_kitlga_${PREFIX}_blastout_start_stop); do
    IFS=$'\t' read -r col1 col2 col3 col4 col5 <<< "$line"
    SEQ="$col1"
    if [ "$col4" -lt "$col5" ]; then
        echo "${SEQ} is forward"
        seqkit grep -p ${SEQ} ${FASTA} -o ${PREFIX2}_${SEQ}.fa
    else
        echo "${SEQ} is reverse"
        seqkit grep -p ${SEQ} ${FASTA} | seqkit seq -r -p -t DNA -o ${PREFIX2}_${SEQ}.fa
    fi
done
