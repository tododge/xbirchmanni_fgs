#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: ./extract_blast_results.sh <ref.fa> <query> <flank_size>"
    exit 1
fi

ml biology ncbi-blast+ R
#requires seqkit

FASTA=$1
QUERY=$2
FLANK=$3

if [[ -n $FASTA ]]; then
  PREFIX=$(echo "$FASTA" | cut -d'_' -f1)
  if [[ -z $PREFIX ]]; then
    PREFIX=$(echo "$FASTA" | cut -d'.' -f1)
  fi
else
  echo "Usage: $0 <FASTA>"
  exit 1
fi

makeblastdb -in ${FASTA} -dbtype nucl
blastn -db ${FASTA} -query ${QUERY} -outfmt 6 | awk '$4 > 650' > blastresults_kitlga_${PREFIX}_blastout
echo "blast results for "${PREFIX} " > 650bp:"
cat blastresults_kitlga_${PREFIX}_blastout

rm ${FASTA}.n*

Rscript get_coords_blast.R blastresults_kitlga_${PREFIX}_blastout

IFS=$'\n'
for line in $(cat blastresults_kitlga_${PREFIX}_blastout_start_stop); do
    IFS=$'\t' read -r col1 col2 col3 col4 col5 <<< "$line"
    SEQ="$col1"
    if [ "$col4" -lt "$col5" ]; then
        echo "${SEQ} is forward" "$((col4+${FLANK}))" "$((col5-${FLANK}))"
        seqkit subseq --chr ${SEQ} --region $((col4+${FLANK})):$((col5-${FLANK})) ${FASTA} -o ${PREFIX}_${SEQ}_$((col4+${FLANK}))_$((col5-${FLANK})).fa
    else
        echo "${SEQ} is reverse" "$((col4-${FLANK}))" "$((col5+${FLANK}))"
        seqkit subseq --chr ${SEQ} --region $((col5+${FLANK})):$((col4-${FLANK})) ${FASTA} | seqkit seq -r -p -t DNA -o ${PREFIX}_${SEQ}_$((col5+${FLANK}))_$((col4-${FLANK})).fa
    fi
done
