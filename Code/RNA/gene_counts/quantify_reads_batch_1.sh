#!/bin/bash
ulimit -v 30258917744
INDEX_FOLDER=../../../Data/External/GencodeM12/GencodeM12
SEQ_FOLDER=../../../Data/Raw/seqs/RNA
function join { local IFS="$1"; shift; echo "$*"; }
SALMON_EXEC=/data_genome1/SharedSoftware/Salmon/bin/salmon
GENE_MODEL_FOLDER=../../../Data/External/GencodeM12
OUTPUT_FOLDER=../../../Data/Processed/RNA_counts

cat ../../../Data/metadata/RNA/Batch1/all_samples_aggregated.txt | while read sample_line
do # 1st pass
  [ -z "$sample_line" ] && continue
  set $sample_line
  sample=$(echo $1)
  reads1=$(echo $2)
#  reads2=$(echo $3)

  $SALMON_EXEC quant -i $INDEX_FOLDER -r $reads1 -o ${OUTPUT_FOLDER}/$sample --seqBias -g ${GENE_MODEL_FOLDER}/gencode.vM12.annotation.gtf -l A

done

