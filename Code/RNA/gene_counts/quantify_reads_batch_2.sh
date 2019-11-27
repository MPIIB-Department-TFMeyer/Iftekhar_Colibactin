#!/bin/bash
ulimit -v 30258917744
INDEX_FOLDER=../../../Data/External/GencodeM12/GencodeM12
SEQ_FOLDER=../../../../seqs/RNA/FASTQ/batch2_all
function join { local IFS="$1"; shift; echo "$*"; }
SALMON_EXEC=/data_genome1/SharedSoftware/Salmon/bin/salmon
GENE_MODEL=../../../Data/External/GencodeM12/GencodeM12/gencode.vM12.annotation.gtf
OUTPUT_FOLDER=../../../Data/Processed/RNA_counts

cat ../../../metadata/RNA/batch2/all_samples_aggregated.txt | while read sample_line
do # 1st pass
  [ -z "$sample_line" ] && continue
  set $sample_line
  sample=$(echo $1)
  reads1=$(echo $2 | sed -e 's/,/ /g')
  reads2=$(echo $3 | sed -e 's/,/ /g')

  $SALMON_EXEC quant -i $INDEX_FOLDER -1 $reads1 -2 $reads2 -o ${OUTPUT_FOLDER}/$sample --seqBias -g $GENE_MODEL -l A

done

