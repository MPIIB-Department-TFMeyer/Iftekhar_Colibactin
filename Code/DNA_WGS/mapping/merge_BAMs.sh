#!/bin/bash

SAMTOOLS=/data_genome1/SharedSoftware/samtools_htslib/samtools
METAFILE_FOLDER=../../../metadata/reseq_mouse_2_and_3
function join { local IFS="$1"; shift; echo "$*"; }


cat $METAFILE_FOLDER/all_samples_aggregated_paired.txt | while read sample_line
do
  [ -z "$sample_line" ] && continue
  set $sample_line
  sample=$(echo $1)

  OUTFILE=${sample}_merged_sorted.bam

  $SAMTOOLS merge -@ 4 -O BAM $OUTFILE ../BWA_MEM/${sample}_*.bam
  $SAMTOOLS index -b $OUTFILE

done


