BAM_FOLDER=../../../Data/Processed/BAM_Dedup
pushd $BAM_FOLDER
samtools merge 2702_NI.bam 2702_{A,B,C}_merged_dedup_realign_dedup.bam 
samtools merge 2702_WI.bam 2702_{D,E,F,G,H}_merged_dedup_realign_dedup.bam 
samtools merge 2702_NI.bam 3466_{A,B,C}_merged_dedup_realign_dedup.bam 
samtools merge 2702_WI.bam 3466_{D,E,F,G,H}_merged_dedup_realign_dedup.bam 
popd
