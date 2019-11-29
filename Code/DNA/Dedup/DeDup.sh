BAM_FOLDER=../../../../Data/Processed/BAM_GATK
OUTPUT_FOLDER=../../../../Data/Processed/BAM_Dedup

samples=$(ls ${BAM_FOLDER}/*.bam)

for s in $samples; do 
	samtools rmdup -s $s ${OUTPUT_FOLDER}/$(basename $s .bam)_dedup.bam
	samtools index ${OUTPUT_FOLDER}/$(basename $s .bam)_dedup.bam
done

