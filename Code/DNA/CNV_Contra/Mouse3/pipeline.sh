BAM_FOLDER=../../../../Data/Processed/BAM_Dedup

CONTRA_FOLDER=/data_genome1/SharedSoftware/CONTRA/CONTRA.v2.0.8
REGIONS=../../../../Data/External/AgilentSureSelect/MouseExome_mm10_liftover_Covered_fixed_chroms.bed
REF_FASTA=$(readlink -f /data/References/MusMusculus/Sequences/Genome/GRCm38_Ensembl/Mus_musculus.GRCm38.dna.primary_assembly_fixed_chr.fa)

OUTPUT_FOLDER=../../../../Data/Processed/DNA_Variants_Contra

control="3466_A"
all_chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

#samples=$(ls ${BAM_FOLDER}/*_sorted_dedup.bam)
samples="3466_B 3466_C 3466_D 3466_E 3466_F 3466_G 3466_H"

for sample in $samples; do 
	ANALYSIS_FOLDER=${OUTPUT_FOLDER}/CN_${sample}

	python ${CONTRA_FOLDER}/contra.py -t $REGIONS -s $BAM_FOLDER/${sample}_merged_dedup_realign_dedup.bam -c $BAM_FOLDER/${control}_merged_dedup_realign_dedup.bam -o $ANALYSIS_FOLDER -l

done