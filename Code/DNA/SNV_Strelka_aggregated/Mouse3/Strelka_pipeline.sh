BAM_FOLDER=../../../../Data/Processed/BAM_Dedup

STRELKA_FOLDER=/data_genome1/SharedSoftware/Strelka/v1.0.11
REGIONS=../../../../Data/External/AgilentSureSelect/MouseExome_mm10_liftover_Covered_fixed_chroms_sorted.bed.bgz
REF_FASTA=$(readlink -f /data/References/MusMusculus/Sequences/Genome/GRCm38_Ensembl/Mus_musculus.GRCm38.dna.primary_assembly_fixed_chr.fa)

control="3466_NI"
all_chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

samples="3466_WI"

for sample in $samples; do 
	ANALYSIS_FOLDER=./strelkaAnalysis_${sample}

	perl $STRELKA_FOLDER/bin/configureStrelkaWorkflow.pl \
	       --normal=${BAM_FOLDER}/${control}.bam \
	       --tumor=${BAM_FOLDER}/${sample}.bam\
	       --ref=$REF_FASTA \
	       --output-dir=$ANALYSIS_FOLDER \
	       --config $STRELKA_FOLDER/etc/strelka_config_bwa_targetedseq.ini

	make -j 12 -C $ANALYSIS_FOLDER

done
