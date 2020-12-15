BAM_FOLDER=../GATK_reseq_mouse_2_and_3

STRELKA_FOLDER=/data_genome1/SharedSoftware/Strelka/v1.0.11

#REF_FASTA=/data_genome1/SharedSoftware/GATK/resources/human_g1k_v37_decoy.fasta
REF_FASTA=/data_genome1/References/GenomeIndices/MusMusculus/BWA/GRCm38/GRCm38

# MOUlyuRAAAAAAA MOUlyuRAAABAAA MOUlyuRAAACAAA MOUlyuRAAADAAA MOUlyuRAAAEAAA MOUlyuRAAAFAAA

control="MOUlyuRAAAAAAA"

SAMPLES="MOUlyuRAAABAAA MOUlyuRAAACAAA"

for sample in $SAMPLES; do
        ANALYSIS_FOLDER=./strelkaAnalysis_${sample}
        #[ ! -e $ANALYSIS_FOLDER ] && mkdir $ANALYSIS_FOLDER

        perl $STRELKA_FOLDER/bin/configureStrelkaWorkflow.pl \
               --normal=${BAM_FOLDER}/${control}_merged_dedup_realign.bam \
               --tumor=${BAM_FOLDER}/${sample}_merged_dedup_realign.bam \
               --ref=$REF_FASTA \
               --output-dir=$ANALYSIS_FOLDER \
               --config $STRELKA_FOLDER/etc/strelka_config_bwa_default.ini

        make -j 12 -C $ANALYSIS_FOLDER

done

control="MOUlyuRAAADAAA"

SAMPLES="MOUlyuRAAAEAAA MOUlyuRAAAFAAA"

for sample in $SAMPLES; do
        ANALYSIS_FOLDER=./strelkaAnalysis_${sample}
        #[ ! -e $ANALYSIS_FOLDER ] && mkdir $ANALYSIS_FOLDER

        perl $STRELKA_FOLDER/bin/configureStrelkaWorkflow.pl \
               --normal=${BAM_FOLDER}/${control}_merged_dedup_realign.bam \
               --tumor=${BAM_FOLDER}/${sample}_merged_dedup_realign.bam \
               --ref=$REF_FASTA \
               --output-dir=$ANALYSIS_FOLDER \
               --config $STRELKA_FOLDER/etc/strelka_config_bwa_default.ini

        make -j 12 -C $ANALYSIS_FOLDER

done


