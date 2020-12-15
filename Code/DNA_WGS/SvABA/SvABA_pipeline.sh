BAM_FOLDER=$(readlink -e ../../mapping/reseq_mouse_2_and_3/BWA_MEM_merged_BAMs)

BIN_FOLDER=/data_genome1/SharedSoftware/SvABA/bin
# NovoBreak requires a BWA index, not just a FASTA sequence
BWA_INDEX=/data_genome1/References/GenomeIndices/MusMusculus/BWA/GRCm38/GRCm38
#BWA_INDEX=/data_genome1/References/GenomeIndices/Human/human_g1k_v37_decoy_guideseq/human_g1k_v37_decoy.fasta

control="MOUlyuRAAAAAAA"

SAMPLES="MOUlyuRAAABAAA MOUlyuRAAACAAA"


for sample in $SAMPLES; do 
	[ ! -e $sample ] && mkdir $sample
	control_BAM=$(ls ${BAM_FOLDER}/${control}_merged_sorted.bam)
	target_BAM=$(ls ${BAM_FOLDER}/${sample}_merged_sorted.bam)
	pushd $sample
	PATH=$BIN_FOLDER:$PATH; ${BIN_FOLDER}/svaba run -G $BWA_INDEX -t $target_BAM -n $control_BAM -p 8 -a $sample
	popd
done

control="MOUlyuRAAADAAA"

SAMPLES="MOUlyuRAAAEAAA MOUlyuRAAAFAAA"


for sample in $SAMPLES; do
        [ ! -e $sample ] && mkdir $sample
        control_BAM=$(ls ${BAM_FOLDER}/${control}_merged_sorted.bam)
        target_BAM=$(ls ${BAM_FOLDER}/${sample}_merged_sorted.bam)
        pushd $sample
        PATH=$BIN_FOLDER:$PATH; ${BIN_FOLDER}/svaba run -G $BWA_INDEX -t $target_BAM -n $control_BAM -p 8 -a $sample
        popd
done

