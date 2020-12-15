

if [ $# -lt 2 ]; then
   echo "Usage:"
   echo "DNASeq_calling_pipeline_single_sample BAM_FOLDER SAMPLE"
   exit
fi

BAM_FOLDER=$1
sample=$2

BAM_SUFFIX="_sorted.bam"

GATK_FOLDER=/data_genome1/SharedSoftware/GATK/GATK_v3.7
PICARD_FOLDER=/data_genome1/SharedSoftware/Picard
REF_FASTA=/data_genome2/References/MusMusculus/Sequences/Genome/GRCm38_Ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa
#DBSNP_FILE=/data_genome1/References/MusMusculus/Variation/dbSNP/dbSNP_v146_GRCm38/VCF/dbSNP_all.chr.vcf
DBSNP_FILE=/data/projects/MB208_Colibactin/analysis/Genotype_mouse/129P2_and_C57BL6_calls.vcf
GOLD_INDELS_FILE=/data_genome1/References/MusMusculus/Variation/MGP/mgp.v3.indels.rsIDdbSNPv137_chromorder_1_10.vcf

# put the temp folder where we know there will be enough space available 
TMP_FOLDER=$(pwd)/tmp
if ! [ -e $TMP_FOLDER ]; then mkdir $TMP_FOLDER; fi

#BASH specific
#http://stackoverflow.com/questions/1527049/bash-join-elements-of-an-array
function join { local IFS="$1"; shift; echo "$*"; }

# This is the final BAM
REALIGNED_MERGED_BAM=${sample}_merged_dedup_realign.bam
REALIGNED_MERGED_BAI=${sample}_merged_dedup_realign.bai

if ! [ -e $REALIGNED_MERGED_BAM ]; then

    all_sample_lanes=""
    all_sample_lane_cnt=0


    MERGED_BAM=${sample}_merged.bam
    for sample_lane in ${BAM_FOLDER}/${sample}_*${BAM_SUFFIX}; do

	INPUT_BAM_FILE=$sample_lane
	i=$(basename $sample_lane $BAM_SUFFIX)

	# skip per lane files that already have been generated
	if [ -e ${i}_rg_dedup_realign_recal.bam -o -e $MERGED_BAM ]; then 
            all_sample_lanes+="INPUT="
            all_sample_lanes+=$i
            all_sample_lanes+="_rg_dedup_realign_recal.bam "
            all_sample_lane_cnt+=1
	    continue
	fi

	run_lane=$(echo $i | awk 'BEGIN{FS="_"} {print $3"_"$4;}')
	RG=${run_lane}_${sample}
	LIBRARY=$sample
	PLATFORM="illumina"
	RGPU=$run_lane
	SAMPLE=$sample

	RG_ADDED_SORTED_BAM=${i}_rg.bam
	if ! [ -e $RG_ADDED_SORTED_BAM ]; then
	   samtools view -F 8 -u -O BAM $INPUT_BAM_FILE | samtools addreplacerg -r "@RG\tID:${RG}\tPL:${PLATFORM}\tPU:${RGPU}\tLB:${LIBRARY}\tSM:${SAMPLE}" -@ 4 -O BAM -o $RG_ADDED_SORTED_BAM - 
	   #java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/AddOrReplaceReadGroups.jar I=$INPUT_BAM_FILE O=$RG_ADDED_SORTED_BAM SO=coordinate RGID=$RG RGLB=$LIBRARY RGPL=$PLATFORM RGPU=$RGPU RGSM=$sample TMP_DIR=$TMP_FOLDER
	fi

	DEDUPPED_BAM=${i}_rg_dedup.bam
	DEDUP_METRICS_FILE=${i}.dedup.metrics
        if ! [ -e $DEDUPPED_BAM ]; then
	    java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MarkDuplicates.jar I=$RG_ADDED_SORTED_BAM O=$DEDUPPED_BAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$DEDUP_METRICS_FILE TMP_DIR=$TMP_FOLDER
        fi

	#REORDERED_BAM=${i}_rg_dedup_reordered.bam
	#java -jar ${PICARD_FOLDER}/ReorderSam.jar I=$DEDUPPED_BAM O=$REORDERED_BAM REFERENCE=$REF_FASTA ALLOW_INCOMPLETE_DICT_CONCORDANCE=true
	#samtools index $REORDERED_BAM

	TARGET_INTERVAL_LIST=${i}_target_intervals.list
	REALIGNED_BAM=${i}_rg_dedup_realign.bam
	# Realign around InDels
	if ! [ -e $REALIGNED_BAM ]; then
	    # 1st define target list
	    java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_FASTA -I $DEDUPPED_BAM -nt 8 -known $GOLD_INDELS_FILE -o $TARGET_INTERVAL_LIST
	    # then realign BAM
	    java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_FASTA  -I $DEDUPPED_BAM -targetIntervals $TARGET_INTERVAL_LIST -known $GOLD_INDELS_FILE -o $REALIGNED_BAM 
        fi

	RECAL_TAB_FILE=${i}_recal_data.tab
	POST_RECAL_TAB_FILE=${i}_post_recal_data.tab
	RECAL_PLOTS=${i}_recalibration_plots.pdf
	RECALIBRATED_BAM=${i}_rg_dedup_realign_recal.bam

	if ! [ -e $RECALIBRATED_BAM ]; then
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $REF_FASTA  -I $REALIGNED_BAM  -knownSites $DBSNP_FILE -knownSites $GOLD_INDELS_FILE -o $RECAL_TAB_FILE
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $REF_FASTA  -I $REALIGNED_BAM  -knownSites $DBSNP_FILE -knownSites $GOLD_INDELS_FILE -BQSR $RECAL_TAB_FILE -o $POST_RECAL_TAB_FILE
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF_FASTA   -before $RECAL_TAB_FILE -after $POST_RECAL_TAB_FILE -plots $RECAL_PLOTS
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T PrintReads -nct 8 -I $REALIGNED_BAM -R $REF_FASTA  -BQSR $RECAL_TAB_FILE -o $RECALIBRATED_BAM

	fi 

	if [ -e $RECALIBRATED_BAM ]; then
	    #only remove intermediate files if we reached the last step
	    rm $RG_ADDED_SORTED_BAM $DEDUPPED_BAM $REALIGNED_BAM
	fi

	all_sample_lanes+="INPUT="
	all_sample_lanes+=$i
	all_sample_lanes+="_rg_dedup_realign_recal.bam "
    all_sample_lane_cnt+=1

   done 

   # now merge per-lane BAMs to per-sample BAMs
   #inputfiles=$(join ' INPUT=' $all_sample_lanes)
   echo $all_sample_lanes
   inputfiles=$all_sample_lanes
   
   # If we have more than one read group/lane/run file per sample
   if [ $all_sample_lane_cnt -gt 1 ]; then
       # If merged bam does not exist yet, merge individual files
       if ! [ -e $MERGED_BAM ]; then
	   java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MergeSamFiles.jar $inputfiles OUTPUT=$MERGED_BAM TMP_DIR=$TMP_FOLDER 
       fi
    
       # rerun dedup and realign on per-sample BAMs
       DEDUPPED_MERGED_BAM=${sample}_merged_dedup.bam
       DEDUP_METRICS_FILE=${sample}.merged.dedup.metrics
       if ! [ -e $DEDUPPED_MERGED_BAM ]; then
	   java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MarkDuplicates.jar I=${sample}_merged.bam O=$DEDUPPED_MERGED_BAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$DEDUP_METRICS_FILE TMP_DIR=$TMP_FOLDER 
       fi
    
       TARGET_INTERVAL_LIST=${sample}_target_intervals.list
       # Realign around InDels
       # 1st define target list
       if ! [ -e $REALIGNED_MERGED_BAM ]; then
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_FASTA  -I $DEDUPPED_MERGED_BAM -nt 8 -known $GOLD_INDELS_FILE -o $TARGET_INTERVAL_LIST
	   # then realign BAM
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_FASTA  -I $DEDUPPED_MERGED_BAM -targetIntervals $TARGET_INTERVAL_LIST -known $GOLD_INDELS_FILE -o $REALIGNED_MERGED_BAM 
       fi
    
       if [ -e $REALIGNED_MERGED_BAM ]; then
	   rm $MERGED_BAM $DEDUPPED_MERGED_BAM
       fi
   # single RG/lane/run per sample: Do not run dedup/realignment again
   else
	if ! [ -e $REALIGNED_MERGED_BAM ]; then        
	  mv $(echo $all_sample_lanes | sed -e 's/INPUT=//' | sed -e 's/"//g' | sed -e 's/ //g') $REALIGNED_MERGED_BAM
	  mv $(echo $all_sample_lanes | sed -e 's/INPUT=//' | sed -e 's/"//g' | sed -e 's/ //g' | sed -e 's/.bam/.bai/') $REALIGNED_MERGED_BAI
	fi
   fi
fi
   
   ############################################################################
   # Call variants
   ############################################################################

   OUTPUT_VCF=${sample}.vcf

   if ! [ -e $OUTPUT_VCF ]; then
	   java -Xmx16g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_FASTA  -nct 1 -nt 1 -glm BOTH --dbsnp $DBSNP_FILE -I $REALIGNED_MERGED_BAM -stand_call_conf 30.0 -o $OUTPUT_VCF 
   fi

   RECAL_FILE=${sample}.recal
   TRANCHES_FILE=${sample}_var_recal.tranches
   RSCRIPT_FILE=${sample}_var_recal.plots.R
   RECAL_SNP_FILE=${sample}_recal_snp_raw_indel.vcf

    if ! [ -e $RECAL_SNP_FILE ]; then
    java -Xmx4g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_FASTA -input $OUTPUT_VCF \
	   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP_FILE \
	   -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP \
	   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	   --recal_file $RECAL_FILE \
	   --tranches_file $TRANCHES_FILE \
	   --rscript_file $RSCRIPT_FILE

    java -Xmx6g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_FASTA -input $OUTPUT_VCF \
	   -mode SNP --ts_filter_level 99.0 \
	   --recal_file $RECAL_FILE \
	   --tranches_file $TRANCHES_FILE \
	   -o $RECAL_SNP_FILE

    fi

    RECAL_INDEL_FILE=${sample}_indel.recal
    TRANCHES_INDEL_FILE=${sample}_indel_var_recal.tranches
    RSCRIPT_INDEL_FILE=${sample}_indel_var_recal.plots.R
    RECAL_SNP_INDEL_FILE=${sample}_recal_snp_indel.vcf

    if ! [ -e $RECAL_SNP_INDEL_FILE ]; then
	   java -Xmx4g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_FASTA -input $RECAL_SNP_FILE \
		   -resource:mills,known=true,training=true,truth=true,prior=12.0 $GOLD_INDELS_FILE \
		   -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL \
		   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 \
		   --recal_file $RECAL_INDEL_FILE \
		   --tranches_file $TRANCHES_INDEL_FILE \
		   --rscript_file $RSCRIPT_INDEL_FILE

	   java -Xmx6g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_FASTA -input $RECAL_SNP_FILE \
		   -mode INDEL --ts_filter_level 99.0 \
		   --recal_file $RECAL_INDEL_FILE \
		   --tranches_file $TRANCHES_INDEL_FILE \
		   -o $RECAL_SNP_INDEL_FILE
    fi

#   FILTERED_VCF=${sample}_filtered.vcf
   #java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_FASTA -V $OUTPUT_VCF -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $FILTERED_VCF


