

if [ $# -lt 3 ]; then
   echo "Usage:"
   echo "DNASeq_calling_pipeline_single_sample BAM_FOLDER SAMPLE ROI.bed [MAX_COVERAGE]"
   exit
fi


BAM_FOLDER=$1
sample=$2
ROI=$3
if [ $# -eq 4 ]; then
   MAX_COV=$4
else 
   MAX_COV=6000
fi

echo "MAX_COV: " $MAX_COV

# GATK versions can be retrieved from https://software.broadinstitute.org/gatk/download/archive
GATK_FOLDER=/data_genome1/SharedSoftware/GATK/GATK_v3.4.0
# Picard is available from https://broadinstitute.github.io/picard/
PICARD_FOLDER=/data_genome1/SharedSoftware/Picard

# GRCm38 from Ensembl ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
REF_FASTA=/data_genome2/References/MusMusculus/Sequences/Genome/GRCm38_Ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa

# based on ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/genotype/SC_MOUSE_GENOMES.genotype.vcf.gz ; known SNP for ancestral strains of current mouse strain
DBSNP_FILE=${BAM_FOLDER}/../../External/MouseGenotype/129P2_and_C57BL6_calls.vcf
# from ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/, reordered chromosomes in VCF
GOLD_INDELS_FILE=${BAM_FOLDER}/../../External/MouseGenotype/mgp.v3.indels.rsIDdbSNPv137_chromorder_1_10.vcf

# put the temp folder where we know there will be enough space available 
TMP_FOLDER=$(pwd)/tmp
if ! [ -e $TMP_FOLDER ]; then mkdir $TMP_FOLDER; fi

#BASH specific
#http://stackoverflow.com/questions/1527049/bash-join-elements-of-an-array
function join { local IFS="$1"; shift; echo "$*"; }

	all_sample_lanes=""
        all_sample_lane_cnt=0

    for sample_lane in ${BAM_FOLDER}/${sample}*_sorted.bam; do

	INPUT_BAM_FILE=$sample_lane
	i=$(basename $sample_lane _sorted.bam)

	# skip per lane files that already have been generated
	if [ -e ${i}_rg_dedup_realign_recal.bam ]; then 
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
	   java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/AddOrReplaceReadGroups.jar I=$INPUT_BAM_FILE O=$RG_ADDED_SORTED_BAM SO=coordinate RGID=$RG RGLB=$LIBRARY RGPL=$PLATFORM RGPU=$RGPU RGSM=$sample TMP_DIR=$TMP_FOLDER
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
	    java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T RealignerTargetCreator -L $ROI -R $REF_FASTA -I $DEDUPPED_BAM -nt 8 -known $GOLD_INDELS_FILE -o $TARGET_INTERVAL_LIST
	    # then realign BAM
	    java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T IndelRealigner -L $ROI -R $REF_FASTA  --maxReadsForRealignment 100000 -I $DEDUPPED_BAM -targetIntervals $TARGET_INTERVAL_LIST -known $GOLD_INDELS_FILE -o $REALIGNED_BAM 
        fi

	RECAL_TAB_FILE=${i}_recal_data.tab
	POST_RECAL_TAB_FILE=${i}_post_recal_data.tab
	RECAL_PLOTS=${i}_recalibration_plots.pdf
	RECALIBRATED_BAM=${i}_rg_dedup_realign_recal.bam

	if ! [ -e $RECALIBRATED_BAM ]; then
#	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -dcov $MAX_COV -L $ROI -R $REF_FASTA -I $REALIGNED_BAM  -knownSites $DBSNP_FILE -knownSites $GOLD_INDELS_FILE -o $RECAL_TAB_FILE
#	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -dcov $MAX_COV -L $ROI -R $REF_FASTA -I $REALIGNED_BAM  -knownSites $DBSNP_FILE -knownSites $GOLD_INDELS_FILE -BQSR $RECAL_TAB_FILE -o $POST_RECAL_TAB_FILE
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -dcov $MAX_COV -L $ROI -R $REF_FASTA -I $REALIGNED_BAM  -knownSites $DBSNP_FILE  -o $RECAL_TAB_FILE
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -dcov $MAX_COV -L $ROI -R $REF_FASTA -I $REALIGNED_BAM  -knownSites $DBSNP_FILE  -BQSR $RECAL_TAB_FILE -o $POST_RECAL_TAB_FILE
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T AnalyzeCovariates -dcov $MAX_COV -L $ROI -R $REF_FASTA  -before $RECAL_TAB_FILE -after $POST_RECAL_TAB_FILE -plots $RECAL_PLOTS
	   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T PrintReads -dcov $MAX_COV -nct 8 -I $REALIGNED_BAM -R $REF_FASTA  -BQSR $RECAL_TAB_FILE -o $RECALIBRATED_BAM

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

   REALIGNED_MERGED_BAM=${sample}_merged_dedup_realign.bam

   if [ ! -e $REALIGNED_MERGED_BAM ]; then

       MERGED_BAM=${sample}_merged.bam

       if [ $all_sample_lane_cnt -gt 1 ]; then 

	   if [ ! -e $MERGED_BAM ]; then
	       java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MergeSamFiles.jar $inputfiles OUTPUT=$MERGED_BAM TMP_DIR=$TMP_FOLDER 
	   fi

	   # optional: rerun dedup and realign on per-sample BAMs
	   DEDUPPED_MERGED_BAM=${sample}_merged_dedup.bam
	   DEDUP_METRICS_FILE=${sample}.merged.dedup.metrics
	   if ! [ -e $DEDUPPED_MERGED_BAM ]; then
	       java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MarkDuplicates.jar I=${sample}_merged.bam O=$DEDUPPED_MERGED_BAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$DEDUP_METRICS_FILE TMP_DIR=$TMP_FOLDER 
	   fi

	   TARGET_INTERVAL_LIST=${sample}_target_intervals.list
	   # Realign around InDels
	   # 1st define target list
	   if ! [ -e $REALIGNED_MERGED_BAM ]; then
	       java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T RealignerTargetCreator -dcov $MAX_COV -R $REF_FASTA -I $DEDUPPED_MERGED_BAM -nt 8 -known $GOLD_INDELS_FILE -o $TARGET_INTERVAL_LIST
	       # then realign BAM
	       java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T IndelRealigner -dcov $MAX_COV -R $REF_FASTA  -I $DEDUPPED_MERGED_BAM --maxReadsForRealignment 100000 -targetIntervals $TARGET_INTERVAL_LIST -known $GOLD_INDELS_FILE -o $REALIGNED_MERGED_BAM 
	   fi
       else
	   REALIGNED_MERGED_BAI=${sample}_merged_dedup_realign.bai
	   if ! [ -e $REALIGNED_MERGED_BAM ]; then	
	       cp $(echo $all_sample_lanes | sed -e 's/INPUT=//' | sed -e 's/"//g' | sed -e 's/ //g') $REALIGNED_MERGED_BAM
	       cp $(echo $all_sample_lanes | sed -e 's/INPUT=//' | sed -e 's/"//g' | sed -e 's/ //g' | sed -e 's/.bam/.bai/') $REALIGNED_MERGED_BAI
	   fi
       fi

       if [ -e $REALIGNED_MERGED_BAM ]; then
           rm $MERGED_BAM $DEDUPPED_MERGED_BAM
       fi
   fi


   OUTPUT_VCF=${sample}.vcf

   if ! [ -e $OUTPUT_VCF ]; then
	   java -Xmx16g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T UnifiedGenotyper -L $ROI -R $REF_FASTA -nct 1 -nt 1 -glm BOTH --dbsnp $DBSNP_FILE -I $REALIGNED_MERGED_BAM -stand_call_conf 30.0 -stand_emit_conf 10.0 -o $OUTPUT_VCF 
   fi

   RAW_SNP_VCF=${sample}_raw_snps.vcf
   RAW_INDEL_VCF=${sample}_raw_indels.vcf

   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T SelectVariants -R $REF_FASTA -V $OUTPUT_VCF -selectType SNP -o $RAW_SNP_VCF 
   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T SelectVariants -R $REF_FASTA -V $OUTPUT_VCF -selectType INDEL -o $RAW_INDEL_VCF

   FILTERED_SNP_VCF=${sample}_filtered_snps.vcf
   FILTERED_INDEL_VCF=${sample}_filtered_indels.vcf
   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_FASTA -V $RAW_SNP_VCF -window 35 -cluster 3 --filterName QualThreshold --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -o $FILTERED_SNP_VCF
   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_FASTA -V $RAW_INDEL_VCF -window 35 -cluster 3 --filterName QualThreshold --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" -o $FILTERED_INDEL_VCF

   if [ -e $FILTERED_SNP_VCF ]; then
      rm $RAW_SNP_VCF
   fi

   if [ -e $FILTERED_INDEL_VCF ]; then
      rm $RAW_INDEL_VCF
   fi
