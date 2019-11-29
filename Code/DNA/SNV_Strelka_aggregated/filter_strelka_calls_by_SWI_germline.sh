MOUSE2_calls=./Mouse2/Strelka/called_variants/SnpEff/Colibactin_Strelka_calls_merged_SnpEff.vcf.gz
MOUSE3_calls=./Mouse3/Strelka/called_variants/SnpEff/Colibactin_Strelka_calls_merged_SnpEff.vcf.gz

OUTPUT_FOLDER=../../../Data/Processed/DNA_Variants_Strelka_aggregated

FEEDER_CELL_GERMLINE=../../../Data/Processed/Feeder_genome/Mm_filtered_SNV_InDel.vcf.gz

bcftools view -O z -o ${OUTPUT_FOLDER}/mouse2_calls_filtered.vcf.gz -T ^$FEEDER_CELL_GERMLINE $MOUSE2_calls
bcftools view -O z -o ${OUTPUT_FOLDER}/mouse3_calls_filtered.vcf.gz -T ^$FEEDER_CELL_GERMLINE $MOUSE3_calls

python ./vcf_to_ts_anno_SnpEff_plus_filters.py -o ${OUTPUT_FOLDER}/Variant_annotation_mouse2.txt $MOUSE2_calls
python ./vcf_to_ts_anno_SnpEff_plus_filters.py -o ${OUTPUT_FOLDER}/Variant_annotation_mouse3.txt $MOUSE3_calls


