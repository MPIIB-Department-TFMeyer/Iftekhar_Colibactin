MOUSE_calls=/data_genome1/MB208_WGS/analysis/Strelka_reseq_mouse_1_and_2/called_variants/SnpEff/Colibactin_Strelka_calls_WGS_reseq_1_2_merged_SnpEff.vcf.gz

FEEDER_CELL_GERMLINE=/data_genome2/projects/MB208_Colibactin/analysis/DNA/AdditionalData/3T3_Genome/Mm_filtered_SNV_InDel.vcf.gz

bcftools view -O z -o mouse_1_2_wgs_calls_filtered.vcf.gz -T ^$FEEDER_CELL_GERMLINE $MOUSE_calls

#cat /data_genome2/projects/MB208_Colibactin/analysis/DNA/Mouse2/Strelka/called_variants/SnpEff/Colibactin_Strelka_Exome_calls.txt | cut -f 1-6,12-15 | uniq > Variant_annotation_mouse2.txt
#cat /data_genome2/projects/MB208_Colibactin/analysis/DNA/Mouse3/Strelka/called_variants/SnpEff/Colibactin_Strelka_Exome_calls.txt | cut -f 1-6,12-15 | uniq > Variant_annotation_mouse3.txt
python /data_genome2/projects/MB208_Colibactin/src/vcf_to_ts_anno_SnpEff_plus_filters.py -o Variant_annotation_mouse1_2.txt ../Strelka_reseq_mouse_1_and_2/called_variants/SnpEff/Colibactin_Strelka_calls_WGS_reseq_1_2_merged_SnpEff.vcf.gz

