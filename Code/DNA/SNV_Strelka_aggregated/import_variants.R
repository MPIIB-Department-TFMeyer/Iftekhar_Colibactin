suppressMessages(library(VariantAnnotation))
suppressMessages(library(ggplot2))
library(data.table)

setwd("/data_genome2/projects/MB208_Colibactin/analysis/DNA/combined_results/Strelka_aggregated/")

all_samples_folder = "/data_genome2/projects/MB208_Colibactin/analysis/DNA/combined_results/Strelka_aggregated/"

all_vcf_raw = list()
# annotated variants, all samples
all_vcf_raw[["mouse2"]] = readVcf(paste(all_samples_folder, "mouse2_calls_filtered.vcf.gz", sep=""))
all_vcf_raw[["mouse3"]] = readVcf(paste(all_samples_folder, "mouse3_calls_filtered.vcf.gz", sep=""))

all_samples = list()
all_samples[["mouse2"]] = c(paste("2702_",c("WI"), sep=""))
all_samples[["mouse3"]] = c(paste("3466_",c("WI"), sep=""))

# Annotations from SnpEff are per matching transcript, i.e. there can be more than one annotations per line. 
# This must be processed appropriately, e.g. by the script vcf_to_ts_csv_SnpEff_plus_filters.py
#ann_header_info = as.data.frame(info(header(vcf_raw)))["ANN","Description"]
#aa = regexpr("'[^\\']*'",ann_header_info, perl=T)
#ann_cols = make.names(unlist(strsplit(substr(ann_header_info, aa+1, aa + attr(aa, "match.length")-2), " \\| ")))
#tmp = info(vcf_raw)[["ANN"]]

get_alt_allele_depth <- function(var_allele, a, c, g, t) {
  switch(var_allele, 
         A=a,
         C=c,
         G=g,
         T=t,
         NA
  )
}

DNAstringlist2char = function(x) sapply(lapply(x, as.character), paste, collapse = ",") 
DNAstringlist2maxlen = function(x) sapply(lapply(x, as.character), function(y) max(nchar(y))) 

first = 1
all_variants = c()
for (mouseid in names(all_vcf_raw)) {

  vcf_raw = all_vcf_raw[[mouseid]]
  
  dp = geno(vcf_raw)[["DP"]]
  tir = geno(vcf_raw)[["TIR"]][,,1]
  ac = list(A=geno(vcf_raw)[["AU"]][,,1], C=geno(vcf_raw)[["CU"]][,,1], G=geno(vcf_raw)[["GU"]][,,1], T=geno(vcf_raw)[["TU"]][,,1])
  

  maxlen_alt = DNAstringlist2maxlen(rowRanges(vcf_raw)$ALT)
  maxlen_ref = DNAstringlist2maxlen(rowRanges(vcf_raw)$REF)
  
  vartype = ifelse(maxlen_ref >1 | maxlen_alt > 1, "Indel", "SNV")
  
  for (s in all_samples[[mouseid]]) {
    r = rowRanges(vcf_raw)
    r$sampleID = s
    r$DP = dp[,s]

    r$multiallelic = unlist(lapply(r$ALT, length))>1
    r$alt_allele = unlist(lapply(r$ALT, paste, collapse=","))
    r$AC = mapply(get_alt_allele_depth, r$alt_allele, ac[["A"]][,s], ac[["C"]][,s], ac[["G"]][,s],ac[["T"]][,s])
    indel_ac = tir[,s]
    r$vartype = vartype
    r$AC = ifelse(vartype=="SNV", r$AC, indel_ac)
    
    # remove variants not detected in this sample
    r = r[!is.na(r$AC)]
    
    if (first) {
      all_variants = r
      first=0
    } else {
      all_variants = c(all_variants, r)
    }
  }
}

# Variant Annotation
# We take annotations already preprocessed by a custom python script
all_anno = list()
all_anno[["mouse2"]] = fread("./Variant_annotation_mouse2.txt", sep="\t", header=T, stringsAsFactors = F)
all_anno[["mouse3"]] = fread("./Variant_annotation_mouse3.txt", sep="\t", header=T, stringsAsFactors = F)

impact_order = c("HIGH"=4,"MODERATE"=3,"LOW"=2,"MODIFIER"=1)

all_anno_df = unique(do.call(rbind, all_anno))
all_anno_df[,IMPACT_ORDER := impact_order[Impact]]
all_anno_df[,MAX_IMPACT := max(IMPACT_ORDER), by=c("CHROM","POS","REF","ALT")]

anno_unique_sites = all_anno_df[,.(Gene=list(Gene), Impact=list(Impact), AnnotationType=list(AnnotationType), dbSNP=unique(dbSNP) ),by=c("CHROM","POS","REF","ALT")]

anno_unique_sites_highest_impact = all_anno_df[IMPACT_ORDER==MAX_IMPACT,][,.(Gene=list(Gene), ENSG=list(ENSG), HighestImpact=unique(Impact), AnnotationType=list(AnnotationType), dbSNP=unique(dbSNP) ),by=c("CHROM","POS","REF","ALT")]

save(all_variants, all_samples, anno_unique_sites, anno_unique_sites_highest_impact, file="Colibactin_Exome_Strelka_Calls_both_mice.Rdata")
