library(data.table)
library(VariantAnnotation)
setwd("/data_genome2/projects/MB208_Colibactin/MB208_WGS/analysis/All_mutations_reseq_mouse_2_and_3/Signatures_indel/")
load("../Colibactin_WGS_Strelka_Calls_reseq_2_3.Rdata")

# SIMPLE file format for indel annotation (Alexandrov et al, 2018, biorXiv)
#ACUTE MYELOID LEUKEMIA  TCGA-AB-2803-03B-01W-0728-08    TCGA-LAML       GRCh37  SNV     22      35689684        35689684        T       C       5       SOMATICSNIPER|RADIA|MUTECT|MUSE|VARSCANS

sample_count = table(names(all_variants))
all_variants$sample_count = integer(length(all_variants))
all_variants$sample_count = sample_count[names(all_variants)]
all_variants$af = all_variants$AC / all_variants$DP

d_indel = all_variants[all_variants$vartype=="Indel" & all_variants$DP < 150 & all_variants$DP >= 10 & all_variants$AC >= 3 & all_variants$sample_count == 1 & !all_variants$multiallelic & all_variants$af >= 0.25]

if(any(elementNROWS(d_indel$ALT)>1)) stop("More than one alternative allele in at least one entry")

# Ins: Ref: X  Alt: X<INS>  POS: Pos of X (i.e. insertion after X)
# Del: Ref: XYZ Alt: X POS: Pos of X to Pos of Z

d_indel$ref_len = nchar(d_indel$REF)
d_indel$alt_seq = sapply(lapply(d_indel$ALT, `[`,1), as.character)
d_indel$alt_len = nchar(d_indel$alt_seq)
d_indel$ref_seq = as.character(d_indel$REF)
d_indel$mut_type = ifelse(d_indel$ref_len < d_indel$alt_len, "INS","DEL")
d_indel$start_fixed = ifelse(d_indel$mut_type == "INS", start(d_indel), start(d_indel)+1) 
d_indel$end_fixed = ifelse(d_indel$mut_type == "INS", end(d_indel), end(d_indel))

d_indel_simple = data.table(disease = rep("Colibactin", length(d_indel)), 
                                      sample_id = d_indel$sampleID, 
                                      project = rep("WGS", length(d_indel)), 
                                      genome = rep("GRCm38", length(d_indel)),
                                      mut_type = d_indel$mut_type,
                                      chrom = as.character(seqnames(d_indel)),
                                      start = d_indel$start_fixed,
                                      end = d_indel$end_fixed,
                                      ref = ifelse(d_indel$mut_type=="INS", "-", substr(d_indel$ref_seq, 2, d_indel$ref_len)),
                                      alt = ifelse(d_indel$mut_type=="INS", substr(d_indel$alt_seq, 2, d_indel$alt_len), "-"),
                                      annotation_cnt = rep(0,length(d_indel)),
                                      annotations = rep("", length(d_indel)),
                                      stringsAsFactors = F)

fwrite(d_indel_simple, file="MB208_WGS_indels_clonal.simple", sep="\t")

system("/data_genome1/SharedSoftware/MutSigAlexandrov/make_spectra_indels.py --genome --style dukenus --fasta /data_genome1/References/MusMusculus/Sequences/Genome/GRCm38_Ensembl/GRCm38.fasta ./MB208_WGS_indels_clonal.simple")
system("mv indel_counts.csv indel_counts_clonal.csv")
file.copy("./cache/catalogs/annotated/MB208_WGS_indels_clonal.simple-indels.gz", "./MB208_WGS_indels_clonal.simple-indels.gz", overwrite = T)

annotated_indels = fread("gzip -dc ./MB208_WGS_indels_clonal.simple-indels.gz", header=F, sep="\t")
setnames(annotated_indels, colnames(annotated_indels), c("sample","genome","chr","pos","ref","alt","indel_class","tx","tx_strand","gene","var_1") )

barplot(table(annotated_indels$indel_class), las=2)

# NOTE: indels are not guaranteed to be unique sites
save(annotated_indels, file="MB208_WGS_indels_clonal_annotated.Rdata")

dd = fread("indel_counts_clonal.csv", sep="\t")
mm = t(as.matrix(dd[, grepl("Colibactin",colnames(dd)), with=F]))
labels = apply(dd[,1:4, with=F],1, function(x) paste(x, collapse="_", sep=""))
rownames(mm) = labels
barplot(mm, beside = T, names.arg=labels, las=2, cex.names = 0.5)
mm_prop = sweep(mm, 1, apply(mm, 1, sum), "/")
barplot(mm_prop, beside = T, names.arg=labels, las=2, cex.names = 0.5)
