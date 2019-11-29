library(GenomicFeatures)
# define genomic bins
chrom_info = read.table("/data_genome2/References/MusMusculus/Sequences/Genome/GRCm38_Ensembl/GRCm38.genome",sep="\t",header=F, stringsAsFactors = F)

# Ens_gtf_file = "/data_genome1/References/Human/Annotations/hg19/Gencode/v19/gencode.v19.annotation_fixed_chroms.gtf"
# chrom_info_txdb = chrom_info
# colnames(chrom_info_txdb)=c("chrom","length")
# chrom_info_txdb$is_circular=ifelse(chrom_info_txdb$chrom=="MT",T,F)
# GencodeTxDB = makeTxDbFromGFF(Ens_gtf_file, format="gtf", chrominfo=chrom_info_txdb, dataSource="GENCODEv19 hg19", organism="Homo sapiens")
# saveDb(GencodeTxDB, file="/data_genome1/References/R_transformed/GENCODEv19_hg19_TxDB.db")

basic_gtf_file = "/data_genome2/References/MusMusculus/Annotation/GENCODE/M17/gencode.vM17.basic.annotation_NCBI_chrom_names.gtf"
chrom_info_txdb = chrom_info
colnames(chrom_info_txdb)=c("chrom","length")
chrom_info_txdb$is_circular=ifelse(chrom_info_txdb$chrom=="MT",T,F)
Gencode_basicTxDB = makeTxDbFromGFF(basic_gtf_file, format="gtf", chrominfo=chrom_info_txdb, dataSource="GENCODE basic M17", organism="Mus musculus")
saveDb(Gencode_basicTxDB, file="/data_genome1/References/R_transformed/GENCODE_basic_M17_TxDB.db")

#TxDb object:
  # Db type: TxDb
  # Supporting package: GenomicFeatures
  # Data source: GENCODE basic M17
  # Organism: Mus musculus
  # Taxonomy ID: 10090
  # miRBase build ID: NA
  # Genome: NA
  # transcript_nrow: 77835
  # exon_nrow: 320014
  # cds_nrow: 212838
  # Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2018-05-29 10:29:44 +0200 (Tue, 29 May 2018)
# GenomicFeatures version at creation time: 1.30.3
# RSQLite version at creation time: 2.0
# DBSCHEMAVERSION: 1.2