  rm(list=ls())
  library(circlize)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(data.table)
  setwd("/data_genome1/MB208_WGS/analysis/All_mutations_reseq_mouse_2_and_3/Circos/")
  
  #######################################################################################################
  # CNV data
  #######################################################################################################
  
  load("../../ControlFreeC_reseq_mouse_2_and_3/analysis/All_CNV_segments_paired_recalled.Rdata")
  
  # segment_length_threshold = 1e5
  p_value_threshold = 1e-2
  #num_of_targets_threshold = 2
  
  all_input_data = list()
  for (f in names(all_long_CN_calls)) {
    tmp = all_long_CN_calls[[f]]
    #tmp2 = subset(tmp, WilcoxonRankSumTestPvalue < p_value_threshold | KolmogorovSmirnovPvalue < p_value_threshold)
    tmp2 = subset(tmp, status != "nc")
    if(nrow(tmp2)<1) next
    tmp2$LogRatios = log2(tmp2$copy.number/2)
    tmp2$CBS.Mean = tmp2$LogRatios
    input_df = tmp2[,c("chr","start","end", "CBS.Mean", "LogRatios")]
    input_df$LogRatios = ifelse(is.infinite(input_df$LogRatios) | abs(input_df$LogRatios)>2, sign(input_df$LogRatios) * 2, input_df$LogRatios)
    input_df$chr = paste0("chr", input_df$chr)
    all_input_data[[f]] = input_df
  }
  
  #######################################################################################################
  # SNV data 
  #######################################################################################################
  
  snv_env = new.env()
  load("../Colibactin_WGS_Strelka_Calls_reseq_2_3.Rdata", envir=snv_env)
  snv_indel_variants = snv_env[["all_variants"]]
  snv_indel_variants$AF = snv_indel_variants$AC / snv_indel_variants$DP
  
  sample_count = table(names(snv_indel_variants))
  snv_indel_variants$sample_count = integer(length(snv_indel_variants))
  snv_indel_variants$sample_count = sample_count[names(snv_indel_variants)]
  
  snv_indel_variants = snv_indel_variants[snv_indel_variants$DP < 150 & snv_indel_variants$DP >= 10 & snv_indel_variants$AC >= 3 & snv_indel_variants$sample_count == 1 & !snv_indel_variants$multiallelic & snv_indel_variants$AF >= 0.25]
  
  #######################################################################################################
  # SV data 
  #######################################################################################################
  
  sv_env = new.env()
  sv_calls = load("../SV_Svaba/All_SV_sites.Rdata", envir = sv_env)
  sv_sites = get("all_sv_sites", sv_env)
  sv_sites = sv_sites[sv_sites$QUAL >= 10 & sv_sites$AD >= 5]

  sv_sites_dt = data.table(SV_ID = sv_sites$SV_ID, breapoint_ID = sv_sites$breakpoint_id, Sample = sv_sites$Sample, chr = as.character(seqnames(sv_sites)), pos = start(sv_sites), SPAN = sv_sites$SPAN)
  # fix ID for one SV with multimapping end for one BP
  sv_sites_dt[ , SV_ID:=ifelse(SV_ID=="c_6_70070001_70095001_21C", ifelse(SPAN==44754203, paste0(SV_ID, "a"), SV_ID), SV_ID) ]

  # order 
  setorder(sv_sites_dt, SV_ID, chr, pos)
  sv_sites_paired = sv_sites_dt[ ,.(chr1=chr[1], pos1 = pos[1], chr2=chr[2], pos2 = pos[2]) , by=c("SV_ID","Sample", "SPAN")]
  
  #######################################################################################################
  # Mutated genes with annotations
  #######################################################################################################
  
  load("../MutatedGenes/All_variants_combined_annotated_WGS1_2_mice_recalled_CNV.Rdata")

  M17db = loadDb("/data_genome1/References/R_transformed/GENCODE_basic_M17_TxDB.db")
  tx_anno = fread("/data_genome1/References/MusMusculus/Annotation/GENCODE/M17/gencode.vM17.transcript.anno.txt", sep="\t", header=T, stringsAsFactors = F)
  setkey(tx_anno, transcript)

  #######################################################################################################
  
  
  #######################################################################################################
  # adjust this to your needs
  
  circos_fun <- function(cnv_data, snv_data, sv_data, gene_pos, output_file_name, sample_names) {
    
    #par(mar = c(2, 2, 2, 2), lwd = 0.1, cex = 0.7)
    
    #circos.par("track.height" = 0.97, gap.degree = 0, canvas.xlim=c(-1.5,1.5), canvas.ylim=c(-1.5,1.5))
    
    cairo_pdf(file=paste(output_file_name,"_a.pdf", sep=""), width=6, height=6)
    
    selected_chroms = paste("chr", c(as.character(1:19),"X","Y"), sep="")
    
    track_order = length(sample_names):1
    names(track_order) = sample_names
    
    track_names = names(sample_names)
    names(track_names) = sample_names
    
    ordered_sample_names = names(track_order)[order(track_order)] # from outside to inside
    
    #dummy_data = data.frame(chr=selected_chroms, start=0, end=0, copy.number = 1, value1=0, stringsAsFactors = FALSE)
    
    circos.initializeWithIdeogram(species="mm10", chromosome.index = selected_chroms, track.height = 0.1, labels.cex=1.5, axis.labels.cex=0.5 )
    
    circos.par("track.height" = 0.15)
    circos.par("cell.padding" = c(0, 0, 0, 0))
    
    # all_input_data=list()
    # all_input_data[["TEST"]] = data.frame(chr=c("chr1","chr1"), start=c(1,100e6), end=c(1e6, 110e6), foo=c(3,-3), value=c(-2, 2), stringsAsFactors = FALSE)
    
    my_panel_fun = function(region, value, ...) {
      up = value[,2] > 0
      circos.genomicLines(region[up, ], value[up,], lwd=2, numeric.column = 2, 
                          type="segment", col="red")
      circos.genomicLines(region[!up, ], value[!up,], lwd=2, numeric.column = 2,
                          type="segment", col="blue")
      if (!all(up==FALSE)) {
        x0 = region[up,1]
        y0 = rep(0, length(x0))
        x1 = region[up,2]
        y1 = value[up,2]
        circos.rect(x0, y0, x1, y1, col="red")
      }
      if (!all(up==TRUE)) {
        x0 = region[!up,1]
        y0 = rep(0, length(x0))
        x1 = region[!up,2]
        y1 = value[!up,2]
        circos.rect(x0, y0, x1, y1, col="blue")
      }
      
    }
    
    #######################################################
    # Genes
    posTransform.fun = function(region) {
      return(region)
    }
    
    p_fun1 <- function(region, value, ...) {
      circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5), cex = 1, posTransform = posTransform.fun)
    }
    
    #circos.genomicTrackPlotRegion(subset(gene_pos, class!="P53"), ylim = c(0, 1), panel.fun = p_fun1, track.height = 0.05, bg.border = NA)
    circos.genomicPosTransformLines(gene_pos, posTransform = posTransform.fun, track.height = 0.04, col = ifelse(gene_pos$class=="P53", "red", ifelse(gene_pos$class=="PanCancer","black", "green")))
    
    #######################################################
    # SNV rainfall
    circos.genomicRainfall(snv_data, pch = 16, cex = 0.5, col = c("#FF000080"))
    
    for (n in ordered_sample_names) {
      sector_data = cnv_data[[n]]
      circos.genomicTrack(sector_data, numeric.column = 5, ylim=c(-2,2), stack=FALSE,
                          panel.fun = my_panel_fun)
    }
    
    # cols = rainbow(length(selected_chroms))
    # names(cols) = selected_chroms
    
    for (cc in selected_chroms) {
      for (tt in ordered_sample_names) {
        ti = track_order[tt] + 4 # first two tracks are karyogram + position + 1 track for genes + 1 Rainfall
        #circos.updatePlotRegion(cc, tt, bg.col = cols[g]) 
        sector_data = subset(cnv_data[[tt]], chr==cc)
        
        sector.xlim = get.cell.meta.data("xlim", sector.index = cc, track.index = ti)
        circos.lines(sector.xlim, c(0,0), lty=2, sector.index = cc, track.index = ti)
        
        if (cc == "chr4") 
          circos.yaxis("left", at=c(seq(-2, 2, by=1)), labels.cex = 0.5, 
                       tick.length=1, sector.index = cc, track.index = ti)
        min_y = 1.9
        if (cc == "chr6") {
          a = 2
          circos.text(mean(sector.xlim), min_y , track_names[tt], facing = "bending", niceFacing = T, 
                      cex = 0.7, adj=c(0.50,1.2), sector.index = cc, track.index = ti, col="darkgreen")
        }
      }
    }

    # take the first samples in the list (the innermost in circos plot) 
    sv_data_sel = sv_sites_paired[sv_sites_paired$Sample==sample_names[1]]
    
    sv_data_sel[, col := 1:.N, by="chr1"]
    max_col = max(sv_data_sel$col)
    colors = rainbow(max_col)
    
    for (i in 1:nrow(sv_data_sel)) {
      tmp = sv_data_sel[i]
      #circos.link(sector.index1 = paste0("chr",tmp$chr1), point1 = tmp$pos1, sector.index2=paste0("chr",tmp$chr2), point2 = tmp$pos2, col=colors[tmp$col])
      circos.link(sector.index1 = paste0("chr",tmp$chr1), point1 = tmp$pos1, sector.index2=paste0("chr",tmp$chr2), point2 = tmp$pos2, col="red")
    }
    
    circos.clear()
    
    dev.off()
    
  }
  
  #############################################################################################################
  
  all_samples = list()
  all_samples[["Mouse1"]] = c("WI"="MOUlyuRAAACAAA", "NI" = "MOUlyuRAAABAAA")
  all_samples[["Mouse2"]] = c("WI"="MOUlyuRAAAFAAA", "NI" = "MOUlyuRAAAEAAA")
  
  sample_to_SNV_calls = c("Mouse1"="MOUlyuRAAACAAA", "Mouse2"="MOUlyuRAAAFAAA")
  
  for (m in names(all_samples)) {
    # filter SNV sites
    snv_sites = reduce(snv_indel_variants[snv_indel_variants$sampleID==sample_to_SNV_calls[m]])
    suppressMessages(seqlevelsStyle(snv_sites) <- "UCSC")
    snv_sites_df = data.frame(chr=seqnames(snv_sites), start=start(snv_sites), end=end(snv_sites), stringsAsFactors = F)
    
    # filter mutated genes
    sel_genes_dt = all_variants_combined[mutated_only_WI & (KEGG_WNT | Reactome_WNT | pan_cancer_driver | KEGG_P53 ) & sampleID %in% all_samples[[m]]["WI"] & mutationType2 %in% c("CNV only, expressed","SNV only, MODERATE/HIGH Impact", "SNV/CNV, MODERATE/HIGH Impact") ]
    sel_genes_unique = sel_genes_dt[, .(GeneSymbol=unique(GeneSymbol), Wnt=ifelse(any(KEGG_WNT|Reactome_WNT), T, F), PanCancerDriver=ifelse(any(pan_cancer_driver), T, F), KEGG_P53=ifelse(any(KEGG_P53), T, F) ), by=geneid_fixed]
    #sel_genes_dt = all_variants_combined[mutated_only_WI & (KEGG_WNT | Reactome_WNT | pan_cancer_driver ) & sampleID %in% all_samples[[m]]["WI"] & mutationType2 %in% c("CNV only, expressed","SNV only, MODERATE/HIGH Impact", "SNV/CNV, MODERATE/HIGH Impact") ]
    #sel_genes_unique = sel_genes_dt[, .(GeneSymbol=unique(GeneSymbol), Wnt=ifelse(any(KEGG_WNT|Reactome_WNT), T, F), PanCancerDriver=ifelse(any(pan_cancer_driver), T, F) ), by=geneid_fixed]
    setkey(sel_genes_unique, "geneid_fixed")
    
    sel_genes = sel_genes_unique$geneid_fixed
    
    tmp = genes(M17db)
    tmp$geneid_fixed = sapply(strsplit(tmp$gene_id,"\\."), function(x) x[1])
    gene_pos = tmp[tmp$geneid_fixed %in% sel_genes]
    suppressMessages(seqlevelsStyle(gene_pos) <- "UCSC")
    mid_pos = (start(gene_pos) + end(gene_pos))/2
    gene_pos_df = data.frame(chr=seqnames(gene_pos), start=mid_pos, end=mid_pos, row.names=gene_pos$geneid_fixed )
    gene_pos_df$GeneSymbol = sel_genes_unique[rownames(gene_pos_df)]$GeneSymbol
    #gene_pos_df$PanCancerDriver = sel_genes_unique[rownames(gene_pos_df)]$PanCancerDriver
    gene_pos_df$class = with(gene_pos_df, ifelse(sel_genes_unique[rownames(gene_pos_df)]$KEGG_P53,"P53", ifelse(sel_genes_unique[rownames(gene_pos_df)]$PanCancerDriver, "PanCancer", "Wnt")))
    #gene_pos_df$class = with(gene_pos_df, ifelse(sel_genes_unique[rownames(gene_pos_df)]$PanCancerDriver, "PanCancer", "Wnt"))
    
    circos_fun(cnv_data = all_input_data, snv_data = snv_sites_df, sv_data = sv_sites_paired, gene_pos = gene_pos_df,  output_file_name = paste("Circos",m,format(Sys.time(), "%Y-%m-%d"),"v5", sep="_"), all_samples[[m]] )
    
  }
