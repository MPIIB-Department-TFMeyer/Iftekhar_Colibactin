  rm(list=ls())
  library(circlize)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(data.table)
  
  #######################################################################################################
  # CNV data
  #######################################################################################################
  
  load("../../../Data/Processed/DNA_Variants_Contra/Long_CNV_recalled.Rdata")
  
  # segment_length_threshold = 1e5
  # p_value_threshold = 1e-2
  num_of_targets_threshold = 2
  
  all_input_data = list()
  for (f in names(all_long_CN_calls)) {
    tmp = all_long_CN_calls[[f]]
    tmp2 = subset(tmp,num.mark > num_of_targets_threshold)
    colnames(tmp2) = c("ID","chr","start","end", "num.mark","LogRatios","sample")
    tmp2$CBS.Mean = tmp2$LogRatios
    input_df = tmp2[,c("chr","start","end", "CBS.Mean", "LogRatios")]
    input_df$LogRatios = ifelse(is.infinite(input_df$LogRatios) | abs(input_df$LogRatios)>2, sign(input_df$LogRatios) * 2, input_df$LogRatios)
    all_input_data[[f]] = input_df
  }
  
  #######################################################################################################
  # SNV data 
  #######################################################################################################
  
  snv_env = new.env()
  load("../SNV_Strelka_aggregated/Colibactin_Exome_Strelka_Calls_both_mice.Rdata", envir=snv_env)
  snv_indel_variants = snv_env[["all_variants"]]
  snv_indel_variants = snv_indel_variants[snv_indel_variants$DP < 150]
  
  
  #######################################################################################################
  # Mutated genes with annotations
  #######################################################################################################
  
  load("../../../Results/DNASeq/All_variants_combined_annotated.Rdata")
  
  M17db = loadDb("../../../Data/External/GencodeM17/GENCODE_basic_M17_TxDB.db")
  tx_anno = fread("../../../Data/External/GencodeM17/gencode.vM17.transcript.anno.txt", sep="\t", header=T, stringsAsFactors = F)
  setkey(tx_anno, transcript)
  
  #######################################################################################################
  
  
  #######################################################################################################
  # adjust this to your needs
  
  circos_fun <- function(cnv_data, snv_data, gene_pos, output_file_name, sample_names) {
    
    #par(mar = c(2, 2, 2, 2), lwd = 0.1, cex = 0.7)
    
    #circos.par("track.height" = 0.97, gap.degree = 0, canvas.xlim=c(-1.5,1.5), canvas.ylim=c(-1.5,1.5))
    
    cairo_pdf(file=paste(output_file_name,".pdf", sep=""), width=20, height=20)
    
    selected_chroms = paste("chr", c(as.character(1:19),"X","Y"), sep="")
    
    track_order = length(sample_names):1
    names(track_order) = sample_names
    
    track_names = LETTERS[1:length(sample_names)]
    names(track_names) = sample_names
    
    ordered_sample_names = names(track_order)[order(track_order)] # from outside to inside
    
    #dummy_data = data.frame(chr=selected_chroms, start=0, end=0, copy.number = 1, value1=0, stringsAsFactors = FALSE)
    
    circos.initializeWithIdeogram(species="mm10", chromosome.index = selected_chroms )
    
    circos.par("track.height" = 0.08)
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
    
    circos.genomicTrackPlotRegion(subset(gene_pos, class=="P53"), ylim = c(0, 1), panel.fun = function(region, value, ...) {
      circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5), cex = 1, posTransform = posTransform.fun)
    }, track.height = 0.05, bg.border = NA)
    
    circos.genomicPosTransformLines(gene_pos, posTransform = posTransform.fun, track.height = 0.04, col = ifelse(gene_pos$class=="P53", "red", ifelse(gene_pos$class=="PanCancer","black", "green")))
    
    #######################################################
    # SNV rainfall
    circos.genomicRainfall(snv_data, pch = 16, cex = 2, col = c("#FF000080"))
    
    for (n in ordered_sample_names) {
      sector_data = cnv_data[[n]]
      circos.genomicTrack(sector_data, numeric.column = 5, ylim=c(-2,2), stack=FALSE,
                          panel.fun = my_panel_fun)
    }
    
    # cols = rainbow(length(selected_chroms))
    # names(cols) = selected_chroms
    
    for (cc in selected_chroms) {
      for (tt in ordered_sample_names) {
        ti = track_order[tt] + 5 # first two tracks are karyogram + position + 2 tracks for genes + 1 Rainfall
        #circos.updatePlotRegion(cc, tt, bg.col = cols[g]) 
        sector_data = subset(cnv_data[[tt]], chr==cc)
        
        sector.xlim = get.cell.meta.data("xlim", sector.index = cc, track.index = ti)
        circos.lines(sector.xlim, c(0,0), lty=2, sector.index = cc, track.index = ti)
        
        if (cc == "chr4") 
          circos.yaxis("left", at=c(seq(-2, 2, by=1)), labels.cex = 0.7, 
                       tick.length=1, sector.index = cc, track.index = ti)
        min_y = 2
        if (cc == "chr6") {
          a = 1
          circos.text(mean(sector.xlim), min_y + 0.1, track_names[tt], facing = "bending", niceFacing = T, 
                      cex = 1.3, adj=c(0.5,1), sector.index = cc, track.index = ti, col="red")
        }
      }
    }
    
    circos.clear()
    
    dev.off()
    
  }
  
  #############################################################################################################
  
  output_folder = "../../../Results/DNASeq"
  
  all_samples = list()
  all_samples[["Mouse2"]] = c(paste("2702_",c("B","C","D","E","F","G","H"), sep=""))
  all_samples[["Mouse3"]] = c(paste("3466_",c("B","C","D","E","F","G","H"), sep=""))
  
  sample_to_SNV_calls = c("Mouse2"="2702_WI", "Mouse3"="3466_WI")

  #WI_samples = c(paste0("2702_",c("D","E","F","G","H")),paste0("3466_",c("D","E","F","G","H")) )
    
  #snv_indel_variants$sampleID2 = paste0(substr(snv_indel_variants$sampleID, 1, 5),ifelse(snv_indel_variants$sampleID %in% WI_samples, "WI","NI"))
  
  for (m in names(all_samples)) {
    # filter SNV sites
    snv_sites = reduce(snv_indel_variants[snv_indel_variants$sampleID==sample_to_SNV_calls[m]])
    suppressMessages(seqlevelsStyle(snv_sites) <- "UCSC")
    snv_sites_df = data.frame(chr=seqnames(snv_sites), start=start(snv_sites), end=end(snv_sites), stringsAsFactors = F)
    
    # filter mutated genes
    sel_genes_dt = all_variants_combined[mutated_only_WI & (KEGG_WNT | Reactome_WNT | pan_cancer_driver | KEGG_P53 ) & sampleID %in% all_samples[[m]] ]
    sel_genes_unique = sel_genes_dt[, .(GeneSymbol=unique(GeneSymbol), Wnt=ifelse(any(KEGG_WNT|Reactome_WNT), T, F), PanCancerDriver=ifelse(any(pan_cancer_driver), T, F), KEGG_P53=ifelse(any(KEGG_P53), T, F) ), by=geneid_fixed]
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
    
    circos_fun(cnv_data = all_input_data, snv_data = snv_sites_df, gene_pos = gene_pos_df,  output_file_name = file.path(result_folder, paste("Circos",m,format(Sys.time(), "%Y-%m-%d"),"v5", sep="_")), all_samples[[m]] )
    
  }
