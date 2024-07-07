# store DESeq analysis (for all result names) as dataframes in a list
# dds = DESeq object; differentially expressed genes from the DESEq() function
resNames_list <- function(dds){
  resName <- as.list(resultsNames(dds)) 
  resName <- resName[-1]
  
  resNames_list <- list()
  for (i in 1:length(resName)) {
    
    signif_genes <- as.data.frame(results(dds, name = resName[[i]]))
    signif_genes <- signif_genes[order(signif_genes$pvalue, decreasing = FALSE),] %>%
      as.data.frame(.)
    
    resNames_list[[paste0((resName)[[i]])]] <- signif_genes
  }
  return(resNames_list)
}

# filters top genes from a DESeq analysis using baseMean and log2FC thresholds
# dfList = list; a list of DESEq2 results converted into a dataframe 
# pvalue, L2FC = integer; numerical thresholds for pvalue and log2FoldChange, respectively

threshFilter <- function(dfList, padjust, L2FC){
  filt <- list()
  
  for (i in 1:length(dfList)) {
    
    signif_genes <- dfList[[i]] %>% filter(dfList[[i]]$padj < padjust) %>%
      filter(abs(.$log2FoldChange) > L2FC)
    
    signif_genes <- signif_genes[order(signif_genes$padj, decreasing = F),] %>%
      as.data.frame(.)
    
    filt[[paste0(names(dfList)[i])]] <- signif_genes
  }
  return(filt)
}

# adds a column to a DESeq2 table and tags genes as UP and DOWN regulated 
# according to baseMean and log2FC thresholds
# dfList = list; a list of DESEq2 results converted into a dataframe 
# pvalue, L2FC = integer; numerical thresholds for pvalue and log2FoldChange, respectively
diff_expressed <- function(dfList, padjust, L2FC){
  
  results <- list()
  
  for (i in 1:length(dfList)) {
    # categorise to upregulated and downregulated genes
    y <- dfList[[i]] %>% 
      mutate(diffexpressed = case_when(
        log2FoldChange > L2FC & padj < padjust ~ 'UP',
        log2FoldChange < (L2FC*-1) & padj < padjust ~ 'DOWN',
        padj > padjust ~ 'NO',
        .default = 'NO')) # none of the cases match 
    
    results[[paste0(names(dfList[i]))]] <- y
  }
  
  return(results)
}


# creates a clusterProfiler-ready input for ORA or GSEA from DESeq
# NOTE: this assumes that the dataframe has a column named ENSEMBL containing
#       ensembl ids
#df = dataframe; DESeq2 results
#gene_name = name of the column that contains the gene names 
#padj = numeric; cutoff FDR or padj value (default = 0.05)
#lf2c = numeric; cutoff log2FoldChange (default = 2)

# preVs3hrPost_stat <- preVs3hrPost %>% 
#   dplyr::select(ENSEMBL, stat) %>% 
#   na.omit() %>% 
#   distinct() %>% 
#   group_by(ENSEMBL) %>% 
#   summarize(stat=mean(stat))
# 
# ranks_preVs3hrPost <- deframe(signif_preVs3hrPost_stat)

cProfilr_input <- function(df){
  input <- list()
  
  # we want the log2 fold change 
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- df$ENSEMBL
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # Exctract significant results (padj < 0.05)
  sig_genes_df = subset(df, padj < 0.05)
  
  # From significant results, we want to filter on log2fold change
  genes <- sig_genes_df$log2FoldChange
  
  # Name the vector
  names(genes) <- sig_genes_df$ENSEMBL
  
  # omit NA values
  genes <- na.omit(genes)
  
  # filter on min log2fold change (log2FoldChange > 2)
  genes <- names(genes)[abs(genes) > 1]
  
  input[["signif_genes"]] <- genes
  input[["gene_list"]] <- gene_list
  
  return(input)
}


# set ScoreType for fgsea ("std" if neg and pos values present; "pos" if only pos; "neg" if only neg)
score_test <- function(gsea_data_vec){
  if(min(gsea_data_vec) < 0 & max(gsea_data_vec) > 0){
    type <- "std"
  } else if(min(gsea_data_vec) < 0 & max(gsea_data_vec) < 0){
    type <- "neg"
  } else{
    type <- "pos"
  }
  return(type)
}

# GSEA function 
# source: https://bioinformaticsbreakdown.com/how-to-gsea/
GSEA = function(gene_list, pathway, pval, condition_name) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  fgRes <- fgsea::fgsea(pathways = pathway,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        nPermSimple=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES)) %>%
    mutate(pathway = gsub("GOBP_","", pathway),
           pathway = gsub("_"," ", pathway))
  
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  # message("Collapsing Pathways -----")
  # concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
  #                                     pathways = pathway,
  #                                     stats = gene_list)
  # 
  # fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  # message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  path_info = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  fig_name = paste0(condition_name,":")
  
  colors = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  theme_set(theme_classic())
  #My_Theme = theme(title=element_text(size=15, face='bold'))
  
  g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colors ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y=paste0("Normalized Enrichment Score (NES): ", path_info),
         title = fig_name,
         subtitle = "GO Biological Pathways NES from GSEA") + 
    theme(plot.title=element_text(hjust=0.5),
          plot.subtitle=element_text(hjust=0.5))
  
  g2 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=Enrichment)) +
    coord_flip() +
    labs(x="Pathway", y=paste0("Normalized Enrichment Score (NES): ", path_info),
         title = fig_name,
         subtitle = "GO Biological Pathways NES from GSEA") + 
    theme(plot.title=element_text(hjust=0.5),
          plot.subtitle=element_text(hjust=0.5))
  
  output = list("Results" = fgRes, "Plot1" = g1, "Plot2" = g2)
  return(output)
}

# calculates the z-score of normalised DESeq data
# norm_data = normalised DESeq data usually using vst()
# signif_genes_filt = dataframe; prefiltered genes according to thresholds
z_score <- function(norm_data, signif_genes_filt){
  count_mat <- norm_data[rownames(signif_genes_filt),]
  
  count_mat_z <- t(apply(count_mat, 1, scale))
  colnames(count_mat_z) <- colnames(count_mat)
  
  return(count_mat_z)
}

# creates a complex heatmap
# numkeep = integer; number of rows (genes) to select for the heatmap
# signif_genes_filt = dataframe; prefiltered genes according to thresholds 
# count_mat_z = matrix; z-score of the normalised DESeq data
hmap <- function(numkeep, signif_genes_filt, count_mat_z, hmap_name){
  
  numkeep <- numkeep
  rowskeep <- c(seq(1:numkeep))
  
  # log2foldchange values
  l2_val <- as.matrix(signif_genes_filt[rowskeep,]$log2FoldChange)
  colnames(l2_val) <- "logFC"
  
  # colours for values between b/w/r for min/max of l2 values
  col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c("blue", "white", "red"))
  
  
  # plot heatmaps
  ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                                 height = unit(2, "cm")))
  
  h1 <- Heatmap(count_mat_z[rowskeep,], column_title = hmap_name,
                column_title_gp = gpar(fontsize = 23),
                cluster_rows = F, # cluster only columns to compare samples
                column_labels = colnames(count_mat_z), name="Z-score",
                cluster_columns = T)
  
  h2 <- Heatmap(l2_val, row_labels = rownames(signif_genes_filt)[rowskeep], 
                cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
                #cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                #grid.text(round(l2_val[i, j],2), x, y)}
  )
  
  h<-h1+h2
  
  return(h)
}


# creates a complex heatmap with sample grouping
# numkeep = integer; number of rows (genes) to select for the heatmap
# signif_genes_filt = dataframe; prefiltered genes according to thresholds 
# count_mat_z = matrix; z-score of the normalised DESeq data
hmap_colGrouped <- function(numkeep, signif_genes_filt, count_mat_z, 
                            hmap_name, sampling_group, sampling_group_col){
  
  # create groupings
  sampling_group <- sampling_group
  sampling_group_col <- sampling_group_col
  
  dend1 = cluster_within_group(count_mat_z, sampling_group)
  
  numkeep <- numkeep
  rowskeep <- c(seq(1:numkeep))
  
  # log2foldchange values
  l2_val <- as.matrix(signif_genes_filt[rowskeep,]$log2FoldChange)
  colnames(l2_val) <- "logFC"
  
  # colours for values between b/w/r for min/max of l2 values
  col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c("blue", "white", "red"))
  
  
  # plot heatmaps
  
  h1 <- Heatmap(count_mat_z[rowskeep,], 
                column_title_gp = gpar(fontsize = 23),
                cluster_rows = F, # cluster only columns to compare samples
                column_labels = colnames(count_mat_z), name="Z-score",
                cluster_columns = dend1, column_split = 3,
                column_names_rot = 45,
                top_annotation = HeatmapAnnotation(Sampling = sampling_group, 
                                                   show_annotation_name = T,
                                                  col = list(Sampling = sampling_group_col)))
  
  h2 <- Heatmap(l2_val, row_labels = rownames(signif_genes_filt)[rowskeep], 
                cluster_rows = F, name="logFC", col = col_logFC,
                #cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                #grid.text(round(l2_val[i, j],2), x, y)}
  )
  
  h<-h1+h2
  
  return(h)
}

# saves a pheatmap file
# from https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
# where x is the pheatmap object
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
