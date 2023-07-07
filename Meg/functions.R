# store DESeq analysis (for all result names) as dataframes in a list
# dds = DESeq object; differentially expressed genes from the DESEq() function
resNames_list <- function(dds){
  resName <- as.list(resultsNames(dds)) 
  resName <- resName[-1]
  
  resNames_list <- list()
  for (i in 1:length(resName)) {
    
    signif_genes <- as.data.frame(results(dds, name = resName[[i]]))
    signif_genes <- signif_genes[order(signif_genes$padj, decreasing = FALSE),] %>%
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
