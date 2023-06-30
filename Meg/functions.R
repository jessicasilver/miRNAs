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
