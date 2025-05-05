# Function to convert Ensembl ID to Gene Symbol in a seurat object
convertEnsemblRownames <- function(seu) {
  # obtain the genes in the current object
  genes <- as.data.frame(rownames(seu))
  # convert ENSEMBL to SYMBOL
  master <- mapIds(org.Hs.eg.db, 
                   keys = genes$`rownames(seu)`, 
                   keytype = "ENSEMBL", 
                   column="SYMBOL")
  master <- as.data.frame(master)
  # match ENSEMBL to SYMBOL
  genes$geneID <- master$master
  # remove ENSEMBL genes that do not map to SYMBOLs
  seu <- seu[!is.na(genes$geneID),]
  genes <- genes[!is.na(genes$geneID),]
  # remove duplicate SYMBOLs
  dis <- duplicated(genes$geneID)
  seu <- seu[!dis,]
  genes <- genes[!dis,]
  # rename rows in appropriate slots (note: scale.data slot is empty since we have not done any scaling so do not rename them here)
  rownames(seu@assays$RNA@counts) <- genes$geneID
  rownames(seu@assays$RNA@data) <- genes$geneID
  return(seu)
}