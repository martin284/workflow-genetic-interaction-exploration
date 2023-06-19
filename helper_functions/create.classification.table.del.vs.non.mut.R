create.classification.table.del.vs.non.mut <- function(CCLE_mut,cancer_names,
                                                       gene_ID,DRIVE,Achilles) {
  
  # keep only cell lines with mutations in the gene given by gene_ID
  CCLE_mut_gene_ID <- CCLE_mut[CCLE_mut$Hugo_Symbol==gene_ID,]
  
  # remove silent mutations
  CCLE_mut_gene_ID <- CCLE_mut_gene_ID[
    CCLE_mut_gene_ID$Variant_Classification!="Silent",]
  
  # remove mutation entries with more harmless mutation
  # in genes with multiple mutations
  CCLE_mut_gene_ID <- CCLE_mut_gene_ID[order(CCLE_mut_gene_ID$isDeleterious,
                                             decreasing=T),]
  CCLE_mut_gene_ID <- CCLE_mut_gene_ID[!duplicated(CCLE_mut_gene_ID$DepMap_ID),]
  
  # replace DepMap_IDs by CCLE names
  CCLE_mut_table <- merge(CCLE_mut_gene_ID, cancer_names, by="DepMap_ID",
                          all=TRUE)
  CCLE_mut_table$DepMap_ID <- NULL
  CCLE_mut_table <- 
    CCLE_mut_table[,c("CCLE_ID","isDeleterious")]
  
  # remove cell lines with non-deleterious mutation
  CCLE_mut_table$keep <- TRUE
  CCLE_mut_table$keep <- ifelse(
    CCLE_mut_table$isDeleterious=="True" |
      is.na(CCLE_mut_table$isDeleterious),
      T,F)
  CCLE_mut_table <- CCLE_mut_table[CCLE_mut_table$keep,]
  CCLE_mut_table$keep <- NULL

  # classify all cell lines without any mutations in FUBP1 as non-deleterious
  # and change char True/False to boolean TRUE/FALSE
  CCLE_mut_table$isDeleterious <- ifelse(is.na(CCLE_mut_table$isDeleterious),
                                        F,T)
  
  # keep only cell lines for which dependency scores are available
  CCLE_mut_table_DRIVE <- merge(CCLE_mut_table,DRIVE['CCLE_ID'],by="CCLE_ID")
  CCLE_mut_table_Achilles <- merge(CCLE_mut_table, Achilles['CCLE_ID'],
                                   by="CCLE_ID")
  
  return(list(DRIVE=CCLE_mut_table_DRIVE, Achilles=CCLE_mut_table_Achilles)) 
}
