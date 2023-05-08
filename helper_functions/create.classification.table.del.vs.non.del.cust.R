create_classification_table <- function(CCLE_mut, cancer_names, gene_ID, DRIVE,
                                        Achilles) {
  # keep only relevant columns
  CCLE_mut <- CCLE_mut[,c(1,8,16)]
  cancer_names <- cancer_names[,c(1,2)]
  
  # create data structure for all cell lines and corresponding mutation type
  cell_lines_depmap_ID <- CCLE_mut[,3]
  cell_lines_depmap_ID <- cell_lines_depmap_ID[!duplicated(
    cell_lines_depmap_ID)]
  cell_lines_depmap_ID <- data.frame(cbind(DepMap_ID=cell_lines_depmap_ID,
                                           mutation_type=NA))
  
  # keep only cell lines with mutations in the gene given by gene_ID
  CCLE_mut_gene_ID <- CCLE_mut[CCLE_mut$Hugo_Symbol==gene_ID,][,c(1,2,3)]
  
  # change ins and del to indel as the difference is not important here
  for (i in seq(nrow(CCLE_mut_gene_ID))) {
    if(CCLE_mut_gene_ID[i,]$Variant_Classification=="In_Frame_Ins" |
       CCLE_mut_gene_ID[i,]$Variant_Classification=="In_Frame_Del") {
      CCLE_mut_gene_ID[i,]$Variant_Classification <- "In_Frame_InDel"
    }
    else if (CCLE_mut_gene_ID[i,]$Variant_Classification=="Frame_Shift_Ins" |
             CCLE_mut_gene_ID[i,]$Variant_Classification=="Frame_Shift_Del") {
      CCLE_mut_gene_ID[i,]$Variant_Classification <- "Frame_Shift_InDel"
    }
  }
  
  # factorize mutation types
  mut_levels <- c("Non_Mutated","Silent","Missense_Mutation","In_Frame_InDel",
                  "Nonsense_Mutation","Frame_Shift_InDel","Splice_Site")
  CCLE_mut_gene_ID$Variant_Classification <- factor(
    CCLE_mut_gene_ID$Variant_Classification, ordered = T, levels = mut_levels)
  
  CCLE_mut_gene_ID$Variant_Classification <- 
    sort(CCLE_mut_gene_ID$Variant_Classification,decreasing = T)
  
  # remove mutation entries with more harmless mutation
  # in genes with multiple mutations
  CCLE_mut_gene_ID <- CCLE_mut_gene_ID[!duplicated(CCLE_mut_gene_ID$DepMap_ID),]

  # replace DepMap_IDs by CCLE names
  CCLE_mut_table <- merge(CCLE_mut_gene_ID, cancer_names, by="DepMap_ID",
                          all=TRUE)
  CCLE_mut_table$DepMap_ID <- NULL
  CCLE_mut_table <- CCLE_mut_table[,c("CCLE_ID","Variant_Classification")]
  colnames(CCLE_mut_table) <- c("CCLE_ID","mutation_type")
  
  # decide which mutation types are deleterious
  # default: Nonsense_Mutation, Splice_Site, Frame_Shift_Del, Frame_Shift_Ins
  deleterious_mutations <- c("Nonsense_Mutation", "Splice_Site",
                             "Frame_Shift_InDel","In_Frame_InDel",
                             "Missense_Mutation")
  CCLE_mut_table$isDeleterious <- ifelse(
    CCLE_mut_table$mutation_type %in% deleterious_mutations, TRUE, FALSE)
  CCLE_mut_table$mutation_type[is.na(CCLE_mut_table$mutation_type)] <- 
    "Non_Mutated"
  
  # keep only cell lines for which dependency scores are available
  CCLE_mut_table_DRIVE <- merge(CCLE_mut_table, DRIVE, by="CCLE_ID")
  CCLE_mut_table_DRIVE <- CCLE_mut_table_DRIVE[,c(1,2,3)]
  colnames(CCLE_mut_table_DRIVE) <- c("CCLE_ID",paste0("mut_type_",gene_ID),
                                      paste0("isDeleterious_",gene_ID))
  
  
  CCLE_mut_table_Achilles <- merge(CCLE_mut_table, Achilles, by="CCLE_ID")
  CCLE_mut_table_Achilles <- CCLE_mut_table_Achilles[,c(1,2,3)]
  colnames(CCLE_mut_table_Achilles) <- c("CCLE_ID",paste0("mut_type_",gene_ID),
                                      paste0("isDeleterious_",gene_ID))
  
  return(list(DRIVE=CCLE_mut_table_DRIVE, Achilles=CCLE_mut_table_Achilles)) 
}
