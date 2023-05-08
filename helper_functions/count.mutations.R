count.mutations <- function(dep.scores.DRIVE,dep.scores.Achilles,cancer.names,
                            CCLE.mut,gene1,gene2) {
  
  # DRIVE
  depmapIDs.DRIVE <- merge(dep.scores.DRIVE['CCLE_ID'], cancer.names,
                           by="CCLE_ID")
  CCLE.mut.DRIVE <- subset(CCLE.mut[,c(1,3)],
                           (DepMap_ID %in% depmapIDs.DRIVE$DepMap_ID))[,1]
  mutation.counts.DRIVE <- as.data.frame(table(factor(CCLE.mut.DRIVE)))
  colnames(mutation.counts.DRIVE) <- c("gene","Freq")
  
  mut.gene1.DRIVE <- mutation.counts.DRIVE[mutation.counts.DRIVE$gene==gene1,]
  mut.gene2.DRIVE <- mutation.counts.DRIVE[mutation.counts.DRIVE$gene==gene2,]
  
  # Achilles
  depmapIDs.Achilles <- merge(dep.scores.Achilles['CCLE_ID'], cancer.names,
                           by="CCLE_ID")
  CCLE.mut.Achilles<- subset(CCLE.mut[,c(1,3)],
                           (DepMap_ID %in% depmapIDs.Achilles$DepMap_ID))[,1]
  mutation.counts.Achilles <- as.data.frame(table(factor(CCLE.mut.Achilles)))
  colnames(mutation.counts.Achilles) <- c("gene","Freq")
  
  mut.gene1.Achilles <- mutation.counts.Achilles[
    mutation.counts.Achilles$gene==gene1,]
  mut.gene2.Achilles <- mutation.counts.Achilles[
    mutation.counts.Achilles$gene==gene2,]
  
  mutation.table <- as.data.frame(cbind(
    DRIVE=c(mut.gene1.DRIVE$Freq,mut.gene2.DRIVE$Freq),
    Achilles=c(mut.gene1.Achilles$Freq,mut.gene2.Achilles$Freq)))
  
  rownames(mutation.table) <- c(gene1,gene2)
  
  return(mutation.table)
  
}
