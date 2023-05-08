create.mutation.counts.plot <- function(dataset,dep.scores,cancer.names,
                                        CCLE.mut) {
  
  depmapIDs <- merge(dep.scores['CCLE_ID'], cancer.names,by="CCLE_ID")
  CCLE.mut.subset <- subset(CCLE.mut[,c(1,3)],
                           (DepMap_ID %in% depmapIDs$DepMap_ID))[,1]
  mutation.counts <- as.data.frame(table(factor(CCLE.mut.subset)))
  gene.cut.off <- 150
  mutation.counts[mutation.counts$Freq>gene.cut.off,]$Freq <- 
    gene.cut.off
  colnames(mutation.counts) <- c("gene","Freq")
  
  mutation.plot <- ggplot(data=mutation.counts,aes(x=Freq)) +
    geom_histogram(binwidth = 5) +
    ylab("# genes") +
    xlab("# mutations") +
    ggtitle(paste("Mutation Counts in",dataset,"Cell Lines")) +
    scale_x_continuous(breaks = c(0,50,100,150), labels = c(0,50,100, ">=150"))
    
  return(mutation.plot)
  
}
