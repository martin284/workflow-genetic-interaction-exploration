create.dep.score.boxplots.MSI <- function(dep.scores,target.gene,class.table,
                                          dataset) {
  dep.scores.target.gene <- dep.scores[c("CCLE_ID",target.gene)]
  dep.scores.target.gene <- merge(dep.scores.target.gene,
                                  class.table,by="CCLE_ID")
  dep.scores.target.gene <- na.omit(dep.scores.target.gene)
  colnames(dep.scores.target.gene)[2] <- "dep.score"
  
  # order boxplots
  dep.scores.target.gene$CCLE_MSI <- factor(dep.scores.target.gene$CCLE_MSI,
                                            levels=c("MSS","MSI"))
  
  boxplot.target.gene <- 
    ggplot(dep.scores.target.gene, aes(x=CCLE_MSI,y=dep.score)) +
    geom_boxplot() +
    geom_jitter(shape=21,width=0.3,fill="darkgreen") +
    stat_compare_means(method = "t.test") +
    ylab(paste(target.gene,"Dependency Score")) +
    scale_x_discrete(labels=c("MSS","MSI")) +
    theme(legend.position="bottom") +
    ggtitle(dataset)
  
  return(boxplot.target.gene)
}
