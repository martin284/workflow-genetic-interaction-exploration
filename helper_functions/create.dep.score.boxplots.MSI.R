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
    geom_jitter(shape=21,width=0.3,aes(fill=CCLE_MSI)) +
    stat_compare_means(method = "wilcox.test",label.x=1.7) +
    ylab(paste(target.gene,"Dependency Score")) +
    scale_x_discrete(labels=c("MSS","MSI")) +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15),
          axis.text = element_text(size = 10)) +
    ggtitle(dataset) +
    xlab("")
  
  return(boxplot.target.gene)
}
