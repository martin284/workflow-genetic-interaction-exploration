create.dep.score.boxplots <- function(dep.scores,target.gene,class.table,
                                      dataset,pos.p.value=1.7,
                                      sample.name.1 = "-non-mutated",
                                      sample.name.2 = "-mutated") {
  dep.scores.target.gene <- dep.scores[c("CCLE_ID",target.gene)]
  dep.scores.target.gene <- merge(dep.scores.target.gene,
                                        class.table,by="CCLE_ID")
  dep.scores.target.gene <- na.omit(dep.scores.target.gene)
  colnames(dep.scores.target.gene)[2] <- "dep.score"
  
  boxplot.target.gene <-
    ggplot(dep.scores.target.gene, aes(x=isMutated,y=dep.score)) +
    geom_boxplot() +
    geom_jitter(shape=21,width=0.3,aes(fill=isMutated)) +
    stat_compare_means(method = "wilcox.test",label.x=pos.p.value) +
    ylab(paste(target.gene,"Dependency Score")) +
    scale_x_discrete(labels=c(paste0(target.gene,sample.name.1),
    paste0(target.gene,sample.name.2))) +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15),
          axis.text = element_text(size = 10)) +
    ggtitle(dataset)  +
    xlab("")
  
  return(boxplot.target.gene)
}
