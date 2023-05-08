create.dep.score.boxplots <- function(dep.scores,target.gene,class.table,
                                      dataset) {
  dep.scores.target.gene <- dep.scores[c("CCLE_ID",target.gene)]
  dep.scores.target.gene <- merge(dep.scores.target.gene,
                                        class.table,by="CCLE_ID")
  dep.scores.target.gene <- na.omit(dep.scores.target.gene)
  colnames(dep.scores.target.gene)[2] <- "dep.score"
  
  boxplot.target.gene <-
    ggplot(dep.scores.target.gene, aes(x=isMutated,y=dep.score)) +
    geom_boxplot() +
    geom_jitter(shape=21,width=0.3,aes(fill=isDeleterious)) +
    stat_compare_means(method = "t.test") +
    ylab(paste(target.gene,"Dependency Score")) +
    scale_x_discrete(labels=c(paste0(genotype,"-non-mutated"),
                              paste0(genotype,"-mutated"))) +
    theme(legend.position="bottom") +
    ggtitle(dataset)
  
  return(boxplot.target.gene)
}
