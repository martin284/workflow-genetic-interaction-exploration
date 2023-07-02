create.mutation.counts.plot.mark <- function(dataset,
                                             dep.scores,
                                             cancer.names,
                                             CCLE.mut) {
  
  # keep cell lines that appear in the CCLE mutation dataset
  # and the gene-screening dataset
  depmapIDs <- merge(dep.scores['CCLE_ID'], cancer.names, by = "CCLE_ID")
  CCLE.mut.subset <- subset(CCLE.mut[ , c(1, 3)],
                            (DepMap_ID %in% depmapIDs$DepMap_ID))[ , 1]
  
  # count number of reported mutations per gene 
  mutation.counts <- as.data.frame(table(factor(CCLE.mut.subset)))
  
  # define cut off for better plotting
  gene.cut.off <- 150
  mutation.counts[mutation.counts$Freq>gene.cut.off,]$Freq <- 
    gene.cut.off
  colnames(mutation.counts) <- c("gene","Freq")
  
  interval.size <- 10
  genes.per.mut.interval <- c()
  gene.of.interest <- c()
  
  # count genes that have mutations within the selected intervals
  for (step in seq(0, gene.cut.off - interval.size, interval.size)) {
    
    count.table <- mutation.counts[
      mutation.counts$Freq > step &
        mutation.counts$Freq <= step + interval.size, ]
    
    counts <- nrow(count.table)
    genes.per.mut.interval <- append(genes.per.mut.interval,counts)
    
    if (gene1 %in% unique(count.table$gene) & 
        gene2 %in% unique(count.table$gene)) {
      gene.of.interest <- append(gene.of.interest, paste(gene1, "and", gene2))
    }
    else if (gene1 %in% unique(count.table$gene)) {
      gene.of.interest <- append(gene.of.interest, gene1)
    }
    else if (gene2 %in% unique(count.table$gene)) {
      gene.of.interest <- append(gene.of.interest, gene2)
    }
    else {
      gene.of.interest <- append(gene.of.interest, "No Gene of Interest")
    }
    
  }
  
  # combine lists to dataframe
  genes.per.mut.interval.df <- as.data.frame(cbind(
    gene.counts = genes.per.mut.interval, gene.of.interest = gene.of.interest))
  genes.per.mut.interval.df$gene.counts <- as.numeric(
    genes.per.mut.interval.df$gene.counts)
  
  # create histogram 
  plot <- ggplot(genes.per.mut.interval.df, 
         aes(x = seq(interval.size, gene.cut.off,interval.size),
             y = gene.counts,
             fill = gene.of.interest)) +
    geom_bar(stat = 'identity') +
    xlab("#mutations") +
    ylab("#genes") +
    scale_x_continuous(breaks = c(0,50,100,150),
                       labels = c(0, 50, 100, ">=150")) +
    scale_fill_discrete(name = "") +
    ggtitle(paste("Mutation Counts in", dataset, "Cell Lines")) +
    theme(legend.position = "bottom")
  
  return(plot)
}
