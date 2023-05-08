plot.mean.diffs <- function(res) {
  
  mean.diffs.DRIVE <- res$mean_differences_DRIVE
  mean.diffs.Achilles <- res$mean_differences_Achilles
  
  plot.DRIVE <- ggplot(mean.diffs.DRIVE,aes(x=D_ND)) +
    geom_histogram(bins = 40,color="tomato3",fill="tomato3") +
    ggtitle("Distribution of Mean Differences in DRIVE") +
    xlab("Mean Difference") +
    geom_vline(aes(xintercept=0))
  
  plot.Achilles <- ggplot(mean.diffs.Achilles,aes(x=D_ND)) +
    geom_histogram(bins = 40,color="tomato3",fill="tomato3") +
    ggtitle("Distribution of Mean Differences in Achilles") +
    xlab("Mean Difference") +
    geom_vline(aes(xintercept=0))
  
  
  return(plots=list(DRIVE=plot.DRIVE,Achilles=plot.Achilles))
}
