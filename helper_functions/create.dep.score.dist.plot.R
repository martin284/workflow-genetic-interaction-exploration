create.dep.score.dist.plot <- function(dep.scores, class.table, dataset) {
  
  dep.scores.merged <- merge(
    dep.scores,class.table[c("CCLE_ID","isDeleterious")],by="CCLE_ID")
  
  # split dataset in D and ND
  dep.scores.D <- dep.scores.merged[dep.scores.merged$isDeleterious==T,
                            c(-1,-ncol(dep.scores.merged))]
  dep.scores.ND <- dep.scores.merged[dep.scores.merged$isDeleterious==F,
                            c(-1,-ncol(dep.scores.merged))]
  
  # join all dependency scores of the two subgroups
  dep.scores.D.joined <- na.omit(data.frame(dep.score=unlist(dep.scores.D))) %>%
    mutate(Group = "Deleterious")
  dep.scores.ND.joined <- na.omit(data.frame(dep.score=unlist(dep.scores.ND))) %>%
    mutate(Group = "Non-deleterious")
  
  dep.scores.joined <- rbind(dep.scores.D.joined, dep.scores.ND.joined)
  
  # create plot
  plot <- ggplot(dep.scores.joined, aes(x = dep.score, color= Group)) +
    geom_density()  +
    ggtitle(paste0("Dependency Score Distribution in ",dataset)) +
    theme(legend.position = "bottom")
  
  return(plot)
}
