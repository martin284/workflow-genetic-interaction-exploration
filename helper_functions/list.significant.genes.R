list.significant.genes <- function(significant.genes){
  
  # create Boolean matrix for gene in sets
  sign.genes.total <- unique(unlist(significant.genes))
  genes.set.table <- as.data.frame(cbind(geneID=sign.genes.total))
  
  check.dataset <- function(gene.name,dataset) {
    return(gene.name %in% dataset)
  }
  
  if(length(significant.genes$welch.DRIVE)>0) {
    genes.set.table$welch.DRIVE <- 
      sapply(genes.set.table$geneID, check.dataset,significant.genes$welch.DRIVE)}
  if(length(significant.genes$wilc.DRIVE)>0) {
    genes.set.table$wilc.DRIVE <- 
      sapply(genes.set.table$geneID, check.dataset,significant.genes$wilc.DRIVE)}
  if(length(significant.genes$bayes.DRIVE)>0) {
    genes.set.table$bayes.DRIVE <- 
      sapply(genes.set.table$geneID, check.dataset,significant.genes$bayes.DRIVE)}
  if(length(significant.genes$welch.Achilles)>0) {
    genes.set.table$welch.Achilles <- 
      sapply(genes.set.table$geneID, check.dataset,significant.genes$welch.Achilles)}
  if(length(significant.genes$wilc.Achilles)>0) {
    genes.set.table$wilc.Achilles <- 
      sapply(genes.set.table$geneID, check.dataset,significant.genes$wilc.Achilles)}
  if(length(significant.genes$bayes.Achilles)>0) {
    genes.set.table$bayes.Achilles <- 
      sapply(genes.set.table$geneID, check.dataset,significant.genes$bayes.Achilles)}
  
  
  # count how often genes occur in sets
  genes.set.table$sum <- rowSums(genes.set.table[,-1])
  
  # extract all genes that occur in maximal number of sets
  max.set <- max(genes.set.table$sum)
  final.sign.genes <- filter(genes.set.table,sum==max.set)
  
  return(list(max.set=max.set,final.significant.genes=final.sign.genes))
}
