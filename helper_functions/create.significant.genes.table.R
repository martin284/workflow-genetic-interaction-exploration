create.significant.genes.table <- function(test.results,
                                           dataset, test.method,
                                           significance.level=0.05) {
  
  p.values.str <- paste0("p_values_",test.method,'_',dataset)
  p.values <- test.results[[p.values.str]]
  
  # order genes by adjusted p-value
  p.values.ordered <- p.values[order(p.values$p_adj),]
  
  # filter for adjusted p-values lower than given significance level
  p.values.ordered.filtered <- p.values.ordered[
    p.values.ordered$p_adj<=significance.level,]
  
  return(p.values.ordered.filtered)
}
