find.significant.genes <- function(threshold = 0.05, res) {
  
  # extract significant genes
  sign.genes.welch.DRIVE <- 
    res$p_values_welch_DRIVE %>% filter(p_adj < threshold) %>% pull(gene)
  sign.genes.welch.Achilles <- 
    res$p_values_welch_Achilles %>% filter(p_adj < threshold) %>% pull(gene)
  sign.genes.wilc.DRIVE <- 
    res$p_values_wilc_DRIVE %>% filter(p_adj < threshold) %>% pull(gene)
  sign.genes.wilc.Achilles <- 
    res$p_values_wilc_Achilles %>% filter(p_adj < threshold) %>% pull(gene)
  sign.genes.bayes.DRIVE <- 
    res$p_values_bayes_DRIVE %>% filter(p_adj < threshold) %>% pull(gene)
  sign.genes.bayes.Achilles <- 
    res$p_values_bayes_Achilles %>% filter(p_adj < threshold) %>% pull(gene)
  
  # collect in list
  sign.genes <- list()
  if(length(sign.genes.welch.DRIVE) > 0){
    sign.genes[[ "welch.DRIVE" ]] = sign.genes.welch.DRIVE}
  if(length(sign.genes.wilc.DRIVE) > 0){
    sign.genes[[ "wilc.DRIVE" ]] = sign.genes.wilc.DRIVE}
  if(length(sign.genes.bayes.DRIVE) > 0){
    sign.genes[[ "bayes.DRIVE" ]] = sign.genes.bayes.DRIVE}
  if(length(sign.genes.welch.Achilles) > 0){
    sign.genes[[ "welch.Achilles" ]] = sign.genes.welch.Achilles}
  if(length(sign.genes.wilc.Achilles) > 0){
    sign.genes[[ "wilc.Achilles" ]] = sign.genes.wilc.Achilles}
  if(length(sign.genes.bayes.Achilles) > 0){
    sign.genes[[ "bayes.Achilles" ]] = sign.genes.bayes.Achilles}
  
  return(sign.genes)
}