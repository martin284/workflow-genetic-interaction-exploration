create.upset.plot <- function(significant.genes){
  
  comb.mat = make_comb_mat(list_to_matrix(significant.genes),
                           mode = "intersect")
  
  # define order of sets
  set.order <- c()
  if(length(significant.genes$welch.DRIVE) > 0) {
    set.order <- append(set.order,"welch.DRIVE")}
  if(length(significant.genes$wilc.DRIVE) > 0) {
    set.order <- append(set.order,"wilc.DRIVE")}
  if(length(significant.genes$bayes.DRIVE) > 0) {
    set.order <- append(set.order,"bayes.DRIVE")}
  if(length(significant.genes$welch.Achilles) > 0) {
    set.order <- append(set.order,"welch.Achilles")}
  if(length(significant.genes$wilc.Achilles) > 0) {
    set.order <- append(set.order,"wilc.Achilles")}
  if(length(significant.genes$bayes.Achilles) > 0) {
    set.order <- append(set.order,"bayes.Achilles")}
  
  # plot
  return(UpSet(comb.mat, set_order = set.order))
}