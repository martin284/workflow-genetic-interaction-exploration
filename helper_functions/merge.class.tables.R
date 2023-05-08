merge.class.tables <- function(class.table.gene1, class.table.gene2) {
  class.table.comb <- merge(class.table.gene1,class.table.gene2,by="CCLE_ID",
                            all=T)
  
  class.table.comb$isDeleterious.x <- replace(
    class.table.comb$isDeleterious.x,is.na(class.table.comb$isDeleterious.x),F)
  class.table.comb$isDeleterious.y <- replace(
    class.table.comb$isDeleterious.y,is.na(class.table.comb$isDeleterious.y),F)
  
  class.table.comb$isDeleterious <- 
    class.table.comb$isDeleterious.x |
    class.table.comb$isDeleterious.y
  class.table.comb$isDeleterious.x <- NULL
  class.table.comb$isDeleterious.y <- NULL
  
  return(class.table.comb)
}
