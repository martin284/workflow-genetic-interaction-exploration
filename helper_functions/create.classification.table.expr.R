create.classification.table.expr <- function(sample.share, 
                                             expression.class.gene,
                                             dep.scores.DRIVE,
                                             dep.scores.Achilles) {
  
  # keep only cell lines that appear in expression dataset and DRIVE
  expr.DRIVE <- merge(expression.class.gene, dep.scores.DRIVE["CCLE_ID"],
                      by = "CCLE_ID")
  # order the cell lines for upcoming classification 
  expr.DRIVE <- expr.DRIVE[order(expr.DRIVE$expr.level), ]
  
  # calculate the sample sizes resulting from the selected sample share
  sample.size.DRIVE <- round(sample.share * nrow(expr.DRIVE))
  
  # extract cell lines with highest and lowest expression in selected gene
  cell.lines.DRIVE.low.expr <- expr.DRIVE[1:sample.size.DRIVE, 1]
  cell.lines.DRIVE.high.expr <- 
    expr.DRIVE[(nrow(expr.DRIVE) - sample.size.DRIVE+1):nrow(expr.DRIVE), 1]
  
  # create classification table
  # to reuse code the variable is named isDeleterious again
  # isDeleterious==T stands for low expression
  class.table.DRIVE.low.expr <- 
    as.data.frame(cbind(CCLE_ID = cell.lines.DRIVE.low.expr, isDeleterious = F))
  class.table.DRIVE.high.expr <- 
    as.data.frame(cbind(CCLE_ID = cell.lines.DRIVE.high.expr, 
                        isDeleterious = T))
  class.table.DRIVE <- 
    rbind(class.table.DRIVE.low.expr, class.table.DRIVE.high.expr)
  
  # remove unnecessary objects
  rm(expr.DRIVE, cell.lines.DRIVE.low.expr, cell.lines.DRIVE.high.expr,
     class.table.DRIVE.low.expr, class.table.DRIVE.high.expr)
  
  # repeat the process for Achilles 
  expr.Achilles <- merge(expression.class.gene, dep.scores.Achilles["CCLE_ID"],
                         by="CCLE_ID")
  expr.Achilles <- expr.Achilles[order(expr.Achilles$expr.level), ]
  
  sample.size.Achilles <- round(sample.share * nrow(expr.Achilles))
  
  cell.lines.Achilles.low.expr <- expr.Achilles[1:sample.size.Achilles, 1]
  cell.lines.Achilles.high.expr <- 
    expr.Achilles[
      (nrow(expr.Achilles) - sample.size.Achilles + 1):nrow(expr.Achilles), 1]
  
  class.table.Achilles.low.expr <- 
    as.data.frame(cbind(CCLE_ID=cell.lines.Achilles.low.expr, 
                        isDeleterious = F))
  class.table.Achilles.high.expr <- 
    as.data.frame(cbind(CCLE_ID=cell.lines.Achilles.high.expr,
                        isDeleterious = T))
  
  class.table.Achilles <- 
    rbind(class.table.Achilles.low.expr, class.table.Achilles.high.expr)
  
  return(list(DRIVE = class.table.DRIVE, Achilles = class.table.Achilles))
}
