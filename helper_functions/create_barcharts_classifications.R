create_barcharts_classifications <- function(DRIVE,
                                             Achilles,
                                             gene_ID) {
  
  # count rows of gene_ID-D and gene_ID-ND classified cell lines
  counts_class_DRIVE <- data.frame(table(DRIVE$isDeleterious))
  counts_class_Achilles <- data.frame(table(Achilles$isDeleterious))
  
  barchart_DRIVE <- 
    ggplot(data = counts_class_DRIVE,
           aes(x = Var1, y = Freq, fill = Var1, label = Freq)) +
    theme(legend.position = "none") +
    geom_bar(stat = "identity", width = 0.5) +
    ggtitle("DRIVE") +
    labs(x = "", y = "number of cell lines") +
    scale_x_discrete(labels = c(paste0(gene_ID,"-ND"), paste0(gene_ID, "-D"))) +
    geom_text(size = 3, position = position_dodge(width = 0.9), vjust = -0.25)
  
  barchart_Achilles <-
    ggplot(data = counts_class_Achilles,
           aes(x = Var1, y = Freq, fill = Var1, label = Freq)) +
    theme(legend.position = "none") +
    geom_bar(stat = "identity", width = 0.5) +
    ggtitle("Achilles") +
    labs(x = "", y = "number of cell lines") +
    scale_x_discrete(labels = c(paste0(gene_ID, "-ND"), paste0(gene_ID, "-D")))+
    geom_text(size = 3, position = position_dodge(width = 0.9), vjust = -0.25)
  
  return(list(DRIVE = barchart_DRIVE, Achilles = barchart_Achilles))
}
