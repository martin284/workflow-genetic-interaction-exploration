create_volcano_plot <- function(test_results, dataset, test_method,
                                labeled_genes=NULL,
                                show_sign_line=TRUE){
  mean_diffs_str <- paste0("mean_differences_",dataset)
  mean_diffs <- test_results[[mean_diffs_str]]
  p_adj_log_str <- paste0("p_values_",test_method,'_',dataset)
  p_adj_log <- test_results[[p_adj_log_str]]['log']
  results <- cbind(p_adj_log, mean_diffs)
  results$gene <- rownames(results)

  results$to_label <- ""

  for (gene in labeled_genes) {
    results$to_label[results$gene == gene] <- gene
  }
  
  test_method_UC <- paste(toupper(substr(test_method, 1, 1)),
                          substr(test_method, 2, nchar(test_method)), sep="")
  title <- paste0(test_method_UC," Test with ",dataset)
  
  if (show_sign_line) {
    sign_line <- -log(0.05)
  }
  else {
    sign_line <- NULL
  }
  
  
  result_plot <- ggplot(data=results, aes(x=D_ND, y=log, label=to_label)) +
    geom_point(color="tomato3") +
    geom_point(data = results[results$to_label != "",], color="black") +
    xlab('mean difference') +
    ylab('-log(P_adj)') +
    geom_label_repel(max.overlaps = Inf, nudge_y = 0.5) +
    geom_hline(yintercept=sign_line) +
    ggtitle(title)
  
  return(result_plot)
}
