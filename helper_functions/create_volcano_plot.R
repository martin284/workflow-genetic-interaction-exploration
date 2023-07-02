create_volcano_plot <- function(test_results,
                                dataset,
                                test_method,
                                labeled_genes = NULL,
                                show_sign_line = TRUE,
                                show_y_label = TRUE) {  
  
  # extract mean differences and p-values for the given dataset
  mean_diffs_str <- paste0("mean_differences_", dataset)
  mean_diffs <- test_results[[mean_diffs_str]]
  p_adj_log_str <- paste0("p_values_", test_method, '_', dataset)
  p_adj_log <- test_results[[p_adj_log_str]]['log']
  results <- cbind(p_adj_log, mean_diffs)
  results$gene <- rownames(results)
  
  # define labeled genes  
  results$to_label <- ""
  for (gene in labeled_genes) {
    results$to_label[results$gene == gene] <- gene
  }
  
  # define if line representing significance threshold is shown
  if (show_sign_line) {
    sign_line <- -log10(0.05)
  }
  else {
    sign_line <- NULL
  }
  
  # define if y-label is shown
  if (show_y_label) {
    y_label <- expression(-log[10](P[adj]))
  } else {
    y_label <- ""
  }
  
  # create result plot
  result_plot <- 
    ggplot(data = results, aes(x = D_ND, y = log, label = to_label)) +
    geom_point(color = "tomato3") +
    geom_point(data = results[results$to_label != "",], color = "black") +
    xlab('Mean difference') +
    ylab(y_label) +
    geom_label_repel(max.overlaps = Inf, nudge_y = 0.5) +
    geom_hline(yintercept = sign_line) +
    ggtitle(dataset) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15))
  
  return(result_plot)
}
