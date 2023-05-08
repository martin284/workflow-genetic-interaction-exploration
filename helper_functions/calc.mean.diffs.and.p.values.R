calc.mean.diffs.and.p.values <- function(dep_scores_DRIVE, 
                                         dep_scores_Achilles, 
                                         CCLE_mut_DRIVE, 
                                         CCLE_mut_Achilles) {
  
  ######################### calculate mean differences #########################
  
  # add D or ND classification to each cell line
  dep_scores_DRIVE_merged <- merge(
    dep_scores_DRIVE,
    CCLE_mut_DRIVE[c("CCLE_ID","isDeleterious")],
    by="CCLE_ID")
  
  # split dataset in D and ND
  dep_scores_DRIVE_D_uncleaned <- 
    dep_scores_DRIVE_merged[dep_scores_DRIVE_merged$isDeleterious==T,
                     c(-1,-ncol(dep_scores_DRIVE_merged))]
  dep_scores_DRIVE_ND_uncleaned <- 
    dep_scores_DRIVE_merged[dep_scores_DRIVE_merged$isDeleterious==F,
                     c(-1,-ncol(dep_scores_DRIVE_merged))]
  
  # check if there are at least cell lines with deleteriously mutated gene
  if (nrow(dep_scores_DRIVE_D_uncleaned)<2) {
    stop("Less than two cell lines with deleteriously mutated gene.")
  }
  
  # remove all gene columns where less than 2 non NA values are left
  dep_scores_DRIVE_D <- 
    dep_scores_DRIVE_D_uncleaned[ ,
                          colSums(!is.na(dep_scores_DRIVE_D_uncleaned))>1]
  dep_scores_DRIVE_ND <- 
    dep_scores_DRIVE_ND_uncleaned[ ,
                          colSums(!is.na(dep_scores_DRIVE_D_uncleaned))>1]
  
  # show warning if gene columns were removed
  if (ncol(dep_scores_DRIVE_D)!=ncol(dep_scores_DRIVE_D_uncleaned)) {
    warning(as.character(
      ncol(dep_scores_DRIVE_D_uncleaned)-ncol(dep_scores_DRIVE_D)),
      " from ", 
      as.character(ncol(dep_scores_DRIVE_D_uncleaned)),
      " gene columns have been removed from DRIVE due to low sample size."
      )
  }
  
  # calculate means
  means_DRIVE_D <- data.frame(colMeans(dep_scores_DRIVE_D, na.rm=T))
  means_DRIVE_ND <- data.frame(colMeans(dep_scores_DRIVE_ND, na.rm=T))
  
  # calculate differences D-ND
  mean_differences_DRIVE <- means_DRIVE_D - means_DRIVE_ND
  colnames(mean_differences_DRIVE) <- "D_ND"
  
  # do the same for Achilles without warning
  dep_scores_Achilles_merged <- merge(
    dep_scores_Achilles,
    CCLE_mut_Achilles[c("CCLE_ID","isDeleterious")],
    by="CCLE_ID")
  
  dep_scores_Achilles_D <-
    dep_scores_Achilles_merged[dep_scores_Achilles_merged$isDeleterious==T,
                             c(-1,-ncol(dep_scores_Achilles_merged))]
  dep_scores_Achilles_ND <- 
    dep_scores_Achilles_merged[dep_scores_Achilles_merged$isDeleterious==F,
                             c(-1,-ncol(dep_scores_Achilles_merged))]
  
  means_Achilles_D <- data.frame(colMeans(dep_scores_Achilles_D, na.rm=T))
  means_Achilles_ND <- 
    data.frame(colMeans(dep_scores_Achilles_ND, na.rm=T))
  mean_differences_Achilles <- means_Achilles_D - means_Achilles_ND
  colnames(mean_differences_Achilles) <- "D_ND"
  
  ############### calculate p-values with welch test ###########################
  
  # DRIVE
  test_results_welch_DRIVE <- col_t_welch(dep_scores_DRIVE_D,
                                          dep_scores_DRIVE_ND)
  p_values_welch_DRIVE <- data.frame(
    matrix(NA, nrow = nrow(test_results_welch_DRIVE), ncol = 2))
  colnames(p_values_welch_DRIVE) <- c('gene','p_value')
  p_values_welch_DRIVE$gene <- rownames(test_results_welch_DRIVE)
  p_values_welch_DRIVE$p_value <- test_results_welch_DRIVE$pvalue
  
  # Achilles
  test_results_welch_Achilles <-
    col_t_welch(dep_scores_Achilles_D, dep_scores_Achilles_ND)
  p_values_welch_Achilles <-
    data.frame(matrix(NA, nrow = nrow(test_results_welch_Achilles), ncol = 2))
  colnames(p_values_welch_Achilles) <- c('gene','p_value')
  p_values_welch_Achilles$gene <- rownames(test_results_welch_Achilles)
  p_values_welch_Achilles$p_value <- test_results_welch_Achilles$pvalue
  
  ### transform p-values ###
  # DRIVE
  p_values_welch_DRIVE$p_adj <- p.adjust(p_values_welch_DRIVE$p_value,
                                         method = "BH")
  p_values_welch_DRIVE$log <- -log(p_values_welch_DRIVE$p_adj)
  
  # Achilles
  p_values_welch_Achilles$p_adj <- p.adjust(p_values_welch_Achilles$p_value,
                                         method = "BH")
  p_values_welch_Achilles$log <- -log(p_values_welch_Achilles$p_adj)
  
  ############### calculate p-values with wilcoxon test ########################
  
  # DRIVE
  test_results_wilc_DRIVE <- col_wilcoxon_twosample(dep_scores_DRIVE_D,
                                                    dep_scores_DRIVE_ND)
  p_values_wilc_DRIVE <- data.frame(
    matrix(NA, nrow = nrow(test_results_wilc_DRIVE), ncol = 2))
  colnames(p_values_wilc_DRIVE) <- c('gene','p_value')
  p_values_wilc_DRIVE$gene <- rownames(test_results_wilc_DRIVE)
  p_values_wilc_DRIVE$p_value <- test_results_wilc_DRIVE$pvalue
  
  # Achilles
  test_results_wilc_Achilles <-
    col_wilcoxon_twosample(dep_scores_Achilles_D, dep_scores_Achilles_ND)
  p_values_wilc_Achilles <-
    data.frame(matrix(NA, nrow = nrow(test_results_wilc_Achilles), ncol = 2))
  colnames(p_values_wilc_Achilles) <- c('gene','p_value')
  p_values_wilc_Achilles$gene <- rownames(test_results_wilc_Achilles)
  p_values_wilc_Achilles$p_value <- test_results_wilc_Achilles$pvalue
  
  ### transform p-values ###
  # DRIVE
  p_values_wilc_DRIVE$p_adj <- p.adjust(p_values_wilc_DRIVE$p_value,
                                         method = "BH")
  p_values_wilc_DRIVE$log <- -log(p_values_wilc_DRIVE$p_adj)
  
  # Achilles
  p_values_wilc_Achilles$p_adj <- p.adjust(p_values_wilc_Achilles$p_value,
                                            method = "BH")
  p_values_wilc_Achilles$log <- -log(p_values_wilc_Achilles$p_adj)
  
  ########## calculate p-values with Bayes moderated t-test ####################
  
  # calculate p-values
  vec_Achilles <- dep_scores_Achilles_merged$isDeleterious
  mat_Achilles <- 
    dep_scores_Achilles_merged[,-ncol(dep_scores_Achilles_merged)]
  mat_Achilles$isDeleterious <- NULL
  rownames(mat_Achilles) <- mat_Achilles$CCLE_ID
  mat_Achilles$CCLE_ID <- NULL
  mat_Achilles <- as.matrix(mat_Achilles)
  
  design_Achilles <- model.matrix(~vec_Achilles)
  fit_Achilles <- lmFit(t(mat_Achilles),design=design_Achilles)
  fit_Achilles <- eBayes(fit = fit_Achilles)
  test_results_bayes_Achilles <- topTable(fit = fit_Achilles, n=Inf, coef = 2,
                                          sort.by = "none")
  
  p_values_bayes_Achilles <- data.frame(
    matrix(NA, nrow = nrow(test_results_bayes_Achilles), ncol = 2))
  colnames(p_values_bayes_Achilles) <- c('gene','p_adj')
  p_values_bayes_Achilles$gene <- rownames(test_results_bayes_Achilles)
  p_values_bayes_Achilles$p_adj <- test_results_bayes_Achilles$adj.P.Val
  
  # transform p-values
  p_values_bayes_Achilles$log <- -log(p_values_bayes_Achilles$p_adj)
  
  # do the same with DRIVE
  vec_DRIVE <- dep_scores_DRIVE_merged$isDeleterious
  # consider removed NAs -> unlike Achilles
  mat_DRIVE <- rbind(dep_scores_DRIVE_D, dep_scores_DRIVE_ND)
  mat_DRIVE <- mat_DRIVE[order(as.numeric(rownames(mat_DRIVE))),]
  mat_DRIVE$isDeleterious <- NULL
  rownames(mat_DRIVE) <- mat_DRIVE$CCLE_ID
  mat_DRIVE$CCLE_ID <- NULL
  mat_DRIVE <- as.matrix(mat_DRIVE)
  
  design_DRIVE <- model.matrix(~vec_DRIVE)
  fit_DRIVE <- lmFit(t(mat_DRIVE),design=design_DRIVE)
  fit_DRIVE <- eBayes(fit = fit_DRIVE)
  test_results_bayes_DRIVE <- topTable(fit = fit_DRIVE, n=Inf, coef = 2,
                                       sort.by = "none")
  
  p_values_bayes_DRIVE <- data.frame(
    matrix(NA, nrow = nrow(test_results_bayes_DRIVE), ncol = 2))
  colnames(p_values_bayes_DRIVE) <- c('gene','p_adj')
  p_values_bayes_DRIVE$gene <- rownames(test_results_bayes_DRIVE)
  p_values_bayes_DRIVE$p_adj <- test_results_bayes_DRIVE$adj.P.Val
  
  p_values_bayes_DRIVE$log <- -log(p_values_bayes_DRIVE$p_adj)
  
  ############## return all results ############################################
  
  return(list(mean_differences_DRIVE=mean_differences_DRIVE,
         mean_differences_Achilles=mean_differences_Achilles,
         p_values_welch_DRIVE=p_values_welch_DRIVE,
         p_values_welch_Achilles=p_values_welch_Achilles,
         p_values_wilc_DRIVE=p_values_wilc_DRIVE,
         p_values_wilc_Achilles=p_values_wilc_Achilles,
         p_values_bayes_DRIVE=p_values_bayes_DRIVE,
         p_values_bayes_Achilles=p_values_bayes_Achilles))
}

  
