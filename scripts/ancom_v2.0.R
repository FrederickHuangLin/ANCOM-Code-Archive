library(exactRankTests)
library(nlme)
library(dplyr)

# OTU table should be a matrix/data.frame with each sample in rows and OTUs in columns. 
# OTU table should contains sample identifier with column name “Sample.ID”.
# Similarly, metadata should be a matrix/data.frame containing the sample identifier with column name “Sample.ID”. 

ANCOM = function(otu_data, meta_data, main_var, zero_cut = 0.90, p_adjust_method = "BH", alpha = 0.05,
                 adj_formula = NULL, rand_formula = NULL){
  # Let both the OTU table and metadata in the format of data.frame.
  otu_data = data.frame(otu_data, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  
  # Subset the metadata so that it contains the same group of samples as the OTU table and in the same order.
  meta_data = meta_data[match(otu_data$Sample.ID, meta_data$Sample.ID), ]
  
  # OTU table transformation: Add pseudocount (1) and take logarithm.
  comp_data = log(as.matrix(otu_data[, -1]) + 1)
  rownames(comp_data) = otu_data$Sample.ID
  
  # Filtering: Discard taxa with proportion of zeros > zero_cut.
  zero_prop = colMeans(comp_data == 0, na.rm = TRUE) 
  # colMeans(...) calculates column-wise means of the indicator matrix (TRUE if 0) to get the proportion of zeros.
  comp_data = comp_data[, zero_prop <= zero_cut]
  n_samp = dim(comp_data)[1]
  n_otu = dim(comp_data)[2]
  otu_id = colnames(comp_data)
  
  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    # Whether the main variable of interest has two levels or more?
    if (length(unique(meta_data%>%pull(main_var))) == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun = exactRankTests::wilcox.exact
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun = stats::kruskal.test
    }
    # Formula
    tformula = formula(paste("x ~", main_var, sep = " "))
  }else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun = stats::aov
    # Formula
    tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
  }else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var))
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }
  
  # Calculate the p-value for each pairwise comparison of taxa.
  p_data = matrix(NA, nrow = n_otu, ncol = n_otu)
  colnames(p_data) = otu_id
  rownames(p_data) = otu_id
  for (i in 1:(n_otu - 1)) {
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
    # e.g. the first alr matrix will be the log abundance data (comp_data) recursively subtracted 
    # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
    # the first (i - 1) columns were calculated by previous iterations, and
    # the i^th column contains all zeros.
    alr_data = apply(comp_data, 2, function(x) x - comp_data[, i]) 
    # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
    # Here, we basically want to iteratively subtract each column of the comp_data by its i^th column.
    alr_data = alr_data[, -(1:i), drop = FALSE]
    n_lr = dim(alr_data)[2] # number of log-ratios (lr)
    alr_data = cbind(alr_data, meta_data) # merge with the metadata
    
    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        tfun(tformula, data = data.frame(x, alr_data, check.names = FALSE))$p.value
        }
      ) 
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE), 
                   na.action = na.omit)
        summary(fit)[[1]][main_var, "Pr(>F)"]
        }
      )
    }else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(fixed = tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE),
                   random = formula(rand_formula),
                   na.action = na.omit)
        anova(fit)[main_var, "p-value"]
        }
      ) 
    }
  }
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  
  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adjust_method))
  
  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))
  
  # Organize results
  res = data.frame(otu_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_otu - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_otu - 1).
  res = res%>%mutate(detected_0.9 = ifelse(W > 0.9 * (n_otu -1), TRUE, FALSE),
                     detected_0.8 = ifelse(W > 0.8 * (n_otu -1), TRUE, FALSE),
                     detected_0.7 = ifelse(W > 0.7 * (n_otu -1), TRUE, FALSE),
                     detected_0.6 = ifelse(W > 0.6 * (n_otu -1), TRUE, FALSE))
  
  return(res)
}






