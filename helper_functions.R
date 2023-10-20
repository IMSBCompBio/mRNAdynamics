
library(reshape2)
library(rmutil)

# Define the loss function (a negative log likelihood) for the nuclear compartment
negloss_nuc = function(param,
                       tp,
                       traf_nuc,
                       tot_nuc){
  
  deg_nuc = param[1]           
  if (deg_nuc <= 0) return(Inf) # reject impossible parameters
  
  # Prediction of the nuclear new/total ratio for the given timepoints
  predict_nuc = 1 - exp(- deg_nuc * tp)  
  
  if (any(is.na(predict_nuc))) return(Inf)  # reject impossible outcomes
  if (any(predict_nuc < 0)) return(Inf)
  
  # - residual sum of squares for the comparison with the asin(sqrt()) transformed,
  # labeling bias corrected, observed new/total ratios in nucleus
  res = - sum((asin(sqrt(predict_nuc)) - traf_nuc)^2 * 2 * tot_nuc)
  return(-res)
}

# Define the loss function (a negative log likelihood) for the cytosolic compartment
negloss_cyt = function(param,
                       tp,
                       deg_nuc,
                       traf_cyt,
                       tot_cyt){
  
  deg_cyt = param[1]
  if (deg_cyt <= 0) return(Inf) # reject impossible parameters
  
  # Prediciton of the cytosolic new/total ratio for the given timepoints
  predict_cyt = 1 - ( exp(- deg_cyt * tp) +
                        deg_cyt / (deg_cyt - deg_nuc) * (exp(- deg_nuc * tp) - 
                                                           exp(- deg_cyt * tp)) )
  
  if (any(is.na(predict_cyt))) return(Inf)  # reject impossible outcomes
  if (any(predict_cyt < 0)) return(Inf)
  
  # - residual sum of squares for the comparison with the asin(sqrt()) transformed,
  # labeling bias corrected, observed new/total ratios in nucleus resp. cytosol
  res = - sum((asin(sqrt(predict_cyt)) - traf_cyt)^2 * 2 * tot_cyt)
  return(-res)
}

load_summarytable <- function(summarytable_file) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summarytable_file:  file storing summary table
  
  # *** RETURN ***
  
  # a data frame containing the summary table (rows: genomic regions, column 1: region name, column 2: measurement
  # description, column 3: library size, column 4: total transcript counts, column 5: labeled transcript counts,
  # column 6: average potential conversion positions, column 7: conversion efficiency column 8: newly synthesized ratio)  
  
  # *******************************************************************************************************************
  
  df_summarytable <- read.table(summarytable_file, sep = "\t")        # reading in summary table as data frame
  colnames(df_summarytable) <- c("name", "des", "lib", "tot", "mod", "cp", "ce", "nr")    # column names
  df_summarytable <- transform(df_summarytable, name = as.character(name), lib = as.numeric(lib), 
                               tot = as.numeric(tot), mod = as.numeric(mod), cp = as.numeric(cp), 
                               ce = as.numeric(ce), nr = as.numeric(nr))
  # transforming col 1 (gene names) to strings, 
  # 3, 4, 5 (read counts) to intergers and col 6, 7, 8 to floats
  
  return(df_summarytable)       # returning
}

summarytable_to_rawdata <- function(summarytable_file, timepoints) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summarytable_file:  file storing summary table
  # timepoints:         vector containing one string representation for each time points measured
  
  # *** RETURN ***
  
  # a named list containing 6 matrices storing (in the following order): nuclear modified read counts, nuclear total
  # read counts, nuclear estimated new/total ratios, cytosolic modified read counts, cytosolic total read counts,
  # cytosolic estimated new/total ratios; matrices store time points as columns and genomic region names as rows
  
  # *******************************************************************************************************************
  
  stable <- load_summarytable(summarytable_file)    # loading in summary table
  n_timepoints <- length(timepoints)                # number of time points measured
  n_measurements <- n_timepoints * 2                # number of measurements (time series for both nucleus and cytosol)
  n_regions <- nrow(stable) / n_measurements        # number of genomic regions
  region_names <- unique(stable[, "name"])          # names of genomic regions
  n_matrix_entries <- n_timepoints * n_regions      # number of output matrices' entries
  
  # initializing matrices
  
  mod_nu <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(mod_nu) <- timepoints
  rownames(mod_nu) <- region_names
  mod_cy <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(mod_cy) <- timepoints
  rownames(mod_cy) <- region_names
  tot_nu <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(tot_nu) <- timepoints
  rownames(tot_nu) <- region_names
  tot_cy <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(tot_cy) <- timepoints
  rownames(tot_cy) <- region_names
  ratio_nu <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(ratio_nu) <- timepoints
  rownames(ratio_nu) <- region_names
  ratio_cy <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(ratio_cy) <- timepoints
  rownames(ratio_cy) <- region_names
  
  # filling matrices
  
  for (i in 1:n_regions) {  # iterating through number of genomic regions
    
    nuc_startidx = (i-1) * n_measurements + 1     # start of nuclear measurements' data rows
    nuc_endidx = nuc_startidx + n_timepoints - 1  # end of nuclear measurements' data rows
    cyt_startidx = nuc_startidx + n_timepoints    # start of cytosolic measurements' data rows
    cyt_endidx = cyt_startidx + n_timepoints - 1  # end of cytosolic measurements' data rows
    
    mod_nu[i, ] <- stable[nuc_startidx:nuc_endidx, "mod"]
    mod_cy[i, ] <- stable[cyt_startidx:cyt_endidx, "mod"]
    tot_nu[i, ] <- stable[nuc_startidx:nuc_endidx, "tot"]
    tot_cy[i, ] <- stable[cyt_startidx:cyt_endidx, "tot"]
    ratio_nu[i, ] <- stable[nuc_startidx:nuc_endidx, "nr"]
    ratio_cy[i, ] <- stable[cyt_startidx:cyt_endidx, "nr"]
    
  }
  
  # returning
  
  return(list(mod_nu=mod_nu, tot_nu=tot_nu, ratio_nu=ratio_nu, mod_cy=mod_cy, tot_cy=tot_cy, ratio_cy=ratio_cy))
}

expr_robust_average <- function(summary_table, n_timepoints) {
  
  # *** PARAMETERS ***
  
  # summary_table:      summary table matrix
  # n_timepoints:       number of time points at which measurements were taken
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)
  
  # *** RETURN ***
  
  # a vector containing the average of the total counts distribution, computed for the middle 50% of the distribution
  # (for robustness), for each measurement time point
  
  n_rows <- nrow(summary_table)                                   # number of data rows in summary table
  n_genes <- n_rows / (2 * n_timepoints)                          # number of genes in the summary table
  quant_25 <- ceiling(n_genes/100 * 25)                           # 25% quantile index  
  quant_75 <- floor(n_genes/100 * 75)                             # 75% quantile index
  interval_50_size <- quant_75 - quant_25 + 1                     # 50% interval sample size
  total_collections <- sapply(1:(n_timepoints*2), function(i) {   # collecting total counts for each time point
    summary_table[seq(i, n_rows, by=(n_timepoints*2)), "tot"]
  })
  total_avgs <- sapply(1:ncol(total_collections), function(i) {   # computing averages of time points' total counts
    total_coll <- total_collections[, i]
    total_coll <- sort(total_coll)
    return(sum(total_coll[quant_25:quant_75]) / interval_50_size)
  })
  return(total_avgs)
  
}

expr_level_regression <- function(summary_table, time_series, norm="avg") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summary_table:      summary table file
  # time_series:        vector of time points at which measurements were taken (not including time point 0)
  # norm:               value by which total counts are normalized to compute expression levels; choose 'lib' to use 
  #                     the library size, choose 'avg' to use the average of the total counts distribution, computed
  #                     for the middle 50% of the distribution (for robustness) (default: "avg")
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)
  
  # *** RETURN ***
  
  # a matrix containing each gene's linear regression model (rows: genes, columns (for both nucleus and cytosol, in
  # that order): intercept, slope, slope normalized w.r.t average expression level, intercept p-value, slope p-value, 
  # model p-value, R2 value, coefficient of variation, average expression level)
  
  # *******************************************************************************************************************
  
  summary_table <- load_summarytable(summary_table)   # loading in summary tables
  n_timepoints <- length(time_series)                 # number of measurement time points
  gene_names <- unique(summary_table[, "name"])       # gene names
  
  # computing expression levels (rows: genes, columns: time series' expression levels)
  if (norm == "lib") {          # total/library_size ratios
    elvl <- summary_table[, "tot"] / summary_table[, "lib"]         
  }
  else if (norm == "avg") {     # total/total_avg ratios
    total_avgs <- expr_robust_average(summary_table, n_timepoints)
    elvl <- summary_table[, "tot"] / total_avgs
  }
  
  elvl <- matrix(elvl, ncol = 2 * n_timepoints, byrow = TRUE)     # matrix for expression
  nu_elvl <- elvl[,1:n_timepoints]                                # nuclear expression levels
  cy_elvl <- elvl[,(n_timepoints+1):(2*n_timepoints)]             # cytosolic expression levels
  
  # iterating through all genes, performing linear regression (matrix; columns: genes, rows (for both nucleus and 
  # cytosol, in that order): intercept, slope, intercept p-value, slope p-value, model p-value, R2 value
  lm_matrix <- sapply(1:nrow(elvl), function(i) {
    
    # data frames storing time series and corresponding expression levels of current gene
    gene_nu_elvl <- data.frame(matrix(c(time_series, nu_elvl[i, ]), ncol=2))
    colnames(gene_nu_elvl) <- c("time", "elvl")
    gene_cy_elvl <- data.frame(matrix(c(time_series, cy_elvl[i, ]), ncol=2))
    colnames(gene_cy_elvl) <- c("time", "elvl")
    
    # computing average expression levels
    nu_elvl_avg <- sum(gene_nu_elvl$elvl)/n_timepoints
    cy_elvl_avg <- sum(gene_cy_elvl$elvl)/n_timepoints
    
    # fitting linear model
    lm_nu <- lm(elvl ~ time, data=gene_nu_elvl)
    lm_cy <- lm(elvl ~ time, data=gene_cy_elvl)
    
    # computing model p-values
    lm_summary_nu <- summary(lm_nu)   # contains coefficients and corresponding p-values as well as F-statistic
    lm_summary_cy <- summary(lm_cy)
    lm_nu_cv <- cv(lm_nu)             # model coefficient of variation
    lm_cy_cv <- cv(lm_cy)
    
    fstat_nu <- lm_summary_nu$fstatistic                                        # getting F-statistic
    model_p_nu <- pf(fstat_nu[1], fstat_nu[2], fstat_nu[3], lower=FALSE)[[1]]   # getting model p-value
    fstat_cy <- lm_summary_cy$fstatistic                                    
    model_p_cy <- pf(fstat_cy[1], fstat_cy[2], fstat_cy[3], lower=FALSE)[[1]]
    
    # creating return column
    intercept_nu <- lm_summary_nu$coefficients[1,1]
    slope_nu <- lm_summary_nu$coefficients[2,1]
    slope_nu_norm <- slope_nu / nu_elvl_avg
    intercept_cy <- lm_summary_cy$coefficients[1,1]
    slope_cy <- lm_summary_cy$coefficients[2,1]
    slope_cy_norm <- slope_cy / cy_elvl_avg
    intercept_p_nu <- lm_summary_nu$coefficients[1,4]
    slope_p_nu <- lm_summary_nu$coefficients[2,4]
    intercept_p_cy <- lm_summary_cy$coefficients[1,4]
    slope_p_cy <- lm_summary_cy$coefficients[2,4]
    
    return_col <- 
      c(intercept_nu, slope_nu, slope_nu_norm, intercept_p_nu, slope_p_nu, model_p_nu, 
        lm_summary_nu$r.squared, lm_nu_cv, nu_elvl_avg,
        intercept_cy, slope_cy, slope_cy_norm, intercept_p_cy, slope_p_cy, model_p_cy, 
        lm_summary_cy$r.squared, lm_cy_cv, cy_elvl_avg)
    return(return_col)
  })
  
  # naming rows and columns of output matrix and returning
  
  colnames(lm_matrix) <- gene_names
  rownames(lm_matrix) <- 
    c("intercept_nu", "slope_nu", "slope_nu_norm", "intercept_p_nu", "slope_p_nu", "model_p_nu", 
      "R2_nu", "cv_nu", "avg_elvl_nu",
      "intercept_cy", "slope_cy", "slope_cy_norm", "intercept_p_cy", "slope_p_cy", "model_p_cy", 
      "R2_cy", "cv_cy", "avg_elvl_cy")
  lm_matrix <- t(lm_matrix)
  return(lm_matrix)
}

quantile_overlap <- function(r1, r2) {
  
  common_loci <- intersect(rownames(r1), rownames(r2))
  overlap_quantiles <- sapply(common_loci, function(i) {
    
    # lower borders of 95% quantiles
    
    lower_1_nu <- r1[i, "deg_nuc_q2.5%"]
    lower_2_nu <- r2[i, "deg_nuc_q2.5%"]
    lower_1_cy <- r1[i, "deg_cyt_q2.5%"]
    lower_2_cy <- r2[i, "deg_cyt_q2.5%"]
    
    # upper borders of 95% quantiles
    
    upper_1_nu <- r1[i, "deg_nuc_q97.5%"]
    upper_2_nu <- r2[i, "deg_nuc_q97.5%"]
    upper_1_cy <- r1[i, "deg_cyt_q97.5%"]
    upper_2_cy <- r2[i, "deg_cyt_q97.5%"]
    
    # nuclear half life quantiles
    if ((lower_1_nu >= upper_2_nu) | (lower_2_nu >= upper_1_nu)) {   # checking if quantiles do not overlap
      overlap_lower_nu <- NA
      overlap_upper_nu <- NA
    }
    else {
      if (lower_1_nu >= lower_2_nu) {
        overlap_lower_nu <- lower_1_nu
      }
      if (lower_2_nu >= lower_1_nu) {
        overlap_lower_nu <- lower_2_nu
      }
      if (upper_1_nu <= upper_2_nu) {
        overlap_upper_nu <- upper_1_nu
      }
      if (upper_2_nu <= upper_1_nu) {
        overlap_upper_nu <- upper_2_nu
      }
    }
    
    # cytosolic half life quantiles
    if ((lower_1_cy >= upper_2_cy) | (lower_2_cy >= upper_1_cy)) {   # checking if quantiles do not overlap
      overlap_lower_cy <- NA
      overlap_upper_cy <- NA
    }
    else {
      if (lower_1_cy >= lower_2_cy) {
        overlap_lower_cy <- lower_1_cy
      }
      if (lower_2_cy >= lower_1_cy) {
        overlap_lower_cy <- lower_2_cy
      }
      if (upper_1_cy <= upper_2_cy) {
        overlap_upper_cy <- upper_1_cy
      }
      if (upper_2_cy <= upper_1_cy) {
        overlap_upper_cy <- upper_2_cy
      }
    } 
    
    # returning overlap quantiles
    return(c(overlap_lower_nu, overlap_upper_nu, overlap_lower_cy, overlap_upper_cy))
  })
  
  overlap_quantiles <- t(overlap_quantiles)
  rownames(overlap_quantiles) <- common_loci
  colnames(overlap_quantiles) <- c("lower_nu", "upper_nu", "lower_cy", "upper_cy")
  return(overlap_quantiles)
  
}

estimation_data_table <- function(rawdata_1, rawdata_2, estimation_results_1, estimation_results_2, 
                                  elvl_regression_s1, elvl_regression_s2, timepoints_1, timepoints_2,
                                  min_cov, max_elvl_slope_stringent, max_elvl_slope_lessstringent, max_estdev, 
                                  max_ciq, min_r2_stringent, min_r2_lessstringent) {
  
  common_loci <- 
    intersect(rownames(estimation_results_1), rownames(estimation_results_2))   # loci measured in both samples
  
  # generate parameter estimation table
  
  deg_nuc_avg <- (estimation_results_1[common_loci, "deg_nuc"] + estimation_results_2[common_loci, "deg_nuc"]) / 2
  deg_cyt_avg <- (estimation_results_1[common_loci, "deg_cyt"] + estimation_results_2[common_loci, "deg_cyt"]) / 2
  rel_mu_s1 <- elvl_regression_s1[common_loci, "avg_elvl_nu"] * estimation_results_1[common_loci, "deg_nuc"] 
  rel_mu_s1_median <- median(rel_mu_s1)
  rel_mu_s1 <- rel_mu_s1 / rel_mu_s1_median
  rel_mu_s2 <- elvl_regression_s2[common_loci, "avg_elvl_nu"] * estimation_results_2[common_loci, "deg_nuc"]
  rel_mu_s2_median <- median(rel_mu_s2)
  rel_mu_s2 <- rel_mu_s2 / rel_mu_s2_median
  
  parameter_estimation_table <- cbind(
    estimation_results_1[common_loci, "half_life_nuc"], estimation_results_2[common_loci, "half_life_nuc"],
    (estimation_results_1[common_loci, "half_life_nuc"] + estimation_results_2[common_loci, "half_life_nuc"]) / 2,
    estimation_results_1[common_loci, "half_life_cyt"], estimation_results_2[common_loci, "half_life_cyt"],
    (estimation_results_1[common_loci, "half_life_cyt"] + estimation_results_2[common_loci, "half_life_cyt"]) / 2,
    estimation_results_1[common_loci, "deg_nuc"], estimation_results_2[common_loci, "deg_nuc"], deg_nuc_avg[common_loci],
    estimation_results_1[common_loci, "deg_cyt"], estimation_results_2[common_loci, "deg_cyt"], deg_cyt_avg[common_loci],
    rel_mu_s1[common_loci], rel_mu_s2[common_loci], (rel_mu_s1[common_loci] + rel_mu_s2[common_loci]) / 2,
    abs(estimation_results_1[common_loci, "deg_nuc"] - 
          estimation_results_2[common_loci, "deg_nuc"]) / deg_nuc_avg[common_loci],
    abs(estimation_results_1[common_loci, "deg_cyt"] -
          estimation_results_2[common_loci, "deg_cyt"]) / deg_cyt_avg[common_loci]
  )

  colnames(parameter_estimation_table) <- 
    c("half_life_nuc_s1", "half_life_nuc_s2", "half_life_nuc_avg", "half_life_cyt_s1", "half_life_cyt_s2", 
      "half_life_cyt_avg", "deg_nuc_s1", "deg_nuc_s2", "deg_nuc_avg", "deg_cyt_s1", "deg_cyt_s2", "deg_cyt_avg",
      "rel_mu_s1", "rel_mu_s2", "rel_mu_avg", "estimate_deviation_nuc", "estimate_deviation_cyt")
  
  # generate quantiles table
  
  quantiles_overlap_deg <- quantile_overlap(estimation_results_1, estimation_results_2)
  quantiles_overlap_hf <- log(2) / quantiles_overlap_deg
  
  quantiles_table <- cbind(
    estimation_results_1[common_loci, "halflife_nuc_q2.5%"], estimation_results_1[common_loci, "halflife_nuc_q97.5%"],
    estimation_results_2[common_loci, "halflife_nuc_q2.5%"], estimation_results_2[common_loci, "halflife_nuc_q97.5%"],
    quantiles_overlap_hf[common_loci, 1], quantiles_overlap_hf[common_loci, 2],
    estimation_results_1[common_loci, "halflife_cyt_q2.5%"], estimation_results_1[common_loci, "halflife_cyt_q97.5%"],
    estimation_results_2[common_loci, "halflife_cyt_q2.5%"], estimation_results_2[common_loci, "halflife_cyt_q97.5%"],
    quantiles_overlap_hf[common_loci, 3], quantiles_overlap_hf[common_loci, 4],
    estimation_results_1[common_loci, "deg_nuc_q2.5%"], estimation_results_1[common_loci, "deg_nuc_q97.5%"],
    estimation_results_2[common_loci, "deg_nuc_q2.5%"], estimation_results_2[common_loci, "deg_nuc_q97.5%"],
    quantiles_overlap_deg[common_loci, 1], quantiles_overlap_deg[common_loci, 2],
    estimation_results_1[common_loci, "deg_cyt_q2.5%"], estimation_results_1[common_loci, "deg_cyt_q97.5%"],
    estimation_results_2[common_loci, "deg_cyt_q2.5%"], estimation_results_2[common_loci, "deg_cyt_q97.5%"],
    quantiles_overlap_deg[common_loci, 3], quantiles_overlap_deg[common_loci, 4])
  
  colnames(quantiles_table) <-
    c("half_life_nuc_s1_q2.5%", "half_life_nuc_s1_q97.5%", "half_life_nuc_s2_q2.5%", "half_life_nuc_s2_q97.5%",
      "half_life_nuc_qoverlap_lower", "half_life_nuc_qoverlap_upper", "half_life_cyt_s1_q2.5%", "half_life_cyt_s1_q97.5%",
      "half_life_cyt_s2_q2.5%", "half_life_cyt_s2_q97.5%", "half_life_cyt_qoverlap_lower", "half_life_cyt_qoverlap_upper",
      "deg_nuc_s1_q2.5%", "deg_nuc_s1_q97.5%", "deg_nuc_s2_q2.5%", "deg_nuc_s2_q97.5%", "deg_nuc_qoverlap_lower",
      "deg_nuc_qoverlap_upper", "deg_cyt_s1_q2.5%", "deg_cyt_s1_q97.5%", "deg_cyt_s2_q2.5%", "deg_cyt_s2_q97.5%",
      "deg_cyt_qoverlap_lower", "deg_cyt_qoverlap_upper")
  
  # generate metadata table
  
  negloss_1_nuc <- sapply(rownames(estimation_results_1), function(g){
    negloss_nuc(estimation_results_1[g, "deg_nuc"], timepoints_1, 
                asin(sqrt(rawdata_1$ratio_nu[g, ])), rawdata_1$tot_nu[g, ])
  })
  names(negloss_1_nuc) <- rownames(estimation_results_1)
  negloss_2_nuc <- sapply(rownames(estimation_results_2), function(g){
    negloss_nuc(estimation_results_2[g, "deg_nuc"], timepoints_2, 
                asin(sqrt(rawdata_2$ratio_nu[g, ])), rawdata_2$tot_nu[g, ])
  })
  names(negloss_2_nuc) <- rownames(estimation_results_2)
  negloss_1_cyt <- sapply(rownames(estimation_results_1), function(g){
    negloss_cyt(estimation_results_1[g, "deg_cyt"], timepoints_1, estimation_results_1[g, "deg_nuc"],
                asin(sqrt(rawdata_1$ratio_cy[g, ])), rawdata_1$tot_cy[g, ])
  })
  names(negloss_1_cyt) <- rownames(estimation_results_1)
  negloss_2_cyt <- sapply(rownames(estimation_results_2), function(g){
    negloss_cyt(estimation_results_2[g, "deg_cyt"], timepoints_2, estimation_results_2[g, "deg_nuc"],
                asin(sqrt(rawdata_2$ratio_cy[g, ])), rawdata_2$tot_cy[g, ])
  })
  names(negloss_2_cyt) <- rownames(estimation_results_2)
  
  elvl_nu_avg <- (elvl_regression_s1[common_loci, "avg_elvl_nu"] + elvl_regression_s2[common_loci, "avg_elvl_nu"]) / 2
  elvl_cy_avg <- (elvl_regression_s1[common_loci, "avg_elvl_cy"] + elvl_regression_s2[common_loci, "avg_elvl_cy"]) / 2
  elvl_ratio_s1 <- (elvl_regression_s1[common_loci, "avg_elvl_cy"] / elvl_regression_s1[common_loci, "avg_elvl_nu"])
  elvl_ratio_s2 <- (elvl_regression_s2[common_loci, "avg_elvl_cy"] / elvl_regression_s2[common_loci, "avg_elvl_nu"])
  elvl_ratio_avg <- (elvl_ratio_s1 + elvl_ratio_s2) / 2
  
  metadata_table <- cbind(
    unname(negloss_1_nuc[common_loci]), unname(negloss_2_nuc[common_loci]),
    unname(negloss_1_cyt[common_loci]), unname(negloss_2_cyt[common_loci]),
    estimation_results_1[common_loci, "Rsquared_nuc"], estimation_results_2[common_loci, "Rsquared_nuc"],
    estimation_results_1[common_loci, "Rsquared_cyt"], estimation_results_2[common_loci, "Rsquared_cyt"],
    estimation_results_1[common_loci, "mean_coverage_nuc"], estimation_results_2[common_loci, "mean_coverage_nuc"],
    estimation_results_1[common_loci, "mean_coverage_cyt"], estimation_results_2[common_loci, "mean_coverage_cyt"],
    elvl_regression_s1[common_loci, "avg_elvl_nu"], elvl_regression_s2[common_loci, "avg_elvl_nu"], elvl_nu_avg[common_loci],
    elvl_regression_s1[common_loci, "avg_elvl_cy"], elvl_regression_s2[common_loci, "avg_elvl_cy"], elvl_cy_avg[common_loci],
    elvl_ratio_s1[common_loci], elvl_ratio_s2[common_loci], elvl_ratio_avg[common_loci],
    elvl_regression_s1[common_loci, "slope_nu_norm"], elvl_regression_s2[common_loci, "slope_nu_norm"], 
    (elvl_regression_s1[common_loci, "slope_nu_norm"] + elvl_regression_s2[common_loci, "slope_nu_norm"]) / 2,
    elvl_regression_s1[common_loci, "slope_cy_norm"], elvl_regression_s2[common_loci, "slope_cy_norm"],
    (elvl_regression_s1[common_loci, "slope_cy_norm"] + elvl_regression_s2[common_loci, "slope_cy_norm"]) / 2
  )
  
  colnames(metadata_table) <-
    c("negloss_nuc_s1", "negloss_nuc_s2", "negloss_cyt_s1", "negloss_cyt_s2", "r2_nuc_s1", "r2_nuc_s2", "r2_cyt_s1",
      "r2_cyt_s2", "mean_cov_nuc_s1", "mean_cov_nuc_s2", "mean_cov_cyt_s1", "mean_cov_cyt_s2", "mean_elvl_nuc_s1",
      "mean_elvl_nuc_s2", "mean_elvl_nuc_avg", "mean_elvl_cyt_s1", "mean_elvl_cyt_s2", "mean_elvl_cyt_avg",          
      "mean_elvl_ratio_cyt_nuc_s1", "mean_elvl_ratio_cyt_nuc_s2", "mean_elvl_ratio_cyt_nuc_avg", 
      "expr_lvl_slope_nuc_s1", "expr_lvl_slope_nuc_s2", "expr_lvl_slope_nuc_avg", 
      "expr_lvl_slope_cyt_s1", "expr_lvl_slope_cyt_s2", "expr_lvl_slope_cyt_avg")
  
  # generate reliability table
  
  reliability_booleans <- cbind(
    as.numeric(metadata_table[common_loci, "mean_cov_nuc_s1"] >= min_cov), 
    as.numeric(metadata_table[common_loci, "mean_cov_nuc_s2"] >= min_cov),
    as.numeric(metadata_table[common_loci, "mean_cov_cyt_s1"] >= min_cov), 
    as.numeric(metadata_table[common_loci, "mean_cov_cyt_s2"] >= min_cov),
    as.numeric(abs(metadata_table[common_loci, "expr_lvl_slope_nuc_avg"]) <= max_elvl_slope_stringent),
    as.numeric(abs(metadata_table[common_loci, "expr_lvl_slope_cyt_avg"]) <= max_elvl_slope_stringent),
    as.numeric(abs(metadata_table[common_loci, "expr_lvl_slope_nuc_avg"]) <= max_elvl_slope_lessstringent),
    as.numeric(abs(metadata_table[common_loci, "expr_lvl_slope_cyt_avg"]) <= max_elvl_slope_lessstringent),
    as.numeric(parameter_estimation_table[common_loci, "estimate_deviation_nuc"] <= max_estdev),
    as.numeric(parameter_estimation_table[common_loci, "estimate_deviation_cyt"] <= max_estdev),
    
    as.numeric(((parameter_estimation_table[common_loci, "deg_nuc_s1"] - quantiles_table[common_loci, "deg_nuc_s1_q2.5%"]) /
                  parameter_estimation_table[common_loci, "deg_nuc_s1"]) <= max_ciq),
    as.numeric(((quantiles_table[common_loci, "deg_nuc_s1_q97.5%"] - parameter_estimation_table[common_loci, "deg_nuc_s1"]) /
                  parameter_estimation_table[common_loci, "deg_nuc_s1"]) <= max_ciq),
    as.numeric(((parameter_estimation_table[common_loci, "deg_nuc_s2"] - quantiles_table[common_loci, "deg_nuc_s2_q2.5%"]) /
                  parameter_estimation_table[common_loci, "deg_nuc_s2"]) <= max_ciq),
    as.numeric(((quantiles_table[common_loci, "deg_nuc_s2_q97.5%"] - parameter_estimation_table[common_loci, "deg_nuc_s2"]) /
                  parameter_estimation_table[common_loci, "deg_nuc_s2"]) <= max_ciq),
    as.numeric(((parameter_estimation_table[common_loci, "deg_cyt_s1"] - quantiles_table[common_loci, "deg_cyt_s1_q2.5%"]) /
                  parameter_estimation_table[common_loci, "deg_cyt_s1"]) <= max_ciq),
    as.numeric(((quantiles_table[common_loci, "deg_cyt_s1_q97.5%"] - parameter_estimation_table[common_loci, "deg_cyt_s1"]) /
                  parameter_estimation_table[common_loci, "deg_cyt_s1"]) <= max_ciq),
    as.numeric(((parameter_estimation_table[common_loci, "deg_cyt_s2"] - quantiles_table[common_loci, "deg_cyt_s2_q2.5%"]) /
                  parameter_estimation_table[common_loci, "deg_cyt_s2"]) <= max_ciq),
    as.numeric(((quantiles_table[common_loci, "deg_cyt_s2_q97.5%"] - parameter_estimation_table[common_loci, "deg_cyt_s2"]) /
                  parameter_estimation_table[common_loci, "deg_cyt_s2"]) <= max_ciq),
    
    as.numeric(estimation_results_1[common_loci, "Rsquared_nuc"] >= min_r2_stringent),
    as.numeric(estimation_results_2[common_loci, "Rsquared_nuc"] >= min_r2_stringent),
    as.numeric(estimation_results_1[common_loci, "Rsquared_nuc"] >= min_r2_lessstringent),
    as.numeric(estimation_results_2[common_loci, "Rsquared_nuc"] >= min_r2_lessstringent),
    as.numeric(estimation_results_1[common_loci, "Rsquared_cyt"] >= min_r2_stringent),
    as.numeric(estimation_results_2[common_loci, "Rsquared_cyt"] >= min_r2_stringent),
    as.numeric(estimation_results_1[common_loci, "Rsquared_cyt"] >= min_r2_lessstringent),
    as.numeric(estimation_results_2[common_loci, "Rsquared_cyt"] >= min_r2_lessstringent)
  )
  
  rownames(reliability_booleans) <- common_loci
  colnames(reliability_booleans) <- 
    c("relab_cov_nuc_s1", "relab_cov_nuc_s2", "relab_cov_cyt_s1", "relab_cov_cyt_s2", 
      "relab_elvlslope_nuc", "relab_less_elvlslope_nuc", "relab_elvlslope_cyt", "relab_less_elvlslope_cyt", 
      "relab_estdev_nuc", "relab_estdev_cyt", 
      "relab_ciq_nuc_s1_lower", "relab_ciq_nuc_s1_upper", "relab_ciq_nuc_s2_lower", "relab_ciq_nuc_s2_upper",
      "relab_ciq_cyt_s1_lower", "relab_ciq_cyt_s1_upper", "relab_ciq_cyt_s2_lower", "relab_ciq_cyt_s2_upper",
      "relab_r2_nuc_s1", "relab_r2_nuc_s2", "relab_less_r2_nuc_s1", "relab_less_r2_nuc_s2",
      "relab_r2_cyt_s1", "relab_r2_cyt_s2", "relab_less_r2_cyt_s1", "relab_less_r2_cyt_s2")
  
  reliability_scores <- cbind(
    rowSums(reliability_booleans[common_loci, 
                                 c("relab_cov_nuc_s1", "relab_cov_nuc_s2", "relab_elvlslope_nuc", 
                                   "relab_estdev_nuc", 
                                   "relab_ciq_nuc_s1_lower", "relab_ciq_nuc_s1_upper",
                                   "relab_ciq_nuc_s2_lower", "relab_ciq_nuc_s2_upper",
                                    "relab_r2_nuc_s1", "relab_r2_nuc_s2")]),
    rowSums(reliability_booleans[common_loci, 
                                 c("relab_cov_nuc_s1", "relab_cov_nuc_s2", "relab_less_elvlslope_nuc",
                                   "relab_estdev_nuc",
                                   "relab_ciq_nuc_s1_lower", "relab_ciq_nuc_s1_upper",
                                   "relab_ciq_nuc_s2_lower", "relab_ciq_nuc_s2_upper",
                                   "relab_less_r2_nuc_s1", "relab_less_r2_nuc_s2")]),
    rowSums(reliability_booleans[common_loci, 
                                 c("relab_cov_nuc_s1", "relab_cov_nuc_s2", 
                                   "relab_cov_cyt_s1", "relab_cov_cyt_s2", 
                                   "relab_elvlslope_nuc", "relab_elvlslope_cyt",
                                   "relab_estdev_nuc", "relab_estdev_cyt", 
                                   "relab_ciq_nuc_s1_lower", "relab_ciq_nuc_s1_upper",
                                   "relab_ciq_nuc_s2_lower", "relab_ciq_nuc_s2_upper",
                                   "relab_ciq_cyt_s1_lower", "relab_ciq_cyt_s1_upper",
                                   "relab_ciq_cyt_s2_lower", "relab_ciq_cyt_s2_upper",
                                   "relab_r2_nuc_s1", "relab_r2_nuc_s2", 
                                   "relab_r2_cyt_s1", "relab_r2_cyt_s2")]),
    rowSums(reliability_booleans[common_loci, 
                                 c("relab_cov_nuc_s1", "relab_cov_nuc_s2", 
                                   "relab_cov_cyt_s1", "relab_cov_cyt_s2", 
                                   "relab_less_elvlslope_nuc", "relab_less_elvlslope_cyt",
                                   "relab_estdev_nuc", "relab_estdev_cyt", 
                                   "relab_ciq_nuc_s1_lower", "relab_ciq_nuc_s1_upper",
                                   "relab_ciq_nuc_s2_lower", "relab_ciq_nuc_s2_upper",
                                   "relab_ciq_cyt_s1_lower", "relab_ciq_cyt_s1_upper",
                                   "relab_ciq_cyt_s2_lower", "relab_ciq_cyt_s2_upper",
                                   "relab_less_r2_nuc_s1", "relab_less_r2_nuc_s2", 
                                   "relab_less_r2_cyt_s1", "relab_less_r2_cyt_s2")]) 
  )
  rownames(reliability_scores) <- common_loci
  colnames(reliability_scores) <- c("relab_score_stringent_nuc", "relab_score_lessstringent_nuc",
                                    "relab_score_stringent_both", "relab_score_lessstringent_both")
  
  reliability_outcomes <- cbind(
    as.numeric(reliability_scores[common_loci, "relab_score_stringent_nuc"] == 10),
    as.numeric(reliability_scores[common_loci, "relab_score_lessstringent_nuc"] == 10),
    as.numeric(reliability_scores[common_loci, "relab_score_stringent_both"] == 20),
    as.numeric(reliability_scores[common_loci, "relab_score_lessstringent_both"] == 20)
  )
  rownames(reliability_outcomes) <- common_loci
  colnames(reliability_outcomes) <- c("reliability_stringent_nuc", "reliability_lessstringent_nuc",
                                      "reliability_stringent_both", "reliability_lessstringent_both")
  
  reliability_table <- cbind(reliability_booleans[common_loci, ], reliability_scores[common_loci, ],
                             reliability_outcomes[common_loci, ])
  
  # all data table
  
  all_data_table = cbind(parameter_estimation_table[common_loci, ], quantiles_table[common_loci, ],
    metadata_table[common_loci, ], reliability_table[common_loci, ])
  
  # returning
  
  return(list(params=parameter_estimation_table, quants=quantiles_table, 
              meta=metadata_table, relab=reliability_table, all=all_data_table))
}

