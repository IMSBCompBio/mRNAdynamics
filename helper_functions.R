# VERSION: 0.1.0

library(ggplot2)
library(gplots)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(scales)

library(foreach)
library(doParallel)

library(deSolve)
library(R2Cuba)
library(DEoptim)
library(rmutil)
library(sjstats)

# defining plotting colours ...........................................................................................

colours_1 = c("yellow4", "yellow3", "gold", "darkgoldenrod3", "darkgoldenrod4",
              "salmon1", "salmon3", "pink1", "pink3", "pink4")
colours_2 = c("dodgerblue4", "dodgerblue2", "skyblue3", "lightblue1",
              "turquoise1", "turquoise3", "turquoise4", "seagreen1", "seagreen3", "seagreen4")
colours_3 = c("grey1", "grey30", "grey60", "mediumpurple4", "mediumpurple2", 
              "magenta4", "magenta2", "hotpink2", "deeppink3", "deeppink4")
colours_4 = c("lightsalmon", "maroon3", "mediumorchid1", "slateblue1", "navy",
              "royalblue1", "lightsteelblue2", "darkseagreen4", "khaki3", "gold")

### Basic Computations ################################################################################################

find_peak <- function(f, f_borders, eval_points = 1000, width = 10) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # f:              function of which the peak is to be found
  # f_borders:      (minimum, maximum) argument to be given to function f (vector)
  
  # eval_points:    optional; number of points the function is to be evaluated at to find the maximum (default: 1000)
  # width:          optional; minimum width the peak should span as percent of the function's domain's range 
  #                 (default: 10)
  
  # *** RETURN ***
  
  # a vector containing the lower and upper peak-spanning x-value
  
  # *******************************************************************************************************************
  
  # sequence of x-values to evaluate the function at
  x_vals <- seq(f_borders[1], f_borders[2], length.out = eval_points)
  # evaluating function
  fv <- Vectorize(f)
  y_vals <- fv(x_vals)
  # finding maximum y-value
  peak_idx <- which.max(y_vals)
  peak_yval <- y_vals[peak_idx]
  
  # defining peak region according to width
  domain <- f_borders[2] - f_borders[1]               # width of the function's domain 
  width <- width/100                                  # transforming min_width percent to ratio
  missing_width <- eval_points * width / 2            # missing width in indices
  lowerindex <- peak_idx - missing_width              # peak region lower index
  upperindex <- peak_idx + missing_width              # peak region upper index
  
  if (lowerindex < 1) {                               # secure that lower peak region does not exceed domain
    lowerindex <- 1
  }                         
  if (upperindex > eval_points) {                     # secure that lower peak region does not exceed domain
    upperindex <- eval_points
  }     
  
  # returning
  return(c(x_vals[lowerindex], x_vals[upperindex]))
}

approx1D <- function(f, lower, upper, ip_points = 100) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # f:            function to be interpolated
  # lower:        lower domain border of the function f
  # upper:        upper domain border of the function f
  # ip_points:    optional; number of points to be used for interpolation (default: 100)
  
  # *** NOTE ***
  
  # data points outside of the domain specified will be interpolated with 0
  
  # *** RETURN ***
  
  # an interpolated function to f
  
  # *******************************************************************************************************************
  
  fv <- Vectorize(f)
  x_vals <- seq(lower, upper, length.out = ip_points)
  f_vals <- fv(x_vals)
  f_interpolation <- approxfun(x_vals, f_vals, yleft = 0, yright = 0)
  return(f_interpolation)
}

approx2D <- function(f, lower, upper, ip_points = 100) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # f:          function to be interpolated; takes an vector as argument, defining (first, second) dimension values
  # lower:      lower domain borders of (first, second) dimension of the function f (vector)
  # upper:      upper domain borders of (first, second) dimension of the function f (vector)
  # ip_points:  optional; number of points to be used for interpolation (per dimension) (default: 100)
  
  # *** NOTE ***
  
  # data points outside of the domain specified will be interpolated with 0
  
  # *** RETURN ***
  
  # an interpolated function to f
  
  # *******************************************************************************************************************
  
  # computing x-values for both dimensions as well as x-axis step size
  
  x1_vals <- seq(lower[1], upper[1], length.out = ip_points)
  x2_vals <- seq(lower[2], upper[2], length.out = ip_points)
  
  x1_stepsize <- x1_vals[2] - x1_vals[1]    # stepsize and offset are needed to correctly locate a request within
  x2_stepsize <- x2_vals[2] - x2_vals[1]    #   the grid (find the corresponding grid indices)
  x1_offset <- lower[1]
  x2_offset <- lower[2]
  d1_maxidx <- length(x1_vals)
  d2_maxidx <- length(x2_vals)
  
  # computing y-values for the function (as matrix: dimension 1 as columns, dimension 2 as rows)
  
  f_vals <- sapply(x1_vals, function(i) {
    print(i)
    sapply(x2_vals,  function(j){f(c(i, j))})
  })
  
  # defining interpolation function:
  # for a data point request, this function computes a weighted average value based on the data points grid available
  
  f_interpolation <- function(request) {
    
    x1_req <- request[1]    # first dimension x-value of request
    x2_req <- request[2]    # second dimension x-value of request
    
    # determining where the request is located within the grid available
    d1_idx_lower <- floor((x1_req - x1_offset) / x1_stepsize) + 1
    d1_idx_upper <- d1_idx_lower + 1
    d2_idx_lower <- floor((x2_req - x2_offset) / x2_stepsize) + 1
    d2_idx_upper <- d2_idx_lower + 1
    # checking if request is outside of specified domain
    if (d1_idx_lower < 1){return(0)}
    if (d1_idx_upper > d1_maxidx){return(0)}
    if (d2_idx_lower < 1){return(0)}
    if (d2_idx_upper > d2_maxidx){return(0)}
    # getting request surrounding grid x-values
    x1_lower <- x1_vals[d1_idx_lower]
    x1_upper <- x1_vals[d1_idx_upper]
    x2_lower <- x2_vals[d2_idx_lower]
    x2_upper <- x2_vals[d2_idx_upper]
    # getting request surrounding grid y-values (rows: dim2, columns: dim1)
    y_low_low <- f_vals[d2_idx_lower, d1_idx_lower]
    y_low_up <- f_vals[d2_idx_upper, d1_idx_lower]
    y_up_low <- f_vals[d2_idx_lower, d1_idx_upper]
    y_up_up <- f_vals[d2_idx_upper, d1_idx_upper]
    # computing weighted average
    request_approximation <- y_low_low * (x1_upper - x1_req) * (x2_upper - x2_req) +
      y_low_up * (x1_upper - x1_req) * (x2_req - x2_lower) + y_up_low * (x1_req - x1_lower) * (x2_upper - x2_req) +
      y_up_up + (x1_req - x1_lower) * (x2_req - x2_lower) * 1
    # returning
    return(request_approximation)
  }
  
  # returning interpolated function
  
  return(f_interpolation)
}

### Plotting ##########################################################################################################

plot_format <- function(ggplot_object, title, x_label, y_label, legend_title = NULL, colour_palette = NULL, 
                        hline = NULL, vline = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # ggplot_object:      ggplot object to be formatted
  # title:              title of the plot
  # x_label:            x-axis description
  # y_label:            y-axis description
  
  # legend_title:       optional; title of the legend (default: NULL)
  # colour_palette:     optional; colour palette to use (default: NULL)
  # hline:              optional; position of a horizontal line to draw (default: NULL)
  # vline:              optional; position of a vertical line to draw (default: NULL)
  
  # *** RETURN ***
  
  # A ggplot object containing all formattings.
  
  # *******************************************************************************************************************
  
  ggplot_object <- ggplot_object +
    scale_color_manual(legend_title, values = colour_palette) +
    theme(panel.background = element_rect(fill = 'white', colour = 'gray')) +
    theme(panel.grid.major = element_line(color = "grey", linetype = "dashed")) +
    theme(panel.grid.minor = element_line(color = "gray85", linetype = "dashed")) +
    theme(legend.key = element_blank()) +
    xlab(x_label) + ylab(y_label) + ggtitle(title)
  if (!is.null(hline)) {ggplot_object <- ggplot_object + geom_hline(yintercept = hline, size=0.5)}
  if (!is.null(vline)) {ggplot_object <- ggplot_object + geom_vline(xintercept = vline, size=0.5)}
  return(ggplot_object)
}

heatmap_format <- function(data_matrix, col_names=NULL, row_names=NULL, plot_title=NULL, 
                           plot_xlab=NULL, plot_ylab=NULL, key_xlab=NULL, hline = NULL, vline = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # data_matrix:    input data matrix to be plotted as heatmap
  
  # col_names:      optional; heatmap column names (vector) (default: NULL)
  # row_names:      optional; heatmap row names (vector) (default: NULL)
  # plot_title:     optional; heatmap title (default: NULL)
  # plot_xlab:      optional; x-axis description (default: NULL)
  # plot_ylab:      optional; y-axis description (default: NULL)
  # key_xlab:       optional; colour code description (default: NULL)
  # hline:          optional; position of a horizontal line to draw (default: NULL)
  # vline:          optional: position of a vertical line to draw (default: NULL)
  
  # *** RETURN ***
  
  # A heatmap object containing all formattings.
  
  # *******************************************************************************************************************
  
  colnames(data_matrix) <- col_names
  rownames(data_matrix) <- row_names
  
  hm <- heatmap.2(data_matrix, Rowv = F, Colv = F, trace = "none", density.info = "none",
                  labRow = NA, labCol = NA, margins = c(2, 2), 
                  symbreaks = F, rowsep=c(hline), colsep = c(vline),
                  key.xlab = key_xlab, xlab = plot_xlab, ylab = plot_ylab, main = plot_title,
                  col=colorRampPalette(rev(brewer.pal(11, "Spectral")))(1000))
  return(hm)
}

plot_dynamics <- function(dynamics_dataframe, rna_species_labels, title, x_label, y_label, legend_title, colours) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # dynamics_dataframe:     matrix containing RNA transcript counts (columns: RNA species, rows: time points)
  # rna_species_labels:     RNA species names (vector)
  # title:                  title of the plot
  # x_label:                x-axis description
  # y_label:                y-label description
  # legend_title:           legend title
  # colours:                colour palette to use (vector)
  
  # *** RETURN ***
  
  # A ggplot object containing the RNA dynamics plot.
  
  # *******************************************************************************************************************
  
  colnames(dynamics_dataframe) <- rna_species_labels
  dynamics_dfm <- melt(dynamics_dataframe, id.vars = "time")
  dynamics_plot <- ggplot(dynamics_dfm) + 
    geom_point(aes(x = dynamics_dfm$time, y = dynamics_dfm$value, color = dynamics_dfm$variable), size = 1.5) 
  dynamics_plot <- plot_format(dynamics_plot, title, x_label, y_label, legend_title, colours) 
  return(dynamics_plot)
}

plot_densities <- function(x_vals, density_functions, 
                           function_labels, title, x_label, y_label, legend_title, colours, p_true = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # x_vals:               x-values to be plotted (vector)
  # density_functions:    density functions to be plotted (list)
  
  # function_labels:      function descriptions
  # title:                plot title
  # x_label:              x-axis description
  # y_label:              y-axis description
  # legend_title:         legend title
  # colours:              colour palette to use (vector)
  
  # p_true:               optional; parameter value which should yield the highest density
  
  # *** RETURN ***
  
  # A ggplot object containing the density functions' plot.
  
  # *******************************************************************************************************************
  
  # computing matrix: columns: density functions, rows: x values
  y_vals <- sapply(1:length(density_functions), function(i) {
    current_density <- density_functions[[i]]
    sapply(x_vals, function(x) {
      current_density(x)
    })
  })
  
  # appending x values to be able to melt data frame
  y_vals <- cbind(y_vals, x_vals)
  # constructing melted data frame
  density_df <- data.frame(y_vals)
  colnames(density_df) <- c(function_labels, "xvals")
  density_dfm <- melt(density_df, id.vars = "xvals")
  
  # plotting
  density_plot <- ggplot(density_dfm) +
    geom_point(aes(x = density_dfm$xvals, y = density_dfm$value, color = density_dfm$variable), size = 1)
  density_plot <- plot_format(density_plot, title, x_label, y_label, legend_title, colours, vline = p_true)
  return(density_plot)
}

plot_clusters <- function(cluster_table, time_series, titles, x_label, y_label, legend_title, colours,
                          plot_lines = FALSE, ylim_nu = NULL, ylim_cy = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # cluster_table:  a file storing cluster center values (in the following format):
  #                 First Column: cluster names
  #                 Second Column: cluster sizes
  #                 Other Columns: time points' ratios
  #                 Rows:
  #                   first row:  time points itself (two times listing of the time series, for nuclear and cytosolic
  #                               ratios)
  #                   other rows: cluster numbers' center values
  
  # time_series:    vector of time points at which measurements were taken
  # titles:         vector; titles to give the plots (nuclear ratios, cytosolic ratios)
  # x_label:        x-axis description
  # y_labels:       y-axis description
  # legend_title:   legend title
  # colours:        colour palette to use for plotting
  
  # plot_lines:     optional; indicates whether plots should display lines instead of points (default: FALSE)
  # ylim_nu:        optional; (minimum, maximum) nuclear ratio to show in plot (default: NULL) 
  # ylim_cy:        optional; (minimum, maximum) cytosolic ratio to show in plot (default: NULL) 
  
  # *** RETURN ***
  
  # a list containing the dynamics cluster plots of the nuclear ratios and the cytosolic ratios
  
  # *******************************************************************************************************************
  
  cluster_data <- read.table(cluster_table, sep="\t", header = TRUE)  # reading in clustering data (rows: clusters, 
                                                                      # columns: time point ratios)
  n_clusters <- nrow(cluster_data)          # number of clusters in the clustering data
  cluster_names <- cluster_data[, 1]        # cluster names
  cluster_sizes <- cluster_data[, 2]        # cluster sizes
  
  n_timepoints <- length(time_series)       # number of time points
  
  # getting data frame for each compartment's ratios (rows: time points, first col: time points, other cols: clusters)
                                              
  nu_ratio_df <- data.frame(                # nuclear ratios
    sapply(1:n_clusters, function(i) {
      sapply(1:n_timepoints, function(t) {
        cluster_data[i, t + 2]
      })
    })
  )
  nu_ratio_df <- cbind(time = time_series, nu_ratio_df)
  
  cy_ratio_df <- data.frame (               # cytosolic ratios
    sapply(1:n_clusters, function(i) {
      sapply((n_timepoints+1):(2*n_timepoints), function(t) {
        cluster_data[i, t + 2]
      })
    })
  )
  cy_ratio_df <- cbind(time = time_series, cy_ratio_df)
  
  # setting data frame column names (time and cluster names + cluster sizes)
  
  col_names <- c("time", sapply(1:n_clusters, function(i) {
    paste(cluster_names[i], " (#", cluster_sizes[i], ")")
  }))
  colnames(nu_ratio_df) <- col_names
  colnames(cy_ratio_df) <- col_names
  
  # melting data frames (for plotting multiple data rows)
  
  nu_ratio_df <- melt(nu_ratio_df, id.vars = "time")
  cy_ratio_df <- melt(cy_ratio_df, id.vars = "time")
  
  # plotting dynamics clusters
  
  if(!plot_lines) {
    nu_ratio_plot <- ggplot(nu_ratio_df) +
      geom_point(aes(x = nu_ratio_df$time, y = nu_ratio_df$value, color = nu_ratio_df$variable), size = 1.5) +
      scale_y_continuous(limits = ylim_nu)
    nu_ratio_plot <- plot_format(nu_ratio_plot, titles[1], x_label, y_label, legend_title, colours)
    
    cy_ratio_plot <- ggplot(cy_ratio_df) +
      geom_point(aes(x = cy_ratio_df$time, y = cy_ratio_df$value, color = cy_ratio_df$variable), size = 1.5) +
      scale_y_continuous(limits = ylim_cyl)
    cy_ratio_plot <- plot_format(cy_ratio_plot, titles[2], x_label, y_label, legend_title, colours)
  }
  else {
    nu_ratio_plot <- ggplot(nu_ratio_df) +
      geom_line(aes(x = nu_ratio_df$time, y = nu_ratio_df$value, color = nu_ratio_df$variable), size = 1) +
      scale_y_continuous(limits = ylim_nu)
    nu_ratio_plot <- plot_format(nu_ratio_plot, titles[1], x_label, y_label, legend_title, colours)
    
    cy_ratio_plot <- ggplot(cy_ratio_df) +
      geom_line(aes(x = cy_ratio_df$time, y = cy_ratio_df$value, color = cy_ratio_df$variable), size = 1) +
      scale_y_continuous(limits = ylim_cy)
    cy_ratio_plot <- plot_format(cy_ratio_plot, titles[2], x_label, y_label, legend_title, colours)
  }
  
  # returning
  
  return(list(nu_ratio_plot, cy_ratio_plot))
}

plot_signals <- function(original_table_s1, original_table_s2, norm_factors_s1, norm_factors_s2, time_series, genes,
                         outdir, file_prefix, title_prefix = "") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # original_table_s1:    data frame containing measurement data of sample 1
  # original_table_s2:    data frame containing measurement data of sample 2
  # norm_factors_s1:      vector of normalization factors used to compute expression levels for sample 1, for each 
  #                       time point (not including time point 0), for nucleus and cytosol
  # norm_factors_s2:      vector of normalization factors used to compute expression levels for sample 2, for each 
  #                       time point (not including time point 0), for nucleus and cytosol
  # time_series:          vector of time points at which measurements were taken (not including time point 0)
  # genes:                EITHER number of genomic regions to select randomly and plot time curves for,
  #                       OR vector of genomic region names for which to plot time curves
  
  # outdir:         output directory to store plots to
  # file_prefix:    file name prefix of the plots, will be extended by the genes' names
  # title_prefix:   optional; text to set as prefix for the plot titles (default: "")
  
  # *** NOTE ***
  
  # data frames storing measurement data must be sorted by gene names and time points (nuclear time series, 
  # cytosolic time series) and have the following format:
  # rows: genomic regions
  # col 1: region name, col 2: measurement description, col 3: library size, 
  # col 4: total transcript counts, col 5: labeled transcript counts, col 6: average potential conversion positions, 
  # col 7: conversion efficiency, col 8: newly synthesized ratio
  
  # *** RETURN ***
  
  # vector of plotted genes' names
  
  # *******************************************************************************************************************
  
  n_timepoints <- length(time_series)                             # number of time series measurements
  time_series <- c(0, time_series)                                # adding time point 0 to time series
  
  # selecting gene names for which to plot time curves
  gene_names <- c()
  if(typeof(genes) == "double") {             # --- case: select genes randomly ---
    # all gene names stored in both sample 1 and 2
    gene_names_collection <- intersect(rownames(estimation_table_s1), rownames(estimation_table_s2))     
    # selecting genes for which time curves are to be plotted (randomly)
    gene_names <- gene_names_collection[sample.int(length(gene_names_collection), genes)] 
  }
  else if(typeof(genes) == "character") {     # --- case: gene names are pre-defined ---
    gene_names <- genes
  }
  
  for(g in gene_names) {  # iterating through gene names 
    
    subtable_s1 <- original_table_s1[which(original_table_s1[, "name"] == g), ]   # loading summary table 1
    subtable_s2 <- original_table_s2[which(original_table_s2[, "name"] == g), ]   # loading summary table 2
    
    # *** new/total ratios ***
    
    new_total_measured_nu_s1 <-     # sample 1, nucleus (NA as missing value for time point 0)
      c(NA, subtable_s1[1:n_timepoints, "nr"])                      
    new_total_measured_cy_s1 <-     # sample 1, cytosol (NA as missing value for time point 0)
      c(NA, subtable_s1[(n_timepoints+1):(2*n_timepoints), "nr"])
    new_total_measured_nu_s2 <-     # sample 2, nucleus (NA as missing value for time point 0)
      c(NA, subtable_s2[1:n_timepoints, "nr"])                      
    new_total_measured_cy_s2 <-     # sample 2, cytosol (NA as missing value for time point 0)
      c(NA, subtable_s2[(n_timepoints+1):(2*n_timepoints), "nr"])
    
    # *** measured expression levels ***
    
    elvl_nu_s1 <-     # sample 1, nucleus 
      c(subtable_s1[1:n_timepoints, "tot"] / norm_factors_s1[1:7])
    elvl_cy_s1 <-     # sample 1, cytosol
      c(subtable_s1[(n_timepoints+1):(2*n_timepoints), "tot"] / norm_factors_s1[8:14])
    elvl_nu_s2 <-     # sample 2, nucleus
      c(subtable_s2[1:n_timepoints, "tot"] / norm_factors_s2[1:7])
    elvl_cy_s2 <-     # sample 2, cytosol
      c(subtable_s2[(n_timepoints+1):(2*n_timepoints), "tot"] / norm_factors_s2[8:14])
    
    # normalizing expression level ratios w.r.t. average (NA as missing value for time point 0)
    
    elvl_nu_s1 <- c(NA, elvl_nu_s1 / (sum(elvl_nu_s1) / n_timepoints))
    elvl_nu_s2 <- c(NA, elvl_nu_s2 / (sum(elvl_nu_s2) / n_timepoints))
    elvl_cy_s1 <- c(NA, elvl_cy_s1 / (sum(elvl_cy_s1) / n_timepoints))
    elvl_cy_s2 <- c(NA, elvl_cy_s2 / (sum(elvl_cy_s2) / n_timepoints))

    # plotting ........................................................................................................
    
    # building data frames
    
    new_total_nu_df <- data.frame(cbind(time_series, new_total_measured_nu_s1, new_total_measured_nu_s2))
    new_total_cy_df <- data.frame(cbind(time_series, new_total_measured_cy_s1, new_total_measured_cy_s2))
    elvl_nu_df <- data.frame(cbind(time_series, elvl_nu_s1, elvl_nu_s2))
    elvl_cy_df <- data.frame(cbind(time_series, elvl_cy_s1, elvl_cy_s2))
    
    colnames(new_total_nu_df) <- c("time", "1", "2")
    colnames(new_total_cy_df) <- c("time", "1", "2")
    colnames(elvl_nu_df) <- c("time", "1", "2")
    colnames(elvl_cy_df) <- c("time", "1", "2")
    
    new_total_nu_df <- melt(new_total_nu_df, id.vars = "time")
    new_total_cy_df <- melt(new_total_cy_df, id.vars = "time")
    elvl_nu_df <- melt(elvl_nu_df, id.vars = "time")
    elvl_cy_df <- melt(elvl_cy_df, id.vars = "time")
    
    # generating & formatting & storing plot
      
    new_total_nu_plot <- ggplot(new_total_nu_df, aes(x = time, y = value, colour = variable)) + geom_point() + geom_line()
    new_total_cy_plot <- ggplot(new_total_cy_df, aes(x = time, y = value, colour = variable)) + geom_point() + geom_line()
    elvl_nu_plot <- ggplot(elvl_nu_df, aes(x = time, y = value, colour = variable)) + geom_point() + geom_line()
    elvl_cy_plot <- ggplot(elvl_cy_df, aes(x = time, y = value, colour = variable)) + geom_point() + geom_line()
    
    new_total_nu_plot <- plot_format(new_total_nu_plot, paste("new/total ratio nucleus\n", g, sep=""), 
                                     "time [min]", "ratio", legend_title = "sample", 
                                     colour_palette = c("dodgerblue4", "darkgoldenrod3"))
    new_total_cy_plot <- plot_format(new_total_cy_plot, paste("new/total ratio cytosol\n", g, sep=""), 
                                     "time [min]", "ratio", legend_title = "sample", 
                                     colour_palette = c("dodgerblue4", "darkgoldenrod3"))
    elvl_nu_plot <- plot_format(elvl_nu_plot, paste("normalized expression level nucleus\n", g, sep=""), 
                                "time [min]", "ratio", legend_title = "sample", 
                                colour_palette = c("dodgerblue4", "darkgoldenrod3"))
    elvl_cy_plot <- plot_format(elvl_cy_plot, paste("normalized expression level cytosol\n", g, sep=""), 
                                "time [min]", "ratio", legend_title = "sample", 
                                colour_palette = c("dodgerblue4", "darkgoldenrod3"))
      
    pdf(paste(outdir, file_prefix, g, ".pdf", sep=""))
    grid.arrange(new_total_nu_plot, new_total_cy_plot, elvl_nu_plot, elvl_cy_plot, nrow=2)
    dev.off()
  }
  
  # returning
  return(gene_names)
}

halflife_estimation_timecurves <- function(original_table_s1, original_table_s2,
                                           estimation_table_s1, estimation_table_s2, time_series, genes,
                                           conveff_s1, conveff_s2, fn_error_nu, fn_error_cy, fp_error_nu, fp_error_cy,
                                           outdir, file_prefix, title_prefix = "",
                                           lab_measured = FALSE, lab_inferred = FALSE,
                                           new_predicted = FALSE, new_recorded = FALSE, new_inferred = FALSE,
                                           split_samples = FALSE) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # original_table_s1:    data frame containing measurement data of sample 1
  # original_table_s2:    data frame containing measurement data of sample 2
  # estimation_table_s1:  data frame containing estimation data of sample 1
  # estimation_table_s2:  data frame containing estimation data of sample 2
  # time_series:          vector of time points at which measurements were taken (not including time point 0)
  # genes:                EITHER number of genomic regions to select randomly and plot time curves for,
  #                       OR vector of genomic region names for which to plot time curves
 
  # conveff_s1:         vector of sample 1 conversion efficiencies (at single measurement time points, ordered by
  #                     nuclear and cytosolic time series); probability by which a single base is labeled
  # conveff_s2:         vector of sample 2 conversion efficiencies (at single measurement time points, ordered by
  #                     nuclear and cytosolic time series); probability by which a single base is labeled
  # fn_error_nu:        nuclear false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for nuclear transcripts)
  # fn_error_cy:        cytosolic false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for cytosolic transcripts)
  # fp_error_nu:        nuclear false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for nuclear transcripts)
  # fp_error_cy:        cytosolic false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for cytosolic transcripts)
  
  # outdir:         output directory to store plots to
  # file_prefix:    file name prefix of the plots, will be extended by the genes' names
  
  # title_prefix:   optional; text to set as prefix for the plot titles (default: "")
  # lab_measured:   set to TRUE to plot labeled/total ratios measured (default: FALSE)
  # lab_inferred:   set to TRUE to plot labeled/total ratios inferred from new/total ratios predicted (default: FALSE)
  # new_predicted:  set to TRUE to plot new/total ratios predicted (default: FALSE)
  # new_recorded:   set to TRUE to plot new/total ratios recorded in summary table (default: FALSE)
  # new_inferred:   set to TRUE to plot new/total ratios inferred from labeled/total ratios measured (default: FALSE)
  # split_samples:  set to TRUE to plot both biological samples in two separate plots (default: FALSE)
  
  # *** NOTE ***
  
  # data frames storing measurement data must be sorted by gene names and time points (nuclear time series, 
  # cytosolic time series) and have the following format:
  # rows: genomic regions
  # col 1: region name, col 2: measurement description, col 3: library size, 
  # col 4: total transcript counts, col 5: labeled transcript counts, col 6: average potential conversion positions, 
  # col 7: conversion efficiency, col 8: newly synthesized ratio
  
  # data frames storing estimation data should be in the following format:
  # rows: genomic regions,
  # col 1: lambda_n + tau estimator, col 2: lambda_c estimator, col 3: sum-of-squares values
  
  # *** RETURN ***
  
  # vector of plotted genes' names
  
  # *******************************************************************************************************************
  
  n_timepoints <- length(time_series)                             # number of time series measurements
  time_series <- c(0, time_series)                                # adding time point 0 to time series
  conveff_nu_s1 <- c(0, conveff_s1[1:n_timepoints])                       # nuclear conversion efficiency, sample 1
  conveff_cy_s1 <- c(0, conveff_s1[(n_timepoints+1):(2*n_timepoints)])    # cytosolic conversion efficiency, sample 1
  conveff_nu_s2 <- c(0, conveff_s2[1:n_timepoints])                       # nuclear conversion efficiency, sample 2
  conveff_cy_s2 <- c(0, conveff_s2[(n_timepoints+1):(2*n_timepoints)])    # cytosolic conversion efficiency, sample 2
  
  # selecting gene names for which to plot time curves
  gene_names <- c()
  if(typeof(genes) == "double") {             # --- case: select genes randomly ---
    # all gene names stored in both sample 1 and 2
    gene_names_collection <- intersect(rownames(estimation_table_s1), rownames(estimation_table_s2))     
    # selecting genes for which time curves are to be plotted (randomly)
    gene_names <- gene_names_collection[sample.int(length(gene_names_collection), genes)] 
  }
  else if(typeof(genes) == "character") {     # --- case: gene names are pre-defined ---
    gene_names <- genes
  }
  
  for(g in gene_names) {  # iterating through gene names 
    
    # original measurements ...........................................................................................
    
    subtable_s1 <- original_table_s1[which(original_table_s1[, "name"] == g), ]
    subtable_s2 <- original_table_s2[which(original_table_s2[, "name"] == g), ]
    
    # *** measured labeled/total ratios ***
    
    lab_total_measured_nu_s1 <-     # sample 1, nucleus (NA as missing value for time point 0)
      c(NA, subtable_s1[1:n_timepoints, "mod"] / subtable_s1[1:n_timepoints, "tot"])                      
    lab_total_measured_cy_s1 <-     # sample 1, cytosol (NA as missing value for time point 0)
      c(NA, subtable_s1[(n_timepoints+1):(2*n_timepoints), "mod"] / subtable_s1[(n_timepoints+1):(2*n_timepoints), "tot"])   # (NA as missing value for time point 0)
    lab_total_measured_nu_s2 <-     # sample 2, nucleus (NA as missing value for time point 0)
      c(NA, subtable_s2[1:n_timepoints, "mod"] / subtable_s2[1:n_timepoints, "tot"])                      
    lab_total_measured_cy_s2 <-     # sample 2, cytosol (NA as missing value for time point 0)
      c(NA, subtable_s2[(n_timepoints+1):(2*n_timepoints), "mod"] / subtable_s2[(n_timepoints+1):(2*n_timepoints), "tot"])
    
    # *** recorded new/total ratios ***
    
    new_total_recorded_nu_s1 <-     # sample 1, nucleus (NA as missing value for time point 0)
      c(NA, subtable_s1[1:n_timepoints, "nr"])
    new_total_recorded_cy_s1 <-     # sample 1, cytosol (NA as missing value for time point 0)
      c(NA, subtable_s1[(n_timepoints+1):(2*n_timepoints), "nr"])
    new_total_recorded_nu_s2 <-     # sample 2, nucleus (NA as missing value for time point 0)
      c(NA, subtable_s2[1:n_timepoints, "nr"])
    new_total_recorded_cy_s2 <-     # sample 2, cytosol (NA as missing value for time point 0)
      c(NA, subtable_s2[(n_timepoints+1):(2*n_timepoints), "nr"])
    
    # estimated values ................................................................................................
    
    # weighted average number of potential conversion positions
    convpos_s1 <- sum(original_table_s1[which(original_table_s1[, "name"] == g), "cp"] * 
                        original_table_s1[which(original_table_s1[, "name"] == g), "tot"]) / 
      sum(original_table_s1[which(original_table_s1[, "name"] == g), "tot"])      # sample 1 
    convpos_s2 <- sum(original_table_s2[which(original_table_s2[, "name"] == g), "cp"] * 
                        original_table_s2[which(original_table_s2[, "name"] == g), "tot"]) / 
      sum(original_table_s2[which(original_table_s2[, "name"] == g), "tot"])      # sample 2
    
    # true positive and true negative rates
    p_mod_nu_s1 <- 1 -                                # nuclear true positive rate, sample 1
      (1 - conveff_nu_s1 * (1 - fn_error_nu)) ** convpos_s1    
    p_mod_cy_s1 <- 1 -                                # cytosolic true positive rate, sample 1
      (1 - conveff_cy_s1 * (1 - fn_error_cy)) ** convpos_s1    
    p_un_nu_s1 <- (1 - fp_error_nu) ** convpos_s1     # nuclear true negative rate, sample 1
    p_un_cy_s1 <- (1 - fp_error_cy) ** convpos_s1     # cytosolic true negative rate, sample 1
    p_mod_nu_s2 <- 1 -                                # nuclear true positive rate, sample 2
      (1 - conveff_nu_s2 * (1 - fn_error_nu)) ** convpos_s2    
    p_mod_cy_s2 <- 1 -                                # cytosolic true positive rate, sample 2
      (1 - conveff_cy_s2 * (1 - fn_error_cy)) ** convpos_s2    
    p_un_nu_s2 <- (1 - fp_error_nu) ** convpos_s2     # nuclear true negative rate, sample 2
    p_un_cy_s2 <- (1 - fp_error_cy) ** convpos_s2     # cytosolic true negative rate, sample 2
    
    # sample 1 lambda_n + tau and lambda_c estimators
    lnt_est_s1 <- estimation_table_s1[g, "deg_nuc"]
    lc_est_s1 <- estimation_table_s1[g, "deg_cyt"]
    
    # sample 2 lambda_n + tau and lambda_c estimators
    lnt_est_s2 <- estimation_table_s2[g, "deg_nuc"]
    lc_est_s2 <- estimation_table_s2[g, "deg_cyt"]

    # *** predicted new/total ratios ***
    
    lnt_split_s1 <- lnt_est_s1 / 2  # splitting lambda_n + tau estimator into two halfs (one for each variable)
    lnt_split_s2 <- lnt_est_s2 / 2  
    
    new_total_predicted_nu_s1 <- ratio_n_sol(time_series, lnt_split_s1, lnt_split_s1)           
    new_total_predicted_cy_s1 <- ratio_c_sol(time_series, lnt_split_s1, lnt_split_s1, lc_est_s1)
    new_total_predicted_nu_s2 <- ratio_n_sol(time_series, lnt_split_s2, lnt_split_s2)           
    new_total_predicted_cy_s2 <- ratio_c_sol(time_series, lnt_split_s2, lnt_split_s2, lc_est_s2)

    # *** inferred labeled/total ratios ***
    
    lab_total_inferred_nu_s1 <- 
      new_total_predicted_nu_s1 * p_mod_nu_s1 + (1 - new_total_predicted_nu_s1) * (1 - p_un_nu_s1)
    lab_total_inferred_cy_s1 <- 
      new_total_predicted_cy_s1 * p_mod_cy_s1 + (1 - new_total_predicted_cy_s1) * (1 - p_un_cy_s1)
    lab_total_inferred_nu_s2 <- 
      new_total_predicted_nu_s2 * p_mod_nu_s2 + (1 - new_total_predicted_nu_s2) * (1 - p_un_nu_s2)
    lab_total_inferred_cy_s2 <- 
      new_total_predicted_cy_s2 * p_mod_cy_s2 + (1 - new_total_predicted_cy_s2) * (1 - p_un_cy_s2)
    
    # *** inferred new/total ratios ***
    
    new_total_inferred_nu_s1 <- 
      (lab_total_measured_nu_s1 + p_un_nu_s1 - 1) / (p_mod_nu_s1 + p_un_nu_s1 - 1)
    new_total_inferred_cy_s1 <- 
      (lab_total_measured_cy_s1 + p_un_cy_s1 - 1) / (p_mod_cy_s1 + p_un_cy_s1 - 1)
    new_total_inferred_nu_s2 <- 
      (lab_total_measured_nu_s2 + p_un_nu_s2 - 1) / (p_mod_nu_s2 + p_un_nu_s2 - 1)
    new_total_inferred_cy_s2 <- 
      (lab_total_measured_cy_s2 + p_un_cy_s2 - 1) / (p_mod_cy_s2 + p_un_cy_s2 - 1)
    
    # plotting ........................................................................................................
    
    # initializing data & colour of data to display in the plots
    
    data_rows <- 21
    data_vector <- 
      c(time_series, 
        lab_total_measured_nu_s1, new_total_predicted_nu_s1, lab_total_inferred_nu_s1, new_total_inferred_nu_s1, new_total_recorded_nu_s1,
        lab_total_measured_cy_s1, new_total_predicted_cy_s1, lab_total_inferred_cy_s1, new_total_inferred_cy_s1, new_total_recorded_cy_s1,
        lab_total_measured_nu_s2, new_total_predicted_nu_s2, lab_total_inferred_nu_s2, new_total_inferred_nu_s2, new_total_recorded_nu_s2,
        lab_total_measured_cy_s2, new_total_predicted_cy_s2, lab_total_inferred_cy_s2, new_total_inferred_cy_s2, new_total_recorded_cy_s2)
    data_names <- 
      c("time",
        "labeled/total measured (nu, 1)", "new/total predicted (nu, 1)", "labeled/total inferred (nu, 1)", "new/total inferred (nu, 1)", "new/total corrected (nu, 1)",
        "labeled/total measured (cy, 1)", "new/total predicted (cy, 1)", "labeled/total inferred (cy, 1)", "new/total inferred (cy, 1)", "new/total corrected (cy, 1)",
        "labeled/total measured (nu, 2)", "new/total predicted (nu, 2)", "labeled/total inferred (nu, 2)", "new/total inferred (nu, 2)", "new/total corrected (nu, 2)",
        "labeled/total measured (cy, 2)", "new/total predicted (cy, 2)", "labeled/total inferred (cy, 2)", "new/total inferred (cy, 2)", "new/total corrected (cy, 2)")
    colour_palette <- c("navy", "dodgerblue4", "dodgerblue2", "skyblue3", "lightblue1",
                        "darkgreen", "limegreen", "lightgreen", "yellow4", "yellow3",
                        "salmon3", "salmon1", "pink4", "pink3", "pink1",
                        "grey1", "grey32", "grey48", "grey65", "grey80")
    line_type <- rep(c(1, 3, 1, 1, 1), 4)
    file_suffix <- c("lm", "li", "np", "nr", "ni")
    
    data_matrix <- matrix(data_vector, ncol=data_rows)    # data matrix containing all data for both samples
    df <- data.frame(data_matrix)                         # corresponding data frame
    colnames(df) <- data_names                            # setting column (data row) names
    
    # filtering data for data rows to be plotted
    
    if (!lab_measured) {
      df[, data_names[c(2, 7, 12, 17)]] <- NULL
      colour_palette[c(1, 6, 11, 16)] <- NA
      line_type[c(1, 6, 11, 16)] <- NA
      file_suffix[1] <- NA
    }
    if (!new_predicted) {
      df[, data_names[c(3, 8, 13, 18)]] <- NULL
      colour_palette[c(2, 7, 12, 17)] <- NA
      line_type[c(2, 7, 12, 17)] <- NA
      file_suffix[3] <- NA
    }
    if (!lab_inferred) {
      df[, data_names[c(4, 9, 14, 19)]] <- NULL
      colour_palette[c(3, 8, 13, 18)] <- NA
      line_type[c(3, 8, 13, 18)] <- NA
      file_suffix[2] <- NA
    }
    if (!new_inferred) {
      df[, data_names[c(5, 10, 15, 20)]] <- NULL
      colour_palette[c(4, 9, 14, 19)] <- NA
      line_type[c(4, 9, 14, 19)] <- NA
      file_suffix[5] <- NA
    }
    if (!new_recorded) {
      df[, data_names[c(6, 11, 16, 21)]] <- NULL
      colour_palette[c(5, 10, 15, 20)] <- NA
      line_type[c(5, 10, 15, 20)] <- NA
      file_suffix[4] <- NA
    }
    
    colour_palette <- colour_palette[!is.na(colour_palette)]    # filtered data row colours
    line_type <- line_type[!is.na(line_type)]                   # filtered line types
    file_suffix <- 
      paste(file_suffix[!is.na(file_suffix)], collapse = "_")   # file suffix indicating data rows plotted
    
    # generating & formatting & storing plot
    
    if (isTRUE(split_samples)) {
      df1 <- df[, seq(1, ((ncol(df) - 1) / 2 + 1))]
      df2 <- df[, c(1, seq(((ncol(df) - 1) / 2 + 2), ncol(df)))]
      df1m <- melt(df1, id.vars = "time")   # melted data frame (for plotting multiple data series)
      df2m <- melt(df2, id.vars = "time")   # melted data frame (for plotting multiple data series)
      
      dynamics_plot_1 <- 
        ggplot(df1m, aes(x = df1m$time, y = df1m$value, color = df1m$variable, linetype=df1m$variable)) + 
        geom_point(size = 1.5) + geom_line(show.legend = FALSE) + scale_linetype_manual(values=line_type)
      dynamics_plot_2 <- 
        ggplot(df2m, aes(x = df2m$time, y = df2m$value, color = df2m$variable, linetype=df2m$variable)) + 
        geom_point(size = 1.5) + geom_line(show.legend = FALSE) + scale_linetype_manual(values=line_type)
      dynamics_plot_1 <- plot_format(dynamics_plot_1, paste(title_prefix, ", s1, ", g), 
                                     "time [min]", "ratio w.r.t total read counts", legend_title = "time course", 
                                     colour_palette = colour_palette)
      dynamics_plot_2 <- plot_format(dynamics_plot_2, paste(title_prefix, ", s2, ", g), 
                                     "time [min]", "ratio w.r.t total read counts", legend_title = "time course", 
                                     colour_palette = colour_palette)
      
      pdf(paste(outdir, file_prefix, g, "_s1_", file_suffix, ".pdf", sep=""))
      print(dynamics_plot_1)
      dev.off()
      pdf(paste(outdir, file_prefix, g, "_s2_", file_suffix, ".pdf", sep=""))
      print(dynamics_plot_2)
      dev.off()
    }
    
    else {
      dfm <- melt(df, id.vars = "time")   # melted data frame (for plotting multiple data series)
      
      dynamics_plot <- ggplot(dfm, aes(x = dfm$time, y = dfm$value, color = dfm$variable, linetype=dfm$variable)) + 
        geom_point(size = 1.5) + geom_line(show.legend = FALSE)
      dynamics_plot <- plot_format(dynamics_plot, paste(title_prefix, g, sep=""), 
                                   "time [min]", "ratio w.r.t total read counts", legend_title = "time course", 
                                   colour_palette = colour_palette) + scale_linetype_manual(values=line_type)
      
      pdf(paste(outdir, file_prefix, g, "_", file_suffix, ".pdf", sep=""))
      print(dynamics_plot)
      dev.off()
    }
    
  }
  
  # returning
  return(gene_names)
}

halflife_estimation_comparison_timecurves <- function(original_table, estimation_tables, time_series, genes,
                                                      outdir, file_prefix, colour_palette,
                                                      estim_des = NULL, title_prefix = "") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # original_table:       data frame containing measurement data
  # estimation_tables:    list of data frames containing estimation data
  # time_series:          vector of time points at which measurements were taken (not including time point 0)
  # genes:                EITHER number of genomic regions to select randomly and plot time curves for,
  #                       OR vector of genomic region names for which to plot time curves
  
  # outdir:           output directory to store plots to
  # file_prefix:      file name prefix of the plots, will be extended by the genes' names
  # colour_palette:   vector of colours to use for plotting
  # estim_des:        optional; vector with estimation table names (used for plotting) (default: tables will be
  #                   numbered in ascending order)
  # title_prefix:     optional; text to set as prefix for the plot titles (default: "")
  
  # *** NOTE ***
  
  # data frames storing measurement data must be sorted by gene names and time points (nuclear time series, 
  # cytosolic time series) and have the following format:
  # rows: genomic regions
  # col 1: region name, col 2: measurement description, col 3: library size, 
  # col 4: total transcript counts, col 5: labeled transcript counts, col 6: average potential conversion positions, 
  # col 7: conversion efficiency, col 8: newly synthesized ratio
  
  # data frames storing estimation data should be in the following format:
  # rows: genomic regions,
  # col 1: lambda_n + tau estimator, col 2: lambda_c estimator, col 3: sum-of-squares values
  
  # *** RETURN ***
  
  # vector of plotted genes' names
  
  # *******************************************************************************************************************
  
  n_timepoints <- length(time_series)         # number of time series measurements
  time_series <- c(0, time_series)            # adding time point 0 to time series
  n_estims <- length(estimation_tables)       # number fo estimations to compare
  
  # selecting gene names for which to plot time curves
  gene_names <- c()
  if(typeof(genes) == "double") {             # --- case: select genes randomly ---
    # all gene names stored in both sample 1 and 2
    gene_names_collection <- rownames(estimation_tables[[1]])     
    # selecting genes for which time curves are to be plotted (randomly)
    gene_names <- gene_names_collection[sample.int(length(gene_names_collection), genes)] 
  }
  else if(typeof(genes) == "character") {     # --- case: gene names are pre-defined ---
    gene_names <- genes
  }
  
  for(g in gene_names) {  # iterating through gene names 
    
    # original measurements ...........................................................................................
    
    subtable_original <- original_table[which(original_table[, "name"] == g), ]
    
    # *** measured labeled/total ratios ***
    
    lab_total_measured_nu <-      # nucleus (NA as missing value for time point 0)
      c(NA, subtable_original[1:n_timepoints, "mod"] / subtable_original[1:n_timepoints, "tot"])                      
    lab_total_measured_cy <-      # cytosol (NA as missing value for time point 0)
      c(NA, subtable_original[(n_timepoints+1):(2*n_timepoints), "mod"] / 
          subtable_original[(n_timepoints+1):(2*n_timepoints), "tot"])
    
    # estimated values ................................................................................................
    
    # lambda_n + tau estimators
    lnt_est <- sapply(1:n_estims, function(i) {
      e_table <- estimation_tables[[i]]
      return(e_table[which(rownames(e_table) == g), "lntau"])
    })
    
    # lambda_c estimators
    lc_est <- sapply(1:n_estims, function(i) {
      e_table <- estimation_tables[[i]]
      return(e_table[which(rownames(e_table) == g), "lc"])
    })
    
    # *** predicted new/total ratios ***
    # matrices (rows: time points, columns: estimation approaches)
    
    lnt_split_s1 <- lnt_est_s1 / 2  # splitting lambda_n + tau estimator into two halfs (one for each variable)
    lnt_split_s2 <- lnt_est_s2 / 2  
    
    new_total_predicted_nu <- sapply(1:n_estims, function(i) {
      ratio_n_sol(time_series, lnt_est[i] / 2, lnt_est[i] / 2)
    })
    new_total_predicted_cy <- sapply(1:n_estims, function(i) {
      ratio_c_sol(time_series, lnt_est[i] / 2, lnt_est[i] / 2, lc_est[i])
    })
    
    # plotting ........................................................................................................
    
    # initializing data to display in the plots
    
    data_rows <- (n_estims + 1) * 2                             # number of data rows
    if (is.null(estim_des)) {estim_des <- paste(1:n_estims)}    # estimation table names (ordered numbers if not defined)
    
    data_vector <- c(time_series, lab_total_measured_nu, lab_total_measured_cy, 
                     new_total_predicted_nu, new_total_predicted_cy)
    data_names <- 
      c("time", "labeled/total measured (nu)", "labeled/total measured (cy)", "new/total predicted (cy, 1)", 
        estim_des)
    
    data_matrix <- matrix(data_vector, ncol=data_rows)    # data matrix containing all original data and estimations
    df <- data.frame(data_matrix)                         # corresponding data frame
    colnames(df) <- data_names                            # setting column (data row) names
    dfm <- melt(df, id.vars = "time")                     # melted data frame (for plotting multiple data series)
    
    # generating & formatting & storing plot
      
    dynamics_plot <- ggplot(dfm, aes(x = dfm$time, y = dfm$value, color = dfm$variable)) + 
      geom_point(size = 1.5) + geom_line(show.legend = FALSE)
    dynamics_plot <- plot_format(dynamics_plot, paste(title_prefix, g, sep=""), 
                                 "time [min]", "ratio w.r.t total read counts", legend_title = "time course", 
                                 colour_palette = colour_palette)
    
    pdf(paste(outdir, file_prefix, g, "_", file_suffix, ".pdf", sep=""))
    print(dynamics_plot)
    dev.off()
  }

  # returning
  return(gene_names)
  
}

halflife_estimation_Rsquared_histograms <- function(Rsquared_table, nu_timecurve_hist, cy_timecurve_hist,
                                                    bin_number = 30, log_scale = FALSE, 
                                                    title_prefix = "", colour = "steelblue") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # Rsquared_table:       data frame storing R-squared values of parameter estimations; rows: genes, column1: R-squared
  #                       nuclear time curve, column2: R-squared cytosolic time curve
  # nu_halflife_hist:     file to store the histogram of nuclear time curves' R-squared values
  # cy_halflife_hist:     file to store the histogram of cytosolic time curves' R-squared values
  
  # bin_number:     number of bins to use for the histograms (default: 30)
  # log_scale:      set to TRUE to if R-squared values should be displayed in log10-scale (default: FALSE)
  # title_prefix:   a prefix to set for the histogram plot titles (default: "")
  # colour:         colour to use for histogram bars (default: 'steelbue')
  
  # *** RETURN ***
  
  # a list containing the nuclear and cytosolic time curves' R-squared histogram plot objects
  
  # *******************************************************************************************************************
  
  # initializing
  
  xlab_nu <- "nuclear time curves' R-squared"       # initializing nuclear x-axis label
  xlab_cy <- "cytosolic time curves' R-squared"     # initializing cytosolic x-axis label
  
  if (log_scale == TRUE) {                  # ajusting halflife table and x-axis labels if log-scale is applied
    Rsquared_table[, "nu_timecurve"] <- 
      log10(Rsquared_table[, "nu_timecurve"])                 # log10 scaling of nuclear time curves' R-squared
    Rsquared_table[, "cy_timecurve"] <- 
      log10(Rsquared_table[, "cy_timecurve"])                 # log10 scaling of cytosolic time curves' R-squared
    xlab_nu <- "log10(nuclear time curves' R-squared)"        # adjusted nuclear x-axis label
    xlab_cy <- "log10(cytosolic time curves' R-squared)"      # adjusted cytosolic x-axis label
  }
  
  # generating and formatting plots
  
  nu_hist <- ggplot(data=Rsquared_table, aes(Rsquared_table$nu_timecurve)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # nuclear time curves' R-squared histogram
  cy_hist <- ggplot(data=Rsquared_table, aes(Rsquared_table$cy_timecurve)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # cytosolic time curves' R-squared histogram
  
  nu_hist <- plot_format(nu_hist, paste(title_prefix, "histogram of nuclear time curves' R-squared values", sep=""), 
                         xlab_nu, "gene count")
  cy_hist <- plot_format(cy_hist, paste(title_prefix, "histogram of cytosolic time curves' R-squared values", sep=""), 
                         xlab_cy, "gene count")
  
  # storing plots and returning
  
  pdf(nu_timecurve_hist)
  print(nu_hist)
  dev.off()
  pdf(cy_timecurve_hist)
  print(cy_hist)
  dev.off()
  
  return(list(nu_hist, cy_hist))
}

halflife_estimation_halflife_histograms <- function(halflife_table, nu_halflife_hist, cy_halflife_hist,
                                                    ratio_halflife_hist, bin_number = 30, log_scale = TRUE, 
                                                    title_prefix = "", colour = "steelblue") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # halflife_table:       data frame storing halflife estimations; rows: genes, column1: nuclear halflife (based on
  #                       lambda_n + tau), column2: cytosolic halflife (based on lambda_c), column3: sum of squares
  # nu_halflife_hist:     file to store the histogram of nuclear halflifes
  # cy_halflife_hist:     file to store the histogram of cytosolic halflifes
  # ratio_halflife_hist:  file to store the histogram of parameter ratios
  
  # bin_number:     number of bins to use for the histograms (default: 30)
  # log_scale:      set to TRUE to if halflife values should be displayed in log10-scale (default: TRUE)
  # title_prefix:   a prefix to set for the histogram plot titles (default: "")
  # colour:         colour to use for histogram bars (default: 'steelbue')
  
  # *** RETURN ***
  
  # a list containing the nuclear and cytosolic halflifes' histogram plot objects as well as parameter ratios'
  # histogram plot
  
  # *******************************************************************************************************************
  
  # initializing
  
  xlab_nu <- "nuclear halflife (based on lambda_n + tau)"         # initializing nuclear x-axis label
  xlab_cy <- "cytosolic halflife (based on lambda_c)"             # initializing cytosolic x-axis label
  xlab_ratio <- "parameter ratio [lambda_c / (lambda_n + tau)]"   # initializing ratio x-axis label
  
  halflife_table <- 
    cbind(halflife_table, "ratio" = c(halflife_table[, 1] / halflife_table[, 2]))   # computing parameter ratios
  colnames(halflife_table) <- c("lntau_halflife", "lc_halflife", "sos", "ratio")
  
  if (log_scale == TRUE) {                  # ajusting halflife table and x-axis labels if log-scale is applied
    halflife_table[, 1] <- 
      log10(halflife_table[, 1])                              # log10 scaling of nuclear halflife
    halflife_table[, 2] <- 
      log10(halflife_table[, 2])                              # log10 scaling of cytosolic halflife
    halflife_table[, "ratio"] <-
      log10(halflife_table[, "ratio"])
    xlab_nu <- "log10 nuclear halflife (based on lambda_n + tau)"         # adjusted nuclear x-axis label
    xlab_cy <- "log10 cytosolic halflife (based on lambda_c)"             # adjusted cytosolic x-axis label
    xlab_ratio <- "log10 parameter ratio [lambda_c / (lambda_n + tau)]"   # adjusted ratio x-axis label
  }
  
  # generating and formatting plots
  
  nu_hist <- ggplot(data=halflife_table, aes(halflife_table$lntau_halflife)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # nuclear halflifes' histogram
  cy_hist <- ggplot(data=halflife_table, aes(halflife_table$lc_halflife)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # cytosolic halflifes' histogram
  ratio_hist <- ggplot(data=halflife_table, aes(halflife_table$ratio)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # parameter ratios' histogram 
  
  nu_hist <- plot_format(nu_hist, paste(title_prefix, "histogram of nuclear halflifes", sep=""), 
                         xlab_nu, "gene count")
  cy_hist <- plot_format(cy_hist, paste(title_prefix, "histogram of cytosolic halflifes", sep=""), 
                         xlab_cy, "gene count")
  ratio_hist <- 
    plot_format(ratio_hist, 
                paste(title_prefix, "histogram of parameter ratios (lambda_c by (lambda_n + tau))", sep=""),
                xlab_ratio, "gene count")
  
  # storing plots and returning
  
  pdf(nu_halflife_hist)
  print(nu_hist)
  dev.off()
  pdf(cy_halflife_hist)
  print(cy_hist)
  dev.off()
  pdf(ratio_halflife_hist)
  print(ratio_hist)
  dev.off()
  
  return(list(nu_hist, cy_hist, ratio_hist))
}

halflife_estimation_halflife_scatter <- function(halflife_table, scatter_plot, colour_table = "", colour_by = "",
                                                 title_prefix = "", legend_title = "", colour = c("steelblue")) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # halflife_table:     data frame storing halflife estimations; rows: genes, column1: nuclear halflife (based on
  #                     lambda_n + tau), column2: cytosolic halflife (based on lambda_c), column3: sum of squares
  # scatter_plot:       file to store the scatter plot
  
  # colour_table:       data frame storing additional information (columns) on genes (rows) by which data points are 
  #                     to be coloured (default: "")
  # colour_by:          name of the 'colour_table' column by which data points (genes) are to be coloured (default: "") 
  
  # title_prefix:   a prefix to set for the scatter plot titles (default: "")
  # legend_title:   legend title (to set for data point colouring) (default: "")
  # colour:         colour palette to use for data points (default: c('steelbue'))
  
  # *** RETURN ***
  
  # the scatter plot object
  
  # *******************************************************************************************************************
  
  # initializing
  
  xlab <- "nuclear halflife (based on lambda_n + tau)"      # initializing x-axis label
  ylab <- "cytosolic halflife (based on lambda_c)"          # initializing y-axis label
  colnames(halflife_table) <- c("lntau_halflife", "lc_halflife")
  
  # initializing data points colouring
  
  colour_column <- data.frame(rep(0, nrow(halflife_table)))   # initializing colour column (all data points equal)
  rownames(colour_column) <- row.names(halflife_table)        # initializing colour column's row names
  show_legend <- FALSE                                        # variable indicating if legend is shown in plot
  if (colour_table != "") {                                   # adjusting colour column, if given
    colour_column <- data.frame(colour_table[, colour_by])      # extracting colour column from colour table
    rownames(colour_column) <- row.names(colour_table)          # setting colour column's row names
    show_legend <- TRUE
  }
  
  colnames(colour_column) <- c("colour_column")                             # setting colour column's name
  halflife_table <- merge(halflife_table, colour_column, by="row.names")    # adding to halflife data frame
  halflife_table <- halflife_table[order(halflife_table$colour_column), ]   # sorting by colouring
  
  # generating and formatting plot

  nu_cy_scatter <- ggplot(data=halflife_table, 
                          aes(x = halflife_table$lntau_halflife, y = halflife_table$lc_halflife)) + 
    geom_point(aes(colour = as.factor(halflife_table$colour_column)), 
               size = 1.5, show.legend = show_legend)                       # scatter plot
  nu_cy_scatter <- plot_format(nu_cy_scatter, paste(title_prefix, "halflife estimators for all genes", sep=""), 
                               xlab, ylab, legend_title = legend_title, colour_palette = colour)
  nu_cy_scatter <- nu_cy_scatter + scale_x_log10() + scale_y_log10()
  
  # storing plot and returning
  
  pdf(scatter_plot)
  print(nu_cy_scatter)
  dev.off()
  
  return(nu_cy_scatter)
}

halflife_estimation_mu_tau_ln_histograms <- function(mu_tau_ln_table, mu_hist, tau_hist, ln_hist, bin_number = 30,
                                                     log_scale = TRUE, title_prefix = "", colour = "steelblue") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mu_tau_ln_table:  data frame storing mu, tau and lambda_n normalized relative estimates; rows: genes, column1: 
  #                   mu estimate, column2: tau estimate, column3: mu stdev, column4: tau stdev
  # mu_hist:          file to store the histogram of mu estimates
  # tau_hist:         file to store the histogram of tau estimates
  # ln_hist:          file to store the histogram of lambda_n estimates
  
  # bin_number:     number of bins to use for the histograms (default: 30)
  # log_scale:      set to TRUE to if mu and tau estimates should be displayed in log10-scale (default: TRUE)
  # title_prefix:   a prefix to set for the histogram plot titles (default: "")
  # colour:         colour to use for histogram bars (default: 'steelbue')
  
  # *** RETURN ***
  
  # a list containing the mu, tau and lambda_n normalized relative estimates' histogram plot objects
  
  # *******************************************************************************************************************
  
  # initializing
  
  mu_label <- "mu normalized, relative estimate"        # initializing mu histogram's x-axis label
  tau_label <- "tau normalized, relative estimate"      # initializing tau histogram's x-axis label
  ln_label <- "lambda_n relative estimate"              # initializing lambda_n histogram's x-axis label
  
  if (log_scale == TRUE) {                  # ajusting estimation table and axis labels if log-scale is applied
    
    valid_indices <- intersect(             # valid indices (0 values yield -Inf and are invalid)
      intersect(which(mu_tau_ln_table[, "mu_rel"] != 0), which(mu_tau_ln_table[, "tau_rel"] != 0)),
      which(mu_tau_ln_table[, "ln_rel"] != 0))
    mu_tau_ln_table <- data.frame(             # log10 scaling of valid mu and tau estimates
      matrix(c(log10(mu_tau_ln_table[valid_indices, "mu_rel"]), log10(mu_tau_ln_table[valid_indices, "tau_rel"]),
               log10(mu_tau_ln_table[valid_indices, "ln_rel"])), ncol = 3)
    )
    colnames(mu_tau_ln_table) <- c("mu_rel", "tau_rel", "ln_rel")
    mu_label <- "log10 mu normalized, relative estimate"        # adjusted mu histogram's x-axis label
    tau_label <- "log10 tau normalized, relative estimate"      # adjusted tau histogram's x-axis label
    ln_label <- "log10 lambda_n relative estimate"              # adjusted lambda_n histogram's x-axis label
  }
  
  # generating and formatting plots
  
  mu_hist_plot <- ggplot(data=mu_tau_ln_table, aes(mu_tau_ln_table$mu_rel)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # mu normalized relative estimates' histogram
  tau_hist_plot <- ggplot(data=mu_tau_ln_table, aes(mu_tau_ln_table$tau_rel)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # tau normalized relative estimates' histogram
  ln_hist_plot <- ggplot(data=mu_tau_ln_table, aes(mu_tau_ln_table$ln_rel)) + 
    geom_histogram(bins=bin_number, col="black", fill=colour)    # lambda_n relative estimates' histogram
  
  mu_hist_plot <- plot_format(mu_hist_plot, paste(title_prefix, 
                                                  "histogram of mu normalized, relative estimates", sep=""), 
                         mu_label, "gene count")
  tau_hist_plot <- plot_format(tau_hist_plot, paste(title_prefix, 
                                                    "histogram of tau normalized, relative estimates", sep=""), 
                         tau_label, "gene count")
  ln_hist_plot <- plot_format(ln_hist_plot, paste(title_prefix, 
                                                  "histogram of lambda_n relative estimates", sep=""), 
                              ln_label, "gene count")
  
  # storing plots and returning
  
  pdf(mu_hist)
  print(mu_hist_plot)
  dev.off()
  pdf(tau_hist)
  print(tau_hist_plot)
  dev.off()
  pdf(ln_hist)
  print(ln_hist_plot)
  dev.off()
  
  return(list(mu_hist_plot, tau_hist_plot, ln_hist_plot))
}

halflife_estimation_consistency <- function(estimation_table_s1, estimation_table_s2, heatmap_plot, box_plot, hist_plot,
                                            reldiff_min = 0, reldiff_max = 1, heatmap_step = 0.01, hist_bins = 30, 
                                            title_prefix = "", 
                                            boxplot_colours = c("mediumseagreen", "seagreen1", "dodgerblue4", 
                                                                "dodgerblue", "skyblue2"), 
                                            hist_colour = "steelblue") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # estimation_table_s1:  data frame storing parameter estimations of sample 1; rows: genes, col1: lambda_n + tau, 
  #                       col2: lambda_c, col3: normalized relative mu, col4: normalized relative tau, col5: relative 
  #                       lambda_n
  # estimation_table_s2:  data frame storing parameter estimations of sample 2; rows: genes, col1: lambda_n + tau, 
  #                       col2: lambda_c, col3: normalized relative mu, col4: normalized relative tau, col5: relative 
  #                       lambda_n
  # heatmap_plot:         file to store the heatmap plot
  # box_plot:             file to store the box plot
  # hist_plot:            file prefix to store relative differences' histogram (will be extended by the corresponding
  #                       parameter names)
  
  # reldiff_min:    minimum relative difference to be displayed (all values below will automatically be converted
  #                 to this minimum value) (default: 0)
  # reldiff_max:    maximum relative difference to be displayed (all values above will automatically be converted
  #                 to this maximum value) (default: 1)
  # heatmap_step:   step size to use for relative differences' count matrix, serving as heatmap input data 
  #                 (default: 0.01) 
  # hist_bins:      number of bins to use for histograms (default: 30)
  
  # title_prefix:       a prefix to set for the scatter plot titles (default: "")
  # boxplot_colours:    colours to use for boxplot (default: "meadiumseagreen", "seagreen1")
  # hist_colour:        colour to use for histograms (default: "steelblue")
  
  # *** RETURN ***
  
  # a list containing the heatmap plot object, boxplot object, and the histogram plot objects (in the following order:
  # lmabda_n + tau, lambda_c, mu, tau, lambda_n)
  
  # *******************************************************************************************************************
  
  # initializing (computing halflife means and relative differences) ..................................................
  
  # all genes stored in both sample 1 and 2
  
  gene_names <- 
    intersect(rownames(estimation_table_s1), rownames(estimation_table_s2))   
  
  # means and relative differences
  
  lntau_means <-                                        # lambda_n + tau means
    (estimation_table_s1[gene_names, "lntau"] + estimation_table_s2[gene_names, "lntau"]) / 2   
  lc_means <-                                           # lambda_c means
    (estimation_table_s1[gene_names, "lc"] + estimation_table_s2[gene_names, "lc"]) / 2
  mu_means <-                                           # mu (relative) means
    (estimation_table_s1[gene_names, "mu_rel"] + estimation_table_s2[gene_names, "mu_rel"]) / 2
  tau_means <-                                          # tau (relative) means
    (estimation_table_s1[gene_names, "tau_rel"] + estimation_table_s2[gene_names, "tau_rel"]) / 2
  ln_means <-                                           # lambda_n (relative) means
    (estimation_table_s1[gene_names, "ln_rel"] + estimation_table_s2[gene_names, "ln_rel"]) / 2
  
  lntau_reldiff <-                        # lambda_n + tau relative differences
    abs(estimation_table_s1[gene_names, "lntau"] - estimation_table_s2[gene_names, "lntau"]) / 
    lntau_means
  lc_reldiff <-                           # lambda_c relative differences
    abs(estimation_table_s1[gene_names, "lc"] - estimation_table_s2[gene_names, "lc"]) / 
    lc_means
  mu_reldiff <-                           # mu (relative) relative differences
    abs(estimation_table_s1[gene_names, "mu_rel"] - estimation_table_s2[gene_names, "mu_rel"]) /
    mu_means
  tau_reldiff <-                          # tau (relative) relative differences
    abs(estimation_table_s1[gene_names, "tau_rel"] - estimation_table_s2[gene_names, "tau_rel"]) /
    tau_means
  ln_reldiff <-                           # ln (relative) relative differences
    abs(estimation_table_s1[gene_names, "ln_rel"] - estimation_table_s2[gene_names, "ln_rel"]) /
    ln_means

  # filtering for valid / invalid values
  
  valid_indices <- intersect(                 # NA and Inf values are invalid (only select valid indices)
    intersect(
      intersect(
    which((!is.na(lntau_reldiff)) & (abs(lntau_reldiff != Inf))), 
                which((!is.na(lc_reldiff)) & (abs(lc_reldiff != Inf)))),
       intersect(which((!is.na(mu_reldiff)) & (abs(mu_reldiff != Inf))), 
                 which((!is.na(tau_reldiff)) & (abs(tau_reldiff != Inf))))
     ), which((!is.na(ln_reldiff)) & (abs(ln_reldiff != Inf)))
  )
  
  lntau_reldiff <- lntau_reldiff[valid_indices]     # filter valid values (lambda_n + tau)
  lc_reldiff <- lc_reldiff[valid_indices]           # filter valid values (lambda_c)
  mu_reldiff <- mu_reldiff[valid_indices]           # filter valid values (mu, relative)
  tau_reldiff <- tau_reldiff[valid_indices]         # filter valid values (tau, relative)
  ln_reldiff <- ln_reldiff[valid_indices]           # filter valid values (ln, relative)
  
  # setting all values below minimum / above maximum relative difference border to border's value
  
  n_genes <- length(valid_indices)            # number of all remaining genes
  
  lntau_reldiff[which(lntau_reldiff < reldiff_min)] <- reldiff_min
  lc_reldiff[which(lc_reldiff < reldiff_min)] <- reldiff_min
  mu_reldiff[which(mu_reldiff < reldiff_min)] <- reldiff_min
  tau_reldiff[which(tau_reldiff < reldiff_min)] <- reldiff_min
  ln_reldiff[which(ln_reldiff < reldiff_min)] <- reldiff_min
  
  lntau_reldiff[which(lntau_reldiff > reldiff_max)] <- reldiff_max
  lc_reldiff[which(lc_reldiff > reldiff_max)] <- reldiff_max
  mu_reldiff[which(mu_reldiff > reldiff_max)] <- reldiff_max
  tau_reldiff[which(tau_reldiff > reldiff_max)] <- reldiff_max
  ln_reldiff[which(ln_reldiff > reldiff_max)] <- reldiff_max
  
  # creating heatmap of halflife-determining parameters' relative differences
  # (lambda_n + tau and lambda_c) .....................................................................................
  
  # creating heatmap count matrix
  
  bin_number <- round((reldiff_max - reldiff_min) / heatmap_step) + 1   # count matrix' number of bin
  reldiff_bins <- seq(reldiff_min, reldiff_max, by = heatmap_step)      # count matrix' bins to use

  heatmap_countmatrix <-                      # initialize count matrix for relative differences 
    matrix(rep(0, bin_number ** 2),           # (rows: lambda_c relative differences,  
           ncol = bin_number)                 # columns: lambda_n + tau relative differences)
  
  for(i in 1:n_genes) {                       # iterating though genes' relative differences
    lnt_diff <- lntau_reldiff[i]                                          # lambda_n + tau rel. diff.  
    lc_diff <- lc_reldiff[i]                                              # lambda_c rel. diff.
    lnt_diff_idx <- round((lnt_diff - reldiff_min) / heatmap_step) + 1    # determine count matrix index (lambda_n + tau)
    lc_diff_idx <- round((lc_diff - reldiff_min) / heatmap_step) + 1      # determine count matrix index (lambda_c)

    heatmap_countmatrix[lc_diff_idx, lnt_diff_idx] <-       # increment lambda_n + tau & lambda_c relative 
      heatmap_countmatrix[lc_diff_idx, lnt_diff_idx] + 1    # difference's count matrix entry
  }

  # generating and storing heatmap
  
  row_and_colnames <- c(paste("<=", reldiff_min, sep=""), 
                        paste(reldiff_bins[2:(bin_number-1)]), paste(">=", reldiff_max, sep=""))
  pdf(heatmap_plot)
  heatmap_reldiff <- 
    heatmap_format(heatmap_countmatrix, col_names = row_and_colnames, row_names = row_and_colnames,
                   plot_title = paste(title_prefix, "lambda_n + tau and lambda_c relative difference\n",
                                      "between sample 1 and 2 (of all genes)", sep=""), 
                   plot_xlab = "lambda_n + tau relative difference", 
                   plot_ylab = "lambda_c relative difference", key_xlab = "gene count")
  dev.off()
  
  # creating boxplot of halflife relative differences .................................................................
  
  # creating a boxplot data frame (two columns: relative difference values, parameter indicator)
  
  boxplot_df <- data.frame(
    matrix(c(lntau_reldiff, lc_reldiff, mu_reldiff, tau_reldiff, ln_reldiff,
             rep(1, n_genes), rep(2, n_genes), rep(3, n_genes), rep(4, n_genes), rep(5, n_genes)
             ), ncol = 2)
  )                                                           # creating data frame
  colnames(boxplot_df) <- c("difference", "parameter")        # setting column names
  boxplot_df$parameter <- 
    factor(boxplot_df$parameter,                              # transforming parameter indicators to data type 'factor' 
           labels = c("lambda_n + tau", "lambda_c",           # and adding parameter indicator labels
                      "mu\n(normalized,\nrelative)", "tau\n(normalized,\nrelative)", "lambda_n\n(relative)"))    

  # creating boxplot
  
  boxplot_reldiff <- ggplot(boxplot_df, 
                            aes(boxplot_df$parameter, boxplot_df$difference, fill=parameter)) +
    geom_boxplot() + scale_fill_manual(values = boxplot_colours)
  boxplot_reldiff <- plot_format(boxplot_reldiff, 
                                 paste(title_prefix, "parameter estimators' relative differences between samples"),
                                 "parameter", "relative difference between samples")
  
  # storing boxplot
  
  pdf(box_plot)
  print(boxplot_reldiff)
  dev.off()
  
  # creating histograms of relative differences .......................................................................
  
  reldiff_df <- data.frame(                   # creating relative differences' data frame
    matrix(c(lntau_reldiff, lc_reldiff, mu_reldiff, tau_reldiff, ln_reldiff), ncol = 5)
  )                                            
  colnames(reldiff_df) <- c("lntau", "lc", "mu_rel", "tau_rel", "ln_rel")   # setting column names of data frame
  
  lntau_hist <- ggplot(data=reldiff_df, aes(reldiff_df$lntau)) + 
    geom_histogram(bins=hist_bins, col="black", fill=hist_colour)    # lambda_n + tau rel. diff. histogram
  lc_hist <- ggplot(data=reldiff_df, aes(reldiff_df$lc)) + 
    geom_histogram(bins=hist_bins, col="black", fill=hist_colour)    # lambda_c rel. diff. histogram
  mu_hist <- ggplot(data=reldiff_df, aes(reldiff_df$mu_rel)) +
    geom_histogram(bins=hist_bins, col="black", fill=hist_colour)    # mu (relative) rel. diff. histogram
  tau_hist <- ggplot(data=reldiff_df, aes(reldiff_df$tau_rel)) +
    geom_histogram(bins=hist_bins, col="black", fill=hist_colour)    # tau (relative) rel. diff. histogram
  ln_hist <- ggplot(data=reldiff_df, aes(reldiff_df$ln_rel)) +
    geom_histogram(bins=hist_bins, col="black", fill=hist_colour)    # lambda_n (relative) rel. diff. histogram
  
  lntau_hist <- plot_format(lntau_hist, 
                         paste(title_prefix, "histogram of lambda_n + tau relative differences\nbetween samples", sep=""), 
                         "lambda_n + tau relative difference", "gene count")      
  lc_hist <- plot_format(lc_hist, 
                         paste(title_prefix, "histogram of lambda_c relative differences\nbetween samples", sep=""), 
                         "lambda_c relative difference", "gene count")            
  mu_hist <- plot_format(mu_hist,
                         paste(title_prefix, "histogram of mu relative differences\nbetween samples", sep=""),
                         "mu (normalized, relative) relative difference", "gene count")
  tau_hist <- plot_format(tau_hist,
                         paste(title_prefix, "histogram of tau relative differences\nbetween samples", sep=""),
                         "tau (normalized, relative) relative difference", "gene count")
  ln_hist <- plot_format(ln_hist,
                         paste(title_prefix, "histogram of lambda_n relative differences\nbetween samples", sep=""),
                         "lambda_n (relative) relative difference", "gene count")
  
  pdf(paste(hist_plot, "lntau.pdf", sep=""))
  print(lntau_hist)
  dev.off()
  pdf(paste(hist_plot, "lc.pdf", sep=""))
  print(lc_hist)
  dev.off()
  pdf(paste(hist_plot, "mu_relative.pdf", sep=""))
  print(mu_hist)
  dev.off()
  pdf(paste(hist_plot, "tau_relative.pdf", sep=""))
  print(tau_hist)
  dev.off()
  pdf(paste(hist_plot, "ln_relative.pdf", sep=""))
  print(ln_hist)
  dev.off()
  
  # returning .........................................................................................................
  
  return(list(heatmap_reldiff, boxplot_reldiff, lntau_hist, lc_hist, mu_hist, tau_hist, ln_hist))
}

### Data Handling #####################################################################################################

boxplot_dataframe <- function(params_differences_matrix, param_labels) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # params_differences_matrix:  matrix containing differences between real and estimated parameter values 
  #                             (columns: parameters, rows: repetitions)
  # param_labels:               parameter descriptions (vector)
  
  # *** RETURN ***
  
  # a data frame suitable for constructing box plots
  
  # *******************************************************************************************************************
  
  # restructuring matrix to to data frame format (two columns: difference values, parameter indicator)
  parameters <- c(1, 2, 3, 4)
  params_differences_df <- matrix(c(
    params_differences_matrix[, 1], params_differences_matrix[, 2], 
    params_differences_matrix[, 3], params_differences_matrix[, 4],
    c(sapply(1:4, function(p){ sapply(1:nrow(params_differences_matrix), function(i){parameters[p]})}))
  ), ncol = 2)
  # setting column names
  colnames(params_differences_df) <- c("difference", "parameter")
  # transforming to data frame
  params_differences_df <- data.frame(params_differences_df)
  # transforming parameter indication to data type 'factor'
  params_differences_df$parameter <- factor(                             
    params_differences_df$parameter, labels = param_labels)
  # returning
  return(params_differences_df)
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
  elvl_ratio_s1 <- (elvl_regression_s1[common_loci, "avg_elvl_cy"] / elvl_regression_s1[common_loci, "avg_elvl_nu"]) / 2
  elvl_ratio_s2 <- (elvl_regression_s1[common_loci, "avg_elvl_cy"] / elvl_regression_s2[common_loci, "avg_elvl_nu"]) / 2
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

parameter_labels <- function(sampling) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # sampling:   RNA dynamics parameter sampling; choose 'raw' for raw parameters, choose 'halflife' for halflife
  #             representation, choose "entities" for parameter combinations as occurring in the ODE
  
  # *** RETURN ***
  
  # a vector containing the RNA dynamics parameter labels
  
  # *******************************************************************************************************************
  
  mu_label <- ""
  tau_label <- ""
  lambda_n_label <- ""
  lambda_c_label <- ""
  if (sampling == "raw") {
    mu_label <- "RNA synthesis rate\n"
    tau_label <- "RNA transport rate\n"
    lambda_n_label <- "nuclear RNA degradation rate\n"
    lambda_c_label <- "cytosolic RNA degradation rate\n"
  } else if (sampling == "halflife") {
    mu_label <- "RNA synthesis rate\n"
    tau_label <- "RNA transport half life\n"
    lambda_n_label <- "nuclear RNA half life\n"
    lambda_c_label <- "cytosolic RNA half life\n"
  } else if (sampling == "entities") {
    mu_label <- "RNA nuclear steady-state\n"
    tau_label <- "RNA cytosolic steady-state proportionality factor\n"
    lambda_n_label <- "exponent 1 (w.r.t.: tau + lambda_n)\n"
    lambda_c_label <- "exponent 2 (w.r.t.: lambda_c)\n"
  }
  return(c(mu_label, tau_label, lambda_n_label, lambda_c_label))
}

grid_fill_statics <- function(mu_values, tau_values, ln_values, lc_values, sampling) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mu_values:    values to choose for synthesis rate (sampling == 'raw') (vector), 
  #               or values to choose for nuclear steady-state (sampling == 'deriv') (vector),
  #               or "static" (arbitrary, constant value for mu)
  # tau_values:   values to choose for transport rate (sampling == 'raw') (vector),
  #               or values to choose for cytosolic steady-state (sampling == 'deriv') (vector),
  #               or "static" (arbitrary, constant value for tau)
  # ln_values:    values to choose for nuclear degradation rate (sampling == 'raw') (vector),
  #               or values to choose for (lambda_n + tau) (sampling == 'deriv') (vector),
  #               or "static" (arbitrary, constant value for ln)
  # lc_values:    values to choose for cytosolic degradation rate (sampling == 'raw' and sampling == 'deriv') (vector),
  #               or "static" (arbitrary, constant value for lc)
  # sampling:     RNA dynamics parameter sampling; choose 'raw' for the RNA dynamics parameters itself, choose 'deriv'
  #               for the derived parameter set
  
  # *** RETURN ***
  
  # A list containing each RNA dynamics parameter's values (as vector) for grid creation.
  
  # This function is meant for usage within the function "grid_theta". It fills values for parameters that are meant 
  # to be static within the grid. In case of 'deriv' sampling it takes care for the re-computed raw parameter lambda_n 
  # not to become negative.
  
  # *******************************************************************************************************************
  
  # raw parameter sampling (parameters can be set to any arbitrary value) .............................................
  
  if (sampling == "raw") {
    
    if (lc_values == "static") {
      lc_values <- c(0.3)
    }
    if (tau_values == "static") {
      tau_values <- c(0.2)
    }
    if (ln_values == "static") {
      ln_values <- c(0.25)
    }
    if (mu_values == "static") {
      mu_values <- c(100)
    }
  }
  
  # derived parameters sampling .......................................................................................
  
  if (sampling == "deriv") {
    
    # counting the number of static parameters (to know how many can be chosen arbitrarily; the last static parameter
    # remaining has to be computed s.t. lambda_n cannot get negative)
    # determining max / min parameter values (needed for computing statics s.t. lambda_n remains positive)
    
    params <- list(mu_values, tau_values, ln_values, lc_values)
    statics <- 0
    params_extremes <- c(-1, -1, -1, -1)
    for(p in 1:4) {
      if (params[[p]] == "static") {statics <- statics + 1}
      else if (p != 3) {params_extremes[p] <- max(params[[p]])}
      else {params_extremes[[p]] <- min(params[[p]])}
    }
    
    # setting static lambda_c ("lc")
    
    if (lc_values == "static") {
      # lambda_c can be chosen arbitrarily
      if (statics > 1) {lc_values <- c(0.3)}
      # lambda_c has to be computed
      else {lc_values = c(0.9 * params_extremes[3] * params_extremes[1] / params_extremes[2])}
      params_extremes[4] <- lc_values[1]  # storing lambda_c's max parameter value (needed for following computations)
      statics <- statics - 1              # static lambda_c is set now (one less remaining static parameter)
    }
    
    # setting static cytosolic steady-state ("tau")
    
    if (tau_values == "static") {
      # cytosolic steady-state can be chosen arbitrarily
      if (statics > 1) {tau_values <- c(50)}
      # cytosolic steady-state has to be computed
      else {tau_values <- c(0.9 * params_extremes[1] * params_extremes[3] / params_extremes[4])}
      params_extremes[2] <- tau_values[1] # storing cytosolic steady state's max parameter value (needed for following 
      #                                     computations)
      statics <- statics - 1              # static cytosolic steady-state is set now (one less remaining static 
      #                                     parameter)
    }
    
    # setting static lambda_n + tau ("ln")
    
    if (ln_values == "static") {
      # lambda_n + tau can be chosen arbitrarily
      if (statics > 1) {ln_values <- c(0.5)}
      # lambda_n + tau has to be computed
      else {ln_values <- c(1.1 * params_extremes[2] / params_extremes[1] * params_extremes[4])}
      params_extremes[3] <- lc_values[1]  # storing lambda_c + tau's max parameter value (for following computations)
      statics <- statics - 1              # static lambda_c + tau is set now (one less remaining static parameter)
    }
    
    # seeting static nuclear steady-state ("mu")
    
    if (mu_values == "static") {
      # lambda_n + tau has to be computed (all other parameters are already set)
      mu_values <- c(1.1 * params_extremes[2] / params_extremes[3] * params_extremes[4])
    }
  }
  
  # returning .........................................................................................................
  
  return(list(mu_values, tau_values, ln_values, lc_values))
}

grid_theta <- function(mu_values, tau_values, ln_values, lc_values, sampling = "raw") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mu_values:  values to choose for synthesis rate (sampling == 'raw'; default) (vector), 
  #             or values to choose for nuclear steady-state (sampling == 'deriv') (vector),
  #             or "static" (arbitrary, constant value for mu)
  # tau_values: values to choose for transport rate (sampling == 'raw'; default) (vector),
  #             or values to choose for cytosolic steady-state (sampling == 'deriv') (vector),
  #             or "static" (arbitrary, constant value for tau)
  # ln_values:  values to choose for nuclear degradation rate (sampling == 'raw'; default) (vector),
  #             or values to choose for (lambda_n + tau) (sampling == 'deriv') (vector),
  #             or "static" (arbitrary, constant value for ln)
  # lc_values:  values to choose for cytosolic degradation rate (sampling == 'raw' and sampling == 'deriv'; default)
  #             (vector), or "static" (arbitrary, constant value for lc)
  # sampling:   optional; RNA dynamics parameter sampling; choose 'raw' for the RNA dynamics parameters itself, 
  #             choose 'deriv' for the derived parameter set (default: raw)
  
  # *** RETURN ***
  
  # a matrix containing (in the following order, as columns): a vector of mu grid values, a vector of tau grid values,
  # a vector of lambda_n grid values, a vector of lambda_c grid values
  
  # This function generates a RNA dynamics parameter grid to be applied to the multi RNA dynamics simulator or the
  # parameter estimation testing routine.
  
  # *******************************************************************************************************************
  
  filled_statics <- grid_fill_statics(mu_values, tau_values, ln_values, lc_values, sampling)  # filling static values
  mu_values <- filled_statics[[1]]
  tau_values <- filled_statics[[2]]
  ln_values <- filled_statics[[3]]
  lc_values <- filled_statics[[4]]
  
  mu_length = length(mu_values)
  tau_length = length(tau_values)
  ln_length = length(ln_values)
  lc_length = length(lc_values)
  
  if (sampling == "raw") {              # parameter computation for raw sampling
    
    get_mu <- function(m, n_dummu) {return(m)}
    get_tau <- function(t, m_dummu, c_dummu) {return(t)}
    get_ln <- function(n, t_dummu) {return(n)}
    get_lc <- function(c) {return(c)}
    
  } else if (sampling == "deriv") {     # parameter computation for deriv sampling
    
    get_mu <- function(m_mod, n_mod) {return(m_mod * n_mod)}
    get_tau <- function(t_mod, m_mod, c) {return(t_mod / m_mod * c)}
    get_ln <- function(n_mod, t) {return(n_mod - t)}
    get_lc <- function(c) {return(c)}
    
  } else {
    
    stop("No valid argument given for parameter 'sampling'. Please choose from 'raw', or 'deriv'.")
    
  }
  
  grid_idx <- 0    # variable tracking current index within the grid
  
  lc_grid <- c(sapply(1:mu_length, function(m) {      # filling lambda_c within grid
    sapply(1:tau_length, function(t) {
      sapply(1:ln_length, function(n) {
        sapply(1:lc_length, function(c) {
          get_lc(lc_values[c])
        })
      })
    })
  }))
  tau_grid <- c(sapply(1:mu_length, function(m) {     # filling tau within grid
    sapply(1:tau_length, function(t) {
      sapply(1:ln_length, function(n) {
        sapply(1:lc_length, function(c) {
          grid_idx <- grid_idx + 1
          get_tau(tau_values[t], mu_values[m], lc_grid[grid_idx])
        })
      })
    })
  }))
  grid_idx <- 0    # resetting current grid index
  ln_grid <- c(sapply(1:mu_length, function(m) {      # filling lambda_n within grid
    sapply(1:tau_length, function(t) {
      sapply(1:ln_length, function(n) {
        sapply(1:lc_length, function(c) {
          grid_idx <- grid_idx + 1
          get_ln(ln_values[n], tau_grid[grid_idx])
        })
      })
    })
  }))
  grid_idx <- 0    # resetting current grid index
  mu_grid <- c(sapply(1:mu_length, function(m) {      # filling mu within grid
    sapply(1:tau_length, function(t) {
      sapply(1:ln_length, function(n) {
        sapply(1:lc_length, function(c) {
          grid_idx <- grid_idx + 1
          get_mu(mu_values[m], ln_values[n])
        })
      })
    })
  }))
  
  return(matrix(c(mu_grid, tau_grid, ln_grid, lc_grid), ncol = 4))
}

transform_params <- function(target, borders_1, borders_2, borders_3, lc_borders) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # target:       transformation target; 'derived' for transforming raw to derived parameter space, or 'raw' for
  #               transforming derived to raw parameter space
  # borders_1:    lower and upper parameter space borders for mu (target == 'derived') or lambda_n + tau
  #               (target == 'raw') (vector)
  # borders_2:    lower and upper parameter space borders for tau (target == 'derived') or steady_n (target == 'raw')
  #               (vector)
  # borders_3:    lower and upper parameter space borders for lambda_n (target == 'derived') steady_c 
  #               (target == 'raw') (vector)
  # lc_borders:   vector; lower and upper parameter space borders for lambda_c
  
  # *** RETURN ***
  
  # a vector containing (in the following order, depending on transformation mode): lambda_n + tau or mu lower 
  # border, lambda_n + tau or mu upper border, steady_n or tau lower border, steady_n or tau upper border,
  # steady_c or lambda_n lower border, steady_c or lambda_n upper border
  
  # This function computes the borders of RNA dynamics derived parameter set's space, based on the borders of raw
  # parameter set's space, or the other way round, respectively.
  
  # NOTE: within the 4-dim space of the derived parameter set, some regions may be undefined (in contrast to the raw
  #       parameters, derived parameters are not independent of each other); that is the case where:
  #       1) the underlying raw parameters exceed their natural borders (e.g., lambda_n becomes negative)
  #       2) the underlying raw parameters exceed user-defined borders
  
  # *******************************************************************************************************************
  
  # transformation: raw -> derived
  
  if (target == 'derived') {
    
    mu_min <- borders_1[1]
    mu_max <- borders_1[2]
    tau_min <- borders_2[1]
    tau_max <- borders_2[2]
    ln_min <- borders_3[1]
    ln_max <- borders_3[2]
    
    lntau_lower <- ln_min + tau_min
    lntau_upper <- ln_max + tau_max
    steadyn_lower  <- mu_min / lntau_upper
    steadyn_upper <- mu_max / lntau_lower
    steadyc_lower <- mu_min / lc_borders[2] * tau_min / (tau_min + ln_max)
    steadyc_upper <- mu_max / lc_borders[1] * tau_max / (tau_max + ln_min)
    
    return(c(lntau_lower, lntau_upper, steadyn_lower, steadyn_upper, steadyc_lower, steadyc_upper))
  }
  
  # transformation: derived -> raw (don't question where the tau-formulas come from, they are the result of four terms
  # mixed by subtractions and substitutions of variables ^^)
  
  else if (target == 'raw') {
    
    lntau_min <- borders_1[1]
    lntau_max <- borders_1[2]
    lntau_total = lntau_min + lntau_max
    sn_min <- borders_2[1]
    sn_max <- borders_2[2]
    sc_min <- borders_3[1]
    sc_max <- borders_3[2]
    lc_min <- lc_borders[1]
    lc_max <- lc_borders[2]
    
    mu_lower <- sn_min * lntau_max
    mu_upper <- sn_max * lntau_min
    tau_lower <- (sc_max * lntau_total - mu_upper / lc_min * lntau_max) /
      (mu_lower / lc_max * sc_max / sc_min - mu_lower * mu_upper / (lc_max * lc_min * sc_min) + mu_upper / lc_min)
    tau_upper <- (sc_min * lntau_total - mu_lower / lc_max * lntau_min) /
      (mu_upper / lc_min * sc_min / sc_max - mu_upper * mu_lower / (lc_min * lc_max * sc_max) + mu_lower / lc_max)
    ln_lower <- lntau_min - tau_lower
    ln_upper <- lntau_max - tau_upper
    
    return(c(mu_lower, mu_upper, tau_lower, tau_upper, ln_lower, ln_upper))
  }
  
  # in case no valid transformation target has ben parsed, NaN is returned
  
  else {return(NaN)}
}

sum_likelihoods <- function(nucleus_probs, cytosol_probs, repetitions, grid_size) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # nucleus_probs:      likelihoods based on nuclear transcript counts; list containing one matrix per repetition, 
  #                     each matrix storing parameter setting likelihoods (rows: time, columns: parameter setting) 
  # cytosol_probs:      likelihoods based on cytosol transcript counts; list containing one matrix per repetition,
  #                     each matrix storing parameter setting likelihoods (rows: time, columns: parameter setting)
  # repetitions:        number of repetitions
  # grid_size:          number of parameter settings
  
  # *** RETURN ***
  
  # a list containing (in the following order):
  #   a list of one vector per repetition, each vector storing the parameter settings' likelihoods summed over time
  #     (for nuclear transcript counts)
  #   a list of one vector per repetition, each vector storing the parameter settings' likelihoods summed over time
  #     (for cytosol transcript counts)
  #   a vector storing the parameter settings' likelihoods summed over time and repetitions (for nuclear transcript 
  #     counts)
  #   a vector storing the parameter settings' likelihoods summed over time and repetitions (for cytosol transcript 
  #     counts)
  #   a list of one vector per repetition, each vector storing the parameter seetings' likelihoods summed over time
  #     and over the two compartments nucleus and cytosol
  #   a vector storing the parameter settings' likelihoods summed over time, repetitions and the two compartments
  #     nucleus and cytosol
  
  # *******************************************************************************************************************
  
  # summing over time points (leaving out t = 0)
  probs_n <- sapply(1:repetitions, function(i) { list(
    sapply(1:grid_size, function(g) {        
      sum(nucleus_probs[[i]][-c(1),  g])})  
  )})
  probs_c <- sapply(1:repetitions, function(i) { list(
    sapply(1:grid_size, function(g) {
      sum(cytosol_probs[[i]][-c(1), g])})
  )})
  
  # restructuring probabilities as matrix (columns: grid parameter sets, rows: repetitions)
  summed_probs_n <- sapply(1:grid_size, function(g) {   
    sapply(1:repetitions, function(i) {probs_n[[i]][g]})
  })
  summed_probs_c <- sapply(1:grid_size, function(g) {
    sapply(1:repetitions, function(i) {probs_c[[i]][g]})
  })
  
  # summing over repetitions (already done if repetitions == 1)
  if (repetitions > 1) {
    summed_probs_n <- colSums(summed_probs_n)
    summed_probs_c <- colSums(summed_probs_c)
  }
  
  # summing over compartments
  overall_resolved <- sapply(1:repetitions, function(i) { list(probs_n[[i]] + probs_c[[i]])})
  overall_summed <- summed_probs_n + summed_probs_c
  
  return(list(probs_n, probs_c, summed_probs_n, summed_probs_c, overall_resolved, overall_summed))
}

### Data Processing ###################################################################################################

summarytable_intersection <- function(summary_tables, direction = "description") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summary_tables:   list of summary tables (data frames; rows: genomic regions, col1: region name,
  #                   col2: data description, col3: library size, col4: total counts, col5: labeled counts,
  #                   col6: conversion efficiency, col7: newly synthesized ratio)
  # direction:        optional; filtering process direction (either "description" or "region", see this function's 
  #                   explanation for more details) (default: "description")
  
  # *** RETURN ***
  
  # a list containing two vectors: the first vector storing genomic regions filtered, the second vector storing 
  # data descriptions filtered
  
  # This function filters genomic regions and data descriptions intersecting in a set of input summary tables.
  # The filtering process is specified by parameter 'direction':
  
  # direction == "description":
  # Data descriptions present in all of the input summary tables are selected. Then, genomic regions occurring for
  # all of the selected data descriptions in all of the tables are selected aswell.
  
  # directions == "region":
  # Genomic regions present in all of the input summary tables are selected. Then, data descriptions occurring for
  # all of the selected regions in all of the tables are selected aswell.
  
  # *******************************************************************************************************************
  
  stab_number = length(summary_tables)      # number of tables
  regions_filtered <- NULL                  # initialize filtered regions
  descriptions_filtered <- NULL             # initialize filtered descriptions
  
  # direction == "description"
  
  if (direction == "description") {
    
    # tables' data descriptions (list of vectors, one vector for each matrix)
    stab_descriptions <- sapply(1:stab_number, function(i) { list (
      summary_matrices[[i]][, "des"]
    )})                                       
    
    # filtering data descriptions (intersection)
    descriptions_filtered <- Reduce(intersect, stab_descriptions)   
    
    # for each matrix and each data description, getting genomic regions
    regions_filtered <- c(sapply(descriptions_filtered, function(d) {     # iterating through descriptions
      sapply(1:stab_number, function(i) {                                   # iterating through matrices
        subset(summary_tables[[i]], summary_tables[[i]][, "des"] == d)[, "name"]  # all description d entries' regions
      })
    }))
    
    # filtering genomic regions (intersection)
    regions_filtered <- Reduce(intersect, regions_filtered)
  }
  
  # direction == "region"
  
  if (direction == "region") {
    
    # matrices' genomic regions (list of vectors, one vector for each matrix)
    stab_regions <- sapply(1:stab_number, function(i) { list (
      summary_table[[i]][, "name"]
    )})                                       
    
    # filtering genomic regions (intersection)
    regions_filtered <- Reduce(intersect, stab_regions)   
    
    # for each matrix and each genomic region, getting data descriptions
    descriptions_filtered <- c(sapply(regions_filtered, function(r) {       # iterating through regions
      sapply(1:stab_number, function(i) {                                     # iterating through matrices
        subset(summary_tables[[i]], summary_tables[[i]][, "name"] == r)[, "des"]  # all region r entries' descriptions
      })
    }))
    
    # filtering data descriptions (intersection)
    descriptions_filtered <- Reduce(intersect, descriptions_filtered)  
  }
  
  # returning
  
  return(list(regions_filtered, descriptions_filtered))
}

summarytable_ratios <- function(mode, summary_table, n_timepoints) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mode:               mode to be executed (type of ratios to be computed); use "counts" to compute modified/total
  #                     and total/library transcript count ratios, use "compartment" to compute nucleus/cytosol and 
  #                     inverse transcript counts
  # summary_table:      summary table (data frame; rows: genomic regions, col1: region name, col2: data description, 
  #                     col3: library size, col4: total read counts, col5: labeled read counts, 
  #                     col6: conversion efficiency, col7: newly synthesized ratio)
  # n_timepoints:       number of time points measured
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)
  
  # *** RETURN ***
  
  # a matrix containing transcript count ratios (rows: genomic regions, columns (first half): modified/total ratios or
  # cytosol/nucleus total ratios, columns (second half): total/library_size ratios or nucleus/cytosol total ratios)
  
  # This function computes either modified/total and total/library_size transcript count ratios or cytosol/nucleus and
  # nucleus/cytosol total count ratios for an input summary table and stores them in a matrix. 
  
  # *******************************************************************************************************************
  
  gene_names <- unique(summary_table[, "name"])
  n_genes <- length(gene_names)
  data_des <- unique(summary_table[, "des"])
  
  if (mode == "counts") {

    ratios_1 <- summary_table[, "mod"] / summary_table[, "tot"]               # mod/total ratios
    ratios_1 <- matrix(ratios_1, ncol = 2 * n_timepoints, byrow = TRUE)       # matrix for mod/total ratios
    colnames(ratios_1) <- 
      sapply(data_des, function(d) {paste(d, "mod/total", sep = "_")})        # matrix' column names
    
    ratios_2 <- summary_table[, "tot"] / summary_table[, "lib"]               # total/lib_size ratios
    ratios_2 <- matrix(ratios_2, ncol = 2 * n_timepoints, byrow = TRUE)       # matrix for total/lib_size ratios
    colnames(ratios_2) <- 
      sapply(data_des, function(d) {paste(d, "total/libsize", sep = "_")})    # matrix' column names
  }
  
  else if (mode == "compartment") {
    
    nu_indices <- sapply(1:n_timepoints, function(i) {            # indices of nuclear measurements
      (c(1:n_genes) - 1) * 2 * n_timepoints + i
    })
    print(nu_indices)
    cy_indices <- sapply(1:n_timepoints, function(i) {            # indices of cytosolic measurements
      (c(1:n_genes) - 1) * 2 * n_timepoints + n_timepoints + i
    })
    print(cy_indices)
    ratios_1 <- matrix(
      summary_table[cy_indices, "tot"] / summary_table[nu_indices, "tot"], 
      ncol = n_timepoints)                                  # matrix for cy/nu ratios
    ratios_2 <- matrix(
      summary_table[nu_indices, "tot"] / summary_table[cy_indices, "tot"],
      ncol = n_timepoints)                                  # matrix for nu/cy ratios
    colnames(ratios_1) <- sapply(data_des[1:n_timepoints], function(d) {
      paste(d, "cy_nu", sep = "_")})                                      # matrix' column names
    colnames(ratios_2) <- sapply(data_des[(n_timepoints+1):(n_timepoints*2)], function(d) {
      paste(d, "nu_cy", sep = "_")})                                      # matrix' column names
    print(ratios_1)
  }
  
  all_mat <- cbind(ratios_1, ratios_2)        # matrix for all ratios
  all_df <- (all_mat)                         # data frame for all ratios
  rownames(all_mat) <- gene_names             # data frame row names (gene names)
  
  return(all_mat)   # returning ratios' data frame
}

estimation_ratios <- function(lntau_est, lc_est, convpos, time_series,
                              conveff, fn_error_nu, fn_error_cy, fp_error_nu, fp_error_cy) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # lntau_est:          lambda_n + tau estimator
  # lc_est:             lambda_c estimator
  # convpos:            average number of potential conversion positions in the transcript
  # time_series:        time series for which estimated modified/total ratios should be computed 
  
  # conveff:            vector of conversion efficiencies (at single measurement time points); probability by which 
  #                     a single base is labeled
  # fn_error_nu:        nuclear false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate)
  # fn_error_cy:        cytosolic false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate)
  # fp_error_nu:        nuclear false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment
  # fp_error_cy:        cytosolic false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment
  
  # *** RETURN ***
  
  # a list containing the nuclear and cytosolic vector of estimated modified/total transcript count ratios (resolved
  # by time)
  
  # *******************************************************************************************************************
  
  p_mod_nu <- 1 - 
    (1 - conveff * (1 - fn_error_nu)) ** convpos    # modification probability of nuclear newly synthesized reads
  p_mod_cy <- 1 - 
    (1 - conveff * (1 - fn_error_cy)) ** convpos    # modification probability of cytosolic newly synthesized reads
  p_un_nu <- (1 - fp_error_nu) ** convpos           # non-modfication probability of nuclear pre-existing reads
  p_un_cy <- (1 - fp_error_cy) ** convpos           # non-modfication probability of cytosolic pre-existing reads
  
  # computing new/total transcript count ratios according to parameter estimators
  lnt_split <- lntau_est / 2  # splitting lambda_n + tau estimator into two halfs (one for each variable)
  n_ratio_new_total <- ratio_n_sol(time_series, lnt_split, lnt_split)           
  c_ratio_new_total <- ratio_c_sol(time_series, lnt_split, lnt_split, lc_est)

  # computing modified/total transcript count ratios according to modification probabilities
  n_ratio_mod_total <- n_ratio_new_total * p_mod_nu + (1 - n_ratio_new_total) * (1 - p_un_nu)
  c_ratio_mod_total <- c_ratio_new_total * p_mod_cy + (1 - c_ratio_new_total) * (1 - p_un_cy)

  # returning
  return(list(n_ratio_mod_total, c_ratio_mod_total))
  
}

maxlike_marginalise_discrete <- function(overall_probs, mu_grid, tau_grid, lambda_n_grid, lambda_c_grid) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # overall_probs:  vector containing likelihoods of grid parameter settings
  # x_grid:         RNA dynamics parameter x's values used for grid search (raw values, not containing any duplicates)
  
  # *** RETURN ***
  
  # a list containing the matrices (in the following order): mu-tau marginal likelihoods, mu-lambda_n marginal
  # likelihoods, mu-lambda_c marginal likelihoods, tau-lambda_n marginal likelihoods, tau-lambda_c marginal
  # likelihoods, lambda_n-lambda_c marginal likelihoods
  
  # This function marginalizes the RNA parameters' likelihoods over a grid of likelihood values (discrete values) in a
  # pairwise fashion. For each pairwise parameter combination, marginal likelihoods are computed (for the two 
  # considered parameters, over the two not considered parameters).
  
  # *******************************************************************************************************************
  
  # initializations ...................................................................................................
  
  # computing number of grid values for each parameter
  mu_values = length(mu_grid)
  tau_values = length(tau_grid)
  lambda_n_values = length(lambda_n_grid)
  lambda_c_values = length(lambda_c_grid)
  # reminder: grid structure
  # (mu * mu_values, tau * tau_values, lambda_n * lambda_n_values, lambda_c * lambda_c_values)
  tau_subsamples = lambda_n_values * lambda_c_values    # number of subsamples under a given mu-tau combination
  mu_subsamples = tau_subsamples * tau_values           # number of subsamples under a given mu value
  tau_ln_combinations = tau_values * lambda_n_values    # number of tau-lambda_n combinations
  
  # mu and tau ........................................................................................................
  
  # summing over lambda_n and lambda_c (marginalisation over mu and tau)
  probs_m_t <- c(sapply(0:(mu_values - 1), function(m) {          # iterate through mu values
    sapply(0:(tau_values - 1), function(t) {                        # iterate through tau values
      index_start = (t * tau_subsamples + 1) + (m * mu_subsamples)
      index_end = index_start + tau_subsamples - 1
      sum(overall_probs[index_start:index_end])      # summing all grid likelihoods with same mu-tau combination
    })
  }))
  # bringing vector of probabilities into matrix shape (mu -> tau)
  probs_mt_matrix <- matrix(probs_m_t, ncol = mu_values)
  
  # mu and lambda_n ...................................................................................................
  
  # summing over tau and lambda_c (marginalisation over mu and lambda_n)
  probs_m_n <- c(sapply(0:(mu_values - 1), function(m) {
    sapply(0:(lambda_n_values - 1), function(n) {
      sum(sapply(0:(tau_values - 1), function(t) {  # mu-lambda_n combinations repeat with tau and lambda_c
        idx_start = (m * mu_subsamples + 1) + (n * lambda_c_values) + 
          (t * tau_subsamples)                      # starting index of grid subarea with same mu-lambda_n combination
        idx_end = idx_start + lambda_c_values - 1   # ending index
        sum(overall_probs[idx_start:idx_end])       # summing over this area
      }))
    })
  }))
  # bringing vector of probabilities into matrix shape (mu -> lambda_n)
  probs_mn_matrix <- matrix(probs_m_n, ncol = mu_values)
  
  # mu and lambda_c ...................................................................................................
  
  # summing over tau and lambda_n (marginalisation over mu and lambda_c)
  probs_m_c <- c(sapply(0:(mu_values - 1), function(m) {
    sapply(0:(lambda_c_values - 1), function(c) {
      sum(                                            # summing over all grid indices with same mu-lambda_c values
        c(sapply(0:(tau_ln_combinations - 1), function(b) {  # (lambda_c repeats with each tau-lambda_n combination)
          overall_probs[(m * mu_subsamples + 1) + c + (b * lambda_c_values)]
        }))
      )
    })
  }))
  # bringing vector of probabilities into matrix shape (mu -> lambda_c)
  probs_mc_matrix <- matrix(probs_m_c, ncol = mu_values)
  
  # tau and lambda_n ..................................................................................................
  
  # summing over mu and lambda_c (marginalisation over tau and lambda_n)
  probs_t_n <- c(sapply(0:(tau_values - 1), function(t) {
    sapply(0:(lambda_n_values - 1), function(n) {
      sum(                                            # summing over all grid areas with same tau and lambda_n values
        c(sapply(0:(mu_values - 1), function(m) {       # (areas repeat with mu values)
          idx_start = 1 + (t * tau_subsamples) + (n * lambda_c_values) + (m * mu_subsamples)
          overall_probs[idx_start:(idx_start + lambda_c_values - 1)]
        }))
      )
    })
  }))
  # bringing vector of probabilities into matrix shape (tau -> lambda_n)
  probs_tn_matrix <- matrix(probs_t_n, ncol = tau_values)
  
  # tau and lambda_c ..................................................................................................
  
  # summing over mu and lambda_n (marginalisation over tau and lambda_c)
  probs_t_c <- c(sapply(0:(tau_values - 1), function(t) {
    c(sapply(0:(lambda_c_values - 1), function(c) {
      seed_index = 1 + (t * tau_subsamples) + c   # initial index where a certain tau-lambda_c combination is located
      sum(c(                                      # summing over repeats of initial index (tau-lambda_c combination)
        sapply(0:(mu_values - 1), function(m) {            # repeats w.r.t. mu
          sapply(0:(lambda_n_values - 1), function(n) {      # repeats w.r.t. lambda_n
            overall_probs[seed_index + m * mu_subsamples + n * lambda_c_values]
          })
        })
      ))
    }))
  }))
  # bringing vector of probabilities into matrix shape (tau -> lambda_n)
  probs_tc_matrix <- matrix(probs_t_c, ncol = tau_values)
  
  # lambda_n and lambda_c .............................................................................................
  
  # summing over mu and tau (marginalisation over lambda_n and lambda_c)
  probs_n_c <- c(sapply(1:tau_subsamples, function(i) {  # initial index of a lambda_n-lambda_c combination
    # within the grid
    sum(c(                                         # summing over repeats of the lambda_n-lambda_c combinations
      sapply(0:(mu_values - 1), function(m) {        # repeats w.r.t. mu
        sapply(0:(tau_values - 1), function(t) {     # repeats w.r.t. tau
          overall_probs[i + m * mu_subsamples + t * tau_subsamples]
        })
      })
    ))
  }))
  # bringing vector of probabilities into matrix shape (lambda_n -> lambda_c)
  probs_nc_matrix <- matrix(probs_n_c, ncol = lambda_n_values)
  
  # returning .........................................................................................................
  
  return(list(probs_mt_matrix, probs_mn_matrix, probs_mc_matrix, probs_tn_matrix, probs_tc_matrix, probs_nc_matrix))
}

maxlike_marginalise_continuous <- function(density_function, param_index, 
                                           d1_borders, d2_borders, d3_borders, d4_borders) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # density_function:   the function to be marginalized; function argument should be a vector of length 4 holding the
  #                     values of the 4 dynamics parameters (4 dimensions), respectively
  # param_index:        density function's input vector index of the parameter (number of dimension) for which the 
  #                     marginal distribution is to be computed
  # d1_border:          first parameter (first dimension) lower and upper border (vector)
  # d2_border:          second prameter (second dimension) lower and upper border (vector)
  # d3_border:          third parameter (third dimension) lower and upper border (vector)
  # d4_border:          fourth parameter (fourth dimension) lower and upper border (vector)
  
  # *** NOTE ***
  
  # the dimension borders for the parameter to compute the marginal distribution for are not needed and ignored
  
  # *** RETURN ***
  
  # the marginal distribution function of the chosen parameter (dimension)
  
  # This function constructs a function which computes the marginal density of one RNA dynamics parameter, given the
  # joint probability density of all parameters. The constructed function marginalizes the joint density for the chosen 
  # parameter (chosen dimension), over the other three parameters (dimensions) using multi-dimensional integration.
  
  # *******************************************************************************************************************
  
  # storing dimension borders of parameters
  
  borders = list(d1_borders, d2_borders, d3_borders, d4_borders)  # all borders
  borders_marginalisation = borders[-param_index]                 # borders of parameters to be marginalised over
  lower_1 = borders_marginalisation[[1]][1]
  upper_1 = borders_marginalisation[[1]][2]
  lower_2 = borders_marginalisation[[2]][1]
  upper_2 = borders_marginalisation[[2]][2]
  lower_3 = borders_marginalisation[[3]][1]
  upper_3 = borders_marginalisation[[3]][2]
  
  # defining a function computing the marginal density at a specific value of the chosen parameter
  
  marginal_density <- function(param_val) {
    # to return the marginal density at one specific value, a function is needed which can be integrated over the
    # remaining three parameters (i.e. the 3D density distribution at that specific value)
    density_3d <- function(args) {    # this function takes a vector (args) of length 3 (the parameters to be 
      # marginalised over)
      args_4d = append(args, c(param_val), param_index - 1)   # including fixed value of the chosen parameter
      return(density_function(args_4d))                       # returning density at the 4 parameter values
    }
    
    # integrating over the 3 dimensions
    density_value <- cuhre(3, 1, density_3d, lower = c(lower_1, lower_2, lower_3), 
                           upper = c(upper_1, upper_2, upper_3), max.eval = 1e8, flags = list(verbose = 0))
    # returning marginalisation at the specific parameter value
    return(density_value[[5]])
  }
  
  # returning marginal distribution
  
  return(marginal_density)
  
}

eval_simulation <- function(sim_data, mu, tau, lambda_n, lambda_c,
                            mu_borders, tau_borders, ln_borders, lc_borders, repetitions,
                            plot_margins = FALSE, plot_derived = FALSE, plot_boxplot = FALSE) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # sim_data:     simulation data; a nested list containing the count data list, the density functions' list and the
  #               re-estimated parameters' list
  # mu:           RNA synthesis rate of the transcript
  # tau:          RNA transport rate of the transcript
  # lambda_n:     nuclear degradation rate of the transcript
  # lambda_c:     cytosolic degradation rate of the transcript
  # mu_borders:   (minimum, maximum) mu value used for re-estimation (vector)
  # tau_borders:  (minimum, maximum) tau value for re-estimation (vector)
  # ln_borders:   (minimum, maximum) lambda_n value for re-estimation (vector)
  # lc_borders:   (minimum, maximum) lambda_c value for re-estimation (vector)
  # repetitions:  number of simulation repetitions
  
  # plot_margins:   optional; parameter indicating whether marginal likelihood distributions are to be plotted 
  #                 (default: FALSE)
  # plot_derived:   optional; parameter indicating whether derived parameters' probability distributions are to be
  #                 plotted (default: FALSE)
  # plot_boxplot:   parameter indicating whether boxplots showing parameter estimation accuracy (estimated / true) are
  #                 to be plotted (default: FALSE)
  
  # *** RETURN ***
  
  # a list containing the following plot objects (in the following order): 
  # mu likelihood density, tau likelihood density, lambda_n likelihood density, lambda_c likelihood density, 
  # lambda_n + tau probability distribution, lambda_n + tau and lambda_c joint probability distribution,
  # steady_n probability distribution, steady_c probability distribution, parameter estimation accuracy / error
  # plot objects not having been created are returned as NULL
  
  # *******************************************************************************************************************
  
  count_data <- sim_data[[1]]   # count data
  densities <- sim_data[[2]]    # probability / density functions
  optima <- sim_data[[3]]       # re-estimated parameters
  
  repetition_labels <- sapply(1:repetitions, function(i) {  # labels for different repetitions
    paste("#r", i, sep = "")
  })
  # derived parameters' borders
  dparams_borders <- transform_params("derived", mu_borders, tau_borders, ln_borders, lc_borders)
  
  # initializing all plot objects with NULL
  mu_density_plot <- NULL
  tau_density_plot <- NULL
  ln_density_plot <- NULL
  lc_density_plot <- NULL
  lntau_prob_plot <- NULL
  lnt_lc_heatmaps <- NULL
  sn_prob_plot <- NULL
  sc_prob_plot <- NULL
  param_differences_boxplot <- NULL
  plot_colours <- rep(colours_3, ceiling(repetitions / length(colours_3)))
  
  # plotting marginal likelihood distributions ........................................................................
  
  if (plot_margins) {
    
    # constructing marginal likelihood densities for the single parameters
    
    mu_densities <- sapply(1:repetitions, function(i) {
      joint_density <- densities[[i]][[1]]
      list(maxlike_marginalise_continuous(joint_density, 1, mu_borders, tau_borders, ln_borders, lc_borders))
    })
    tau_densities <- sapply(1:repetitions, function(i) {
      joint_density <- densities[[i]][[1]]
      list(maxlike_marginalise_continuous(joint_density, 2, mu_borders, tau_borders, ln_borders, lc_borders))
    })
    ln_densities <- sapply(1:repetitions, function(i) {
      joint_density <- densities[[i]][[1]]
      list(maxlike_marginalise_continuous(joint_density, 3, mu_borders, tau_borders, ln_borders, lc_borders))
    })
    lc_densities <- sapply(1:repetitions, function(i) {
      joint_density <- densities[[i]][[1]]
      list(maxlike_marginalise_continuous(joint_density, 4, mu_borders, tau_borders, ln_borders, lc_borders))
    })
    
    # mu marginal density distribution
    mu_vals <- seq(from = mu_borders[1], to = mu_borders[2], length = 1000)
    mu_density_plot <- plot_densities(mu_vals, mu_densities, repetition_labels, "mu probability density",
                                      "mu", "probability densiy", "repetitions", plot_colours, p_true = mu)
    print(mu_density_plot)
    
    # tau marginal density distribution
    tau_vals <- seq(from = tau_borders[1], to = tau_borders[2], length = 1000)
    tau_density_plot <- plot_densities(tau_vals, tau_densities, repetition_labels, "tau probability density",
                                       "tau", "probability density", "repetitions", plot_colours, p_true = tau)
    print(tau_density_plot)
    
    # lambda_n marginal density distribution
    ln_vals <- seq(from = ln_borders[1], to = ln_borders[2], length = 1000)
    ln_density_plot <- plot_densities(ln_vals, ln_densities, repetition_labels, "lambda_n probability density",
                                      "lambda_n", "probability density", "repetitions", plot_colours, 
                                      p_true = lambda_n)
    print(ln_density_plot)
    
    # lambda_c marginal density distribution
    lc_vals <- seq(from = lc_borders[1], to = lc_borders[2], length = 1000)
    lc_density_plot <- plot_densities(lc_vals, lc_densities, repetition_labels, "lambda_c probability density",
                                      "lambda_c", "probability density", "repetitions", plot_colours, 
                                      p_true = lambda_c)
    print(lc_density_plot)
    
  }
  
  # plotting derived parameters' probability distributions ............................................................
  
  if (plot_derived) {
    
    # lambda_n + tau probability distribution
    lntau_vals <- seq(from = dparams_borders[1], to = dparams_borders[2], length = 1000)
    lntau_probabilities <- sapply(1:repetitions, function(i) {densities[[i]][[2]]})
    lntau_prob_plot <- plot_densities(lntau_vals, lntau_probabilities, repetition_labels,
                                      "lambda_n + tau probability distribution", "lambda_n + tau",
                                      "probability", "repetitions", plot_colours, p_true = (lambda_n + tau))
    print(lntau_prob_plot)
    
    # lambda_n + tau and lambda_c 2D probability distribution
    # (plotting one heatmap for each repetition)
    lc_vals <- seq(from = lc_borders[1], to = lc_borders[2], length = 1000)
    lnt_lc_probabilities <- sapply(1:repetitions, function(i) {densities[[i]][[3]]})
    lnt_lc_heatmaps <- sapply(1:repetitions, function(i) {          # one heatmap for each repetition
      current_prob <- lnt_lc_probabilities[[i]]
      distribution_matrix <- sapply(lntau_vals, function(lnt) {     # creating matrix of 2D distribution values
        sapply(lc_vals, function(lc) {                              # for heatmap plotting
          current_prob(c(lnt, lc))
        })
      })
      lnt_lc_hm <-
        heatmap_format(distribution_matrix, col_names = c(rep(NA, 1000)), row_names = c(rep(NA, 1000)),
                       plot_title = paste("lambda_n + tau and lambda_c\njoint probability distribution",
                                          repetition_labels[i]),
                       plot_xlab = "lambda_n + tau", plot_ylab = "lambda_c", key_xlab = "likelihood",
                       vline = which(lntau_vals > (lambda_n + tau))[1], hline = which(lc_vals > lambda_c)[1])
      print(lnt_lc_hm)
      return(lnt_lc_hm)
    })
    
    # steady_n probability distribution
    sn_vals <- seq(from = dparams_borders[3], to = dparams_borders[4], length = 1000)
    sn_probabilities <- sapply(1:repetitions, function(i) {densities[[i]][[4]]})
    sn_prob_plot <- plot_densities(sn_vals, sn_probabilities, repetition_labels, "steady_n probability distribution",
                                   "steady_n", "probability", "repetitions", plot_colours, 
                                   p_true = (mu / (lambda_n + tau)))
    print(sn_prob_plot)
    
    # steady_c probability distribution
    sc_vals <- seq(from = dparams_borders[5], to = dparams_borders[6], length = 1000)
    sc_probabilities <- sapply(1:repetitions, function(i) {densities[[i]][[5]]})
    sc_prob_plot <- plot_densities(sc_vals, sc_probabilities, repetition_labels, "steady_c probability distribution",
                                   "steady_c", "probability", "repetitions", plot_colours,
                                   p_true = (mu * tau / ((lambda_n + tau) * lambda_c)))
    print(sc_prob_plot)
    
  }
  
  # box plot of parameter estimation accuracy .........................................................................
  
  if (plot_boxplot) {
    
    # creating matrix with parameter differences
    p_true <- c(mu, tau, lambda_n, lambda_c)
    param_differences_matrix <- sapply(1:4, function(p) {
      sapply(1:repetitions, function(i) {
        p_estimated = optima[[i]][p]
        p_estimated / p_true[p]
      })
    })
    
    # getting box plot data frame
    param_differences_df <- boxplot_dataframe(param_differences_matrix, c("mu", "tau", "lambda_n", "lambda_c"))
    
    # plotting box plot
    param_differences_boxplot <- ggplot(param_differences_df,
                                        aes(param_differences_df$parameter, param_differences_df$difference,
                                            fill = parameter)) +
      geom_boxplot() + scale_fill_manual(values = c("royalblue3", "mediumseagreen", "seagreen1", "salmon"))
    param_differences_boxplot <- plot_format(param_differences_boxplot, "Parameter Estimation\n",
                                             "parameter", "estimation error (estimated/true)")
    print(param_differences_boxplot)
    
  }
  
  # returning .........................................................................................................
  
  return(list(mu_density_plot, tau_density_plot, ln_density_plot, lc_density_plot,
              lntau_prob_plot, lnt_lc_heatmaps, sn_prob_plot, sc_prob_plot, param_differences_boxplot))
  
}

### Data Analysis #####################################################################################################

sample_variance <- function(summary_table_s1, summary_table_s2, mode = "raw") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summary_table_s1:     data frame containing measurement data of sample 1
  # summary_table_s2:     data frame containing measurement data of sample 2
  # mode:                 mode for variance computation; choose "raw" for computing raw variance, choose "raw_norm" for
  #                       computing raw variance normalized by expected variance (using a beta distribution for average
  #                       mod/total signal to account for measurement imprecision, assuming a binomial distribution for
  #                       the single measurements' data),
  #                       choose "transform" for computing variance on arcsin(sqrt) transformed data, normalized by 
  #                       1/(4*n) (arcsin(sqrt) is a common variance-stabilizing transform for ratios and also suitable
  #                       for binomials, assuming an approximated Gaussian distribution of the transformed data which 
  #                       comes down to a variance of 1/(4*n)) (default: raw)
  
  # *** NOTE ***
  
  # data frames storing measurement data must be sorted by gene names and time points (nuclear time series, 
  # cytosolic time series) and have the following format:
  # rows: genomic regions
  # col 1: region name, col 2: measurement description, col 3: library size, 
  # col 4: total transcript counts, col 5: labeled transcript counts, col 6: average potential conversion positions, 
  # col 7: conversion efficiency, col 8: newly synthesized ratio
  
  # *** RETURN ***
  
  # a matrix storing both variances and average signals for the single gene's measurements (rows (first half): time
  # point measurements' variances, rows (second half): time point measurements' average signals, columns: gene names)
  
  # *******************************************************************************************************************
  
  gene_names_s1 <- unique(summary_table_s1[, "name"])     # getting gene names recorded for sample 1
  gene_names_s2 <- unique(summary_table_s2[, "name"])     # getting gene names recorded for sample 2
  
  overlapping_genes <- intersect(gene_names_s1, gene_names_s2)  # gene names recorded in both samples
  
  # iterating through genes, computing variances
  
  variance_vs_signal <- sapply(1:length(overlapping_genes), function(i) {  
    
    print(paste("gene ", i, " from ", length(overlapping_genes), sep=""))
    g <- overlapping_genes[i]
    
    subtable_s1 <- summary_table_s1[which(summary_table_s1[, "name"] == g), ]   # gene's sample 1 subtable
    subtable_s2 <- summary_table_s2[which(summary_table_s2[, "name"] == g), ]   # gene's sample 2 subtable
    total_s1 <- subtable_s1[, "tot"]      # gene's sample 1 total counts
    total_s2 <- subtable_s2[, "tot"]      # gene's sample 2 total counts
    mod_s1 <- subtable_s1[, "mod"]        # gene's sample 1 modified counts
    mod_s2 <- subtable_s2[, "mod"]        # gene's sample 2 modified counts
    
    if (mode == "raw") {              # computing raw variance
      
      vars <- (mod_s1 / total_s1 - mod_s2 / total_s2) ** 2
      avg_signal <- (mod_s1 / total_s1 * total_s1 + mod_s2 / total_s2 * total_s2) / (total_s1 + total_s2)
        
    }
    
    else if (mode == "raw_norm") {    # computing raw, normalized variance
      
      alpha <- 1 + mod_s1 + mod_s2                # alpha parameter for Beta-distribution on true mod/total ratio p
      beta <- 2 + total_s1 + total_s2 - alpha     # beta parameter for Beta-distribution on true mod/total ratio p

      # variance measured from data
      measured_variance <- (mod_s1 / total_s1 - mod_s2 / total_s2) ** 2   
      
      # variance expected for fixed p
      point_expected_variance <- function(p, a = alpha, b = beta) {       
        p_prob <- dbeta(p, a, b)    # probability density for mod/total ratio of p
        var <- (total_s1 * (total_s1 - 1) * p ** 2 + total_s1 * p) / total_s1 ** 2 +
        (total_s2 * (total_s2 - 1) * p ** 2 + total_s2 * p) / total_s2 ** 2 -
        2 * p ** 2                  # variance expected for mod/total ratio of p
        return(p_prob * var)        # probability of p * expected variance for p
      }
      
      # variance expected over parameter space of p
      full_expected_variance <- sapply(1:length(alpha), function(i) {
        var_results <- integrate(point_expected_variance, lower=0, upper=1, a = alpha[i], b = beta[i],
                                 subdivisions = 1e3, rel.tol = 1e-15, stop.on.error = FALSE)                 
        return(var_results$value)
      })
     
      vars <- measured_variance / full_expected_variance                  # normalized variance
      avg_signal <- (mod_s1 / total_s1 * total_s1 + mod_s2 / total_s2 * total_s2) / (total_s1 + total_s2)
      
    }
    
    else if (mode == "transform") {   # computing transformed data's variances
      
      transform_s1 <- asin(sqrt(mod_s1 / total_s1))
      transform_s2 <- asin(sqrt(mod_s2 / total_s2))
      vars <- (transform_s1 - transform_s2) ** 2 / (1 / (4 * total_s1) + 1 / (4 * total_s2))
      avg_signal <- (transform_s1 * total_s1 + transform_s2 * total_s2) / (total_s1 + total_s2)
    
    }
    
    return(c(vars, avg_signal))
    
  })

  # setting return matrix rownames and column names, returning
  
  data_des <- summary_table_s1[which(summary_table_s1[, "name"] == gene_names_s1[1]), "des"]  # measurement descriptions
  var_des <- sapply(data_des, function(i) {paste(i, "var", sep="_")})   # variance rows' descriptions
  sig_des <- sapply(data_des, function(i) {paste(i, "sig", sep="_")})   # signal rows' descriptions
  rownames(variance_vs_signal) <- c(var_des, sig_des)
  colnames(variance_vs_signal) <- overlapping_genes
  return(variance_vs_signal)
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

expr_level_regression <- function(summary_table, time_series, norm="lib") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summary_table:      summary table file
  # time_series:        vector of time points at which measurements were taken (not including time point 0)
  # norm:               value by which total counts are normalized to compute expression levels; choose 'lib' to use 
  #                     the library size, choose 'avg' to use the average of the total counts distribution, computed
  #                     for the middle 50% of the distribution (for robustness) (default: "lib")
  
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

cluster_count_ratios <- function(mode, summary_table, n_timepoints,
                                 n_clusters, outfile = NULL, iter_max = 10, center_sets = 1) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mode:                   mode to be executed (type of count data to be clustered); use "mod_total" to cluster by
  #                         modified/total transcript count ratios, use "total_lib" to cluster by total/library_size
  #                         transcript count ratios
  # summary_table:          summary table file
  # n_timepoints:           number of time points measured
  # n_clusters:             number of clusters to generate
  
  # outfile:                optional; output file to write results to (default: NULL)
  # iter_max:               optional; maximum number of clustering iterations (default: 10)
  # center_sets:            optional; number of random cluster center sets to choose from as initialization 
  #                         (default: 1)
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)
  
  # *** RETURN ***
  
  # a list containing (in the following order):
  #   data table on which the clustering was performed (storing modified/total or total/library_size ratios)
  #       rows: genomic regions, 
  #       columns: data_description 
  #   a vector containing the RNA species' cluster numbers 
  #   a matrix of cluster centers (rows: cluster numbers, columns: time points)
  #   total sum-of-squares
  #   vector of within-cluster sum-of-squares
  #   the number of RNA species within each cluster
  
  # This function clusters RNA species by their modified/total transcript counts or total/library_size counts, using 
  # k-means clustering.
  
  # *******************************************************************************************************************
  
  # loading in summary tables
  
  summary_table <- load_summarytable(summary_table)
  
  # computing transcript count ratios (on which clustering should be performed)
  
  ratios <- summarytable_ratios("counts", summary_table, n_timepoints)
  if (mode == "mod_total") {
    ratios <- ratios[,1:(2*n_timepoints)]
  }
  else if (mode == "total_lib") {
    ratios <- ratios[,(2*n_timepoints+1):(n_timepoints*4)]
  }
  else {
    ratios <- NULL
  }
    
  # clustering
  
  cluster_results <- kmeans(ratios, n_clusters, iter.max = iter_max, nstart = center_sets)
  cluster_results <- list(cluster_results[1], cluster_results[2], cluster_results[3], cluster_results[4], 
                          cluster_results[7])
  
  # writing to output file and returning
  
  if (!is.null(outfile)) {
    sink(outfile)
    write.table(ratios, sep = "\t")
    for (r in cluster_results) {
      write.table(r, sep = "\t", append = TRUE)
    }
    sink()
  }
  
  return(cluster_results)
}

cluster_curve_types <- function(mode, summary_table, n_timepoints,
                                n_clusters, outfile = NULL, iter_max = 10, center_sets = 1) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mode:                   mode to be executed (type of count data to be clustered); use "mod_total" to cluster by
  #                         normalized modified/total transcript count ratios, use "total_lib" to cluster by normalized
  #                         total/library_size transcript count ratios
  # summary_table:          summary table file
  # n_timepoints:           number of time points measured
  # n_clusters:             number of clusters to generate
  
  # outfile:                optional; output file to write results to (default: NULL)
  # iter_max:               optional; maximum number of clustering iterations (default: 10)
  # center_sets:            optional; number of random cluster center sets to choose from as initialization 
  #                         (default: 1)
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)
  
  # *** RETURN ***
  
  # a list containing (in the following order):
  #   data table on which the clustering was performed (storing modified/total and total/library_size ratios)
  #       rows: genomic regions, 
  #       columns: data_description
  #   a vector containing the RNA species' cluster numbers 
  #   a matrix of cluster centers (rows: cluster numbers, columns: time points)
  #   total sum-of-squares
  #   vector of within-cluster sum-of-squares
  #   the number of RNA species within each cluster
  
  # This function clusters RNA species by their normalized modified/total transcript counts or normalized 
  # total/library_size counts, using k-means clustering.
  
  # *******************************************************************************************************************
  
  # loading in summary tables
  
  summary_table <- load_summarytable(summary_table)
  
  # computing transcript count ratios (on which clustering should be performed)

  ratios <- summarytable_ratios("counts", summary_table, n_timepoints)

  if (mode == "mod_total") {
    ratios <- ratios[,1:(n_timepoints*2)]
  }
  else if (mode == "total_lib") {
    ratios <- ratios[,(n_timepoints*2+1):(n_timepoints*4)]
  }
  else {
    ratios <- NULL
  }
  
  # normalizing transcript count ratios (by their mean) and log-transforming (foldchanges up and down come to same absolute scale)
  
  mean_ratios <- sapply(1:nrow(ratios), function(i) {
    c(
      (sum(ratios[i, 1:n_timepoints]) / n_timepoints),
      (sum(ratios[i, (n_timepoints+1):(n_timepoints*2)]) / n_timepoints)
    )
  })
  ratios_norm <- cbind(
    ratios[, 1:n_timepoints] / mean_ratios[1, ], 
    ratios[, (n_timepoints+1):(n_timepoints*2)] / mean_ratios[2, ]
  )
  ratios_norm[, 1:(n_timepoints*2)] <- log10(ratios_norm[, 1:(n_timepoints*2)])

  # clustering
  
  cluster_results <- kmeans(ratios_norm, n_clusters, iter.max = iter_max, nstart = center_sets)
  cluster_results <- list(cluster_results[1], cluster_results[2], cluster_results[3], cluster_results[4], 
                          cluster_results[7])
  
  # writing to output file and returning
  
  if (!is.null(outfile)) {
    sink(outfile)
    write.table(ratios_norm, sep = "\t")
    for (r in cluster_results) {
      write.table(r, sep = "\t", append = TRUE)
    }
    sink()
  }
  
  return(cluster_results)
}

halflife_landscape <- function(mode, gene, summary_table, lntau_grid, lc_grid, 
                               conveff_nu = NULL, conveff_cy = NULL, fn_error_nu = NULL, fn_error_cy = NULL,
                               fp_error_nu = NULL, fp_error_cy = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mode:                   half life landscape to be computed; choose "raw_lsq" for least squares on raw data, 
  #                         choose "trans_lsq" for least squares on variance-stabilizing transformed data, choose
  #                         "binom" for a Bayesian approach assuming underlying binomial distributions 
  # gene:                   gene name for which the half life landscape should be computed
  # summary_table:          summary table file
  # time_series:            time points at which measurements were taken
  # lntau_grid:             vector of lntau values for which half life landscape is to be computed
  # lc_grid:                vector of lc values for which half life landscape is to be computed
  
  # conveff_nu:       needed for mode == "binom": nuclear conversion efficiency, resolved by time (default: NULL)
  # conveff_cy:       needed for mode == "binom": cytosolic conversion efficiency, resolved by time (default: NULL)
  # fn_error_nu:      needed for mode == "binom": nuclear false-negative error; probability by which a single, 
  #                   labeled U position (i.e. 4sU converted to C) does not show up as C position in the alignment 
  #                   (can be set to the C > non-C conversion rate) (default: NULL)
  # fn_error_cy:      needed for mode == "binom": cytosolic false-negative error; probability by which a single, 
  #                   labeled U position (i.e. 4sU converted to C) does not show up as C position in the alignment 
  #                   (can be set to the C > non-C conversion rate) (default: NULL)
  # fp_error_nu:      needed for mode == "binom": nuclear false-positive error; probability by which a single, 
  #                   unlabeled U position shows up as C (i.e. labeled) position in the alignment (default: NULL)
  # fp_error_cy:      needed for mode == "binom": cytosolic false-positive error; probability by which a single, 
  #                   unlabeled U position shows up as C (i.e. labeled) position in the alignment (default: NULL)
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)
  
  # *** RETURN ***
  
  # a matrix storing the half life landscape (sum of squares values for least squares approaches, likelihood values
  # for the Bayesian approaches), rows: lambda_c values (lambda_c grid), columns: lambda_n + tau values (lambda_n +
  # tau grid)
  
  # *******************************************************************************************************************

  # getting gene's measurement data
  
  n_timepoints <- length(time_series)                                 # number of measurement time points
  stable <- load_summarytable(summary_table)                          # loading in summary table
  gene_table <- stable[which(stable[, "name"] == gene), ]             # getting gene's sub-table of summary table
  nu_table <- gene_table[1:n_timepoints, ]                            # gene's nuclear data
  cy_table <- gene_table[(n_timepoints + 1):(n_timepoints * 2), ]     # gene's cytosolic data
  
  # least squares on raw data
  
  if (mode == "raw_lsq") {
    
    nu_newtot <- nu_table[, "nr"]     # nuclear corrected new/total ratios
    cy_newtot <- cy_table[, "nr"]     # cytosolic corrected new/total ratios
    nu_tot <- nu_table[, "tot"]       # nuclear total counts
    cy_tot <- cy_table[, "tot"]       # cytosolic total counts
    
    hf_landscape <- sapply(lntau_grid, function(lnt) {
      sapply(lc_grid, function(lc) {
        halflife_sumofsquares_newtot(c(lnt, lc), time_series, nu_newtot, cy_newtot, nu_tot, cy_tot, var_stab = FALSE)
      })
    })
    
  }
  
  # least squares on transformed data
  
  if (mode == "trans_lsq") {
  
    nu_newtot <- nu_table[, "nr"]     # nuclear corrected new/total ratios
    cy_newtot <- cy_table[, "nr"]     # cytosolic corrected new/total ratios
    nu_tot <- nu_table[, "tot"]       # nuclear total counts
    cy_tot <- cy_table[, "tot"]       # cytosolic total counts
    
    hf_landscape <- sapply(lntau_grid, function(lnt) {
      sapply(lc_grid, function(lc) {
        halflife_sumofsquares_newtot(c(lnt, lc), time_series, nu_newtot, cy_newtot, nu_tot, cy_tot, var_stab = TRUE)
      })
    })
  
  }
  
  # Bayesian approach assuming binomial distributions
  
  if (mode == "binom") {
    
    nu_total <- nu_table[, "tot"]     # nuclear total counts
    nu_mod <- nu_table[, "mod"]       # nuclear modified counts
    cy_total <- cy_table[, "tot"]     # cytosolic total counts
    cy_mod <- cy_table[, "mod"]       # cytosolic modified counts
    
    convpos <- sum(c(nu_table[, "cp"] * nu_total, cy_table[, "cp"] * cy_total)) / 
      sum(c(nu_total, cy_total))                          # weighted average number of potential conversion positions
    p_mod_nu <- 1 - 
      (1 - conveff_nu * (1 - fn_error_nu)) ** convpos     # modification probability of nuclear newly synthesized reads
    p_mod_cy <- 1 - 
      (1 - conveff_cy * (1 - fn_error_cy)) ** convpos     # modification probability of cytosolic newly synthesized reads
    p_un_nu <- (1 - fp_error_nu) ** convpos               # non-modfication probability of nuclear pre-existing reads
    p_un_cy <- (1 - fp_error_cy) ** convpos               # non-modfication probability of cytosolic pre-existing reads 
    
    n_gridvals <- length(lntau_grid)    # number of values in the parameter grid
    joint_probability <-                # getting (lambda_n + tau) and lambda_c joint likelihood function
      halflife_estimation_maxlike(nu_total, cy_total, nu_mod, cy_mod, 
                                  p_un_nu = rep(p_un_nu, n_timepoints), p_un_cy = rep(p_un_cy, n_timepoints), 
                                  p_mod_nu = p_mod_nu, p_mod_cy = p_mod_cy,
                                  lntau_borders = c(lntau_grid[1], lntau_grid[n_gridvals]), 
                                  lc_borders = c(lc_grid[1], lc_grid[n_gridvals]), 
                                  time_series = time_series, mode = "binom", use_lookup = TRUE, ip_points = 1000)[[6]]
    
    hf_landscape <- sapply(lntau_grid, function(lnt) {
      sapply(lc_grid, function(lc) {
        joint_probability(c(lnt, lc))
      })
    })
    
  }
  
  # returning half life landscape
  
  return(hf_landscape)
  
}

halflife_Rsquared <- function(original_table, estimation_table, time_series, mode,
                              conveff = NULL, fn_error_nu = NULL, fn_error_cy = NULL, fp_error_nu = NULL, fp_error_cy = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # original_table:     data frame containing measurement data
  # estimation_table:   data frame containing parameter estimators
  # time_series:        vector of time points at which measurements were taken (not including time point 0)
  # mode:               mode of estimation procedure; choose "infer" for labeled/total ratios inferred from the
  #                     predicted new/total ratios, choose "estim" for back-estimated new/total ratios
  
  # *** the following parameters are mandatory for 'mode' == "infer", but redundant for 'mode' == "estim" ***
  
  # conveff:            vector of conversion efficiencies (at single measurement time points, ordered by nuclear and
  #                     cytosolic time series); probability by which a single base is labeled
  # fn_error_nu:        nuclear false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for nuclear transcripts)
  # fn_error_nu:        cytosolic false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for cytosolic transcripts)
  # fp_error_nu:        nuclear false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for nuclear transcripts)
  # fp_error_cy:        cytosolic false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for cytosolic transcripts)
  
  # *** NOTE ***
  
  # data frame storing measurement data should be sorted by gene names and time points and have the following format:
  # rows: genomic regions
  # col 1: region name, col 2: measurement description, col 3: library size, 
  # col 4: total transcript counts, col 5: labeled transcript counts, col 6: average potential conversion positions, 
  # col 7: conversion efficiency, col 8: newly synthesized ratio
  
  # data frame storing parameter estimatiors should be in the following format:
  # rows: genomic regions,
  # col 1: lambda_n + tau estimator, col 2: lambda_c estimator, col 3: sum-of-squares values
  
  # *** RETURN ***
  
  # a data frame containing two columns storing normalized sum-of-squares for genes' parameter estimations 
  # (column 1: for nuclear time curve, column 2: for cytosolic time curve)
  
  # *******************************************************************************************************************
  
  # number of time series measurements
  n_timepoints <- length(time_series)
  # number of genes within data tables
  n_genes <- nrow(estimation_table)
  # names of genes within data tables
  gene_names <- rownames(estimation_table)
  # initialize matrix storing R-squared values (rows: genes, col1: nuclear time curve, col2: cytosolic time curve)
  r_squared <- matrix(rep(0, n_genes * 2), ncol = 2)
  
  for(i in 1:n_genes) {       # iterating through genes
    
    g <- gene_names[i]          # getting gene name
    
    # original measurements ...........................................................................................
    
    gene_original <- 
      original_table[which(original_table[, "name"] == g), ]                # current gene's total sub-table
    gene_original_nu <- gene_original[1:n_timepoints,]                      # current gene's nucleus sub-table
    gene_original_cy <- gene_original[(n_timepoints+1):(n_timepoints*2), ]  # current gene's cytosol sub-table
    
    # predicted ratios ................................................................................................
    
    # lambda_n + tau and lambda_c estimators
    lnt_estimation <- estimation_table[which(rownames(estimation_table) == g), "lntau"]
    lnt_split <- lnt_estimation / 2       # splitting lambda_n + tau estimator into two halfs (for both variables)
    lc_estimation <- estimation_table[which(rownames(estimation_table) == g), "lc"]
    
    # computing new/total transcript count ratios according to parameter estimators
    ratio_new_total_nu <- ratio_n_sol(time_series, lnt_split, lnt_split)           
    ratio_new_total_cy <- ratio_c_sol(time_series, lnt_split, lnt_split, lc_estimation)
    
    if (mode == "infer") {
      
      # measured ratios ...............................................................................................
      
      gene_original_nu_ratios <-                                              # current gene's nuclear mod/total ratio
        c(gene_original_nu[, "mod"] / gene_original_nu[, "tot"]) * 100        # (multiply by 100 for %)
      gene_original_cy_ratios <-                                              # current gene's cytosolic mod_total ratio
        c(gene_original_cy[, "mod"] / gene_original_cy[, "tot"]) * 100        # (multiply by 100 for %)
      
      # inferred ratios .................................................................................................
      
      convpos <- sum(c(gene_original_nu[, "cp"], gene_original_cy[, "cp"])) / 
        (n_timepoints * 2)                              # average number of potential conversion positions
      
      p_mod_nu <- 1 - 
        (1 - conveff * (1 - fn_error_nu)) ** convpos    # modification probability of nuclear newly synthesized reads
      p_mod_cy <- 1 - 
        (1 - conveff * (1 - fn_error_cy)) ** convpos    # modification probability of cytosolic newly synthesized reads
      p_un_nu <- (1 - fp_error_nu) ** convpos           # non-modfication probability of nuclear pre-existing reads
      p_un_cy <- (1 - fp_error_cy) ** convpos           # non-modfication probability of cytosolic pre-existing reads
      
      # computing modified/total transcript count ratios according to modification probabilities (multiply by 100 for %)
      ratio_mod_total_nu <- (ratio_new_total_nu * p_mod_nu + (1 - ratio_new_total_nu) * (1 - p_un_nu)) * 100
      ratio_mod_total_cy <- (ratio_new_total_cy * p_mod_cy + (1 - ratio_new_total_cy) * (1 - p_un_cy)) * 100
      
      # residual sum of squares .......................................................................................
      
      sum_of_residues_nu <- sum((gene_original_nu_ratios - ratio_mod_total_nu) ** 2)
      sum_of_residues_cy <- sum((gene_original_cy_ratios - ratio_mod_total_cy) ** 2)
      
    }
    
    else if (mode == "estim") {
      
      # measured ratios ...............................................................................................
      
      gene_original_nu_ratios <- gene_original_nu[, "nr"] * 100     # gene's recorded nuclear new/total ratio
      gene_original_cy_ratios <- gene_original_cy[, "nr"] * 100     # gene's recorded cytosolic new/total ratio
      
      # residual sum of squares .......................................................................................
      
      sum_of_residues_nu <- sum(((gene_original_nu_ratios) - (ratio_new_total_nu * 100)) ** 2)  # x100 for %
      sum_of_residues_cy <- sum(((gene_original_cy_ratios) - (ratio_new_total_cy * 100)) ** 2)  # x100 for %
    }
    
    # computing normalized sum-of-squares .............................................................................
    
    # total sum of squares:
    data_avg_nu <- sum(gene_original_nu_ratios) / n_timepoints
    data_avg_cy <- sum(gene_original_cy_ratios) / n_timepoints
    sum_total_nu <- sum((gene_original_nu_ratios - data_avg_nu) ** 2)
    sum_total_cy <- sum((gene_original_cy_ratios - data_avg_cy) ** 2)
    
    # R-squared:
    r_squared_nu <- 1 - sum_of_residues_nu/sum_total_nu
    r_squared_cy <- 1 - sum_of_residues_cy/sum_total_cy
    r_squared[i, 1] <- r_squared_nu
    r_squared[i, 2] <- r_squared_cy
  }
  
  # formatting r-squared data frame and returning
  
  r_squared_df <- data.frame(r_squared)
  rownames(r_squared_df) <- gene_names
  colnames(r_squared_df) <- c("nu_timecurve", "cy_timecurve")
  return(r_squared_df)
}

cy_nu_ratio <- function(data_table, gene_list = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # data_table:     matrix or data frame containing nuclear and cytosolic half life parameters as columns "deg_nuc_avg"
  #                 and "deg_cyt_avg", nuclear and cytosolic normalized expression levels as "mean_elvl_nuc_avg" and
  #                 "mean_elvl_cyt_avg" and genomic regions as rows named by the corresponding IDs
  # gene_list:      optional; vector of genomic region IDs to be used for calculation (default: all genomic regions
  #                 recorded in the input data table)
  
  # *** RETURN ***
  
  # a list containing a named vector storing the single genes' cyt/nuc ratio estimates, a single value being the point 
  # estimate for the cyt/nuc ratio, and another named vector storing the single gene's lambda_n estimates (computed 
  # under the use of the cyt/nuc ratio point estimate)
  
  # This function computes an upper limit estimate of the cellular cytosolic/nuclear transcript count ratio, using the 
  # half life parameter estimators and normalized expression levels.
  
  # *******************************************************************************************************************
  
  if(is.null(gene_list)) {    # if no genomic region selection is specified, use all in the input data table
    gene_list <- rownames(data_table)
  }
  
  # icomputing the cyt/nuc ratio upper limit estimates
  ratio_estimates <- 
    data_table[gene_list, "mean_elvl_nuc_avg"] * data_table[gene_list, "deg_nuc_avg"] / 
    (data_table[gene_list, "mean_elvl_cyt_avg"] * data_table[gene_list, "deg_cyt_avg"])
  
  # computing the median cyt/nuc ratio estimate
  final_ratio_estimate <- median(ratio_estimates)
  
  # plugging in the final estimate for the cyt/nuc ratio upper limit, computing the corresponding lambda_n values
  lambda_n_values <- 
    (data_table[gene_list, "mean_elvl_nuc_avg"] * data_table[gene_list, "deg_nuc_avg"] - 
    final_ratio_estimate * data_table[gene_list, "mean_elvl_cyt_avg"] * data_table[gene_list, "deg_cyt_avg"]) /
    data_table[gene_list, "mean_elvl_nuc_avg"]
  
  # returning
  return(list(ratio_estimates, final_ratio_estimate, lambda_n_values))
}

### Sampling ##########################################################################################################

sample_nnew <- function(n_total, p) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_total:    total number of counts
  # p:          ratio of newly synthesized compared to total RNA transcripts
  
  # *** RETURN ***
  
  # the number of newly synthesized RNA transcript counts
  
  # This functions samples the counts of newly synthesized RNA transcripts from a binomial distribution.
  
  # *******************************************************************************************************************
  
  return(rbinom(1, n_total, p))
}

sample_nmod <- function(n_total, n_new, p_un, p_mod) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_total:    total transcript counts
  # n_new:      newly synthesized transcript counts
  # p_un:       probability that a pre-existing transcript remains unmodified
  # p_mod:      probability that a newly synhesized transcript is modified
  
  # *** RETURN ***
  
  # the number of modified RNA transcript counts
  
  # This function samples the number of modified transcript counts from a folding of two binomial distributions
  # describing the number of modified newly synthesized and modified pre-existing counts, in dependence of the number
  # of newly synthesized and pre-existing transcript counts.
  
  # *******************************************************************************************************************
  
  n_pre = n_total - n_new
  prob_distribution <- sapply(0:n_total, function(n_mod) {
    p_nmod_nnew(n_mod, n_pre, n_new, p_un, p_mod)
  })
  
  return(sample(0:n_total, 1, prob = prob_distribution))
}

sample_transcripts <- function(n_total, c_total, q_nucleus, q_cytosol, p_un, p_mod, number = 1) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_total:      total number of nuclear transcript counts, resolved by time points (vector)
  # c_total:      total number of cytosolic transcript counts, resolved by time points (vector)
  # q_nucleus:    probability of observing a newly synthesized transcript within the nucleus (ratio of newly 
  #               synthesized transcripts within the nucleus), resolved by time points (vector)
  # q_cytosol:    probability of observing a newly synthesized transcript within the cytosol (ratio of newly 
  #               synthesized transcripts within the cytosol), resolved by time points (vector)
  # p_un:         probability that a pre-existing transcript remains unmodified, resolved by time points (vector)
  # p_mod:        probability that a newly synhesized transcript is modified, resolved by time points (vector)
  # number:       optional; number of samples to be generated (default: 1)
  
  # *** RETURN ***
  
  # a list cointaining two matrices, the first matrix storing nucleus and the second cytosol modified transcript counts
  # (columns: repetitions, rows: time points)
  
  # This function samples modified transcript counts for both the nucleus and the cytosol, resolved by time, in
  # dependence of RNA dynamics parameters (as represented by new/total ratios) and global parameters.
  
  # *******************************************************************************************************************
  
  time_points = length(q_nucleus)
  
  # sampling modified transcripts counts within nucleus
  
  sample_mod_n <- sapply(1:number, function(n) {  # sampling as often as specified by number
    
    sample_new_n <- sapply(1:time_points, function(t){    # sampling newly synthesized transcript counts
      sample_nnew(n_total[t], q_nucleus[t])})
    sapply(1:time_points, function(t){                    # sampling modified transcript counts
      sample_nmod(n_total[t], sample_new_n[t], p_un[t], p_mod[t])})
  })
  
  # sampling modified transcripts counts within cytosol
  
  sample_mod_c <- sapply(1:number, function(n) {  # sampling as often as specified by number
    
    sample_new_c <- sapply(1:time_points, function(t){   # sampling newly synthesized transcript counts
      sample_nnew(c_total[t], q_cytosol[t])})
    sapply(1:time_points, function(t){                   # sampling modified transcript counts
      sample_nmod(c_total[t], sample_new_c[t], p_un[t], p_mod[t])})
  })
  
  # returning: list cointaining two matrices, with columns referring repetitions and rows to time points
  
  return(list(sample_mod_n, sample_mod_c))
}

sample_count_data <- function(mu, tau, lambda_n, lambda_c, time_series,
                              rna_total_n, rna_total_c, library_size_n, library_size_c, 
                              p_un, p_mod, number = 1) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mu:                 synthesis rate
  # tau:                transport rate
  # lambda_n:           nuclear degradation rate
  # lambda_c:           cytosolic degradation rate
  # time_series:        time points at which measurements were taken (vector)
  
  # rna_total_n:      total amount of nuclear RNA transcripts [molecules / cell] (summed over all transcript species)
  # rna_total_c:      total amount of cytosolic RNA transcripts [molecules / cell] (summed over all transcript species)
  # library_size_n:   total amount of nuclear read counts (summed over all transcript species), 
  #                   resolved by time (vector)
  # library_size_c:   total amount of cytosolic read counts (summed over all transcript species),
  #                   resolved by time (vector)
  
  # p_un:       probability that a pre-existing transcript remains unmodified, resolved by time (vector)
  # p_mod:      probability that a newly synhesized transcript is modified, resolved by time (vector)
  # number:     optional; number of samples to be generated (default: 1)
  
  # *** RETURN ***
  
  # a list containing the following matrices (rows: time points, columns: repetitions) (in the following order): 
  #   matrix of nuclear transcript counts,
  #   matrix of modified nuclear transcript counts,
  #   matrix of numbers of cytosolic transcript counts,
  #   matrix of modified cytosolic transcript counts
  
  # This function simulates observed count data for a transcript species g.
  
  # *******************************************************************************************************************
  
  time_points = length(time_series)        # number of count samples to generate
  
  # simulating RNA dynamics
  
  raw_data <- multi_rna_dynamics_sol(c(mu), c(tau), c(lambda_n), c(lambda_c), time_series)
  raw_q_nucleus <- raw_data[[5]][,2]  # RNA ratio := q
  raw_q_cytosol <- raw_data[[6]][,2]  # RNA ratio := q
  
  # sampling total nuclear and cytosolic transcript counts
  
  steady_n_raw = mu / (lambda_n + tau)                        # nuclear steady-state value
  steady_c_raw = steady_n_raw * tau / lambda_c                # cytosolic steady-state value
  n_counts = sapply(1:time_points, function(t) {              # nuclear transcript counts based on relative abundance
    round(library_size_n[t] * steady_n_raw / rna_total_n[t], digits = 0)   
  })
  c_counts = sapply(1:time_points, function(t) {              # cytosolic transcript counts based on relative abundance
    round(library_size_c[t] * steady_c_raw / rna_total_c[t], digits = 0)
  })
  n_count_samples <- sapply(1:number, function(i) {           # sampling nuclear transcript counts (rows: time points,
    sapply(1:time_points, function(t) {                       # columns: repetitions)
      rpois(1, n_counts[t])
    })})
  c_count_samples <- sapply(1:number, function(i) {           # sampling cytosolic transcript counts (rows: time points,
    sapply(1:time_points, function(t) {                       # columns: repetitions)
      rpois(1, c_counts[t])     
    })})
  
  # sampling modified nuclear and cytosolic transcript counts
  
  mod_counts <- sample_transcripts(n_counts, c_counts, raw_q_nucleus, raw_q_cytosol, p_un, p_mod, 
                                   number = number)
  n_mod_samples <- mod_counts[[1]]
  c_mod_samples <- mod_counts[[2]]
  
  # returning
  
  return(list(n_count_samples, n_mod_samples, c_count_samples, c_mod_samples))
}

### RNA Dynamics Model Simulation #####################################################################################

rna_dynamics <- function(time, state, params) {
  
  # *******************************************************************************************************************
  
  #   function (time sequence (t), state (y), parameters (params))
  #   -> return list(c(changes), global variables): [(changes), global1, global2, ...]
  
  # *******************************************************************************************************************
  
  with(as.list(c(state, params)), {
    drna_pre_n = - (export + deg_n) * rna_pre_n;
    drna_pre_c = export * rna_pre_n - deg_c * rna_pre_c
    drna_new_n = syn - (export + deg_n) * rna_new_n
    drna_new_c = export * rna_new_n - deg_c * rna_new_c
    return(list(c(drna_pre_n, drna_pre_c, drna_new_n, drna_new_c)))
  })
}

pre_n_sol <- function(time, pre_n_start, export, deg_n) {
  result = pre_n_start  * exp(- (deg_n + export) * time)
  return(result)
}

pre_c_sol <- function(time, pre_n_start, pre_c_start, export, deg_n, deg_c) {
  result = exp(- deg_c * time) * pre_c_start + 
    export / (deg_c - deg_n - export) * (
      pre_n_start * (exp(- (deg_n + export) * time) - exp(- deg_c * time)))
  return(result)
}

new_n_sol <- function(time, new_n_start, syn, export, deg_n) {
  result = syn / (deg_n + export) * (1 - exp(- (deg_n + export) * time)) + 
    new_n_start * exp(- (deg_n + export) * time)
  return(result)
}

new_c_sol <- function(time, new_n_start, new_c_start, syn, export, deg_n, deg_c) {
  result = syn / (deg_n + export) * export / deg_c + exp(- deg_c * time) * new_c_start + 
    export / (deg_c - deg_n - export) * (
      new_n_start * (exp(- (deg_n + export) * time) - exp(- deg_c * time)) +
        syn * (exp(- deg_c * time) / deg_c - exp(- (deg_n + export) * time) / (deg_n + export))
    )
  return(result)
}

ratio_n_sol <- function(time, export, deg_n) {
  result = 1 - exp(- (deg_n + export) * time)
  return(result)
}
  
ratio_c_sol <- function(time, export, deg_n, deg_c) {
  result = 1 - (exp(- deg_c * time) + 
                  deg_c / (deg_c - deg_n - export) * (exp(- (deg_n + export) * time) - exp(- deg_c * time)))
  return(result)
}

multi_rna_dynamics_sim <- function(syn_rate, exp_rate, deg_rate_n, deg_rate_c, time_series,
                                   rna_pre_n = NULL, rna_pre_c = NULL, rna_new_n = NULL, rna_new_c = NULL,
                                   rna_species_labels = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # syn_rate:     vector; RNA synthesis rates (per RNA species)
  # exp_rate:     vector; RNA export rate (per RNA species)
  # deg_rate_n:   vector; RNA degradation rate (nucleus) (per RNA species)
  # deg_rate_c:   vector; RNA degradation rate (cytosol) (per RNA species)
  # time_series:  vector; time points at which measurements were taken
  
  # rna_pre_n:    optional; vector; pre-existing RNA transcripts within the nucleus (per RNA species);
  #               if NULL, steady-state transcripts are used (assuming that no new transcripts are existing at t = 0)
  #               (default: NULL)
  # rna_pre_c:    optional; vector; pre-existing RNA transcripts within the cytosol (per RNA species);
  #               if NULL, steady-state transcripts are used (assuming that no new transcripts are existing at t = 0)
  #               (default: NULL)
  # rna_new_n:    optional; vector; newly synthesized RNA transcripts within the nucleus (per RNA species);
  #               if NULL, it is assumed that no transcripts exist at t = 0 (rna_new_n is set to a vector of 0s)
  #               (default: NULL)
  # rna_new_c:    optional; vector; newly synthesized RNA transcripts within the cytosol (per RNA species);
  #               if NULL, it is assumed that no transcripts exist at t = 0 (rna_new_c is set to a vector of 0s)
  #               (default: NULL)
  # rna_species_labels:   optional; labels to use for RNA species within the plots (default: NULL)
  
  # *** NOTE ***
  
  # As soon as any of the parameters rna_pre_n, rna_pre_c, rna_new_n and rna_new_c evaluates to NULL, steady-state
  # transcripts with the assumption that no new transcripts exist at t = 0 are used; in this case, other values given
  # for these parameters are ignored.
  
  # *** RETURN ***
  
  # a list containing (in the following order): 
  # pre-existing nuclear transcripts' data frame, pre-existing cytosolic transcripts' data frame, 
  # newly synthesized nuclear transcripts' data frame, newly synthesized cytosolic transcripts' data frame, 
  # nuclear transcripts ratio's data frame, cytosolic transcripts ratio's data frame, corresponding plots in the same
  # order as data frames
  
  # This function simulates RNA dynamics for a two-compartment model both using an ODE system simulation. 
  # It creates 6 plots (RNA transcripts over time) if parameter 'plotting' is set to TRUE:
  # nucleus, pre-existing; nucleus, newly synthesized; cytosol, pre-existing; cytosol, newly synthesized;
  # ratio, nucleus new and nucleus total; ratio, cytosol new and cytosol total
  
  # *******************************************************************************************************************
  
  # *******************************************************************************************************************
  # INITIALIZATIONS
  # *******************************************************************************************************************
  
  n = length(syn_rate)                        # sample size
  n_plus_1 = n + 1                            # sample size + 1
  
  # initializing ratio data frames' column names
  
  ratio_n_colnames = vector(mode = "character", length = n)
  ratio_c_colnames = vector(mode = "character", length = n)
  
  # initializing output data frames (one per RNA entity; rows refer to time points, columns to RNA species)
  
  pre_n_df = data.frame("time" = time_series)
  pre_c_df = data.frame("time" = time_series)
  new_n_df = data.frame("time" = time_series)
  new_c_df = data.frame("time" = time_series)
  
  # *******************************************************************************************************************
  # RNA DYNAMICS SIMULATIONS
  # *******************************************************************************************************************
  
  for (i in 1:n) {    # iterating through RNA species
    
    # defining RNA species' column names (to distinguish between outputs from the different RNA species)
    
    pre_n_string <- paste("rna_pre_n", rna_species_labels[i], sep = "_")
    pre_c_string <- paste("rna_pre_c", rna_species_labels[i], sep = "_")
    new_n_string <- paste("rna_new_n", rna_species_labels[i], sep = "_")
    new_c_string <- paste("rna_new_c", rna_species_labels[i], sep = "_")
    ratio_n_colnames[i] <- paste("rna_n", rna_species_labels[i], sep = "_")
    ratio_c_colnames[i] <- paste("rna_c", rna_species_labels[i], sep = "_")
    
    # initializing state variables
    
    state <- c(rna_pre_n = 0, 
               rna_pre_c = 0,
               rna_new_n = 0.,
               rna_new_c = 0.)
    
    # in case any state parameter evaluates to NULL, assigning steady-state transcript numbers 
    
    if (is.null(rna_pre_n) | is.null(rna_pre_c) | is.null(rna_new_n) | is.null(rna_pre_c))  {
      state["rna_pre_n"] = syn_rate[i] / (exp_rate[i] + deg_rate_n[i])
      state["rna_pre_c"] = (syn_rate[i] / (exp_rate[i] + deg_rate_n[i])) * exp_rate[i]/ deg_rate_c[i]   
      state["rna_new_n"] = 0
      state["rna_new_c"] = 0
    }
    
    # else, assigning state parameter values
    
    else {
      state["rna_pre_n"] = rna_pre_n[i]
      state["rna_pre_c"] = rna_pre_c[i]
      state["rna_new_n"] = rna_new_n[i]
      state["rna_new_c"] = rna_new_c[i]
    }
    
    # defining parameter values
    
    params <- c(syn = syn_rate[i],
                export = exp_rate[i],
                deg_n = deg_rate_n[i],
                deg_c = deg_rate_c[i])
    
    # simulating ODE system, transforming output to data frame object
    
    out <- ode(state, time_series, rna_dynamics, params)
    out <- data.frame(out)
    colnames(out) <- c("time", pre_n_string, pre_c_string, new_n_string, new_c_string)
    
    # appending output to output data frame
    
    pre_n_df[, pre_n_string] <- out[, pre_n_string]
    pre_c_df[, pre_c_string] <- out[, pre_c_string]
    new_n_df[, new_n_string] <- out[, new_n_string]
    new_c_df[, new_c_string] <- out[, new_c_string]
    
  }
  
  # *******************************************************************************************************************
  # COMPUTING RATIOS
  # *******************************************************************************************************************
  
  ratio_n_df = data.frame("time" = new_n_df$time, 
                          new_n_df[, 2:n_plus_1] / (new_n_df[, 2:n_plus_1] + pre_n_df[, 2:n_plus_1]))
  colnames(ratio_n_df) <- c("time", ratio_n_colnames)
  ratio_c_df = data.frame("time" = new_c_df$time, 
                          new_c_df[, 2:n_plus_1] / (new_c_df[, 2:n_plus_1] + pre_c_df[, 2:n_plus_1]))
  colnames(ratio_c_df) <- c("time", ratio_c_colnames)
  
  # *******************************************************************************************************************
  # PLOTTING
  # *******************************************************************************************************************
  
  if (is.null(rna_species_labels)) {
    rna_species_labels = c("time", sapply(1:n, function(l){as.character(l)}))
  }
  else {
    rna_species_labels = c("time", rna_species_labels)
  }
  
  pre_n_plot <- plot_dynamics(pre_n_df, rna_species_labels, "nucleus, pre-existing", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  pre_c_plot <- plot_dynamics(pre_c_df, rna_species_labels, "cytosol, pre-existing", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  new_n_plot <- plot_dynamics(new_n_df, rna_species_labels, "nucleus, newly synthesized", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  new_c_plot <- plot_dynamics(new_c_df, rna_species_labels, "cytosol, newly synthesized", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  ratio_n_plot <- plot_dynamics(ratio_n_df, rna_species_labels, "nucleus, ratio (new / total)", "time", "RNA ratio", 
                                "RNA_species", colours_3)
  ratio_c_plot <- plot_dynamics(ratio_c_df, rna_species_labels, "cytosol, ratio (new / total)", "time", "RNA ratio", 
                                "RNA_species", colours_3)
  
  # *******************************************************************************************************************
  # RETURNING
  # *******************************************************************************************************************
  
  return(list(pre_n_df, pre_c_df, new_n_df, new_c_df, ratio_n_df, ratio_c_df, 
              pre_n_plot, pre_c_plot, new_n_plot, new_c_plot, ratio_n_plot, ratio_c_plot))
}

multi_rna_dynamics_sol <- function(syn_rate, exp_rate, deg_rate_n, deg_rate_c, time_series,
                                   rna_pre_n = NULL, rna_pre_c = NULL, rna_new_n = NULL, rna_new_c = NULL,
                                   rna_species_labels = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # syn_rate:     vector; RNA synthesis rates (per RNA species)
  # exp_rate:     vector; RNA export rate (per RNA species)
  # deg_rate_n:   vector; RNA degradation rate (nucleus) (per RNA species)
  # deg_rate_c:   vector; RNA degradation rate (cytosol) (per RNA species)
  # time_series:  vector; time points at which measurements were taken
  
  # rna_pre_n:    optional; vector; pre-existing RNA transcripts within the nucleus (per RNA species);
  #               if NULL, steady-state transcripts are used (assuming that no new transcripts are existing at t = 0)
  #               (default: NULL)
  # rna_pre_c:    optional; vector; pre-existing RNA transcripts within the cytosol (per RNA species);
  #               if NULL, steady-state transcripts are used (assuming that no new transcripts are existing at t = 0)
  #               (default: NULL)
  # rna_new_n:    optional; vector; newly synthesized RNA transcripts within the nucleus (per RNA species);
  #               if NULL, it is assumed that no transcripts exist at t = 0 (rna_new_n is set to a vector of 0s)
  #               (default: NULL)
  # rna_new_c:    optional; vector; newly synthesized RNA transcripts within the cytosol (per RNA species);
  #               if NULL, it is assumed that no transcripts exist at t = 0 (rna_new_c is set to a vector of 0s)
  #               (default: NULL)
  # rna_species_labels:   optional; labels to use for RNA species within the plots (default: NULL)
  
  # *** NOTE ***
  
  # As soon as any of the parameters rna_pre_n, rna_pre_c, rna_new_n and rna_new_c evaluates to NULL, steady-state
  # transcripts with the assumption that no new transcripts exist at t = 0 are used; in this case, other values given
  # for these parameters are ignored.
  
  # *** RETURN ***
  
  # a list containing (in the following order): 
  # pre-existing nuclear transcripts' data frame, pre-existing cytosolic transcripts' data frame, 
  # newly synthesized nuclear transcripts' data frame, newly synthesized cytosolic transcripts' data frame, 
  # nuclear transcripts ratio's data frame, cytosolic transcripts ratio's data frame, corresponding plots in the same
  # order as data frames
  
  # This function simulates RNA dynamics for a two-compartment model both using an ODE system solution. 
  # It creates 6 plots (RNA transcripts over time) if parameter 'plotting' is set to TRUE:
  # nucleus, pre-existing; nucleus, newly synthesized; cytosol, pre-existing; cytosol, newly synthesized;
  # ratio, nucleus new and nucleus total; ratio, cytosol new and cytosol total
  
  # *******************************************************************************************************************
  
  # *******************************************************************************************************************
  # INITIALIZATIONS
  # *******************************************************************************************************************
  
  n = length(syn_rate)                        # sample size
  n_plus_1 = n + 1                            # sample size + 1
  
  # initializing ratio data frames' column name vectors
  
  ratio_n_colnames = vector(mode = "character", length = n)
  ratio_c_colnames = vector(mode = "character", length = n)
  
  # initializing output data frames (one per RNA entity; rows refer to time points, columns to RNA species)
  
  pre_n_df = data.frame("time" = time_series)
  pre_c_df = data.frame("time" = time_series)
  new_n_df = data.frame("time" = time_series)
  new_c_df = data.frame("time" = time_series)
  
  # *******************************************************************************************************************
  # RNA DYNAMICS SIMULATIONS
  # *******************************************************************************************************************
  
  for (i in 1:n) {    # iterating through RNA species
    
    # defining RNA species' column names (to distinguish between outputs from the different RNA species)
    
    pre_n_string <- paste("rna_pre_n", rna_species_labels[i], sep = "_")
    pre_c_string <- paste("rna_pre_c", rna_species_labels[i], sep = "_")
    new_n_string <- paste("rna_new_n", rna_species_labels[i], sep = "_")
    new_c_string <- paste("rna_new_c", rna_species_labels[i], sep = "_")
    ratio_n_colnames[i] <- paste("rna_n", rna_species_labels[i], sep = "_")
    ratio_c_colnames[i] <- paste("rna_c", rna_species_labels[i], sep = "_")
    
    # initializing state variables
    
    state <- c(rna_pre_n = 0, 
               rna_pre_c = 0,
               rna_new_n = 0.,
               rna_new_c = 0.)
    
    # in case any state parameter evaluates to NULL, assigning steady-state transcript numbers 
    
    if (is.null(rna_pre_n) | is.null(rna_pre_c) | is.null(rna_new_n) | is.null(rna_pre_c))  {
      state["rna_pre_n"] = syn_rate[i] / (exp_rate[i] + deg_rate_n[i])
      state["rna_pre_c"] = (syn_rate[i] / (exp_rate[i] + deg_rate_n[i])) * exp_rate[i]/ deg_rate_c[i]   
      state["rna_new_n"] = 0
      state["rna_new_c"] = 0
    }
    
    # else, assigning state parameter values
    
    else {
      state["rna_pre_n"] = rna_pre_n[i]
      state["rna_pre_c"] = rna_pre_c[i]
      state["rna_new_n"] = rna_new_n[i]
      state["rna_new_c"] = rna_new_c[i]
    }
    
    # defining parameter values
    
    params <- c(syn = syn_rate[i],
                export = exp_rate[i],
                deg_n = deg_rate_n[i],
                deg_c = deg_rate_c[i])
    
    # computing analytical solutions
    
    pre_n_df[, pre_n_string] = pre_n_sol(time_series, state["rna_pre_n"], params["export"], params["deg_n"])
    pre_c_df[, pre_c_string] = pre_c_sol(time_series, state["rna_pre_n"], state["rna_pre_c"], params["export"],
                                         params["deg_n"], params["deg_c"])
    new_n_df[, new_n_string] = new_n_sol(time_series, state["rna_new_n"], params["syn"], params["export"],
                                         params["deg_n"])
    new_c_df[, new_c_string] = new_c_sol(time_series, state["rna_new_n"], state["rna_new_c"], params["syn"],
                                         params["export"], params["deg_n"], params["deg_c"])
  }
  
  # *******************************************************************************************************************
  # COMPUTING RATIOS
  # *******************************************************************************************************************
  
  ratio_n_df = data.frame("time" = new_n_df$time, 
                          new_n_df[, 2:n_plus_1] / (new_n_df[, 2:n_plus_1] + pre_n_df[, 2:n_plus_1]))
  colnames(ratio_n_df) <- c("time", ratio_n_colnames)
  ratio_c_df = data.frame("time" = new_c_df$time, 
                          new_c_df[, 2:n_plus_1] / (new_c_df[, 2:n_plus_1] + pre_c_df[, 2:n_plus_1]))
  colnames(ratio_c_df) <- c("time", ratio_c_colnames)
  
  # *******************************************************************************************************************
  # PLOTTING
  # *******************************************************************************************************************
  
  if (is.null(rna_species_labels)) {
    rna_species_labels = c("time", sapply(1:n, function(l){as.character(l)}))
  }
  else {
    rna_species_labels = c("time", rna_species_labels)
  }
  
  pre_n_plot <- plot_dynamics(pre_n_df, rna_species_labels, "nucleus, pre-existing", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  pre_c_plot <- plot_dynamics(pre_c_df, rna_species_labels, "cytosol, pre-existing", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  new_n_plot <- plot_dynamics(new_n_df, rna_species_labels, "nucleus, newly synthesized", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  new_c_plot <- plot_dynamics(new_c_df, rna_species_labels, "cytosol, newly synthesized", "time", "RNA abundance", 
                              "RNA_species", colours_3)
  ratio_n_plot <- plot_dynamics(ratio_n_df, rna_species_labels, "nucleus, ratio (new / total)", "time", "RNA ratio", 
                                "RNA_species", colours_3)
  ratio_c_plot <- plot_dynamics(ratio_c_df, rna_species_labels, "cytosol, ratio (new / total)", "time", "RNA ratio", 
                                "RNA_species", colours_3)
  
  # *******************************************************************************************************************
  # RETURNING
  # *******************************************************************************************************************
  
  return(list(pre_n_df, pre_c_df, new_n_df, new_c_df, ratio_n_df, ratio_c_df, 
              pre_n_plot, pre_c_plot, new_n_plot, new_c_plot, ratio_n_plot, ratio_c_plot))
}

### RNA Dynamics Parameter Estimation #################################################################################

traf <- asin(sqrt(ratio))

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

get_grid_maxlike <- function(summed_probs, overall_probs, mu_grid_extended, tau_grid_extended,
                             lambda_n_grid_extended, lambda_c_grid_extended) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summed_probs:     list of vectors, each vector containing one repetition's likelihoods of the grid parameter sets
  # overall_probs:    likelihoods of grid parameter sets, summed over the repetitions (vector)
  # x_grid_extended:  RNA dynamics parameter x's values used within the grid (this is NOT the values to be used in
  #                   general, but extended to the grid actually used during search)
  
  # *** RETURN ***
  
  # a list containing (in the following order): 
  # four vectors containing maximum likelihood mu / tau / lambda_n / lambda_c values for each repetition
  # four vectors containing maximum likelihood mu / tau / lambda_n / lambda_c values over all repetitions
  
  # *******************************************************************************************************************
  
  repetitions <- length(summed_probs)     # getting number of repetitions
  
  max_likelihoods <- sapply(1:repetitions, function(i) {    # for each repetition, getting maximum likelihood
    max(unlist(summed_probs[i]))
  })
  max_occurrences <- sapply(1:repetitions, function(i) {    # for each repetition, getting all parameter settings with
    list(which(summed_probs[[i]] %in% c(max_likelihoods)))  #  maximum likelihood
  })
  
  max_mu <- c(sapply(1:repetitions, function(i) {           # maximum likelihood mu-values (for each repetition)
    mu_grid_extended[max_occurrences[[i]]]
  })) 
  max_tau <- c(sapply(1:repetitions, function(i) {          # maximum likelihood tau-values (for each repetition)
    tau_grid_extended[max_occurrences[[i]]]
  })) 
  max_ln <- c(sapply(1:repetitions, function(i) {           # maximum likelihood lambda_n-values (for each repetition)
    lambda_n_grid_extended[max_occurrences[[i]]]
  })) 
  max_lc <- c(sapply(1:repetitions, function(i) {           # maximum likelihood lambda_c-values (for each repetition)
    lambda_c_grid_extended[max_occurrences[[i]]]
  }))
  
  overall_max_likelihood <- max(overall_probs)                                      # overall maximum likelihood
  overall_max_occurrences <- which(overall_probs %in% c(overall_max_likelihood))    # overall max. likelihood indices
  overall_max_mu <- mu_grid_extended[overall_max_occurrences]                       # overall max. likelihood mu's
  overall_max_tau <- tau_grid_extended[overall_max_occurrences]                     # overall max. likelihood tau's
  overall_max_ln <- lambda_n_grid_extended[overall_max_occurrences]                 # overall max. likelihood lambda_n
  overall_max_lc <- lambda_c_grid_extended[overall_max_occurrences]                 # overall max. likelihood lambda_c
  
  return(list(max_mu, max_tau, max_ln, max_lc, overall_max_mu, overall_max_tau, overall_max_ln, overall_max_lc))
}

p_nmod_nnew <- function(n_mod, n_pre, n_new, p_un, p_mod, mode = "binom", betabin_od = 2) {
  
  # *******************************************************************************************************************
  # n_mod:    modified transcript counts
  # n_pre:    pre-existing transcript counts
  # n_new:    newly synthesized transcript counts
  # p_un:     probability that a pre-existing transcript remains unmodified
  # p_mod:    probability that a newly synhesized transcript is modified
  
  # mode:         set to "binom" for an underlying binomial distribution, set to "betabinom" for an underlying beta-
  #               binomial distribution (default: "binom")
  # betabin_od:   needed for 'mode' == 'betabinom'; beta-binomial distribution overdispersion parameter (alpha + beta)
  #               (default: 2)
  
  # This function computes the probability of the modified transcript counts observed, in dependence of both newly
  # synthesized and pre-existing transcripts counts and global parameters.
  
  # *******************************************************************************************************************
  
  n_new_range <- 0:n_new  # n_new values which are summed over
  if (mode == "binom") {
    result <- sum(dbinom(n_new_range, n_new, p_mod) * dbinom(n_new_range + n_pre - n_mod, n_pre, p_un)) 
  }
  else if (mode == "betabinom") {
    result <- sum(dbetabinom(n_new_range, n_new, p_mod, betabin_od) * 
                    dbetabinom(n_new_range + n_pre - n_mod, n_pre, p_un, betabin_od))
  }
  
  return(result)
}

p_nmod_nnew_matrix <- function(mod_samples, n_total, repetitions, time_points, p_un, p_mod,
                               mode = "binom", betabin_od = 2) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # mod_samples:    matrix containing modified transcript counts (cols: repetitions, rows: time points)
  # n_total:        matrix containing total number of transcript counts (cols: repetitions, rows: time points)
  # repeptitions:   number of data samples (within mod_samples and n_total)
  # time_points:    number of time points (within mod_samples and n_total)
  # p_un:           vector; probability that a pre-existing transcript remains unmodified (resolved by time)
  # p_mod:          vector; probability that a newly synhesized transcript is modified (resolved by time)
  
  # mode:         set to "binom" for an underlying binomial distribution, set to "betabinom" for an underlying beta-
  #               binomial distribution (default: "binom")
  # betabin_od:   needed for 'mode' == 'betabinom'; beta-binomial distribution overdispersion parameter (alpha + beta)
  #               (default: 2)
  
  # *** RETURN ***
  
  # a nested list containing the modified transcript counts' probabilities for each repetition (first list index:
  # repetitions, second list index: time points, third vector index: newly synthesized counts' probabilities)
  
  # *******************************************************************************************************************
  
  p_nmod_nnew_results <- 
    sapply(1:repetitions, function(i) { list(   # for each repetition
      sapply(1:time_points, function(t) {          # for each time point
        total_ti <- n_total[t, i]
        mod_ti <- mod_samples[t, i]
        p_un_t <- p_un[t]
        p_mod_t <- p_mod[t]
        list(
          sapply(0:total_ti, function(nnew) {         # for each possible newly synthesized count,
            p_nmod_nnew(mod_ti, total_ti - nnew, nnew, p_un_t, p_mod_t, 
                        mode = mode, betabin_od = betabin_od)   # compute probability
          })   
        )
      })             
  )})
  return(p_nmod_nnew_results)
}

p_nnew_ntotal <- function(n_new, n_total, p_new, mode = "binom", betabin_od = 2) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_new:        newly synthesized transcript counts (a single value or a vector of values)
  # n_total:      total number of transcript counts
  # p_new:        ratio of newly synthesized transcripts and the total amount of transcripts (according to RNA
  #               dynamics parameters)
  
  # mode:         set to "binom" for an underlying binomial distribution, set to "betabinom" for an underlying beta-
  #               binomial distribution (default: "binom")
  # betabin_od:   needed for 'mode' == 'betabinom'; beta-binomial distribution overdispersion parameter (alpha + beta)
  #               (default: 2)
  
  # *** RETURN ***
  
  # probablity of observing the newly synthesized transcript counts (a single value or a vector of values)
  
  # This function computes the probability of newly synthesized transcript counts, in dependence of both total
  # transcript counts and RNA dynamics parameters.
  
  # *******************************************************************************************************************
  
  if (mode == "binom") {return(dbinom(n_new, n_total, p_new))}
  else if (mode == "betabinom") {return(dbetabinom(n_new, n_total, p_new, betabin_od))}
  else {return(NULL)}
  
}

p_nmod_theta <- function(n_mod, n_total, p_new, p_un, p_mod, p_mod_pre = NULL,
                         mode = "binom", betabin_od = 2) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_mod:        modified transcript counts
  # n_total:      total number of transcript counts
  # p_new:        ratio of newly synthesized transcripts and the total amount of transcripts (according to RNA
  #               dynamics parameters)
  # p_un:         probability that a pre-existing transcript remains unmodified
  # p_mod:        probability that a newly synhesized transcript is modified
  
  # p_mod_pre:    optional; vector of pre-computed modified transcript counts' probabilities, given newly synthesized
  #               transcript counts, resolved by newly synthesized transcript counts (default: NULL)
  # mode:         set to "binom" for an underlying binomial distribution, set to "betabinom" for an underlying beta-
  #               binomial distribution (default: "binom")
  # betabin_od:   needed for 'mode' == 'betabinom'; beta-binomial distribution overdispersion parameter (alpha + beta)
  #               (default: 2)
  
  # *** RETURN ***
  
  # probablity of observing the modified transcript counts
  
  # This function computes the probability of the modified transcript counts observed, in dependence of both RNA
  # dynamics parameters and global parameters.
  
  # *******************************************************************************************************************
  
  # n_new values which are summed over
  n_new_range <- (0:n_total)    
  # vector containing the first factor of the sum (probability of modified counts given newly synthesized counts)
  if (is.null(p_mod_pre)) {     
    p_mod_pre <- sapply(n_new_range, function(x){
      p_nmod_nnew(n_mod, n_total - x, x, p_un, p_mod, mode = mode, betabin_od = betabin_od)
    })
  }
  # vector containing second factor of sum (probability of newly synthesized counts given new/total ratio)
  p_new_pre <- p_nnew_ntotal(n_new_range, n_total, p_new, mode = mode, betabin_od = betabin_od)

  # summing over (vector containing first factor of the sum * vector containing the second factor of the sum)
  result <- sum(p_mod_pre * p_new_pre) 
  return(result)
}

prob_thetaderiv_1 <- function(n_total, n_mod, lambda_n_tau, time_series, p_un, p_mod, 
                              p_mod_pre = NULL, mode = "binom", betabin_od = 2, log_scale = FALSE) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_total:          nuclear total transcript counts observed, resolved by time (vector)
  # n_mod:            nuclear modified transcript counts observed, resolved by time (vector)
  # lambda_n_tau:     (lambda_n + tau) value suggested
  # time_series:      time points at which measurements were taken (vector)
  # p_un:             probability that a pre-existing transcript remains unmodified, resolved by time (vector)
  # p_mod:            probability that a newly synhesized transcript is modified, resolved by time (vector)
  
  # p_mod_pre:        optional; list of pre-computed nuclear modified transcript counts' probability vectors, given   
  #                   newly synthesized transcript counts (list index: time points, vector index: nnew values)
  #                   (default: NULL)
  # mode:             set to "binom" for an underlying binomial distribution, set to "betabinom" for an underlying 
  #                   beta-binomial distribution (default: "binom")
  # betabin_od:       needed for 'mode' == 'betabinom'; beta-binomial distribution overdispersion parameter
  #                   (alpha + beta) (default: 2)
  # log_scale:        optional; set to TRUE to log-scale probabilities (default: FALSE)
  
  # *** RETURN ***
  
  # conditional probability of the suggested parameter value, given the observed data (summed over time)
  
  # *******************************************************************************************************************
  
  time_points <- length(time_series)            # number of time points
  lnt_split <- lambda_n_tau / 2                 # splitting lambda_n + tau into lambda_n and tau
  
  # vector of nuclear ratios (new / total) resolved by time points
  sim_ratio_nucleus <- ratio_n_sol(time_series, lnt_split, lnt_split)       
  
  # pre-computing matrix of modified transcript counts' probabilities, given newly synthesized transcript counts, 
  # resolved by time points
  if (is.null(p_mod_pre)) {                     
    p_mod_pre <- p_nmod_nnew_matrix(matrix(n_mod, ncol = 1), matrix(n_total, ncol = 1), 
                                    1, time_points, p_un, p_mod, mode = mode, betabin_od = betabin_od)[[1]]
  }
  
  # conditional probability, for each time point
  probs <- sapply(1:time_points, function(t) {
    p_nmod_theta(n_mod[t], n_total[t], sim_ratio_nucleus[t], p_un[t], p_mod[t], p_mod_pre[[t]],
                 mode = mode, betabin_od = betabin_od)
  })

  # log scale: returning summed probability (over all time points)
  if (log_scale) {
    probs <- log10(probs)
    probs[which(probs == -Inf)] <- .Machine$integer.max * -1
    return(sum(probs + 1))
  }
  # normal scale: returning mutiplied probability (over all time points)
  else {return(prod(probs))}
}

prob_thetaderiv_2 <- function(comp_mode, lambda_n_tau, lambda_c, time_series, lookup_table = NULL, 
                              c_total = NULL, c_mod = NULL, p_un = NULL, p_mod = NULL, p_mod_pre = NULL,
                              dist_mode = "binom", betabin_od = 2, log_scale = FALSE) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # comp_mode:        computation of conditional probability; choose 'denovo' for computing the probability (requires:
  #                   c_total, c_mod, p_un, p_mod), choose 'lookup' to retrieve probability from a loookup table
  #                   (requires: lookup_table)
  # lambda_n_tau:     (lambda_n + tau) value suggested
  # lambda_c:         lambda_c value suggested
  # time_series:      time points at which measurements were taken (vector)
  
  # lookup_table:     required in mode == 'lookup'; function for conditional probability of lambda_n + tau and 
  #                   lambda_c, interpolated over new/total transcript counts ratio, resolved by time (list)
  
  # c_total:          required in comp_mode == 'denovo'; cytosolic total transcript counts observed, resolved by time 
  #                   (vector)
  # c_mod:            required in comp_mode == 'denovo'; cytosolic modified transcript counts observed, resolved by time
  #                   (vector)
  # p_un:             required in comp_mode == 'denovo'; probability that a pre-existing transcript remains unmodified,
  #                   resolved by time (vector)
  # p_mod:            required in comp_mode == 'denovo'; probability that a newly synhesized transcript is modified,
  #                   resolved by time (vector)
  
  # p_mod_pre:        optional in comp_mode == 'denovo'; list of pre-computed nuclear modified transcript counts'    
  #                   probability vectors, given newly synthesized transcript counts (list index: time points, vector 
  #                   index: nnew values) (default: NULL)
  # dist_mode:        optional in comp_mode == 'denovo'; set to "binom" for an underlying binomial distribution, 
  #                   set to "betabinom" for an underlying beta-binomial distribution (default: "binom")
  # betabin_od:       optional in 'dist_mode' == 'betabinom'; beta-binomial distribution overdispersion parameter 
  #                   (alpha + beta) (default: 2)
  # log_scale:        optional; set to TRUE to log-scale probabilities (comp_mode == 'denovo') or indicate that
  #                   lookup table returns log-scaled probabilities (comp_mode == "lookup") (default: FALSE)
  
  # *** RETURN ***
  
  # conditional probability of the suggested parameter values, given the observed data (summed over time)
  
  # This function returns the conditional probability of derived parameter lambda_c. It either computes it de novo
  # (setting parameter mode = 'denovo') or returns the value from a lookup table (setting parameter mode = 'lookup').
  # Using the lookup table is especially helpful to reduce runtime (summing over all possible n_new is avoided).
  
  # *******************************************************************************************************************

  time_points = length(time_series)     # number of time points
  lnt_split <- lambda_n_tau / 2         # splitting lambda_n + tau in lambda_n and tau
  
  # vector of cytosolic ratios (new / total) resolved by time points
  sim_ratio_cytosol <- ratio_c_sol(time_series, lnt_split, lnt_split, lambda_c)           
  p <- NULL                             # conditional probability (initialized with NULL)

  # MODE = LOOKUP conditional probability
  
  if (comp_mode == "lookup") {
    probs <- sapply(1:time_points, function(t) {    # conditional probabilities resolved by time
      lookup_table[[t]](sim_ratio_cytosol[t])                                
    })
    p <- NULL
    if (log_scale) {p <- sum(probs)}                # log scale: conditional probabilities summed over time
    else {p <- prod(probs)}                         # normal scale: conditional probabilities multiplied over time
  }
  
  # MODE = DENOVO conditional probability
  
  else if (comp_mode == "denovo") {
    # pre-computing matrix of modified transcript counts' probabilities, given newly synthesized transcript counts, 
    # resolved by time points
    if (is.null(p_mod_pre)) {                     
      p_mod_pre <- p_nmod_nnew_matrix(matrix(c_mod, ncol = 1), matrix(c_total, ncol = 1), 
                                      1, time_points, p_un, p_mod, mode = dist_mode, betabin_od = betabin_od)[[1]]
    }
    # conditional probability, for each time point
    probs <- sapply(1:time_points, function(t){   
      p_nmod_theta(c_mod[t], c_total[t], sim_ratio_cytosol[t], p_un[t], p_mod[t], p_mod_pre[[t]],
                   mode = dist_mode, betabin_od = betabin_od)   
    })
    if (log_scale) {
      probs <- log10(probs)
      probs[which(probs == -Inf)] = .Machine$integer.max * -1
      p <- sum(probs)}              # log scale: conditional probabilities summed over time
    else {p <- prod(probs)}         # normal scale: conditional probabilities multiplied over time
  }

  # returning conditional probability
  return(p)
}

prob_thetaderiv_3 <- function(n_total, n_rna, n_libsize, steady_n) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_total:    nuclear total transcript counts observed, resolved by time (vector)
  # n_rna:      nuclear total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # n_libsize:  nuclear total amount of read counts (summed over all transcript species), resolved by time (vector)
  # steady_n:   steady_n parameter value suggested
  
  # *** RETURN ***
  
  # conditional probability of the suggested parameter value, given the observed data (summed over time)
  
  # *******************************************************************************************************************
  
  time_points <- length(n_total)  # number of time points
  
  # nuclear transcript counts according to suggested parameter value
  
  sim_n_total <- sapply(n_libsize, function(l) {
    round(steady_n * l / n_rna, digits = 0)
  })
  
  # probability of the observed count data
  
  probs <- sapply(1:time_points, function(t) {
    dpois(sim_n_total[t], n_total[t])
  })
  
  # returning
  
  return(sum(probs))
}

prob_thetaderiv_4 <- function(c_total, c_rna, c_libsize, steady_c) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # c_total:    cytosolic total transcript counts observed, resolved by time (vector)
  # c_rna:      cytosolic total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # c_libsize:  cytosolic total amount of read counts (summed over all transcript species), resolved by time (vector)
  # steady_c:   steady_c parameter value suggested
  
  # *** RETURN ***
  
  # conditional probability of the suggested parameter value, given the observed data (summed over time)
  
  # *******************************************************************************************************************
  
  time_points <- length(c_total)    # number of time points
  
  # cytosolic transcript counts according to suggested parameter value
  
  sim_c_total <- sapply(c_libsize, function(l) {
    round(steady_c * l / c_rna, digits = 0)
  })
  
  # probability of the observed count data
  
  probs <- sapply(1:time_points, function(t) {
    dpois(sim_c_total[t], c_total[t])
  })
  
  # returning
  
  return(sum(probs))
}

probability_interpolation <- function(lntau_borders, lc_borders, steadyn_borders, steadyc_borders,
                                      n_total, n_mod, n_rna, n_libsize, c_total, c_mod, c_rna, c_libsize,
                                      time_series, p_un, p_mod, new_by_total_resolution = 1000, 
                                      ip_points = 100, search_points = 1000, width = c(10, 10, 10)) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # lntau_borders:    vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:       vector with (minimum, maximum) lambda_c value 
  # steadyn_borders:  vector with (minimum, maximum) steady_n value 
  # steadyc_borders:  vector with (minimum, maximum) steady_c value 
  
  # n_total:      nuclear total counts observed for the transcript species, resolved by time (vector)
  # n_mod:        nuclear modified counts observed for the transcrip species, resolved by time (vector)
  # n_rna:        nuclear total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # n_libsize:    nuclear total amount of read counts (summed over all transcript species), resolved by time (vector)
  # c_total:      cytosolic total counts observed for the transcript species, resolved by time (vector)
  # c_mod:        cytosolic modified counts observed for the transcript species, resolved by time (vector)
  # c_rna:        cytosolic total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # c_libsize:    cytosolic total amount of read counts (summed over all transcript species), resolved by time (vector)
  
  # time_series:      time points at which measurements were taken (vector)
  # p_un:             probability that a pre-existing transcript remains unmodified, resolved by time (vector)
  # p_mod:            probability that a newly synhesized transcript is modified, resolved by time (vector)
  
  # new_by_total_resolution:    interpolation parameter; new/total transcript counts ratio interval resolution used for 
  #                             interpolation of the conditional probability 2 (default: 1000)
  
  # ip_points:          optional; number of points used for interpolation (default: 100)
  # search_points:      optional; peak-search parameter; number of points the likelihood functions is evaluated at to 
  #                     find the peak (default: 1000)
  # width:              optional; peak-search parameter; vector; width (as percent of a function's domain's range) of
  #                     a peak-spanning region considered for interpolation for likelihood_1, likelihood_3, likelihood_4
  #                     (default: c(10, 10, 10))

  # *** RETURN ***
  
  # a list containing (in the following order): an interpolated function for likelihood_1 (1-dim), 
  # likelihood_2 (2-dim, input parameters as vector), likelihood_3 (1-dim), likelihood_4 (1-dim),
  # a vector with peak region bordering x-values for likelihood_1, likelihood_3, likelihood_4
  
  # This function performs interpolation on the derived parameters' likelihood functions, given observed count data
  # for one transcript species. For the one-dimensional likelihood functions, it determines the peak region first and
  # interpolates on the peak region only (all values beyond are interpolated with 0).
  
  # *******************************************************************************************************************

  time_points <- length(time_series)          # number of time points
  
  # pre-computing conditional probabilities and cond. prob. 2 lookup table --------------------------------------------
  
  # list of matrices for nuclear n_mod probabilities, given n_new (resolved by time points)                 
  p_mod_pre_n <- p_nmod_nnew_matrix(matrix(n_mod, ncol = 1), matrix(n_total, ncol = 1), 
                                    1, time_points, p_un, p_mod)[[1]]
  # list of matrices for cytosolic c_mod probabilities, given c_new (resolved by time points)                 
  p_mod_pre_c <- p_nmod_nnew_matrix(matrix(c_mod, ncol = 1), matrix(c_total, ncol = 1), 
                                    1, time_points, p_un, p_mod)[[1]]
  
  # list of vectors for c_new probabilities, given c_total (resolved by time points) (needed for prob 2 lookup table)
  p_new_pre_c <- p_nnew_ntotal_matrix(matrix(c_total, ncol = 1), 1, time_points)[[1]]
  # list of lookup tables for conditional probability 2 (resolved by time points)
  new_by_total_ratios <- seq(0, 1, 1/new_by_total_resolution)     # new/total ratios to evaluate cond. prob. 2
  lookup_table_p2 <- sapply(1:time_points, function(t) { list(    # lookup tables for conditional probability 2
    sapply(new_by_total_ratios, function(q) {    
      p_nmod_theta(c_mod[t], c_total[t], q, p_un[t], p_mod[t], p_mod_pre = p_mod_pre_c[t], p_new_pre = p_new_pre_c[t])
    })
  )})

  # conditional probability functions with plugged-in in count data ---------------------------------------------------
  
  prob1 <- function(lambda_n_tau) {
    return(prob_thetaderiv_1(n_total, n_mod, lambda_n_tau[1], time_series, p_un, p_mod, p_mod_pre_n))
  }
  prob2 <- function(params) {
    lambda_n_tau = params[1]
    lambda_c = params[2]
    result = prob_thetaderiv_2("lookup", lambda_n_tau, lambda_c, time_series, lookup_table = lookup_table_p2)
    return(result)
  }
  prob3 <- function(steady_n) {
    return(prob_thetaderiv_3(n_total, n_rna, n_libsize, steady_n[1]))
  }
  prob4 <- function(steady_c) {
    return(prob_thetaderiv_4(c_total, c_rna, c_libsize, steady_c[1]))
  }
  
  # finding peak of functions prob1, prob3 and prob4 (interpolation has to be performed only within the region where 
  # the peak is located, this improves interpolation and reduces running time) ----------------------------------------
  
  lntau_borders_adj <- find_peak(prob1, lntau_borders, search_points, width[1])
  steadyn_borders_adj <- find_peak(prob3, steadyn_borders, search_points, width[2])
  steadyc_borders_adj <- find_peak(prob4, steadyc_borders, search_points, width[3])
  
  # interpolating conditional probability functions -------------------------------------------------------------------
  
  print("approximating prob1")
  prob1_interpol <- approx1D(prob1, lntau_borders_adj[1], lntau_borders_adj[2], ip_points)
  print("approximating prob2")
  prob2_interpol <- approx2D(prob2, c(lntau_borders[1], lc_borders[1]), c(lntau_borders[2], lc_borders[2]), ip_points)
  print("approximating prob3")
  prob3_interpol <- approx1D(prob3, steadyn_borders_adj[1], steadyn_borders_adj[2], ip_points)
  print("approximating prob4")
  prob4_interpol <- approx1D(prob4, steadyc_borders_adj[1], steadyc_borders_adj[2], ip_points)
  
  # returning interpolated functions ----------------------------------------------------------------------------------
  
  return(list(prob1_interpol, prob2_interpol, prob3_interpol, prob4_interpol, 
              lntau_borders_adj, steadyn_borders_adj, steadyc_borders_adj))
}

probability_transformation <- function(lntau_probability, lc_probability, steadyn_probability, steadyc_probability,
                                       lntau_borders, lc_borders, steadyn_borders, steadyc_borders,
                                       steadyn_borders_adj = NULL, steadyc_borders_adj = NULL) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # lntau_probability:    conditional probability function of lambda_n + tau (1-dim)
  # lc_probability:       conditional probability function of lambda_c and lambda_n + tau (2-dim)
  # steadyn_probability:  conditional probability function of nuclear steady-state (1-dim)
  # steadyc_probability:  conditional probability function of cytosolic steady-state (1-dim)
  
  # lntau_borders:    vector; lower and upper border of lambda_n + tau
  # lc_borders:       vector; lower and upper border of lambda_c
  # steadyn_borders:  vector; lower and upper border of steady_n
  # steadyc_borders:  vector; lower and upper border of steady_c
  
  # steadyn_borders_adj:    vector; alternative lower and upper border of steady_n to be used for integration 
  #                         (default: NULL)
  # steadyc_borders_adj:    vector; alternative lower and upper border of steady_c to be used for integration 
  #                         (default: NULL)
  
  # *** NOTE ***
  
  # Except for lc_probability, conditional probability functions must be vectorized (as argument, a vector can be
  # given and for each vector element the function computes the result). lc_likelihood takes a vector containing the
  # lambda_n + tau and lambda_c value (in this order).
  
  # *** RETURN ***
  
  # the raw parameter set density function (input parameters as vector in the following order: 
  # mu, tau, lambda_n, lambda_c)
  
  # This function transforms the single conditional probability functions based on the derived parameter set into the 
  # raw parameter set's density function.
  
  # *******************************************************************************************************************
  
  # conditional probability function integrals to determine normalization factor

  if (!steadyn_borders_adj) {             # if no alternative borders have been specified, plug in regular borders
    steadyn_borders_adj <- steadyn_borders
  }
  if (!steadyc_borders_adj) {
    steadyc_borders_adj <- steadyc_borders
  }

  lntau_lc_combined <- function(args) {   # combining lntau_probability (1D) and lc_probability (2D) for integration
    lntau = args[1]                                 #   since they both share parameter lntau
    lc = args[2]
    value = lntau_probability(lntau) * lc_probability(c(lntau, lc))
    return(value)
  }
  
  lntau_lc_nf <- cuhre(2, 1, lntau_lc_combined,
                       lower = c(lntau_borders[1], lc_borders[1]), upper = c(lntau_borders[2], lc_borders[2]),
                       rel.tol = 1e-3, abs.tol = 1e-12, flags = list(verbose = 0), max.eval = 10**9)[[5]][1]
  print(lntau_lc_nf)                                                              # lntau + lc integral (normalization factor)
  sn_nf <- integrate(steadyn_probability, steadyn_borders_adj[1], steadyn_borders_adj[2], 
                     subdivisions = 1e6, rel.tol = 1e-3)[[1]]   # steadyn integral (normalization factor)
  print(sn_nf)
  sc_nf <- integrate(steadyc_probability, steadyc_borders_adj[1], steadyc_borders_adj[2], 
                     subdivisions = 1e6, rel.tol = 1e-3)[[1]]   # steadyc integral (normalization factor)
  print(sc_nf)
  nf <- lntau_lc_nf * sn_nf * sc_nf     # combined normalization factor

  # retrieving raw parameter borders from derived parameter borders
  
  raw_borders <- transform_params("raw", lntau_borders, steadyn_borders, steadyc_borders, lc_borders)

  # setting up raw parameter set's density function (normalized & transformed according to 'change of variable')
  
  raw_density <- function(args) {
    
    mu = args[1]
    tau = args[2]
    lambda_n = args[3]
    lambda_c = args[4]

    lntau = lambda_n + tau
    steadyn <- mu / lntau
    steadyc <- steadyn * tau / lambda_c

    density_value =
      lntau_probability(lntau) * lc_probability(c(lntau, lambda_c)) *
      steadyn_probability(steadyn) * steadyc_probability(steadyc) / nf *
      1 / (steadyn / (lambda_c * lntau))
      
    return(density_value)
  }
  
  # returning
  
  return(raw_density)
}

parameter_probability <- function(n_rna, c_rna, n_libsize, c_libsize, n_total, c_total, n_mod, c_mod, p_un, p_mod,
                                  mu_borders, tau_borders, ln_borders, lc_borders, time_series,
                                  new_by_total_resolution = 1000, ip_points = 100, search_points = 1000, 
                                  width = c(10, 10, 10)) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_rna:        nuclear total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # c_rna:        cytosolic total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # n_libsize:    nuclear total amount of read counts (summed over all transcript species), resolved by time (vector)
  # c_libsize:    cytosolic total amount of read counts (summed over all transcript species), resolved by time (vector)
  # n_total:      nuclear total counts observed for the transcript species, resolved by time (vector)
  # c_total:      cytosolic total counts observed for the transcript species, resolved by time (vector)
  # n_mod:        vector; nuclear modified counts observed for the transcrip species, resolved by time
  # c_mod:        vector; cytosolic modified counts observed for the transcript species, resolved by time
  # p_un:         vector; probability that a pre-existing transcript remains unmodified, resolved by time
  # p_mod:        vector; probability that a newly synhesized transcript is modified, resolved by time
  
  # mu_borders:   vector with (minimum, maximum) lambda_n + tau value 
  # tau_borders:  vector with (minimum, maximum) lambda_c value 
  # ln_borders:   vector with (minimum, maximum) steady_n value 
  # lc_borders:   vector with (minimum, maximum) steady_c value
  # time_series:  time points at which measurements were taken (vector)
  
  # new_by_total_resolution:    optional; new/total transcript counts ratio interval resolution used for computation 
  #                             of the conditional probability 2 lookup table (given total counts, resolved by 
  #                             RNA dynamics parameters i.e. new/total ratios) (default: 1000)
  # ip_points:                  optional; number of points used for interpolations (default: 100)
  # search_points:              optional; peak-search parameter; number of points the likelihood functions is evaluated 
  #                             at to find the peak (default: 1000)
  # width:                      optional; peak-search parameter; vector; width (as percent of a function's domain's
  #                             range) of a peak-spanning region considered for interpolation for likelihood_1, 
  #                             likelihood_3, likelihood_4 (default: c(10, 10, 10))
  
  # *** RETURN ***
  
  # a list containing (in the following order): raw parameters' joint probability density, lambda_n + tau conditional
  # probability, lambda_c and lambda_n + tau joint conditional probability, steady_n conditional probability,
  # steady_c conditional probability
  
  # This function computes RNA dynamics parameters' probability densities (using a Bayesian approach), given observed 
  # data, for one transcript species.
  
  # *******************************************************************************************************************
  
  n_time_points <- length(time_series)
  
  # retrieving derived parameter space from raw parameter space
  
  derived_borders <- transform_params("derived", mu_borders, tau_borders, ln_borders, lc_borders)
  lntau_borders <- derived_borders[1:2]
  sn_borders <- derived_borders[3:4]
  sc_borders <- derived_borders[5:6]
  
  # computing derived parameters' conditional probability functions
  
  derived_probs <- probability_interpolation(lntau_borders, lc_borders, sn_borders, sc_borders,
                                             n_total, n_mod, n_rna, n_libsize, c_total, c_mod, c_rna, c_libsize,
                                             time_series, p_un, p_mod, new_by_total_resolution, ip_points, 
                                             search_points, width)
  
  lntau_probability <- Vectorize(derived_probs[[1]])
  lc_probability <- derived_probs[[2]]                # this is a 2D function, no Vectorize needed
  sn_probability <- Vectorize(derived_probs[[3]])
  sc_probability <- Vectorize(derived_probs[[4]])
  
  sn_borders_peak <- derived_probs[[6]]     # getting functions' peak region borders (over which interpolation was
  sc_borders_peak <- derived_probs[[7]]     #   performed and integration is to be performed, too)
  
  # computing joint raw parameters' density functions
  
  joint_density <- probability_transformation(lntau_probability, lc_probability, sn_probability, sc_probability, 
                                              lntau_borders, lc_borders, sn_borders, sc_borders,
                                              sn_borders_peak, sc_borders_peak)
  
  # returning
  
  return(list(joint_density, lntau_probability, lc_probability, sn_probability, sc_probability))
}

halflife_grid <- function(lntau_borders, lc_borders, time_series, p_un_nu, p_un_cy, p_mod_nu, p_mod_cy, 
                          grid_resolution = 100) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # lntau_borders:  vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:     vector with (minimum, maximum) lambda_c value
  # time_series:    vector of time points at which measurements were taken
  # p_un_nu:        vector; probability that a nuclear, pre-existing transcript remains unmodified, resolved by time
  # p_un_cy:        vector; probability that a cytosolic, pre-existing transcript remains unmodified, resolved by time
  # p_mod_nu:       vector; probability that a nuclear, newly synhesized transcript is modified, resolved by time
  # p_mod_cy:       vector; probability that a cytosolic, newly synhesized transcript is modified, resolved by time
  
  # grid_resolution:    optional; grid search resolution to use per parameter (default: 100)
  
  # *** RETURN ***
  
  # a list of matrices (rows: lambda_n + tau, columns: lambda_c) containing the modified/total ratios in percent as 
  # predicted by the RNA dynamics model; the first half of the list contains nuclear ratios (ordered by time), the 
  # second half cytosolic ratios (ordered by time)
  
  # *******************************************************************************************************************
  
  # creating parameter grid (using a logarithmic scale)
  
  lntau_grid <- sapply(1:grid_resolution, function(i) {                                             # lntau grid
    exp(seq(log(lntau_borders[1]), log(lntau_borders[2]), length.out = grid_resolution))
  })
  lc_grid <- sapply(
    exp(seq(log(lc_borders[1]), log(lc_borders[2]), length.out = grid_resolution)), function(i) {   # lc grid
      rep(i, grid_resolution)
    })
  gridsize <- grid_resolution ** 2                                                                  # grid size
  
  # simulating RNA dynamics through grid
  # matrix; column names: grid indices, first half of rows: nucleus mod/total ratios (ordered by time), second half of 
  # rows: cytosol mod/total ratios (ordered by time)
  
  grid_rna_sim <- sapply(1:gridsize, function(i) {   # iterating through grid
    
    lntau_val <- lntau_grid[i]    # getting grid lambda_n + tau value
    lntau_split <- lntau_val/2    # splitting lambda_n + tau into two halfs (one for tau and one for lambda_n)
    lc_val <- lc_grid[i]          # getting grid lambda_c value
    
    n_ratio_new_total <- ratio_n_sol(time_series, lntau_split, lntau_split)           # nuclear new/total ratio
    c_ratio_new_total <- ratio_c_sol(time_series, lntau_split, lntau_split, lc_val)   # cytosol new/total ratio
    
    n_ratio_mod_total <- n_ratio_new_total * p_mod_nu + 
      (1 - n_ratio_new_total) * (1 - p_un_nu)                                         # nuclear mod/total ratio
    c_ratio_mod_total <- c_ratio_new_total * p_mod_cy + 
      (1 - c_ratio_new_total) * (1 - p_un_cy)                                         # cytosol mod/total ratio
    
    return(c(n_ratio_mod_total, c_ratio_mod_total) * 100)     # returning mod/total ratios in [%]
    
  })
  
  # restructuring data to a list of matrices (list indices: nu and cy time series, rows: lambda_n + tau values,
  # columns: lambda_c values)
  
  list_rna_sim <- sapply(1:nrow(grid_rna_sim), function(i) { list(    # iterating through rows (nu and cy time_series)
    sapply(0:(grid_resolution - 1), function(j) { 
      grid_rna_sim[i, ][(j * grid_resolution + 1):(j * grid_resolution + grid_resolution)]
    })
  )})
  
  # returning
  
  return(list_rna_sim)
}

halflife_sumofsquares_labtot <- function(halflife_pars, nu_measured, cy_measured, time_series, 
                                         p_un_nu, p_un_cy, p_mod_nu, p_mod_cy) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # halflife_pars:  vector with (lambda_n + tau, lambda_c) parameter suggestions 
  # nu_measured:    vector of nuclear modified/total transcript count ratios measured, resolved by time  
  # cy_measured:    vector of cytosolic modified/total transcript count ratios measured, resolved by time  
  # time_series:    vector of time points at which measurements were taken
  # p_un_nu:        vector; probability that a nuclear, pre-existing transcript remains unmodified, resolved by time
  # p_un_cy:        vector; probability that a cytosolic, pre-existing transcript remains unmodified, resolved by time
  # p_mod_nu:       vector; probability that a nuclear, newly synhesized transcript is modified, resolved by time
  # p_mod_cy:       vector; probability that a cytosolic, newly synhesized transcript is modified, resolved by time
  
  # *** RETURN ***
  
  # sum of squares of the estimated vs. measured modified/total count ratios (for the given parameter suggestion);
  # this function infers labeled/total counts observed for a given new/total prediction of a parameter choice, using
  # the error estimates p_mod and p_un
  
  # *******************************************************************************************************************

  lntau <- halflife_pars[1]     # lambda_n + tau suggestion
  lc_val <- halflife_pars[2]    # lambda_c suggestion
  if (lntau == lc_val) {        # lntau == lc yields zero division problem in function 'ratio_c_sol'
    lc_val <- lc_val + 0.00001
  }   
  lntau_split <- halflife_pars[1] / 2   # lambda_n and tau suggestion (splitting into two halfs for the two variables)
  
  # computing estimated values
  n_ratio_new_total <- ratio_n_sol(time_series, lntau_split, lntau_split)           # nuclear new/total ratio
  c_ratio_new_total <- ratio_c_sol(time_series, lntau_split, lntau_split, lc_val)   # cytosol new/total ratio
  n_ratio_mod_total <- (n_ratio_new_total * p_mod_nu + 
                          (1 - n_ratio_new_total) * (1 - p_un_nu)) * 100            # nuclear mod/total ratio in [%]
  c_ratio_mod_total <- (c_ratio_new_total * p_mod_cy + 
                          (1 - c_ratio_new_total) * (1 - p_un_cy)) * 100            # cytosol mod/total ratio in [%]
  
  # computing sum of squares between estimation and measurement
  sum_of_squares <- sum((nu_measured * 100 - n_ratio_mod_total) ** 2, 
                        (cy_measured * 100 - c_ratio_mod_total) ** 2)   # computing sum of squares (measured vs est.)
  return(sum_of_squares)                                                # returning
}

halflife_sumofsquares_newtot <- function(halflife_pars, time_series, new_tot_nu, new_tot_cy, 
                                         total_nu, total_cy, var_stab = TRUE) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # halflife_pars:  vector with (lambda_n + tau, lambda_c) parameter suggestions 
  # time_series:    vector of time points at which measurements were taken
  # new_tot_nu:     vector of nuclear new/total transcript count ratios suggestions, resolved by time  
  # new_tot_cy:     vector of cytosolic new/total transcript count ratios suggestions, resolved by time 
  # total_nu:       vector of nuclear total counts, resolved by time; needed for variance stabililzig transformation,
  #                 an arbitrary value can be specified in case that no transformation is performed
  # total_cy:       vector of cytosolic total counts, resolved by time; needed for variance stabililzig transformation,
  #                 an arbitrary value can be specified in case that no transformation is performed
  # var_stab:       set to TRUE to perform variance-stabilizing transformation of the data points (default: TRUE)
  
  # *** RETURN ***
  
  # sum of squares of the input suggestion vs. ODE predicted modified/total count ratios (for the given parameter 
  # suggestion); this function is based on a pre-defined suggestion of new/total ratios (which can, for example, be
  # the estimates given by the EM algorithm for fixing the conversion efficiency at each time point)
  
  # *******************************************************************************************************************
  
  lntau <- halflife_pars[1]     # lambda_n + tau suggestion
  lc_val <- halflife_pars[2]    # lambda_c suggestion
  if (lntau == lc_val) {        # lntau == lc yields zero division problem in function 'ratio_c_sol'
    indi <- 1
    lc_val <- lc_val + 0.00001
  }   
  lntau_split <- halflife_pars[1] / 2   # lambda_n and tau suggestion (splitting into two halfs for the two variables)
  
  # computing estimated values
  n_ratio_new_total <- 
    ratio_n_sol(time_series, lntau_split, lntau_split)            # nuclear new/total ratio (in %)
  c_ratio_new_total <- 
    ratio_c_sol(time_series, lntau_split, lntau_split, lc_val)    # cytosol new/total ratio (in %)
  n_ratio_new_total[which(n_ratio_new_total < 0)] <- 0            # accounting for computational accuracy problems
  c_ratio_new_total[which(c_ratio_new_total < 0)] <- 0            # accounting for computational accuracy problems
  
  # computing sum of squares between suggestion and ODE prediction
  if (!var_stab) {    # without variance-stabilizing transformation           
    sum_of_squares <- sum((new_tot_nu * 100 - n_ratio_new_total * 100) ** 2,  # x100 to avoid very small values
                          (new_tot_cy * 100 - c_ratio_new_total * 100) ** 2) 
  }
  else {              # with variance-stabilizing transformation
    sum_of_squares <- sum((asin(sqrt(new_tot_nu)) - asin(sqrt(n_ratio_new_total))) ** 2 * 4 * total_nu, 
                          (asin(sqrt(new_tot_cy)) - asin(sqrt(c_ratio_new_total))) ** 2 * 4 * total_cy) 
  }

  if (sum_of_squares == -Inf) {sum_of_squares <- Inf}       # adaptation for DEoptim
  else if (is.na(sum_of_squares)) {sum_of_squares <- Inf}   # adaptation for DEoptim
  return(sum_of_squares)    # returning
}

### Wrappers ##########################################################################################################

wrapper_mu_tau_ln_relative <- function(summary_table_file, estimation_table_file, time_series) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summary_table_file:     file storing summary table (rows: genomic regions, col1: region name, col2: measurement
  #                         description, col3: library size, col4: total transcript counts, col5: labeled transcript  
  #                         counts, col6: average potential conversion positions, col7: conversion efficiency, 
  #                         col8: newly synthesized ratio)
  # estimation_table_file:  file storing lambda_n + tau and lambda_c estimatiors (rows: genomic regions, 
  #                         col1: lambda_n + tau, col2: lambda_c, col3: sum of squares)
  # time_series:            vector of time points at which measurements were taken
  
  # *** NOTE ***
  
  # The summary table file must be sorted by gene name, compartment, and time measurement (in that order). 
  
  # *** RETURN ***
  
  # a matrix storing single genes' relative estimates (rows: genes, col1: mu relative estimate, col2: tau relative 
  # estimate, col3: lambda_n relative estimate, col4: mu estimate's standard deviation, col5: tau estimate's standard 
  # deviation, col6: lambda_n estimate's standard deviation,)
  
  # This function estimates relative synthesis (mu), transport (tau), and nuclear degradation (lambda_n) rates, 
  # based on the estimators of the halflife parameters lambda_n + tau and lambda_c. Estimates of mu and tau are
  # normalized w.r.t. median amongst all genes (normalization is not performed for lambda_n since noise can lead to
  # negative values here).
  
  # *******************************************************************************************************************
  
  # initializing ......................................................................................................
  
  s_table <- load_summarytable(summary_table_file)                  # loading in summary table
  est_table <- read.table(estimation_table_file,                    # loading in estimation table
                          header = TRUE, sep = "\t", quote = "")    
  gene_names <- rownames(est_table)                                 # all gene names
  n_genes <- length(gene_names)                                     # number of genes
  n_timepoints <- length(time_series)                               # number of time points (time series measurements)
  
  # iterating through genes, computing relative mu and tau estimates
  # columns: genes, rows: mu, tau and lambda_n estimates (resolved by time) ...........................................
  
  raw_relative_estimates <- sapply(0:(n_genes-1), function(i) {   # iterating through genes
    
    # getting gene's measurement information
    
    gene_name <- gene_names[i + 1]                    # getting gene's name
    start_idx_nu <- i * 2 * n_timepoints + 1          # row index within summary table where gene's nu entries start
    end_idx_nu <- start_idx_nu + n_timepoints - 1     # row index within summary table where gene's nu entries end
    start_idx_cy <- end_idx_nu + 1                    # row index within summary table where gene's cy entries start
    end_idx_cy <- start_idx_cy + n_timepoints - 1     # row index within summary table where gene's cy entries end
    
    submatrix_nu <- s_table[start_idx_nu:end_idx_nu, ]      # getting gene's submatrix with nuclear entries
    submatrix_cy <- s_table[start_idx_cy:end_idx_cy, ]      # getting gene's submatrix with cytosolic entries
    
    # getting gene's estimation information
    
    lnt_lc_estimators <- est_table[gene_name, ]
    
    # computing relative estimates (from relative expression and halflife estimators)
    
    nu_expression_levels <- submatrix_nu[, "tot"] / submatrix_nu[, "lib"] * 10**5     # expression levels (times
    cy_expression_levels <- submatrix_cy[, "tot"] / submatrix_cy[, "lib"] * 10**5     # 10^5 to avoid small values)
    mu_rel_estimates <- nu_expression_levels * lnt_lc_estimators[1, 1]                # mu relative estimates
    tau_rel_estimates <- cy_expression_levels * lnt_lc_estimators[1, 2]               # tau relative estimates
    ln_rel_estimates <- mu_rel_estimates / nu_expression_levels - tau_rel_estimates   # ln relative estimates
    
    return(c(mu_rel_estimates, tau_rel_estimates, ln_rel_estimates))
  })
  
  # normalizing relative estimates w.r.t. median, 
  # computing average and stdev for each gene .........................................................................
  
  # computing genes' relative estimates' averages (columns: genes, rows: mu and tau average relative estimate)
  raw_estimate_means <- sapply(1:ncol(raw_relative_estimates), function(i) {
    mu_average <- sum(raw_relative_estimates[1:n_timepoints, i]) / n_timepoints
    tau_average <- sum(raw_relative_estimates[(n_timepoints + 1):(n_timepoints * 2), i]) / n_timepoints
    ln_average <- sum(raw_relative_estimates[(n_timepoints * 2 + 1):(n_timepoints * 3), i]) / n_timepoints
    return(c(mu_average, tau_average, ln_average))
  })
  # computing genes' relative estimates' stdevs (columns: genes, rows: mu and tau average relative estimate)
  raw_estimate_stdevs <- sapply(1:ncol(raw_relative_estimates), function(i) {
    mu_stdev <- sd(raw_relative_estimates[1:n_timepoints, i])
    tau_stdev <- sd(raw_relative_estimates[(n_timepoints + 1):(n_timepoints * 2), i])
    ln_stdev <- sd(raw_relative_estimates[(n_timepoints * 2 + 1):(n_timepoints * 3), i])
    return(c(mu_stdev, tau_stdev, ln_stdev))
  })
  
  # computing median of genes' average relative estimates
  mu_median <- median(raw_estimate_means[1, ])    # mu median of genes' average relative estimates
  tau_median <- median(raw_estimate_means[2, ])   # tau median of genes' average relative estimates
  
  # normalizing with median
  mu_rel_normalized <- raw_estimate_means[1, ] / mu_median
  mu_stdev_normalized <- raw_estimate_stdevs[1, ] * sqrt(n_timepoints/ mu_median)
  tau_rel_normalized <- raw_estimate_means[2, ] / tau_median
  tau_stdev_normalized <- raw_estimate_stdevs[2, ] * sqrt(n_timepoints/ tau_median)
  
  # returning matrix with normalized relative estimates and standard deviations .......................................
  
  rel_est_mat <- matrix(c(mu_rel_normalized, tau_rel_normalized, raw_estimate_means[3, ],   # generating matrix
                          mu_stdev_normalized, tau_stdev_normalized, raw_estimate_stdevs[3, ]), ncol = 6)   
  rownames(rel_est_mat) <- gene_names           # rows: genes, columns: mu, tau, lambda_n rel. est. and rel. est. stdev
  colnames(rel_est_mat) <- c("mu_rel", "tau_rel", "ln_rel", "mu_stdev", "tau_stdev", "ln_stdev")    
  return(rel_est_mat)
}

wrapper_tau_ln_absolute <- function(summary_table_file, estimation_table_file, time_series, cy_nu_ratio) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summary_table_file:     file storing summary table (rows: genomic regions, col1: region name, col2: measurement
  #                         description, col3: library size, col4: total transcript counts, col5: labeled transcript  
  #                         counts, col6: average potential conversion positions, col7: conversion efficiency, 
  #                         col8: newly synthesized ratio)
  # estimation_table_file:  file storing lambda_n + tau and lambda_c estimatiors (rows: genomic regions, 
  #                         col1: lambda_n + tau, col2: lambda_c, col3: sum of squares)
  # time_series:            vector of time points at which measurements were taken
  # cy_nu_ratio:            ratio between cytosolic and nuclear cellular RNA (holding for the bulk of RNA species 
  #                         evaluated)
  
  # *** NOTE ***
  
  # The summary table file must be sorted by gene name, compartment (nu, cy), and time measurement (in that order). 
  
  # *** RETURN ***
  
  # a matrix storing single genes' absolute estimates (rows: genes, col1: tau relative estimate, col2: lambda_n 
  # estimate, col3: tau estimate's standard deviation, col4: lambda_n estimate's standard deviation)
  
  # This function estimates absolute transport (tau) and nuclear degradation (lambda_n) rates, based on the 
  # estimators of the halflife parameters lambda_n + tau and lambda_c and the cellular ratio between cytosolic and
  # nuclear total RNA.
  
  # *******************************************************************************************************************
  
  # initializing ......................................................................................................
  
  s_table <- load_summarytable(summary_table_file)                  # loading in summary table
  est_table <- read.table(estimation_table_file,                    # loading in estimation table
                          header = TRUE, sep = "\t", quote = "")    
  gene_names <- rownames(est_table)                                 # all gene names
  n_genes <- length(gene_names)                                     # number of genes
  n_timepoints <- length(time_series)                               # number of time points (time series measurements)
  
  # computing ratio between library size scaling factors ..............................................................
  
  libsizes_nu <- stable[1:n_timepoints, "lib"]                      # nuclear measurements' library sizes
  libsizes_cy <- stable[(n_timepoints+1):(n_timepoints*2), "lib"]   # cytosolic measurements' library sizes
  
  lib_scaling_ratio <- libsizes_cy / (libsizes_nu * cy_nu_ratio)    # library size scaling factors' ratio (cy by nu)
  
  # iterating through genes, computing tau, lambda_n and relative mu estimates
  # columns: genes, rows: mu, tau and lambda_n estimates (resolved by time) ...........................................
  
  raw_estimates <- sapply(0:(n_genes-1), function(i) {    # iterating through genes
    
    # getting gene's measurement information
    
    gene_name <- gene_names[i + 1]                    # getting gene's name
    start_idx_nu <- i * 2 * n_timepoints + 1          # row index within summary table where gene's nu entries start
    end_idx_nu <- start_idx_nu + n_timepoints - 1     # row index within summary table where gene's nu entries end
    start_idx_cy <- end_idx_nu + 1                    # row index within summary table where gene's cy entries start
    end_idx_cy <- start_idx_cy + n_timepoints - 1     # row index within summary table where gene's cy entries end
    
    submatrix_nu <- s_table[start_idx_nu:end_idx_nu, ]      # getting gene's submatrix with nuclear entries
    submatrix_cy <- s_table[start_idx_cy:end_idx_cy, ]      # getting gene's submatrix with cytosolic entries
    
    # getting gene's estimation information
    
    lnt_lc_estimators <- est_table[gene_name, ]
    
    # computing parameter estimates
    
    tau_estimates <- 
      (submatrix_cy[, "tot"] * lnt_lc_estimators[1, 2]) / 
      (submatrix_nu[, "tot"] * lib_scaling_ratio)             # tau estimates
    ln_estimates <- lnt_lc_estimators[1, 1] - tau_estimates   # lambda_n estimates
    
    all <- c(tau_estimates, ln_estimates)                     # all estimates in one vector
    return(all)
  })
  
  # normalizing mu relative estimates w.r.t. median, 
  # computing average and stdev for each gene .........................................................................
  
  # computing genes' estimates' averages (columns: genes, rows: mu, tau, ln average estimate)
  raw_estimate_means <- sapply(1:ncol(raw_estimates), function(i) {
    tau_average <- sum(raw_estimates[1:n_timepoints, i]) / n_timepoints
    ln_average <- sum(raw_estimates[(n_timepoints + 1):(n_timepoints * 2), i]) / n_timepoints
    return(c(tau_average, ln_average))
  })
  # computing genes' estimates' stdevs (columns: genes, rows: mu, tau, ln average estimate)
  raw_estimate_stdevs <- sapply(1:ncol(raw_estimates), function(i) {
    tau_stdev <- sd(raw_estimates[1:n_timepoints, i])
    ln_stdev <- sd(raw_estimates[(n_timepoints + 1):(n_timepoints * 2), i])
    return(c(tau_stdev, ln_stdev))
  })
  
  # returning matrix with average estimates and standard deviations .......................................
  
  abs_est_mat <- 
    matrix(c(raw_estimate_means[1, ], raw_estimate_means[2, ],   
             raw_estimate_stdevs[1, ], raw_estimate_stdevs[2, ]), 
           ncol = 4)                          # generating matrix
  rownames(abs_est_mat) <- gene_names         # rows: genes, columns: mu, tau, lambda_n rel. est. and rel. est. stdev
  colnames(abs_est_mat) <- c("tau_abs", "ln_abs", "tau_stdev", "ln_stdev")    
  return(abs_est_mat)
}

halflife_estimation_maxlike <- function(n_total, c_total, n_mod, c_mod, p_un_nu, p_un_cy, p_mod_nu, p_mod_cy,
                                        lntau_borders, lc_borders, time_series, 
                                        mode = "binom", betabin_od = 2, use_lookup = TRUE, ip_points = 1000) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_total:      vector; nuclear total counts observed for the transcript species, resolved by time
  # c_total:      vector; cytosolic total counts observed for the transcript species, resolved by time
  # n_mod:        vector; nuclear modified counts observed for the transcrip species, resolved by time
  # c_mod:        vector; cytosolic modified counts observed for the transcript species, resolved by time
  # p_un_nu:      vector; probability that a nuclear pre-existing transcript remains unmodified, resolved by time
  # p_un_cy:      vector; probability that a cytosolic pre-existing transcript remains unmodified, resolved by time
  # p_mod_nu:     vector; probability that a nuclear newly synhesized transcript is modified, resolved by time
  # p_mod_cy:     vector; probability that a cytosolic newly synhesized transcript is modified, resolved by time
  
  # lntau_borders:  vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:     vector with (minimum, maximum) lambda_c value
  # time_series:    vector of time points at which measurements were taken
  
  # mode:             set to "binom" for an underlying binomial distribution, set to "betabinom" for an underlying 
  #                   beta-binomial distribution (default: "binom")
  # betabin_od:       needed for 'mode' == 'betabinom'; beta-binomial distribution overdispersion parameter
  #                   (alpha + beta) (default: 2)
  # use_lookup:       optional; set to TRUE to use a lookup table for conditional probability 2, interpolating over
  #                   a grid of new/total ratios (given the total counts, resolved by time i.e. one interpolated
  #                   function for each time point) (default: TRUE)
  # ip_points:        optional; number of points used for interpolations (default: 1000)

  # *** NOTE ***
  
  # In case a lookup table is computed for conditional probability 2, 2D interpolation on this conditional probability
  # is omitted.
  
  # *** RETURN ***
  
  # a list containing (in the following order): 
  #     estimation optimization results
  #     lambda_n + tau raw conditional probability
  #     lambda_c and lambda_n + tau raw joint conditional probability
  #     lambda_n + tau interpolated conditional probability 1
  #     lambda_c and lambda_n + tau interpolated conditional probability 2
  #     joint conditional probability 1 and 2 function
  
  # This function performs RNA dynamics parameter estimation on the half-life determining parameters lambda_n + tau
  # amd lambda_c (using a Bayesian approach), given observed data, for one transcript species.
  
  # *******************************************************************************************************************

  time_points <- length(time_series)                                                    # number of time points
  
  # pre-computing conditional probabilities ---------------------------------------------------------------------------
  
  # list of matrices for nuclear n_mod probabilities, given n_new (resolved by time points)                 
  p_mod_pre_n <- p_nmod_nnew_matrix(matrix(n_mod, ncol = 1), matrix(n_total, ncol = 1), 
                                    1, time_points, p_un_nu, p_mod_nu, mode = mode, betabin_od = betabin_od)[[1]]
  # list of matrices for cytosolic c_mod probabilities, given c_new (resolved by time points)                 
  p_mod_pre_c <- p_nmod_nnew_matrix(matrix(c_mod, ncol = 1), matrix(c_total, ncol = 1), 
                                    1, time_points, p_un_cy, p_mod_cy, mode = mode, betabin_od = betabin_od)[[1]]

  # pre-computing lookup table for cond. prob. 2 ----------------------------------------------------------------------
  
  lookup_table <- NULL
  
  if (use_lookup) {
    
    # list of lookup tables for conditional probability 2 (list (time points) of functions (new/total ratios)) 
    lookup_table <- sapply(1:time_points, function(t) { 
        lookup <- function(new_total_ratio) {
          p <- p_nmod_theta(c_mod[t], c_total[t], new_total_ratio, p_un_cy[t], p_mod_cy[t], 
                            p_mod_pre = p_mod_pre_c[[t]], mode = mode, betabin_od = betabin_od)
          if ((log10(p) == -Inf) | is.na(p)) {return(.Machine$integer.max * -1)}
          else {return(log10(p))}
        }
        lookup_interpolated <- approx1D(lookup, 0, 1, ip_points = ip_points)
        return(lookup_interpolated)
    })
  }
  
  # computing parameter likelihood / density distributions ------------------------------------------------------------
  
  # conditional probability functions with plugged-in in count data
  
  prob1 <- function(lambda_n_tau) {
    return(prob_thetaderiv_1(n_total, n_mod, lambda_n_tau[1], time_series, p_un_nu, p_mod_nu, p_mod_pre_n,
                             mode = mode, betabin_od = betabin_od, log_scale = TRUE))
  }
  prob2 <- function(params) {   # not using cond. prob. 2 lookup table
    lambda_n_tau = params[1]
    lambda_c = params[2]
    result = prob_thetaderiv_2("denovo", lambda_n_tau, lambda_c, time_series, lookup_table = lookup_table,
                               mode = mode, betabin_od = betabin_od, log_scale = TRUE)
    return(result)
  }
  if (use_lookup) {             # using cond. prob. lookup table
    prob2 <- function(params) {
      lambda_n_tau = params[1]
      lambda_c = params[2]
      result = prob_thetaderiv_2("lookup", lambda_n_tau, lambda_c, time_series, lookup_table = lookup_table,
                                 log_scale = TRUE)
      return(result)
    }
  }
  
  # interpolating conditional probability functions (cond. prob. 2 is only interpolated if there is no lookup table)
  
  prob1_interpol <- approx1D(prob1, lntau_borders[1], lntau_borders[2], ip_points)
  prob2_interpol <- prob2
  if (!use_lookup) {
    prob2_interpol <- approx2D(prob2, c(lntau_borders[1], lc_borders[1]), c(lntau_borders[2], lc_borders[2]), ip_points)
  }
  
  lntau_probability <- Vectorize(prob1_interpol)        
  lc_probability <- prob2_interpol
  
  # generating joint probability distribution -------------------------------------------------------------------------
  
  lntau_lc_combined <- function(args) {   
    lntau = args[1]                                
    lc = args[2]
    if (lntau == lc) {lntau <- lntau + 10**(-5)}    # lntau == lc yields zero division problem in function 'ratio_c_sol'
    value = lntau_probability(lntau) + lc_probability(c(lntau, lc))
    return(value)
  }
  
  # estimating parameter values from joint probability density (maximizing probability density) -----------------------
  
  joint_density_inv <- function(args) {     # since global optimization only supports minimization, joint density
    lntau_lc_combined(args) * (-1)          #   is inverted s.t. the minimum of inverted function corresponds to the
  }                                         #   maximum of the original function
  optima <- DEoptim(joint_density_inv,
                    lower = c(lntau_borders[1], lc_borders[1]),
                    upper = c(lntau_borders[2], lc_borders[2]),
                    control = list(itermax = 10**4, trace=FALSE))   # finding parameter values
  
  # returning
  
  return(list(optima, prob1, prob2, lntau_probability, lc_probability, lntau_lc_combined))
}

parameter_estimation_maxlike <- function(n_rna, c_rna, n_libsize, c_libsize, n_total, c_total, n_mod, c_mod, p_un, p_mod,
                                         mu_borders, tau_borders, ln_borders, lc_borders, time_series,
                                         new_by_total_resolution = 1000, ip_points = 100, search_points = 1000, 
                                         width = c(10, 10, 10)) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_rna:        nuclear total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # c_rna:        cytosolic total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # n_libsize:    vector; nuclear total amount of read counts (summed over all transcript species), resolved by time
  # c_libsize:    vector; cytosolic total amount of read counts (summed over all transcript species), resolved by time
  # n_total:      vector; nuclear total counts observed for the transcript species, resolved by time
  # c_total:      vector; cytosolic total counts observed for the transcript species, resolved by time
  # n_mod:        vector; nuclear modified counts observed for the transcrip species, resolved by time
  # c_mod:        vector; cytosolic modified counts observed for the transcript species, resolved by time
  # p_un:         vector; probability that a pre-existing transcript remains unmodified, resolved by time
  # p_mod:        vector; probability that a newly synhesized transcript is modified, resolved by time
  
  # mu_borders:   vector with (minimum, maximum) lambda_n + tau value 
  # tau_borders:  vector with (minimum, maximum) lambda_c value 
  # ln_borders:   vector with (minimum, maximum) steady_n value 
  # lc_borders:   vector with (minimum, maximum) steady_c value
  # time_series:  vector of time points at which measurements were taken
  
  # new_by_total_resolution:    optional; new/total transcript counts ratio interval resolution used for computation 
  #                             of the conditional probability 2 lookup table (given total counts, resolved by 
  #                             RNA dynamics parameters i.e. new/total ratios) (default: 1000)
  # ip_points:                  optional; number of points used for interpolations (default: 100)
  # search_points:              optional; peak-search parameter; number of points the likelihood functions is evaluated
  #                             at to find the peak (default: 1000)
  # width:                      optional; peak-search parameter; vector; width (as percent of a function's domain's
  #                             range) of a peak-spanning region considered for interpolation for likelihood_1,
  #                             likelihood_3, likelihood_4 (default: c(10, 10, 10))

  # *** RETURN ***
  
  # a list containing (in the following order):
  #   a list containing (in the following order): raw parameters' probability density, lambda_n + tau conditional
  #     probability, lambda_c and lambda_n + tau joint conditional probability, steady_n conditional probability,
  #     steady_c conditional probability
  #   a list of estimated raw parameter settings (mu, tau, ln, lc)
  
  # This function performs RNA dynamics parameter estimation (using a Bayesian approach), given observed data,
  # for one transcript species.
  
  # *******************************************************************************************************************
  
  # computing parameter likelihood / density distributions
  
  print("Computing Densities 1,2,3 and 4")
  
  densities <- parameter_probability(n_rna, c_rna, n_libsize, c_libsize, n_total, c_total, n_mod, c_mod, p_un, p_mod, 
                                     mu_borders, tau_borders, ln_borders, lc_borders, time_series, 
                                     new_by_total_resolution = new_by_total_resolution, ip_points = ip_points, 
                                     search_points = search_points, width = width)
  
  # re-estimating parameter values from joint probability density (maximizing probability density)
  
  print("Determining Optima")
  
  joint_density <- densities[[1]]      # getting joint density
  joint_density_inv <- function(args) {     # since global optimization only supports minimization, joint density
    joint_density(args) * (-1)              #   is inverted s.t. the minimum of inverted function corresponds to the
  }                                         #   maximum of the original function
  optima <- DEoptim(joint_density_inv, lower = c(mu_borders[1], tau_borders[1], ln_borders[1], lc_borders[1]),
                    upper = c(mu_borders[2], tau_borders[2], ln_borders[2], lc_borders[2]),
                    control = list(itermax = 10**5, trace = FALSE))$optim$bestmem   # finding parameter values 

  # returning
  
  return(list(densities, optima))
}

model_simulation_maxlike <- function(n_rna, c_rna, n_libsize, c_libsize, p_un, p_mod, mu, tau, lambda_n, lambda_c, 
                                     time_series, mu_borders, tau_borders, ln_borders, lc_borders, 
                                     new_by_total_resolution = 1000, ip_points = 100, search_points = 1000, 
                                     width = c(10, 10, 10), repetitions = 1) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_rna:        nuclear total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # c_rna:        cytosolic total amount of RNA transcripts [molecules / cell] (summed over all transcript species)
  # n_libsize:    nuclear total amount of read counts (summed over all transcript species), resolved by time (vector)
  # c_libsize:    cytosolic total amount of read counts (summed over all transcript species), resolved by time (vector)
  # p_un:         probability that a pre-existing transcript remains unmodified
  # p_mod:        probability that a newly synhesized transcript is modified
  
  # mu:           RNA synthesis rate of the transcript
  # tau:          RNA transport rate of the transcript
  # lambda_n:     nuclear degradation rate of the transcript
  # lambda_c:     cytosolic degradation rate of the transcript
  # time_series:  time points at which measurements were taken (vector)
  
  # mu_borders:   vector with (minimum, maximum) mu value for re-estimation
  # tau_borders:  vector with (minimum, maximum) tau value for re-estimation
  # ln_borders:   vector with (minimum, maximum) lambda_n value for re-estimation
  # lc_borders:   vector with (minimum, maximum) lambda_c value for re-estimation
  
  # new_by_total_resolution:    optional; new/total transcript counts ratio interval resolution used for computation 
  #                             of the conditional probability 2 lookup table (given total counts, resolved by 
  #                             RNA dynamics parameters i.e. new/total ratios) (default: 1000)
  # ip_points:                  optional; number of points used for interpolations (default: 100)
  # search_points:              optional; peak-search parameter; number of points the likelihood functions is evaluated
  #                             at to find the peak (default: 1000)
  # width:                      optional; peak-search parameter; vector; width (as percent of a function's domain's
  #                             range) of a peak-spanning region considered for interpolation for likelihood_1,
  #                             likelihood_3, likelihood_4 (default: c(10, 10, 10))
  
  # repetitions:  number of repetitions for simulation (data sampling and following parameter re-estimation)
  #               (default: 1)
  
  # *** RETURN ***
  
  # a list containing (in the following order):
  #   a list containing the count data (vector of total nuclear counts, vector of modified nuclear counts, vector of
  #     total cytosolic counts, vector of modified cytosolic counts) for the transcript, 
  #   density functions nested list (first index: repetitions, second index: density / conditional probability 
  #     functions (4d joint raw parameters, lntau, lntau + lc 2d joint distribuion, steady_n, steady_c)),
  #   a list of vectors (one vector for each repetition) of estimated raw parameter settings (mu, tau, ln, lc)
  
  # Based on the underlying RNA dynamics model, this function samples transcript count data based on the RNA dynamics
  # parameter values given and estimates parameter values back using a Bayesian approach. The model assumes steady-
  # state conditions.
  
  # *******************************************************************************************************************
  
  # *******************************************************************************************************************
  # SIMULATING MODEL
  # *******************************************************************************************************************
  
  time_points <- length(time_series)    # number of time points
  
  # sampling count data
  # (returns: a list of vectors: total nuclear, modified nuclear, total cytosolic, modified cytosolic counts (each 
  # containing #repetitions count samples))
  
  print("Sampling")
  
  count_data <- sample_count_data(mu, tau, lambda_n, lambda_c, time_series, n_rna, c_rna, n_libsize, c_libsize,
                                  p_un, p_mod, number = repetitions)

  n_total_data <- count_data[[1]]
  n_mod_data <- count_data[[2]]
  c_total_data <- count_data[[3]]
  c_mod_data <- count_data[[4]]
    
  # computing parameter likelihood / density distributions (for each data sample repetition)
  # (result of sapply loop is a nested list, first index: repetitions, second index: density / probability functions)

  print("Computing Densities")
  
  densities <- sapply(1:repetitions, function(i) {   # iterating through repetitions
    print(i)
    list(parameter_estimation(n_rna, c_rna, n_libsize, c_libsize, n_total_data[, i], c_total_data[, i],
                              n_mod_data[, i], c_mod_data[, i], p_un, p_mod, 
                              mu_borders, tau_borders, ln_borders, lc_borders, time_series, 
                              p_mod_pre_n = p_mod_pre_n[[i]], p_mod_pre_c = p_mod_pre_c[[i]], 
                              new_by_total_resolution = new_by_total_resolution, ip_points = ip_points,
                              search_points = search_points, peak_border_ratio = peak_border_ratio,
                              min_width = min_width))
                              # computing derived and raw parameter likelihood / density distributions
  })

  # re-estimating parameter values from joint probability density (density maximum) (for each repetition)
  
  print("Determining Optima")
  
  optima <- sapply(1:repetitions, function(i) {
    
    print(i)                                  # printing progress
    joint_density <- densities[[i]][[1]]      # getting joint density of i-th repetition
    joint_density_inv <- function(args) {     # since global optimization only supports minimization, joint density
      joint_density(args) * (-1)              #   is inverted s.t. the minimum of inverted function corresponds to the
    }                                         #   maximum of the original function
    list(DEoptim(joint_density_inv, lower = c(mu_borders[1], tau_borders[1], ln_borders[1], lc_borders[1]),
                 upper = c(mu_borders[2], tau_borders[2], ln_borders[2], lc_borders[2]),
                 control = list(itermax = 10**5, trace = FALSE))$optim$bestmem)   # finding parameter values 
                 # maximizing probability density (minimizing inverted density)
  })
  
  # returning
  
  return(list(count_data, densities, optima))
}

halflife_estimation_leastsquare <- function(n_total, c_total, n_mod, c_mod, p_un_nu, p_un_cy, p_mod_nu, p_mod_cy,
                                            lntau_borders, lc_borders, time_series, grid_resolution = 100) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # n_total:      vector; nuclear total counts observed, resolved by time
  # c_total:      vector; cytosolic total counts observed, resolved by time
  # n_mod:        vector; nuclear modified counts observed, resolved by time
  # c_mod:        vector; cytosolic modified counts observed, resolved by time
  # p_un_nu:      probability that a nuclear pre-existing transcript remains unmodified
  # p_un_cy:      probability that a cytosolic pre-existing transcript remains unmodified
  # p_mod_nu:     probability that a nuclear newly synhesized transcript is modified
  # p_mod_cy:     probability that a cytosolic newly synhesized transcript is modified
  
  # lntau_borders:  vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:     vector with (minimum, maximum) lambda_c value
  # time_series:    vector of time points at which measurements were taken
  
  # grid_resolution:    optional; grid search resolution to use per parameter (default: 100)
  
  # *** RETURN ***
  
  # a list containing (in the following order):
  #   vector; lambda_n + tau values as used within the grid
  #   vector; lambda_c values as used within the grid
  #   list (nuclear and cytosolic time series) of matrices (rows: lambda_n + tau values, columns: lambda_c values),
  #     of simulated modified/total ratios
  #   matrix; grid sum of squares values, rows: lambda_n + tau, columns: lambda_c
  #   vector; grid least squares indices of (lambda_n + tau) and lambda_c
  
  # This function performs RNA dynamics parameter estimation using a simple least-squares approach on the half-life 
  # determining parameters lambda_n + tao amd lambda_c, given observed data.
  
  # *******************************************************************************************************************
  
  # initializations
  
  n_datapoints <- length(time_series) * 2     # number of data points (measurements)
  lnt_vals <- exp(
    seq(log(lntau_borders[1]), log(lntau_borders[2]), length.out = grid_resolution))  # lambda_n + tau grid values
  lc_vals <- exp(
    seq(log(lc_borders[1]), log(lc_borders[2]), length.out = grid_resolution))        # lambda_c grid values
  
  # computing grid of mod/total ratios (list of matrices, list indices: nu and cy time series, rows: lambda_n + tau
  # values, columns: lambda_c values)
  
  mod_by_total_grids <- halflife_grid(lntau_borders, lc_borders, time_series, p_un_nu, p_un_cy, p_mod_nu, p_mod_cy, 
                                      grid_resolution = grid_resolution)
  
  # computing sum of squares through grid, determining least squares grid indices
  # matrix; rows: lambda_n + tau values, columns: lambda_c values
  
  mod_by_total_nu_input <- n_mod / n_total * 100    # input data nucleus mod/total ratios in [%]
  mod_by_total_cy_input <- c_mod / c_total * 100    # input data cytosol mod/total ratios in [%]
  
  grid_sumofsquares <-
    sapply(1:grid_resolution, function(lc) {                            # iterating through lambda_c values
      sapply(1:grid_resolution, function(lnt) {                           # iterating through lambda_n + tau values
        
        sim_data <- sapply(1:n_datapoints, function(i) {                    # simulated mod/total ratios
          mod_by_total_grids[[i]][lnt, lc]
        })
        input_data <- c(mod_by_total_nu_input, mod_by_total_cy_input)       # data mod/total ratios
        
        sim_data_differences <- sim_data - input_data           # getting differences between simulation and data
        sim_data_squares <- sim_data_differences ** 2           # getting squares of differences
        return(sum(sim_data_squares))                           # returning summed squares
        
      })
    })
  
  # determining minimum sum of squares indices
  # vector; lamba_n + tau min index, lambda_c min index
  
  min_idx <- which(grid_sumofsquares == min(grid_sumofsquares))   # getting minimum sum-of-squares index
  min_idx_row <- min_idx %% grid_resolution                       # getting corresponding row index
  if (min_idx_row == 0) { min_idx_row = grid_resolution}          # (when hitting last row)
  min_idx_col <- ceiling(min_idx/grid_resolution)                 # getting corresponding column index
  grid_sumofsquares_minidx <- c(min_idx_row, min_idx_col)         # storing min indices
  
  # returning
  
  return(list(lnt_vals, lc_vals, mod_by_total_grids, grid_sumofsquares, grid_sumofsquares_minidx))
}

wrapper_halflife_maxlike <- function(summarytable_file, time_series, conveff_nu, conveff_cy,
                                     fn_error_nu, fn_error_cy, fp_error_nu, fp_error_cy,
                                     lntau_borders, lc_borders, mode = "binom", betabin_od = 2,
                                     ip_points = 1000, max_iter = 10000, parallel = 1) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summarytable_file:  file storing summary table (rows: genomic regions, col1: region name, col2: measurement
  #                     description, col3: library size, col4: total transcript counts, col5: labeled transcript  
  #                     counts, col6: average potential conversion positions, col7: conversion efficiency, 
  #                     col8: newly synthesized ratio)
  # time_series:        vector of time points at which measurements were taken
  # conveff_nu:         vector of nuclear conversion efficiencies at measurement time points; probability by which a 
  #                     single base is labeled
  # conveff_cy:         vector of cytosolic conversion efficiencies at measurement time points; probability by which a 
  #                     single base is labeled
  
  # fn_errorprogress.socket <- make.socket(port=1)_nu:        nuclear false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for nuclear transcripts)
  # fn_error_nu:        cytosolic false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for cytosolic transcripts)
  # fp_error_nu:        nuclear false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for nuclear transcripts)
  # fp_error_cy:        cytosolic false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for cytosolic transcripts)
  
  # lntau_borders:      vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:         vector with (minimum, maximum) lambda_c value
  
  # mode:         optional; set to "binom" for an underlying binomial distribution, set to "betabinom" for an  
  #               underlying beta-binomial distribution (default: "binom")
  # betabin_od:   needed for 'mode' == 'betabinom'; beta-binomial distribution overdispersion parameter
  #               (alpha + beta) (default: 2)
  # ip_points:    optional; number of points used for interpolations (default: 1000)
  # max_iter:     optional; maximum number of iteration to use for sum-of-squares minimization (default: 1000)
  # parallel:     optional; the number of cores to use for parallelization (default: 1, meaning no parallelization)
  
  # *** NOTE ***
  
  # The summary table file must be sorted by gene name, compartment, and time measurement (in that order). 
  
  # *** RETURN ***
  
  # a matrix storing single genes' parameter estimations (rows: genes, columns: lambda_n + tau and lambda_c 
  # estimations)
  
  # This function performs halflife maximum likelihood estimation for genomic regions given within a summary table file.
  # It is based on a Bayesian approach, assuming Binomial distributions for underlying count data. 
  
  # *******************************************************************************************************************
  
  # initializing ......................................................................................................
  
  stable <- load_summarytable(summarytable_file)              # loading in summary table
  n_timepoints <- length(time_series)                         # number of time points (time series measurements)
  n_genes <- length(stable[, 1]) / (n_timepoints * 2)         # number of genes recorded in summary table                      
  
  # preparing parallelization
  
  parallel_clusters <- makePSOCKcluster(parallel, outfile="")
  registerDoParallel(parallel_clusters)
  
  # iterating through genes of summary table, 
  # performing halflife estimation ....................................................................................
  
  gene_estimations <- 
    foreach (i = 0:(n_genes-1), .combine = rbind, .packages = c("DEoptim"), 
             .export = c("approx1D",  "ratio_n_sol", "ratio_c_sol", 
                         "p_nmod_nnew", "p_nmod_nnew_matrix", "p_nnew_ntotal", "p_nmod_theta",
                         "prob_thetaderiv_1", "prob_thetaderiv_2", "halflife_estimation_maxlike" 
                        ),
             .verbose = TRUE) %dopar% {
    
    print(paste("gene ", i + 1, " from ", n_genes, sep=""))
    
    start_idx_nu <- i * 2 * n_timepoints + 1          # row index within summary table where gene's nu entries start
    end_idx_nu <- start_idx_nu + n_timepoints - 1     # row index within summary table where gene's nu entries end
    start_idx_cy <- end_idx_nu + 1                    # row index within summary table where gene's cy entries start
    end_idx_cy <- start_idx_cy + n_timepoints - 1     # row index within summary table where gene's cy entries end
    
    # getting gene's measurement information
    
    submatrix_nu <- stable[start_idx_nu:end_idx_nu, ]     # getting gene's submatrix with nuclear entries
    submatrix_cy <- stable[start_idx_cy:end_idx_cy, ]     # getting gene's submatrix with cytosol entries
    
    n_total <- submatrix_nu[, "tot"]      # gene's nuclear total counts
    c_total <- submatrix_cy[, "tot"]      # gene's cytosol total counts
    n_mod <- submatrix_nu[, "mod"]        # gene's nuclear labeled counts
    c_mod <- submatrix_cy[, "mod"]        # gene's cytosol labeled counts
    n_ratio <- n_mod / n_total            # gene's nuclear mod/total ratios
    c_ratio <- c_mod / c_total            # gene's cytosolic mod/total ratios
    
    # getting modification probabilities
    
    convpos <- sum(c(submatrix_nu[, 6] * n_total, submatrix_cy[, 6] * c_total)) / 
      sum(c(n_total, c_total))                            # weighted average number of potential conversion positions
    p_mod_nu <- 1 - 
      (1 - conveff_nu * (1 - fn_error_nu)) ** convpos     # modification probability of nuclear newly synthesized reads
    p_mod_cy <- 1 - 
      (1 - conveff_cy * (1 - fn_error_cy)) ** convpos     # modification probability of cytosolic newly synthesized reads
    p_un_nu <- (1 - fp_error_nu) ** convpos               # non-modfication probability of nuclear pre-existing reads
    p_un_cy <- (1 - fp_error_cy) ** convpos               # non-modfication probability of cytosolic pre-existing reads
    
    # performing halflife estimation
    
    estimation_results <- 
      halflife_estimation_maxlike(n_total, c_total, n_mod, c_mod, 
                                  rep(p_un_nu, n_timepoints), rep(p_un_cy, n_timepoints), p_mod_nu, p_mod_cy, 
                                  lntau_borders, lc_borders, time_series, 
                                  mode = mode, betabin_od = betabin_od, ip_points = ip_points)[[1]]
    
    # storing estimation and corresponding sum of squares value
    
    lntau_choice <- estimation_results$optim$bestmem[[1]]                 # lambda_n + tau choice
    lc_choice <- estimation_results$optim$bestmem[[2]]                    # lambda_c choice
    maxval <- estimation_results$optim$bestval[[1]] * -1                  # corresponding maximum likelihood value
    return(c(i, lntau_choice, lc_choice, maxval))                         # returning parameter estimations
    
  }
  
  # finishing and returning estimation results ........................................................................
  
  stopCluster(parallel_clusters)                                    # giving free parallelization cores
  gene_estimations_sorted <- 
    gene_estimations[order(gene_estimations[, 1]),][, 2:4]          # sorting genes like in input summary table
  colnames(gene_estimations_sorted) <- c("lntau", "lc", "maxval")   # assigning column names
  return(gene_estimations)
}

wrapper_halflife_leastsquare_grid <- function(summarytable_file, time_series, conveff,
                                              fn_error_nu, fn_error_cy, fp_error_nu, fp_error_cy,
                                              lntau_borders, lc_borders, grid_resolution = 100) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summarytable_file:  file storing summary table (rows: genomic regions, col1: region name, col2: measurement
  #                     description, col3: library size, col4: total transcript counts, col5: labeled transcript  
  #                     counts, col6: average potential conversion positions, col7: conversion efficiency, 
  #                     col8: newly synthesized ratio)
  # time_series:        vector of time points at which measurements were taken
  # conveff:            conversion efficiency; probability by which a single base is labeled
  
  # fn_error_nu:        nuclear false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for nuclear transcripts)
  # fn_error_nu:        cytosolic false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for cytosolic transcripts)
  # fp_error_nu:        nuclear false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for nuclear transcripts)
  # fp_error_cy:        cytosolic false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for cytosolic transcripts)
  
  # lntau_borders:      vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:         vector with (minimum, maximum) lambda_c value
  
  # grid_resolution:    optional; grid search resolution to use per parameter (default: 100)
  
  # *** NOTE ***
  
  # The summary table file must be sorted by gene name, compartment, and time measurement (in that order). 
  
  # *** RETURN ***
  
  # a matrix storing single genes' parameter estimations (rows: genes, columns: lambda_n + tau and lambda_c 
  # estimations)
  
  # This function performs halflife leastsquares estimation for genomic regions given within a summary table file.
  
  # *******************************************************************************************************************
  
  # initializing
  
  stable <- load_summarytable(summarytable_file)              # loading in summary table
  n_timepoints <- length(time_series)                         # number of time points (time series measurements)
  n_genes <- length(stable[, 1]) / (n_timepoints * 2)         # number of genes recorded in summary table                      
  
  gene_estimations <- matrix(rep(0, 2 * n_genes), ncol = 2)   # matrix storing parameter estimations, 
  rownames(gene_estimations) <- rep("", n_genes)              # rows: genes
  colnames(gene_estimations) <- c("lntau", "lc")              # columns: lntau and lc estimations
  
  # iterating through genes of summary table, performing halflife estimation
  
  for (i in 0:(n_genes-1)) {
    
    start_idx_nu <- i * 2 * n_timepoints + 1          # row index within summary table where gene's nu entries start
    end_idx_nu <- start_idx_nu + n_timepoints - 1     # row index within summary table where gene's nu entries end
    start_idx_cy <- end_idx_nu + 1                    # row index within summary table where gene's cy entries start
    end_idx_cy <- start_idx_cy + n_timepoints - 1     # row index within summary table where gene's cy entries end
    
    # getting gene's measurement information
    
    submatrix_nu <- stable[start_idx_nu:end_idx_nu, ]     # getting gene's submatrix with nuclear entries
    submatrix_cy <- stable[start_idx_cy:end_idx_cy, ]     # getting gene's submatrix with cytosol entries
    
    n_total <- submatrix_nu[, 4]      # gene's nuclear total counts
    c_total <- submatrix_cy[, 4]      # gene's cytosol total counts
    n_mod <- submatrix_nu[, 5]        # gene's nuclear labeled counts
    c_mod <- submatrix_cy[, 5]        # gene's cytosol labeled counts
    
    convpos <- sum(c(submatrix_nu[, 6] * n_total, submatrix_cy[, 6] * c_total)) / 
      sum(c(n_total, c_total))                        # weighted average number of potential conversion positions
    p_mod_nu <- 1 - 
      (1 - conveff * (1 - fn_error_nu)) ** convpos    # modification probability of nuclear newly synthesized reads
    p_mod_cy <- 1 - 
      (1 - conveff * (1 - fn_error_cy)) ** convpos    # modification probability of cytosolic newly synthesized reads
    p_un_nu <- (1 - fp_error_nu) ** convpos           # non-modfication probability of nuclear pre-existing reads
    p_un_cy <- (1 - fp_error_cy) ** convpos           # non-modfication probability of cytosolic pre-existing reads
    
    # performing halflife estimation
    
    estimation_results <- halflife_estimation_leastsquare(n_total, c_total, n_mod, c_mod, 
                                                          p_un_nu, p_un_cy, p_mod_nu, p_mod_cy, 
                                                          lntau_borders, lc_borders, time_series, 
                                                          grid_resolution = grid_resolution)
    lntau_choice <- estimation_results[[1]][estimation_results[[5]][1]]   # lambda_n + tau choice
    lc_choice <- estimation_results[[2]][estimation_results[[5]][2]]      # lambda_c choice
    rownames(gene_estimations)[i + 1] <- submatrix_nu[1, 1]               # storing genomic region name
    gene_estimations[i + 1, ] <- c(lntau_choice, lc_choice)               # storing parameter estimations
    
  }
  
  # returning estimation results
  
  return(gene_estimations)
  
}

wrapper_halflife_leastsquare_optim_labtot <- function(summarytable_file, time_series, conveff_nu, conveff_cy, 
                                                      fn_error_nu, fn_error_cy, fp_error_nu, fp_error_cy,
                                                      lntau_borders, lc_borders, max_iter = 10000) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summarytable_file:  file storing summary table (rows: genomic regions, col1: region name, col2: measurement
  #                     description, col3: library size, col4: total transcript counts, col5: labeled transcript  
  #                     counts, col6: average potential conversion positions, col7: conversion efficiency, 
  #                     col8: newly synthesized ratio)
  # time_series:        vector of time points at which measurements were taken
  # conveff_nu:         vector of nuclear conversion efficiencies at measurement time points; probability by which a 
  #                     single base is labeled
  # conveff_cy:         vector of cytosolic conversion efficiencies at measurement time points; probability by which a 
  #                     single base is labeled
  
  # fn_error_nu:        nuclear false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for nuclear transcripts)
  # fn_error_nu:        cytosolic false-negative error; probability by which a single, labeled U position (i.e. 4sU
  #                     converted to C) does not show up as C position in the alignment (can be set to the C > non-C 
  #                     conversion rate) (for cytosolic transcripts)
  # fp_error_nu:        nuclear false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for nuclear transcripts)
  # fp_error_cy:        cytosolic false-positive error; probability by which a single, unlabeled U position shows up as 
  #                     C (i.e. labeled) position in the alignment (for cytosolic transcripts)
  
  # lntau_borders:      vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:         vector with (minimum, maximum) lambda_c value
  
  # max_iter:     optional; maximum number of iteration to use for sum-of-squares minimization (default: 1000)
  
  # *** NOTE ***
  
  # The summary table file must be sorted by gene name, compartment, and time measurement (in that order). 
  
  # *** RETURN ***
  
  # a matrix storing single genes' parameter estimations (rows: genes, columns: lambda_n + tau and lambda_c 
  # estimations)
  
  # This function performs halflife leastsquares estimation for genomic regions given within a summary table file.
  # It is based on comparing labeled/total ratios, inferring the labeled/total ratios observed for the new/total ratios
  # of the ODE prediction using the error estimated, comparing to the measured labeled/total ratios. 
  
  # *******************************************************************************************************************
  
  # initializing ......................................................................................................
  
  stable <- load_summarytable(summarytable_file)              # loading in summary table
  n_timepoints <- length(time_series)                         # number of time points (time series measurements)
  n_genes <- length(stable[, 1]) / (n_timepoints * 2)         # number of genes recorded in summary table                      
  
  gene_estimations <- matrix(rep(0, 3 * n_genes), ncol = 3)   # matrix storing parameter estimations, 
  rownames(gene_estimations) <- rep("", n_genes)              # rows: genes, columns: lntau estimator,
  colnames(gene_estimations) <- c("lntau", "lc", "sos")       # lc estimator, minimum sum of squares
  
  # iterating through genes of summary table, 
  # performing halflife estimation ....................................................................................
  
  for (i in 0:(n_genes-1)) {
    
    start_idx_nu <- i * 2 * n_timepoints + 1          # row index within summary table where gene's nu entries start
    end_idx_nu <- start_idx_nu + n_timepoints - 1     # row index within summary table where gene's nu entries end
    start_idx_cy <- end_idx_nu + 1                    # row index within summary table where gene's cy entries start
    end_idx_cy <- start_idx_cy + n_timepoints - 1     # row index within summary table where gene's cy entries end
    
    # getting gene's measurement information
    
    submatrix_nu <- stable[start_idx_nu:end_idx_nu, ]     # getting gene's submatrix with nuclear entries
    submatrix_cy <- stable[start_idx_cy:end_idx_cy, ]     # getting gene's submatrix with cytosol entries
    
    n_total <- submatrix_nu[, 4]      # gene's nuclear total counts
    c_total <- submatrix_cy[, 4]      # gene's cytosol total counts
    n_mod <- submatrix_nu[, 5]        # gene's nuclear labeled counts
    c_mod <- submatrix_cy[, 5]        # gene's cytosol labeled counts
    n_ratio <- n_mod / n_total        # gene's nuclear mod/total ratios
    c_ratio <- c_mod / c_total        # gene's cytosolic mod/total ratios
    
    # getting modification probabilities
    
    convpos <- sum(c(submatrix_nu[, 6] * n_total, submatrix_cy[, 6] * c_total)) / 
      sum(c(n_total, c_total))                            # weighted average number of potential conversion positions
    p_mod_nu <- 1 - 
      (1 - conveff_nu * (1 - fn_error_nu)) ** convpos     # modification probability of nuclear newly synthesized reads
    p_mod_cy <- 1 - 
      (1 - conveff_cy * (1 - fn_error_cy)) ** convpos     # modification probability of cytosolic newly synthesized reads
    p_un_nu <- (1 - fp_error_nu) ** convpos               # non-modfication probability of nuclear pre-existing reads
    p_un_cy <- (1 - fp_error_cy) ** convpos               # non-modfication probability of cytosolic pre-existing reads
    
    # performing halflife estimation
    
    estimation_results <- optim(c(0.01, 0.01), halflife_sumofsquares_labtot, 
                                nu_measured = n_ratio, cy_measured = c_ratio, time_series = time_series, 
                                p_un_nu = p_un_nu, p_un_cy = p_un_cy, p_mod_nu = p_mod_nu, p_mod_cy = p_mod_cy,
                                method = "L-BFGS-B", lower = c(lntau_borders[1], lc_borders[1]), 
                                upper = c(lntau_borders[2], lc_borders[2]), control = list(maxit = max_iter))
    
    # storing estimation and corresponding sum of squares value
    
    lntau_choice <- estimation_results$par[1]         # lambda_n + tau choice
    lc_choice <- estimation_results$par[2]            # lambda_c choice
    sum_of_squares <- estimation_results$value        # corresponding sum of squares
    rownames(gene_estimations)[i + 1] <- submatrix_nu[1, 1]                   # storing genomic region name
    gene_estimations[i + 1, ] <- c(lntau_choice, lc_choice, sum_of_squares)   # storing parameter estimations
    
  }
  
  # returning estimation results ......................................................................................
  
  return(gene_estimations)
}

wrapper_halflife_leastsquare_optim_newtot <- function(summarytable_file, time_series, lntau_borders, lc_borders, 
                                                      var_stab = TRUE, seed = c(0.00139, 0.00139),
                                                      max_iter = 10000, accuracy = 10**7) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summarytable_file:  file storing summary table (rows: genomic regions, col1: region name, col2: measurement
  #                     description, col3: library size, col4: total transcript counts, col5: labeled transcript  
  #                     counts, col6: average potential conversion positions, col7: conversion efficiency, 
  #                     col8: newly synthesized ratio)
  # time_series:        vector of time points at which measurements were taken
  # lntau_borders:      vector with (minimum, maximum) lambda_n + tau value 
  # lc_borders:         vector with (minimum, maximum) lambda_c value
  
  # var_stab:       optional; set to TRUE to perform variance-stabilizing transformation of the data points 
  #                 (default: TRUE)
  # seed:           optional; starting values for half life parameters in sum-of-squares minimization
  #                 (default: c(0.00139, 0.00139) corresponding to a half life of (499min, 499min))
  # max_iter:       optional; maximum number of iteration to use for sum-of-squares minimization (default: 1000)
  # accuracy:       optional; accuracy (as factor of machine tolerance) for sum-of-squares minimization (smaller values
  #                 yield better accuracy) (default: 1e7)
  
  # *** NOTE ***
  
  # The summary table file must be sorted by gene name, compartment, and time measurement (in that order).
  
  # *** RETURN ***
  
  # a matrix storing single genes' parameter estimations (rows: genes, columns: lambda_n + tau and lambda_c 
  # estimations)
  
  # This function performs halflife leastsquares estimation for genomic regions given within a summary table file.
  # It is based on comparing new/total ratios as defined within the summary table (those are EM-derived estimates
  # using a fixed conversion efficiency at each time point) with the ODE predictions. 
  
  # *******************************************************************************************************************
  
  # initializing ......................................................................................................
  
  stable <- load_summarytable(summarytable_file)              # loading in summary table
  n_timepoints <- length(time_series)                         # number of time points (time series measurements)
  n_genes <- length(stable[, 1]) / (n_timepoints * 2)         # number of genes recorded in summary table   
  
  gene_estimations <- matrix(rep(0, 3 * n_genes), ncol = 3)   # matrix storing parameter estimations, 
  rownames(gene_estimations) <- rep("", n_genes)              # rows: genes, columns: lntau estimator,
  colnames(gene_estimations) <- c("lntau", "lc", "sos")       # lc estimator, minimum sum of squares
  
  # iterating through genes of summary table, 
  # performing halflife estimation ....................................................................................
  
  for (i in 0:(n_genes-1)) {
    
    print(i)    
    
    start_idx_nu <- i * 2 * n_timepoints + 1          # row index within summary table where gene's nu entries start
    end_idx_nu <- start_idx_nu + n_timepoints - 1     # row index within summary table where gene's nu entries end
    start_idx_cy <- end_idx_nu + 1                    # row index within summary table where gene's cy entries start
    end_idx_cy <- start_idx_cy + n_timepoints - 1     # row index within summary table where gene's cy entries end
    
    # getting gene's measurement information
    
    submatrix_nu <- stable[start_idx_nu:end_idx_nu, ]     # getting gene's submatrix with nuclear entries
    submatrix_cy <- stable[start_idx_cy:end_idx_cy, ]     # getting gene's submatrix with cytosol entries
    
    n_total <- submatrix_nu[, "tot"]        # gene's nuclear total counts
    c_total <- submatrix_cy[, "tot"]        # gene's cytosolic total counts
    n_ratio <- submatrix_nu[, "nr"]         # gene's nuclear new/total ratios
    c_ratio <- submatrix_cy[, "nr"]         # gene's cytosolic new/total ratios

    # performing halflife estimation
    
    # estimation_results <- optim(seed, halflife_sumofsquares_newtot, time_series = time_series, 
    #                             new_tot_nu = n_ratio, new_tot_cy = c_ratio, total_nu = n_total, total_cy = c_total,
    #                             var_stab = var_stab,
    #                             method = "L-BFGS-B", 
    #                             lower = c(lntau_borders[1], lc_borders[1]), upper = c(lntau_borders[2], lc_borders[2]), 
    #                             control = list(maxit = max_iter, factr = accuracy))
    
    estimation_results <- DEoptim(halflife_sumofsquares_newtot, 
                                  lower = c(lntau_borders[1], lc_borders[1]), upper = c(lntau_borders[2], lc_borders[2]),
                                  control=list(itermax=max_iter, trace=FALSE, reltol=10**(-12)), 
                                  time_series = time_series, new_tot_nu = n_ratio, new_tot_cy = c_ratio, 
                                  total_nu = n_total, total_cy = c_total, var_stab = var_stab)
    
    # storing estimation and corresponding sum of squares value
    
    # lntau_choice <- estimation_results$par[1]         # lambda_n + tau choice
    # lc_choice <- estimation_results$par[2]            # lambda_c choice
    # sum_of_squares <- estimation_results$value        # corresponding sum of squares
    # rownames(gene_estimations)[i + 1] <- submatrix_nu[1, 1]                   # storing genomic region name
    # gene_estimations[i + 1, ] <- c(lntau_choice, lc_choice, sum_of_squares)   # storing parameter estimations
    
    lntau_choice <- estimation_results$optim$bestmem[1]         # lambda_n + tau choice
    lc_choice <- estimation_results$optim$bestmem[2]            # lambda_c choice
    sum_of_squares <- estimation_results$optim$bestval          # corresponding sum of squares
    rownames(gene_estimations)[i + 1] <- submatrix_nu[1, 1]                   # storing genomic region name
    gene_estimations[i + 1, ] <- c(lntau_choice, lc_choice, sum_of_squares)   # storing parameter estimations
    
  }
  
  # returning estimation results ......................................................................................
  
  return(gene_estimations)
}

#######################################################################################################################
# EXECUTION (MAIN)
#######################################################################################################################





