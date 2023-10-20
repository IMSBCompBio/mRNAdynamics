# load the library for Metropolis-Hastings sampling
library(mcmc)
# load the library for ggplot
library(ggplot2)
# load the library for the heatscatter plot
library(LSD)
# load the library for DEoptim
library(DEoptim)
# load the libraries for parallelization
library(foreach)
library(doParallel)

###
### Estimation of the nuclear an cytosolic half life parameters for all genes
###

# Define the loss function (a negative log likelihood) for later optimization 
negloss = function(params,
                   t1,
                   traf_nuc,
                   tot_nuc,
                   traf_cyt,
                   tot_cyt){
  
  if (any(params <= 0)) return(Inf) # reject impossible parameters
  if (any(params >= 1)) return(Inf)
  
  # Prediction of the nuclear new/total ratio for the given timepoints
  predict_nuc = 1 - exp(- params[1] * t1)  
  
  # Prediciton of the cytosolic new/total ratio for the given timepoints
  predict_cyt = 1 - ( exp(- params[2] * t1) +
                        params[2] / (params[2] - params[1]) * (exp(- params[1] * t1) - 
                                                                 exp(- params[2] * t1)) )
  if (any(is.na(c(predict_nuc, predict_cyt)))) return(Inf)  # reject impossible outcomes
  if (any(c(predict_nuc, predict_cyt) < 0)) return(Inf)
  
  # - residual sum of squares for the comparison with the asin(sqrt()) transformed,
  # labeling bias corrected, observed new/total ratios in nucleus resp. cytosol
  res = - sum( (asin(sqrt(predict_nuc)) - traf_nuc)^2 * 2 * tot_nuc) -
    sum( (asin(sqrt(predict_cyt)) - traf_cyt)^2 * 2 * tot_cyt)
  return(-res)
}

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

# Define the loss function (a negative log likelihood) for the cytosolic apartment
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

### doMCMC takes the nuclear data corresponding to one gene in one experiment 
### and outputs an MCMC sample of desired size
doMCMC_nuc = function(timepoints,
                      transformed_nuc,
                      total_nuc,
                      distinct_samples = 10^4,
                      min_acceptance = 0.15,
                      max_acceptance = 0.45,
                      start_scale = 1,
                      start_param = c(0.008)){
  
  # define the loss functions in a way that is suitable for the metrop function
  loss_nuc = function(param,
                      tp=timepoints,
                      traf_nuc=transformed_nuc,
                      tot_nuc=total_nuc){
    return(-1 * negloss_nuc(param, tp, traf_nuc, tot_nuc))
  }
  
  # find a suitable scale for the MCMC proposal function 
  scale_nuc = start_scale
  res_nuc = metrop(loss_nuc,initial=start_param,nbatch=10^4,scale=scale_nuc)
  
  # increase the scale by a factor of 2 if the acceptance rate is too high
  counter = 0
  repeat{
    if (res_nuc$accept < max_acceptance | counter == 100) break()
    scale_nuc = scale_nuc * 2
    res_nuc = metrop(res_nuc,nbatch=10^4,scale=scale_nuc)
    counter = counter + 1
  }
  # decrease the scale by a factor of 2 if the acceptance rate is too low
  counter = 0
  repeat{
    if (res_nuc$accept > min_acceptance | counter == 100) break()
    scale_nuc = scale_nuc / 2
    res_nuc = metrop(res_nuc,nbatch=10^4,scale=scale_nuc)
    counter = counter + 1
  }
  
  # perform the final MCMC chain which, 
  # in expectation, contains the desired number of distinct samples
  res_nuc = metrop(res_nuc,nbatch=round(distinct_samples/res_nuc$accept),scale=scale_nuc)
  
  # return the final result
  return(res_nuc)
  
}

### doMCMC takes the cytosolic data corresponding to one gene in one experiment 
### and outputs an MCMC sample of desired size
doMCMC_cyt = function(timepoints,
                      deg_nuc,
                      transformed_cyt,
                      total_cyt,
                      distinct_samples = 10^4,
                      min_acceptance = 0.15,
                      max_acceptance = 0.45,
                      start_scale = 1,
                      start_param = c(0.01)){
  
  # define the loss function in a way that is suitable for the metrop function
  loss_cyt = function(param,
                      tp=timepoints,
                      deg_n=deg_nuc,
                      traf_cyt=transformed_cyt,
                      tot_cyt=total_cyt){
    return(-1 * negloss_cyt(param, tp, deg_n, traf_cyt, tot_cyt))
  }
  
  # find a suitable scale for the MCMC proposal function 
  scale_cyt = start_scale
  res_cyt = metrop(loss_cyt,initial=start_param,nbatch=10^4,scale=scale_cyt)
  
  # increase the scale by a factor of 2 if the acceptance rate is too high
  counter = 0
  repeat{
    if (res_cyt$accept < max_acceptance | counter == 100) break()
    scale_cyt = scale_cyt * 2
    res_cyt = metrop(res_cyt,nbatch=10^4,scale=scale_cyt)
    counter = counter + 1
  }
  # decrease the scale by a factor of 2 if the acceptance rate is too low
  counter = 0
  repeat{
    if (res_cyt$accept > min_acceptance | counter == 100) break()
    scale_cyt = scale_cyt / 2
    res_cyt = metrop(res_cyt,nbatch=10^4,scale=scale_cyt)
    counter = counter + 1
  }
  
  # perform the final MCMC chain which, 
  # in expectation, contains the desired number of distinct samples
  res_cyt = metrop(res_cyt,nbatch=round(distinct_samples/res_cyt$accept),scale=scale_cyt)
  
  # return the final result
  return(res_cyt)
  
}

### doMCMC takes the data corresponding to one gene in one experiment 
### and outputs an MCMC sample of desired size
doMCMC_2d = function(timepoints,
                  transformed_nuc,
                  total_nuc,
                  transformed_cyt,
                  total_cyt,
                  distinct_samples = 10^5,
                  min_acceptance = 0.15,
                  max_acceptance = 0.45,
                  start_scale = 1,
                  start_params = c(0.008,0.01)){
  
  # define the loss function in a way that is suitable for the metrop function
  # params is a 2-dim parameter vector containing (deg_nuc = lambda_nuc+tau, lambda_cyt)
  loss = function(params,
                  t1=timepoints,
                  traf_nuc=transformed_nuc,
                  tot_nuc=total_nuc,
                  traf_cyt=transformed_cyt,
                  tot_cyt=total_cyt){
    if (any(params <= 0)) return(-Inf) # reject impossible parameters
    if (any(params >= 1)) return(-Inf)
    
    # prediction of the nuclear new/total ratio for the given timepoints
    predict_nuc = 1 - exp(- params[1] * t1)  
    
    # prediciton of the cytosolic new/total ratio for the given timepoints
    predict_cyt = 1 - ( exp(- params[2] * t1) +
                          params[2] / (params[2] - params[1]) * (exp(- params[1] * t1) - 
                                                                   exp(- params[2] * t1)) ) 
    if (any(is.na(c(predict_nuc, predict_cyt)))) return(-Inf) # reject impossible outcomes
    if (any(c(predict_nuc, predict_cyt) < 0)) return(-Inf)
    
    # - residual sum of squares for the comparison with the asin(sqrt()) transformed,
    # labeling bias corrected, observed new/total ratios in nucleus resp. cytosol
    res = - sum( (asin(sqrt(predict_nuc)) - traf_nuc)^2 * 2 * tot_nuc) -
      sum( (asin(sqrt(predict_cyt)) - traf_cyt)^2 * 2 * tot_cyt)
    return(res)
  }
  
  # find a suitable scale for the MCMC proposal function 
  scale = start_scale
  res = metrop(loss,initial=start_params,nbatch=10^4,scale=scale)
  
  # increase the scale by a factor of 2 if the acceptance rate is too high
  repeat{
    if (res$accept < max_acceptance) break()
    scale = scale * 2
    res = metrop(res,nbatch=10^4,scale=scale)    
  }
  # decrease the scale by a factor of 2 if the acceptance rate is too low
  repeat{
    if (res$accept > min_acceptance) break()
    scale = scale / 2
    res = metrop(res,nbatch=10^4,scale=scale)
  }
  
  # perform the final MCMC chain which, 
  # in expectation, contains the desired number of distinct samples
  res = metrop(res,nbatch=round(distinct_samples/res$accept),scale=scale)
  
  # return the final result
  return(res)
  
}

### fitTimeseries takes 4 time series and performs an optimization for the maximum likelihood parameter estimates,
### first for the nuclear and subsequently for the cytosolic degradation rate
fitTimeseries = function(
  timepoints,
  transformed_nuc,
  total_nuc,
  transformed_cyt,
  total_cyt,
  include_samples = F, # if TRUE, the whole MCMC chain will also be returned
  quantiles = c(0.025,0.25,0.5,0.75,0.975), # quantiles to be reported
  distinct_samples = 10^4,
  min_acceptance = 0.15,
  max_acceptance = 0.45,
  start_scale = 1,
  start_params = c(0.008,0.01)){
  

  # Estimate the maximum likelihood distributions with MCMC sampling and
  # degradation rates with optim
  res_nuc = doMCMC_nuc(timepoints=timepoints,
                       transformed_nuc=transformed_nuc,
                       total_nuc=total_nuc,
                       distinct_samples = distinct_samples,
                       min_acceptance = min_acceptance,
                       max_acceptance = max_acceptance,
                       start_scale = start_scale,
                       start_param = c(start_params[1]))
  ML_deg_nuc = median(res_nuc$batch)
  deg_nuc = optim(c(ML_deg_nuc),negloss_nuc,
                  tp=timepoints,traf_nuc=transformed_nuc,tot_nuc=total_nuc,
                  method="Brent", lower=0, upper=1)$par
  names(deg_nuc) = c("deg_nuc")
  
  res_cyt = doMCMC_cyt(timepoints=timepoints,
                       deg_nuc = deg_nuc,
                       transformed_cyt=transformed_cyt,
                       total_cyt=total_cyt,  
                       distinct_samples = distinct_samples,
                       min_acceptance = min_acceptance,
                       max_acceptance = max_acceptance,
                       start_scale = start_scale,
                       start_param = c(start_params[2]))
  ML_deg_cyt = median(res_cyt$batch)
  deg_cyt = optim(c(ML_deg_cyt),negloss_cyt,
                  tp=timepoints,deg_nuc=deg_nuc,traf_cyt=transformed_cyt,tot_cyt=total_cyt,
                  method="Brent", lower=0, upper=1)$par
  names(deg_cyt) = c("deg_cyt")
  
  # Calculate degradation rate quantiles
  quantiles_deg_nuc = quantile(res_nuc$batch[,1],probs=quantiles)
  quantile_names = names(quantiles_deg_nuc)
  names(quantiles_deg_nuc) = paste("deg_nuc_q",quantile_names,sep="")
  
  quantiles_deg_cyt = quantile(res_cyt$batch[,1],probs=quantiles)
  names(quantiles_deg_cyt) = paste("deg_cyt_q",quantile_names,sep="")
  
  # Calculate the half lifes and half life quantiles (order needs to be reversed!)
  half_life_nuc <- log(2)/deg_nuc
  half_life_cyt <- log(2)/deg_cyt
  names(half_life_nuc) <- c("half_life_nuc")
  names(half_life_cyt) <- c("half_life_cyt")
  
  quantiles_halflife_nuc = rev(log(2)/quantiles_deg_nuc)
  names(quantiles_halflife_nuc) = paste("halflife_nuc_q",quantile_names,sep="")
  
  quantiles_halflife_cyt = rev(log(2)/quantiles_deg_cyt)
  names(quantiles_halflife_cyt) = paste("halflife_cyt_q",quantile_names,sep="")
  
  # Calculate the nuclear and the cytosolic Rsquared value (= % of explained variance)
  predict_nuc = 1 - exp(- deg_nuc * timepoints)  
  rsquared_nuc = 1 - sum( (asin(sqrt(predict_nuc)) - transformed_nuc)^2 ) / var(transformed_nuc) 
  names(rsquared_nuc) = "Rsquared_nuc"
  
  predict_cyt = 1 - ( exp(- deg_cyt * timepoints) +
                        deg_cyt / (deg_cyt - deg_nuc) * (exp(- deg_nuc * timepoints) - 
                                                                 exp(- deg_cyt * timepoints)) ) 
  rsquared_cyt = 1 - sum( (asin(sqrt(predict_cyt)) - transformed_cyt)^2 ) / var(transformed_cyt) 
  names(rsquared_cyt) = "Rsquared_cyt"
  
  # Calculate the mean coverage across the total nuc / total cyt time series
  mean_coverage_nuc = mean(total_nuc)
  names(mean_coverage_nuc) = "mean_coverage_nuc"
  mean_coverage_cyt = mean(total_cyt)
  names(mean_coverage_cyt) = "mean_coverage_cyt"
  
  # Calculate the coefficient of variation of the total nuc / total cyt time series
  cv_nuc = sd(total_nuc)/mean(total_nuc)
  names(cv_nuc) = "cv_nuc"
  cv_cyt = sd(total_cyt)/mean(total_cyt)
  names(cv_cyt) = "cv_cyt"
  
  # Compose the list of output variables
  final = list(deg_nuc = deg_nuc,
               deg_cyt = deg_cyt,
               halflife_nuc = half_life_nuc,
               halflife_cyt = half_life_cyt,
               quantiles_deg_nuc = quantiles_deg_nuc,
               quantiles_deg_cyt = quantiles_deg_cyt,
               quantiles_halflife_nuc = quantiles_halflife_nuc,
               quantiles_halflife_cyt = quantiles_halflife_cyt,
               mean_coverage_nuc = mean_coverage_nuc,
               mean_coverage_cyt = mean_coverage_cyt,
               cv_nuc = cv_nuc,
               cv_cyt = cv_cyt,
               rsquared_nuc = rsquared_nuc,
               rsquared_cyt = rsquared_cyt,
               nsamples_nuc = res_nuc$nbatch,
               nsamples_cyt = res_cyt$nbatch,
               acceptance_nuc = res_nuc$accept,
               acceptance_cyt = res_cyt$accept,
               time_nuc = res_nuc$time[1],
               time_cyt = res_cyt$time[1])
  if (include_samples){ final$samples_nuc = res_nuc$batch; final$samples_cyt = res_cyt$batch}
  
  return(final)    
}

### fitGene is a wrapper for fitTimeseries
fitGene = function(gene, # name of the 3'UTR to be fitted
                   timepoints, # vector of effective timepoints
                   rawdata, # data object containing all processed data of one experiment
                   include_samples = F, # if TRUE, the whole MCMC chain will also be returned
                   quantiles = c(0.025,0.25,0.5,0.75,0.975), # quantiles to be reported
                   distinct_samples = 10^4,
                   min_acceptance = 0.15,
                   max_acceptance = 0.45,
                   start_scale = 1){
  
  # fetch the required time series corresponding to <gene> from the raw_data object  
  # obtain the data for nuclear totals and new/total ratios
  total_nuc = rawdata$tot_nu[gene,]
  ratio_nuc = rawdata$ratio_nu[gene,]
  transformed_nuc = asin(sqrt(ratio_nuc))
  
  # obtain the data for cytosolic totals and new/total ratios
  total_cyt = rawdata$tot_cy[gene,]
  ratio_cyt = rawdata$ratio_cy[gene,]
  transformed_cyt = asin(sqrt(ratio_cyt))
  
  # Find the maximum likelihood seed for MCMC chain, for nuclear degradation
  nuc_seed <- DEoptim(negloss_nuc,
                      lower = c(0), upper = c(1),
                      control = list(itermax = 500, trace=FALSE),
                      tp=timepoints, traf_nuc=transformed_nuc, tot_nuc=total_nuc)$optim$bestmem
  
  # Find the maximum likelihood seed for MCMC chain, for cytosolic degradation
  cyt_seed <- DEoptim(negloss_cyt,
                      lower = c(0), upper = c(1),
                      control = list(itermax = 500, trace=FALSE),
                      tp=timepoints, deg_nuc=nuc_seed, traf_cyt=transformed_cyt, tot_cyt=total_cyt)$optim$bestmem
  
  # call fitTimeseries and return its output
  return(
    fitTimeseries( timepoints = timepoints,
                   transformed_nuc = transformed_nuc,
                   total_nuc = total_nuc,
                   transformed_cyt = transformed_cyt,
                   total_cyt = total_cyt,
                   include_samples = include_samples, 
                   quantiles = quantiles,
                   distinct_samples = distinct_samples,
                   min_acceptance = min_acceptance,
                   max_acceptance = max_acceptance,
                   start_scale = start_scale,
                   start_params = c(nuc_seed, cyt_seed))
  )
  
}

### Plot nuclear / cytosolic ratio time series and the maximum likelihood fit
plotGene = function(gene,timepoints,rawdata,results_table){
  ratio_nuc = rawdata$ratio_nu[gene,]
  ratio_cyt = rawdata$ratio_cy[gene,]
  deg_nuc = results_table[gene,"deg_nuc"]
  deg_cyt = results_table[gene,"deg_cyt"]
  predict_nuc =  1 - exp(- deg_nuc * timepoints)  
  predict_cyt =  1 - ( exp(- deg_cyt * timepoints) +
                         deg_cyt / (deg_cyt - deg_nuc) * (exp(- deg_nuc * timepoints) - 
                                                            exp(- deg_cyt * timepoints)) ) 
  plot(timepoints,ratio_nuc,col="blue",lwd=2,main=gene,ylim=c(0,1),type="n")
  abline(h=0:5/5,col="light grey",lty=3)
  points(timepoints,ratio_nuc,type="b",col="blue",lwd=2)
  points(timepoints,ratio_cyt,type="b",col="red",lwd=2)
  points(timepoints,predict_nuc,type="l",col="cyan",lwd=2)
  points(timepoints,predict_cyt,type="l",col="pink",lwd=2)
}

### Plot an MCMC sample of a measurement time series
plotMCMC_1d = function(gene,timepoints,rawdata,results_table,distinct_samples=10^4){
  
  total_nuc = rawdata$tot_nu[gene,]
  ratio_nuc = rawdata$ratio_nu[gene,]
  transformed_nuc = asin(sqrt(ratio_nuc))
  total_cyt = rawdata$tot_cy[gene,]
  ratio_cyt = rawdata$ratio_cy[gene,]
  transformed_cyt = asin(sqrt(ratio_cyt))
  estimates = c(results_table[gene,c("deg_nuc","deg_cyt")])
  
  # Obtain an MCMC sample of desired size, starting at the estimated values
  res_nuc = doMCMC_nuc(timepoints=timepoints,
                       transformed_nuc=transformed_nuc,
                       total_nuc=total_nuc,
                       distinct_samples = distinct_samples,
                       start_param = estimates[1])
  
  res_cyt = doMCMC_cyt(timepoints=timepoints,
                       deg_nuc = estimates[1],
                       transformed_cyt=transformed_cyt,
                       total_cyt=total_cyt,  
                       distinct_samples = distinct_samples,
                       start_param = estimates[2])
  
  # Histograms of the MCMC sample
  res_nuc <- data.frame(matrix(log(2)/res_nuc$batch, ncol=1))
  res_cyt <- data.frame(matrix(log(2)/res_cyt$batch, ncol=1))
  
  hist_nuc <- ggplot(data=res_nuc) + 
    geom_density(aes(x=res_nuc[, 1]), colour="black", fill="grey85") + 
    theme_minimal() + xlab("nuclear half life")
  hist_cyt <- ggplot(data=res_cyt) + 
    geom_density(aes(x=res_cyt[, 1]), colour="black", fill="grey85") + 
    theme_minimal() + xlab("cytosolic half life")
  
  return(list(hist_nuc, hist_cyt))
}

plotMCMC_2d = function(gene,timepoints,rawdata,results_table,outfile_prefix,distinct_samples=10^4){
  
  total_nuc = rawdata$tot_nu[gene,]
  ratio_nuc = rawdata$ratio_nu[gene,]
  transformed_nuc = asin(sqrt(ratio_nuc))
  total_cyt = rawdata$tot_cy[gene,]
  ratio_cyt = rawdata$ratio_cy[gene,]
  transformed_cyt = asin(sqrt(ratio_cyt))
  start_params = c(results_table[gene,c("deg_nuc","deg_cyt")])
  
  # Obtain an MCMC sample of desired size, starting at the ML estimate
  res = doMCMC_2d(timepoints=timepoints,
               transformed_nuc=transformed_nuc,
               transformed_cyt=transformed_cyt,
               total_nuc=total_nuc,
               total_cyt=total_cyt,  
               distinct_samples = distinct_samples,
               start_params = start_params)
  
  xcoord = log(2)/res$batch[,1]
  ycoord = log(2)/res$batch[,2]
  xcoord_df = data.frame(matrix(xcoord, ncol=1))
  ycoord_df = data.frame(matrix(ycoord, ncol=1))
  #xfromto = quantile(xcoord,probs=c(0.01,0.99))
  #yfromto = quantile(ycoord,probs=c(0.01,0.99))
  xfromto=c(51,59)
  yfromto=c(36,52)
  
  # Heatscatterplot of the MCMC sample
  pdf(paste(outfile_prefix, "2d.pdf", sep="_"), useDingbats = FALSE)
  heatscatter(x=xcoord,y=ycoord,
              xlim = xfromto,
              ylim = yfromto,
              xlab="nuclear half life",ylab="cytosolic half life",main=gene, colpal = "matlablike")
  abline(v=results_table[gene,"half_life_nuc"])
  abline(h=results_table[gene,"half_life_cyt"])
  points(results_table[gene,"half_life_nuc"], results_table[gene,"half_life_cyt"], pch=19)
  
  #points(results_table[gene,"half_life_nuc"],results_table[gene,"half_life_cyt"],pch=19,col="green")
  
  # Add a countour plot to the central part of the heatscatterplot
  x_grid = seq(from=xfromto[1],to=xfromto[2],length=250)
  y_grid = seq(from=yfromto[1],to=yfromto[2],length=250)
  z_grid = t(sapply(x_grid,function(x){sapply(y_grid,function(y){
    negloss(log(2)/c(x,y),timepoints,transformed_nuc,total_nuc,transformed_cyt,total_cyt)
  })}))
  contour(x_grid, y_grid, z_grid, nlevels=30,add=T,col="#606060")
  dev.off()
  
  # Density plots for nuclear and cytosolic half lifes
  p_nuc <- ggplot() + 
    ggtitle(gene) +
    geom_density(data=xcoord_df, aes(x=xcoord_df[, 1]), colour="black", fill="grey85") + 
    geom_vline(xintercept = results_table[gene,"half_life_nuc"]) +
    theme_minimal() + xlab("nuclear half life [min]") + scale_x_continuous(limits=c(51,59), breaks=c(52,54,56,58)) +
    theme(axis.title=element_text(size=12), axis.text=element_text(size=12),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  p_cyt <- ggplot(data=ycoord_df) + 
    ggtitle(gene) +
    geom_density(data=xcoord_df, aes(x=ycoord_df[, 1]), colour="black", fill="grey85") + 
    geom_vline(xintercept = results_table[gene,"half_life_cyt"]) +
    theme_minimal() + xlab("cytosolic half life [min]") + scale_x_continuous(limits=c(36,52), breaks=c(40,45,50)) +
    theme(axis.title=element_text(size=12), axis.text=element_text(size=12),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  pdf(paste(outfile_prefix, "1d_nuc.pdf", sep="_"), useDingbats = FALSE)
  print(p_nuc)
  dev.off()
  pdf(paste(outfile_prefix, "1d_cyt.pdf", sep="_"), useDingbats = FALSE)
  print(p_cyt)
  dev.off()
  
  return(list(p_nuc, p_cyt))
}


############
### MAIN ###
############

# load the data
workdir = "~"
setwd(workdir)

filenames = c("rawdata_s1_nonUTRreads_clean.RData","rawdata_s2_nonUTRreads_clean.RData")
#cat("NB: The data file <",filename,"> was not loaded.\n",sep="")
load(filenames[1])
load(filenames[2])

genes1 = rownames(rawdata_s1_nonUTRreads_clean$tot_cy)
cat("3'UTRs measured in experiment 1 : ",length(genes1),"\n")
genes2 = rownames(rawdata_s2_nonUTRreads_clean$tot_cy)
cat("3'UTRs measured in experiment 2 : ",length(genes2),"\n")
common_genes = intersect(genes1,genes2)
cat("3'UTRs measured in both experiments : ",length(common_genes),"\n")

genes1 = setdiff(genes1, common_genes)
genes2 = setdiff(genes2, common_genes)

# define time points at which the samples were taken
measurement_timepoints = c(15, 30, 45, 60, 90, 120, 180)
# calculate the effective time points, due to a delay in 4sU labeling
delay1 = 10.27
delay2 = 10.91
timepoints1 = measurement_timepoints - delay1
timepoints2 = measurement_timepoints - delay2

# preparing parallelization
n_cores = 1
parallel_clusters <- makePSOCKcluster(n_cores, outfile="")
registerDoParallel(parallel_clusters)

# fit 3'UTRs in both experiments
results_both <- foreach(j = 1:length(common_genes), .packages=c("mcmc", "LSD", "DEoptim"),
                        .export = c("fitGene"), .verbose=TRUE) %dopar% {
                          gene = common_genes[j]
                          if (j %% 10 == 0) cat(j," , ")
                          # fits
                          r1 <- fitGene(gene,timepoints1,rawdata_s1,distinct_samples=10^4,include_samples=T)
                          r2 <- fitGene(gene,timepoints2,rawdata_s2,distinct_samples=10^4,include_samples=T)
                          print
                          # cumulative function for MCMC samples
                          cumf1_nuc <- ecdf(r1$samples_nuc)
                          cumf1_cyt <- ecdf(r1$samples_cyt)
                          cumf2_nuc <- ecdf(r2$samples_nuc)
                          cumf2_cyt <- ecdf(r2$samples_cyt)
                          r1$samples_nuc <- NULL
                          r1$samples_cyt <- NULL
                          r2$samples_nuc <- NULL
                          r2$samples_cyt <- NULL
                          # within-sample quantiles
                          within_quantiles1 <- c(cumf1_nuc(r1$deg_nuc),cumf1_cyt(r1$deg_cyt))
                          within_quantiles2 <- c(cumf2_nuc(r2$deg_nuc),cumf2_cyt(r2$deg_cyt))
                          names(within_quantiles1) <- c("nuc_quant", "cyt_quant")
                          names(within_quantiles2) <- c("nuc_quant", "cyt_quant")
                          # between-sample quantiles
                          between_quantiles1 <- c(cumf1_nuc(r2$deg_nuc),cumf1_cyt(r2$deg_cyt))
                          between_quantiles2 <- c(cumf2_nuc(r1$deg_nuc),cumf2_cyt(r1$deg_cyt))
                          names(between_quantiles1) <- c("2in1_nuc_quant", "2in1_cyt_quant")
                          names(between_quantiles2) <- c("1in2_nuc_quant", "1in2_cyt_quant")
                          # returning
                          r1$within_quant = within_quantiles1
                          r1$between_quant = between_quantiles1
                          r2$within_quant = within_quantiles2
                          r2$between_quant = between_quantiles2
                          return(list(r1, r2))
                        }
names(results_both) = common_genes
cat("\n")

# fit 3'UTRs from first experiment
results1 <- foreach(j = 1:length(genes1), .packages=c("mcmc", "LSD", "DEoptim"),
                    .export = c("fitGene"), .verbose=TRUE) %dopar% {
                      gene = genes1[j]
                      if (j %% 10 == 0) cat(j," , ")
                      r1 <- fitGene(gene,timepoints1,rawdata_s1,distinct_samples=10^4)
                      r1$within_quant = c("nuc_quant"=NA, "cyt_quant"=NA)
                      r1$between_quant = c("2in1_nuc_quant"=NA, "2in1_cyt_quant"=NA)
                      return(r1)
                    }
names(results1) = genes1
cat("\n")

# fit 3'UTRs from second experiment
results2 <- foreach(j = 1:length(genes2), .packages=c("mcmc", "LSD", "DEoptim"),
                    .export = c("fitGene", "fitTimeseries"), .verbose=TRUE) %dopar% {
                      gene = genes2[j]
                      if (j %% 10 == 0) cat(j," , ")
                      r2 <- fitGene(gene,timepoints2,rawdata_s2,distinct_samples=10^4)
                      r2$within_quant = c("nuc_quant"=NA, "cyt_quant"=NA)
                      r2$between_quant = c("1in2_nuc_quant"=NA, "1in2_cyt_quant"=NA)
                      return(r2)
                    }
names(results2) = genes2
cat("\n")

# giving free parallelization cores
stopCluster(parallel_clusters) 

# Convert the results_both list into two results tables
results1_table = NULL
results2_table = NULL
for (j in 1:length(results_both)){
  entry_both = results_both[[j]]
  entry_1 = with(entry_both[[1]],c(halflife_nuc, halflife_cyt,
                                   deg_nuc, deg_cyt,
                                   quantiles_halflife_nuc, quantiles_halflife_cyt,
                                   quantiles_deg_nuc, quantiles_deg_cyt,
                                   within_quant, between_quant,
                                   mean_coverage_nuc, mean_coverage_cyt,
                                   cv_nuc,cv_cyt, rsquared_nuc, rsquared_cyt))
  entry_2 = with(entry_both[[2]],c(halflife_nuc, halflife_cyt,
                                   deg_nuc, deg_cyt,
                                   quantiles_halflife_nuc, quantiles_halflife_cyt,
                                   quantiles_deg_nuc, quantiles_deg_cyt,
                                   within_quant, between_quant,
                                   mean_coverage_nuc, mean_coverage_cyt,
                                   cv_nuc,cv_cyt, rsquared_nuc, rsquared_cyt))
  results1_table = rbind(results1_table,entry_1)
  results2_table = rbind(results2_table,entry_2)
}

# Convert the results1 list into a results table
for (j in 1:length(results1)){
  entry = with(results1[[j]],c(halflife_nuc, halflife_cyt,
                               deg_nuc, deg_cyt,
                               quantiles_halflife_nuc, quantiles_halflife_cyt,
                               quantiles_deg_nuc, quantiles_deg_cyt,
                               within_quant, between_quant,
                               mean_coverage_nuc, mean_coverage_cyt,
                               cv_nuc, cv_cyt, rsquared_nuc, rsquared_cyt))
  results1_table = rbind(results1_table,entry)
}
rownames(results1_table) = c(common_genes, genes1)

# Convert the results2 list into a results table
for (j in 1:length(results2)){
  entry = with(results2[[j]],c(halflife_nuc, halflife_cyt,
                               deg_nuc, deg_cyt,
                               quantiles_halflife_nuc, quantiles_halflife_cyt,
                               quantiles_deg_nuc, quantiles_deg_cyt,
                               within_quant, between_quant,
                               mean_coverage_nuc, mean_coverage_cyt,
                               cv_nuc, cv_cyt, rsquared_nuc, rsquared_cyt))
  results2_table = rbind(results2_table,entry)
}
rownames(results2_table) = c(common_genes, genes2)

# Save the results
save(results1,results2,results_both,results1_table,results2_table,file = "Estimation_results.RData")

cat("Done.")
 
# gene = common_genes[1]
# plotGene(gene,timepoints1,rawdata_s1_nonUTRreads_clean,results1_table)
# plotMCMC(gene,timepoints,rawdata_s1_nonUTRreads_clean,results1_table)

# plot(results1_table[,"mean_coverage_nuc"],results2_table[,"mean_coverage_nuc"],main="mean coverage",pch=19,cex=0.5)
# plot(results1_table[,"ML_halflife_nuc"],results2_table[,"ML_halflife_nuc"],pch=19,cex=0.5,
#      xlim=c(0,400), ylim = c(0,400), main="Nuclear half life")
# plot(results1_table[,"ML_halflife_cyt"],results2_table[,"ML_halflife_cyt"],pch=19,cex=0.5,
#      xlim=c(0,400), ylim = c(0,400), main="Cytosolic half life")
# 
# for (gene in common_genes){
#   plotGene(gene,timepoints,rawdata_s2_nonUTRreads_clean,results2_table)
# }
