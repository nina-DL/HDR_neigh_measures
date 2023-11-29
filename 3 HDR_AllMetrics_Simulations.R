# Load Data, Functions, Libraries -----------------------------
library(parallel)

pkgs = c("rgl", "MASS", "VC2copula", "ks", "kdecopula", "kde1d", "VineCopula", "pbapply", "mclust", "KScorrect", "data.table", 
         "gsl", # for Debye function
         "gtools", #for Dirichlet
         "mltools" # for empiricalcdf (M7)
         )

sapply(pkgs, require, character.only = TRUE)

# load simulated scenarios and the true "estimated" HDR on a big ss
load("Sim Results/Simulated Scenarios/Def_scenarios.RData")

# load all functions called in this script
source("2 HDR_functions.R")

# |-- Main function to get performance errors for all scenarios -------------

get_ALL_err_light = function(mvdc = mvd_def_sD3, true_level = true_q_sD3, Dirichlet_alpha = NA, 
                        Nsample = 1e3, coverage_prob = 0.95, plot = T){
  # This function evaluates classification errors of the different measures, denoted by M0-M8.
  # It takes as main input the definition of a scenario of interest and its true alpha-quantile. Specifically:
  # 1. mvdc:: an object of class "mvdc", defining the copula model and the marginals as well as their parameters
  # 2. true_level:: the alpha-quantile of the density of interest following the density quantile approach of Hyndman.
  #                 by default we use alpha = 0.05, targeting a coverage probability of 95%.
  # In case of a scenario other than an mvdc object, e.g., the Dirichlet case, we specify it in a new argument:
  # 3. Dirichlet_alpha:: a vector characterising the parameters of a Dirichlet distrribution. 
  #                      by default this is set to NA.
  # The following additional arguments can be set:
  # 4. Nsample:: sample size, that is, the number of data points to build the HDR on. 
  # 5. coverage_prob:: the desired coverage probability for the HDR. By default, this is set to 0.95.
  # 6. plot:: a boolean object specifying whether to plot the estimated HDR, that is the points classification.
  #           By default this is set to T.

  
  if(any(is.na(Dirichlet_alpha)==T)){
    # Simulate the data from an mvdc object and get the true density
    data_sim = rMvdc(n=Nsample, mvdc = mvdc)
    fxy_true = dMvdc(x = data_sim, mvdc = mvdc)
    margin_family = mvdc@margins
    #main_title = paste("Scenario:", paste(mvdc@margins,collapse=", "), "margins &", class(mvdc@copula)[1])
  }else{
    data_sim = rdirichlet(Nsample, alpha = Dirichlet_alpha)
    fxy_true = ddirichlet(x = data_sim, alpha = Dirichlet_alpha)
    data_sim = data_sim[,1:(length(Dirichlet_alpha)-1)]
    margin_family = c("beta", "beta")
    #main_title = paste0("Scenario: Dirichlet (", paste(Dirichlet_alpha,collapse=", "), ")")
  }
  
  # M0: true density ---> OK
  res_M0 = HDR_class(data = data_sim, true_density = fxy_true, est_density = fxy_true, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  # M1: KDE ---> OK
  # Default: myh = ks::Hpi.diag(data_S2,binned=TRUE) 
  # OPTIMAL bandwidth: Chacón, J.E., Duong, T., Wand, M.P. (2011), Asymptotics for General Multivariate Kernel Density Derivative Estimators, Statistica Sinica, 21, 807–840.
  myh_opt = (4/(dim(data_sim)[2]+4))^(1/(dim(data_sim)[2]+6))*dim(data_sim)[1]^(-1/(dim(data_sim)[2]+6))*abs(cov(data_sim))^(1/2)
  fit_kde = kde(data_sim, eval.points = data_sim, H = myh_opt) # from ks package
  res_M1 = HDR_class(data = data_sim, true_density = fxy_true, est_density = fit_kde$estimate, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  # M2: nonparametric copula ---> OK
  npC_est = NP_cop(data_sim)
  res_M2 = HDR_class(data = data_sim, true_density = fxy_true, est_density = npC_est$dist_mv, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  # M3: parametric copula ---> OK
  parametric_fit = Pfit_margins_cop(data = data_sim, margin_family = margin_family)
  pC_est = P_cop(data = data_sim, margin_family = parametric_fit$margin_family, 
                 margin_params = parametric_fit$margin_params, copula_model = parametric_fit$copula_model) 
  res_M3 = HDR_class(data = data_sim, true_density = fxy_true, est_density = pC_est, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)

  # M4: k-neighborhood distances ---> OK
  # Optimal k: sqrt(myn/2). Rule of thumb based on simulations: sqrt(myn*tau) / sqrt(myn*rho) / sqrt(myn/2)
  myk = round(sqrt(Nsample*.5)) # Rule of thumb based on simulations: sqrt(Nsample*tau) / sqrt(Nsample*rho) 
  dist_M4 = mapply(knn_dist, data_sim[,1], data_sim[,2], MoreArgs = list(k = myk, data = data_sim))
  # pay attention to the measure being a concentration/sparsity measure 
  res_M4 = HDR_class(data = data_sim, true_density = fxy_true, est_density = 1/dist_M4, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  # M5: Weighted sum of CDF distances from the k nearest neighbors (UV CDF estimation) ---> OK
  # Variation of https://core.ac.uk/download/pdf/30276753.pdf
  # Not a Metric (triangular inequality satisfied only in specif ordering between three points)
  myk = 30 # quite robust to k 
  est_cdfx = ecdf(data_sim[,1])
  est_cdfy = ecdf(data_sim[,2])
  dist_M5 = mapply(knn_CDFdist, data_sim[,1], data_sim[,2], MoreArgs = list(k = myk, data = data_sim, cdf_x=est_cdfx, cdf_y=est_cdfy))
  res_M5 = HDR_class(data = data_sim, true_density = fxy_true, est_density = 1/dist_M5, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  # # M: CDF distance in a d-neighborhood region (~PDF in a d-region) with UNI ECDFs 
  # myd = c(sum(abs(range(data_sim[,1])))/sqrt(Nsample), sum(abs(range(data_sim[,2])))/sqrt(Nsample))
  # dist_M = pbmapply(d_CDFdist_uv, data_sim[,1], data_sim[,2], MoreArgs = list(d = myd, data = data_sim, cdf_x = est_cdfx, cdf_y = est_cdfy))
  # res_M = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M, true_level = true_level, 
  #                  coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[6], ")"))
  
  # M6: CDF distance in a d-neighborhood region (~PDF in a d-region) with Mv ECDF ---> OK
  # Optimal d?
  myd = ifelse(any(is.na(Dirichlet_alpha)==T), exp(2.13-0.3*log(Nsample)), 0.10)
  dist_M6 = mapply(d_CDFdist, data_sim[,1], data_sim[,2], MoreArgs = list(data = data_sim, d = rep(myd, 2)))
  res_M6 = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M6, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  # M8: CDF distance in a d-neighborhood region (~PDF in a d-region) with NP Copula ---> OK
  # from kdecopula package: https://cran.r-project.org/web/packages/kdecopula/vignettes/kdecopula.pdf
  # here I use the TLL2nn transformation method (see 2.3. of the above and Fig 7 for a comparison among different methods)
  # Optimal d? ---> fitted nonlinear regression with exponential decay (otherwise myd = 1.5)
  myd = ifelse(any(is.na(Dirichlet_alpha)==T), exp(1.74-0.26*log(Nsample)), exp(-1.22-0.23*log(Nsample)))
  dist_M7 = mapply(d_NPCopCDFdist, data_sim[,1], data_sim[,2], 
                      MoreArgs = list(data = data_sim, d = rep(myd, 2), cdf_x = est_cdfx, cdf_y = est_cdfy, fit_cop = npC_est$fit_cop))
  res_M7 = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M7, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  # M9: CDF distance in a d-neighborhood region with Parametric Copula ---> OK
  # Optimal d? ---> fitted nonlinear regression with exponential decay (otherwise myd = 0.7). For Dirichlet, take small values < 0.05
  myd = ifelse(any(is.na(Dirichlet_alpha)==T), exp(1.60 - 0.41*log(Nsample)), 0.02)
  dist_M8 = mapply(d_PCopCDFdist, data_sim[,1], data_sim[,2], 
                      MoreArgs = list(d = rep(myd, 2), margin_family = parametric_fit$margin_family, 
                                      margin_params = parametric_fit$margin_params, copula_model = parametric_fit$copula_model))
  res_M8 = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M8, true_level = true_level, 
                   coverage_prob = coverage_prob, build_plot = F)
  
  
  #final error matrix
  Error_mat = cbind(c(res_M0$ERR, res_M0$FPR, res_M0$FNR, res_M0$Accuracy, res_M0$F1_score, res_M0$MCC),
                    c(res_M1$ERR, res_M1$FPR, res_M1$FNR, res_M1$Accuracy, res_M1$F1_score, res_M1$MCC),
                    c(res_M2$ERR, res_M2$FPR, res_M2$FNR, res_M2$Accuracy, res_M2$F1_score, res_M2$MCC),
                    c(res_M3$ERR, res_M3$FPR, res_M3$FNR, res_M3$Accuracy, res_M3$F1_score, res_M3$MCC),
                    c(res_M4$ERR, res_M4$FPR, res_M4$FNR, res_M4$Accuracy, res_M4$F1_score, res_M4$MCC),
                    c(res_M5$ERR, res_M5$FPR, res_M5$FNR, res_M5$Accuracy, res_M5$F1_score, res_M5$MCC),
                    #c(res_M$ERR, res_M$FPR, res_M$FNR, res_M$Accuracy, res_M$F1_score, res_M$MCC),
                    c(res_M6$ERR, res_M6$FPR, res_M6$FNR, res_M6$Accuracy, res_M6$F1_score, res_M6$MCC),
                    c(res_M7$ERR, res_M7$FPR, res_M7$FNR, res_M7$Accuracy, res_M7$F1_score, res_M7$MCC),
                    c(res_M8$ERR, res_M8$FPR, res_M8$FNR, res_M8$Accuracy, res_M8$F1_score, res_M8$MCC))

  return(Error_mat)
}


# try function
set.seed(123)
get_ALL_err(mvdc = mvd_def_sA1t, true_level = true_q_sA1t, Dirichlet_alpha = NA, Nsample = 1e2, coverage_prob = 0.95, plot = T)
set.seed(123)
get_ALL_err_light(mvdc = mvd_def_sA1t, true_level = true_q_sA1t, Dirichlet_alpha = NA, Nsample = 1e2, coverage_prob = 0.95, plot = T)

set.seed(123)
get_ALL_err(mvdc = NA, Dirichlet_alpha = c(1,1,2), true_level = true_q_sDir_112, Nsample = 1e3, coverage_prob = 0.95, plot = T)


# Hyper-setup ----------------------
Nsim = 1000
Nsamples = c(50, 100, 250, 500, 1000)

measures = c("M0: Known density","M1: KDE", "M2: KDE-NPCop", "M3: KDE-PCop", 
             "M4: kNN-Euclidean", "M5: kNN-CDFuv", #"M6:eps-CDF-uv", 
             "M6: Ɛ-CDFmv", "M7: Ɛ-CDF-NPCop","M8: Ɛ-CDF-PCop")
metrics = c("Total Error Rate", "FPR", "FNR", "Accuracy", "F1 Score", "MCC")

n_metrics = length(metrics); n_methods = length(measures)

# |-- Get average results for sC1t (extensive simulations; results are saved) ---------------
Res_sD1t = list()
ErrMean = ErrSD = list()
set.seed(123)
for(i in 1:length(Nsamples)){
  #Res_sD1t[[i]] = pbreplicate(Nsim, get_FULL_res(mvdc = mvd_def_sD1t, true_level = true_q_sD1t, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F), simplify=FALSE)
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, varlist = ls())
  Res_sD1t[[i]] <- parLapply(cl, 1:Nsim, \(x) {
    sapply(pkgs, require, character.only = TRUE)
    set.seed(x)
    tryCatch({ get_ALL_err_light(mvdc = mvd_def_sD1t, true_level = true_q_sD1t, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F)}, 
             error = function(e) { matrix(NA, ncol = n_methods, nrow = n_metrics)})
  })
  stopCluster(cl)
  
  #Res_sD1t[[i]] = Filter(function(a) any(!is.na(a)), Res_sD1t[[i]])
  ErrMean[[i]] = apply(array(unlist(Res_sD1t[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), mean, na.rm = T)
  ErrSD[[i]] = apply(array(unlist(Res_sD1t[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), sd, na.rm = T)
  colnames(ErrMean[[i]]) = colnames(ErrSD[[i]]) = measures
  rownames(ErrMean[[i]]) = rownames(ErrSD[[i]]) = metrics
  myf = Res_sD1t[[i]]
  file = paste0("Res_sD1t_n",Nsamples[i],".RData")
  message(paste("MC Completed for SS:", Nsamples[i]))
  save(myf, file = file)
  }

names(ErrMean) = names(ErrSD) = names(Res_sD1t) = Nsamples
save(Res_sD1t, ErrMean, ErrSD, file = "Res_sD1t.RData")


# |-- Get average results for sC1 (extensive simulations; results are saved) ---------------
Res_sD1 = list()
ErrMean = ErrSD = list()
set.seed(123)
for(i in 1:length(Nsamples)){
  #Res_sD1[[i]] = pbreplicate(Nsim, get_FULL_res(mvdc = mvd_def_sD1, true_level = true_q_sD1, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F), simplify=FALSE)
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, varlist = ls())
  Res_sD1[[i]] <- parLapply(cl, 1:Nsim, \(x) {
    sapply(pkgs, require, character.only = TRUE)
    set.seed(x)
    tryCatch({ get_ALL_err_light(mvdc = mvd_def_sD1, true_level = true_q_sD1, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F)}, 
             error = function(e) { matrix(NA, ncol = n_methods, nrow = n_metrics)})
  })
  stopCluster(cl)
  
  #Res_sD1[[i]] = Filter(function(a) any(!is.na(a)), Res_sD1[[i]])
  ErrMean[[i]] = apply(array(unlist(Res_sD1[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), mean, na.rm = T)
  ErrSD[[i]] = apply(array(unlist(Res_sD1[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), sd, na.rm = T)
  colnames(ErrMean[[i]]) = colnames(ErrSD[[i]]) = measures
  rownames(ErrMean[[i]]) = rownames(ErrSD[[i]]) = metrics
  myf = Res_sD1[[i]]
  file = paste0("Res_sD1_n",Nsamples[i],".RData")
  message(paste("MC Completed for SS:", Nsamples[i]))
  save(myf, file = file)
}

names(ErrMean) = names(ErrSD) = names(Res_sD1) = Nsamples
save(Res_sD1, ErrMean, ErrSD, file = "Res_sD1.RData")

# |-- Get average results for sC2 (extensive simulations; results are saved) ---------------
Res_sD2 = list()
ErrMean = ErrSD = list()
set.seed(123)
for(i in 1:length(Nsamples)){
  #Res_sD2[[i]] = pbreplicate(Nsim, get_FULL_res(mvdc = mvd_def_sD2, true_level = true_q_sD2, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F), simplify=FALSE)
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, varlist = ls())
  Res_sD2[[i]] <- parLapply(cl, 1:Nsim, \(x) {
    sapply(pkgs, require, character.only = TRUE)
    set.seed(x)
    tryCatch({ get_ALL_err_light(mvdc = mvd_def_sD2, true_level = true_q_sD2, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F)}, 
             error = function(e) { matrix(NA, ncol = n_methods, nrow = n_metrics)})
  })
  stopCluster(cl)
  
  #Res_sD2[[i]] = Filter(function(a) any(!is.na(a)), Res_sD2[[i]])
  ErrMean[[i]] = apply(array(unlist(Res_sD2[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), mean, na.rm = T)
  ErrSD[[i]] = apply(array(unlist(Res_sD2[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), sd, na.rm = T)
  colnames(ErrMean[[i]]) = colnames(ErrSD[[i]]) = measures
  rownames(ErrMean[[i]]) = rownames(ErrSD[[i]]) = metrics
  myf = Res_sD2[[i]]
  file = paste0("Res_sD2_n",Nsamples[i],".RData")
  message(paste("MC Completed for SS:", Nsamples[i]))
  save(myf, file = file)
}

names(ErrMean) = names(ErrSD) = names(Res_sD2) = Nsamples
save(Res_sD2, ErrMean, ErrSD, file = "Res_sD2.RData")

# |-- Get average results for sC3 (extensive simulations; results are saved) ---------------
Res_sD3 = list()
ErrMean = ErrSD = list()
set.seed(123)
for(i in 1:length(Nsamples)){
  #Res_sD3[[i]] = pbreplicate(Nsim, get_FULL_res(mvdc = mvd_def_sD3, true_level = true_q_sD3, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F), simplify=FALSE)
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, varlist = ls())
  Res_sD3[[i]] <- parLapply(cl, 1:Nsim, \(x) {
    sapply(pkgs, require, character.only = TRUE)
    set.seed(x)
    tryCatch({ get_ALL_err_light(mvdc = mvd_def_sD3, true_level = true_q_sD3, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F)}, 
             error = function(e) { matrix(NA, ncol = n_methods, nrow = n_metrics)})
  })
  stopCluster(cl)
  
  #Res_sD3[[i]] = Filter(function(a) any(!is.na(a)), Res_sD3[[i]])
  ErrMean[[i]] = apply(array(unlist(Res_sD3[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), mean, na.rm = T)
  ErrSD[[i]] = apply(array(unlist(Res_sD3[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), sd, na.rm = T)
  colnames(ErrMean[[i]]) = colnames(ErrSD[[i]]) = measures
  rownames(ErrMean[[i]]) = rownames(ErrSD[[i]]) = metrics
  myf = Res_sD3[[i]]
  file = paste0("Res_sD3_n",Nsamples[i],".RData")
  message(paste("MC Completed for SS:", Nsamples[i]))
  save(myf, file = file)
}

names(ErrMean) = names(ErrSD) = names(Res_sD3) = Nsamples
save(Res_sD3, ErrMean, ErrSD, file = "Res_sD3.RData")


# |-- Get average results for sDir (extensive simulations; results are saved) ---------------
Res_Dir112 = list()
# ERR = ERR_SD = FPR = FPR_SD = FNR = FNR_SD = 
#   Accuracy = Accuracy_SD = F1 = F1_SD = MCC = MCC_SD = matrix(NA, nrow = length(Nsamples), ncol = n_methods)
ErrMean = ErrSD = list()
set.seed(12345)
for(i in 1:length(Nsamples)){
  Res_Dir112[[i]] = pbreplicate(Nsim, get_FULL_res(mvdc = NA, Dirichlet_alpha = c(1,1,2), true_level = true_q_sDir_112, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F), simplify=FALSE)
  ErrMean[[i]] = apply(array(unlist(Res_Dir112[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), mean, na.rm = T)
  ErrSD[[i]] = apply(array(unlist(Res_Dir112[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), sd, na.rm = T)
  colnames(ErrMean[[i]]) = colnames(ErrSD[[i]]) = measures
  rownames(ErrMean[[i]]) = rownames(ErrSD[[i]]) = c("Total Error Rate", "FPR", "FNR", "Accuracy", "F1 Score", "MCC")
  message(paste("MC Completed for SS:", Nsamples[i]))
}

names(ErrMean) = names(ErrSD) = names(Res_Dir112) = Nsamples
save(Res_Dir112, ErrMean, ErrSD, file = "Res_Dir112.RData")


# |-- Get average results for sDir (extensive simulations; results are saved) ---------------
Res_sDir = list()
ErrMean = ErrSD = list()
set.seed(123)
for(i in 1:length(Nsamples)){
  #Res_sDir[[i]] = pbreplicate(Nsim, get_FULL_res(mvdc = mvd_def_sDir, true_level = true_q_sDir, Dirichlet_alpha = NA, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F), simplify=FALSE)
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, varlist = ls())
  Res_sDir[[i]] <- parLapply(cl, 1:Nsim, \(x) {
    sapply(pkgs, require, character.only = TRUE)
    set.seed(x)
    tryCatch({ get_ALL_err_light(mvdc = NA, Dirichlet_alpha = c(1,1,2), true_level = true_q_sDir_112, Nsample = Nsamples[i], coverage_prob = 0.95, plot = F)}, 
             error = function(e) { matrix(NA, ncol = n_methods, nrow = n_metrics)})
  })
  stopCluster(cl)
  
  #Res_sDir[[i]] = Filter(function(a) any(!is.na(a)), Res_sDir[[i]])
  ErrMean[[i]] = apply(array(unlist(Res_sDir[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), mean, na.rm = T)
  ErrSD[[i]] = apply(array(unlist(Res_sDir[[i]]), dim = c(n_metrics,n_methods,Nsim)), c(1,2), sd, na.rm = T)
  colnames(ErrMean[[i]]) = colnames(ErrSD[[i]]) = measures
  rownames(ErrMean[[i]]) = rownames(ErrSD[[i]]) = metrics
  myf = Res_sDir[[i]]
  file = paste0("Res_sDir_n",Nsamples[i],".RData")
  message(paste("MC Completed for SS:", Nsamples[i]))
  save(myf, file = file)
}

names(ErrMean) = names(ErrSD) = names(Res_sDir) = Nsamples
save(Res_sDir, ErrMean, ErrSD, file = "Res_sDir.RData")

