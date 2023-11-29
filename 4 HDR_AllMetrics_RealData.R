# Real Data: Parametric is unfeasible ----------------
pkgs = c("rgl", "MASS", "VC2copula", "ks", "kdecopula", "kde1d", "VineCopula", "pbapply", "mclust", "KScorrect", "data.table", 
         "gsl", # for Debye function
         "gtools", #for Dirichlet
         "mltools" # for empiricalcdf (M7)
)

sapply(pkgs, require, character.only = TRUE)

# load all functions called in this script
source("2 HDR_functions.R")

# Main function with NP measure evaluating all metrics with final plot --------------------

# main function that computes classification errors and returns a plot of the estimated HDR
HDR_NP_class = function(data = df_scaled, est_density = fxy_est, cutoff_level = NA, coverage_prob = 0.95, 
                        build_plot = T, main_title = "Estimated HDR"){
  N = nrow(data)
  
  hdr_mv = data.frame(Order = 1:N, X = data[,1], Y = data[,2], fxy = est_density)
  hdr_mv = hdr_mv[order(hdr_mv$Order),]
  
  # Quantile Method
  if(is.na(cutoff_level)==T){
    cutoff_level = quantile(est_density, probs = (1-coverage_prob), na.rm = T)
  }
  
  hdr_mv$ok = as.integer(hdr_mv$fxy>cutoff_level)
  hdr_mv$col = sapply(hdr_mv$ok, FUN = my_col, col="darkmagenta")
  hdr_mv$pch = sapply(hdr_mv$ok, FUN = my_pch, pch=4)
  
  if(build_plot == T){
    HDR_plot = plot(hdr_mv$X, hdr_mv$Y, xlim=c(min(hdr_mv$X), max(hdr_mv$X)), 
                    ylim=c(min(hdr_mv$Y), max(hdr_mv$Y)), 
                    xlab = expression(x[1]), ylab = expression(x[2]), cex=0.5, pch = hdr_mv$pch, col=hdr_mv$col,
                    main = main_title, bty="n")
  } else {HDR_plot = NA}
  
  return(list(#df = hdr_mv, 
    N = N, HDR_plot = HDR_plot, in_out = hdr_mv$ok))
}

get_NPm_err = function(data = df_scaled, coverage_prob = 0.95, plot = T){
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
  
  Nsample = nrow(data)
    
  # This function evaluates the following 9 measures: M0-M8
  measures = c(expression(M[0]*":KDE"), expression(M[1]*":kNN-Eucl"),
               expression(M[2]*":kNN-CDF"), expression(M[3]*":"*epsilon*"-CDF"),
               expression(M[0]^NPCop*":DE"), expression(M[3]^NPCop*":"*epsilon*"-CDF"))
  
  if(plot == T){
    
    op <- par(mfrow = c(2,3),
              oma = c(2,2,2,1) + 0.1,
              mar = c(2,2,2,1) + 0.1,
              mai = c(0.4, 0.1, 0.4, 0.1),
              bty = 'n')
    
    main_title = "95% HDR on MAGIC Data"
  }
  
  data_sim = data
  # M1: KDE ---> OK
  # Default: myh = ks::Hpi.diag(data_S2,binned=TRUE) 
  # OPTIMAL bandwidth: Chacón, J.E., Duong, T., Wand, M.P. (2011), Asymptotics for General Multivariate Kernel Density Derivative Estimators, Statistica Sinica, 21, 807–840.
  myh_opt = (4/(dim(data_sim)[2]+4))^(1/(dim(data_sim)[2]+6))*dim(data_sim)[1]^(-1/(dim(data_sim)[2]+6))*abs(cov(data_sim))^(1/2)
  fit_kde = kde(data_sim, eval.points = data_sim, H = myh_opt) # from ks package
  res_M1 = HDR_NP_class(data = data_sim, est_density = fit_kde$estimate, coverage_prob = coverage_prob, 
                        build_plot = plot, main_title = measures[1])
  
  # M4: k-neighborhood distances ---> OK
  # Optimal k: sqrt(myn/2). Rule of thumb based on simulations: sqrt(myn*tau) / sqrt(myn*rho) / sqrt(myn/2)
  myk = round(sqrt(Nsample*.5)) # Rule of thumb based on simulations: sqrt(Nsample*tau) / sqrt(Nsample*rho) 
  dist_M4 = pbmapply(knn_dist, data_sim[,1], data_sim[,2], MoreArgs = list(k = myk, data = data_sim))
  # pay attention to the measure being a concentration/sparsity measure 
  res_M4 = HDR_NP_class(data = data_sim, est_density = 1/dist_M4, coverage_prob = coverage_prob, 
                        build_plot = plot, main_title = measures[2])
  
  # M5: Weighted sum of CDF distances from the k nearest neighbors (UV CDF estimation) ---> OK
  # Variation of https://core.ac.uk/download/pdf/30276753.pdf
  # Not a Metric (triangular inequality satisfied only in specif ordering between three points)
  myk = 30 # quite robust to k 
  est_cdfx = ecdf(data_sim[,1])
  est_cdfy = ecdf(data_sim[,2])
  dist_M5 = pbmapply(knn_CDFdist, data_sim[,1], data_sim[,2], MoreArgs = list(k = myk, data = data_sim, cdf_x=est_cdfx, cdf_y=est_cdfy))
  res_M5 = HDR_NP_class(data = data_sim, est_density = 1/dist_M5, coverage_prob = coverage_prob, 
                        build_plot = plot, main_title = measures[3])
  
  # # M: CDF distance in a d-neighborhood region (~PDF in a d-region) with UNI ECDFs 
  # myd = c(sum(abs(range(data_sim[,1])))/sqrt(Nsample), sum(abs(range(data_sim[,2])))/sqrt(Nsample))
  # dist_M = pbmapply(d_CDFdist_uv, data_sim[,1], data_sim[,2], MoreArgs = list(d = myd, data = data_sim, cdf_x = est_cdfx, cdf_y = est_cdfy))
  # res_M = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M, true_level = true_level, 
  #                  coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[6], ")"))
  
  # M6: CDF distance in a d-neighborhood region (~PDF in a d-region) with Mv ECDF ---> OK
  # Optimal d?
  myd = exp(2.13-0.3*log(Nsample))
  dist_M6 = pbmapply(d_CDFdist, data_sim[,1], data_sim[,2], MoreArgs = list(data = data_sim, d = rep(myd, 2)))
  res_M6 = HDR_NP_class(data = data_sim, est_density = dist_M6, coverage_prob = coverage_prob, 
                        build_plot = plot, main_title = measures[4])
  
  
  # M2: nonparametric copula ---> OK
  npC_est = NP_cop(data_sim)
  res_M2 = HDR_NP_class(data = data_sim, est_density = npC_est$dist_mv, coverage_prob = coverage_prob, 
                        build_plot = plot, main_title = measures[5])
  
  # M8: CDF distance in a d-neighborhood region (~PDF in a d-region) with NP Copula ---> OK
  # from kdecopula package: https://cran.r-project.org/web/packages/kdecopula/vignettes/kdecopula.pdf
  # here I use the TLL2nn transformation method (see 2.3. of the above and Fig 7 for a comparison among different methods)
  # Optimal d? ---> fitted nonlinear regression with exponential decay (otherwise myd = 1.5)
  myd = exp(1.74-0.26*log(Nsample))
  dist_M7 = pbmapply(d_NPCopCDFdist, data_sim[,1], data_sim[,2], 
                     MoreArgs = list(data = data_sim, d = rep(myd, 2), cdf_x = est_cdfx, cdf_y = est_cdfy, fit_cop = npC_est$fit_cop))
  res_M7 = HDR_NP_class(data = data_sim, est_density = dist_M7, coverage_prob = coverage_prob, 
                     build_plot = plot, main_title = measures[6])
  
  if(plot == T){title(main_title, line = 0.5, cex.main = 1.5, outer = TRUE)}
  
  #return(NA)
}


# |-- Load data -------------

# https://archive.ics.uci.edu/ml/datasets/MAGIC+Gamma+Telescope
# https://arxiv.org/pdf/1701.00845.pdf
# 1.  fLength:  continuous  # major axis of ellipse [mm]
# 2.  fWidth:   continuous  # minor axis of ellipse [mm] 
# 3.  fSize:    continuous  # 10-log of sum of content of all pixels [in #phot]
# 4.  fConc:    continuous  # ratio of sum of two highest pixels over fSize  [ratio]
# 5.  fConc1:   continuous  # ratio of highest pixel over fSize  [ratio]
# 6.  fAsym:    continuous  # distance from highest pixel to center, projected onto major axis [mm]
# 7.  fM3Long:  continuous  # 3rd root of third moment along major axis  [mm] 
# 8.  fM3Trans: continuous  # 3rd root of third moment along minor axis  [mm]
# 9.  fAlpha:   continuous  # angle of major axis with vector to origin [deg]
# 10.  fDist:    continuous  # distance from origin to center of ellipse [mm]
# 11.  class:    g,h         # gamma (signal), hadron (background)

MAGIC = read.table("Real data/magic04.data", sep = ",")
gMAGIC = MAGIC[MAGIC$V11=="g",]

plot(gMAGIC[,c(5,7)])
plot(gMAGIC[,c(5,8)])
plot(gMAGIC[,c(7,8)])

df = gMAGIC[,c(5,7)]
df_scaled = apply(df, 2, scale)
plot(df_scaled)
myn = dim(df)[1]

get_NPm_err(data = df_scaled, coverage_prob = 0.95, plot = T)

# Notes:
# All sensitive to scale, except KDE_Cop



# MODEL AVERAGING

res_MA = apply(cbind(res_M1$in_out, res_M4$in_out, res_M5$in_out, res_M6$in_out, res_M2$in_out, res_M7$in_out),1,mean, na.rm = T)

N = nrow(data)

hdr_mv = data.frame(Order = 1:N, X = data[,1], Y = data[,2])
hdr_mv = hdr_mv[order(hdr_mv$Order),]

hdr_mv$ok = as.integer(res_MA>.5)
table(hdr_mv$ok)/N

hdr_mv$col = sapply(hdr_mv$ok, FUN = my_col, col="darkmagenta")
hdr_mv$pch = sapply(hdr_mv$ok, FUN = my_pch, pch=4)

dev.off()
plot(hdr_mv$X, hdr_mv$Y, xlim=c(min(hdr_mv$X), max(hdr_mv$X)),
     ylim=c(min(hdr_mv$Y), max(hdr_mv$Y)),
     xlab = expression(x[1]), ylab = expression(x[2]), cex=0.5, pch = hdr_mv$pch, col=hdr_mv$col,
     main = "95% HDR on MAGIC Data (Measure Averaging)", bty="n")


sum(hdr_mv$ok != res_M1$in_out)/N*100
sum(hdr_mv$ok != res_M4$in_out)/N*100
sum(hdr_mv$ok != res_M5$in_out)/N*100
sum(hdr_mv$ok != res_M6$in_out)/N*100
sum(hdr_mv$ok != res_M2$in_out)/N*100
sum(hdr_mv$ok != res_M7$in_out)/N*100
