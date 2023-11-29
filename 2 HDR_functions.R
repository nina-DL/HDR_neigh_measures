# Main functions ---> OK --------------------

# function to color highest-density vs non highest-density points
my_col = function(x, col) ifelse(x==0, col, "cadetblue")
my_pch = function(x, pch) ifelse(x==0, pch, 20)

# main function that computes classification errors and returns a plot of the estimated HDR
HDR_class = function(data = mvt_draws, true_density = fxy_true, est_density = fxy_cop, true_level = NA, 
                     cutoff_level = NA, coverage_prob = 0.95, build_plot = T, mat_plot = NA,
                     main_title = "Estimated HDR"){
  N = nrow(data)
  
  hdr_mv = data.frame(Order = 1:N, X = data[,1], Y = data[,2], fxy = est_density)
  hdr_mv = hdr_mv[order(hdr_mv$Order),]
  
  # Quantile Method
  if(is.na(true_level)==T){
    true_level = quantile(true_density, probs = (1-coverage_prob), na.rm = T)
  }
  
  # Quantile Method
  if(is.na(cutoff_level)==T){
    cutoff_level = quantile(est_density, probs = (1-coverage_prob), na.rm = T)
  }
  
  hdr_mv$ok_true = as.integer(true_density>true_level)
  
  hdr_mv$ok = as.integer(hdr_mv$fxy>cutoff_level)
  hdr_mv$col = sapply(hdr_mv$ok, FUN = my_col, col="red")
  
  # Classification performance measures
  FP = as.numeric(sum(hdr_mv$ok_true>hdr_mv$ok))
  FN = as.numeric(sum(hdr_mv$ok_true<hdr_mv$ok))
  TP = as.numeric(sum((hdr_mv$ok_true==hdr_mv$ok)&(hdr_mv$ok_true==0)))
  TN = as.numeric(sum((hdr_mv$ok_true==hdr_mv$ok)&(hdr_mv$ok_true==1)))
  ER = as.numeric(sum(hdr_mv$ok_true!=hdr_mv$ok))
  
  if(build_plot == T){
    ERR = ER/N
    FPR = FP/(FP+TN)
    FNR = FN/(FN+TP)
    HDR_plot = plot(hdr_mv$X, hdr_mv$Y, xlim=c(min(hdr_mv$X), max(hdr_mv$X)), 
                    ylim=c(min(hdr_mv$Y), max(hdr_mv$Y)), 
                    xlab = expression(x[1]), ylab = expression(x[2]), cex=0.5, pch = 20, col=hdr_mv$col,
                    main = main_title, bty="n",
                    sub = paste0("Total Error = ", bquote(.(round(ERR,3))), 
                                 "; FPR = ", bquote(.(round(FPR,3))),
                                 "; FNR = ", bquote(.(round(FNR,3)))))
    
    HDR_contour = contour(seq(min(data[,1])-1, max(data[,1])+1, length.out = 100), 
                          seq(min(data[,2])-1, max(data[,2])+1, length.out = 100), 
                          mat_plot, levels = true_level, add = T, col = "blue", lwd = 2)
    
  } else {HDR_plot = NA}
  
  return(list(#df = hdr_mv, 
    N = N, TP = TP, TN = TN, FP = FP, FN = FN, 
    ERR = ER/N, 
    FPR = FP/(FP+TN), 
    FNR = FN/(FN+TP),
    Accuracy = 1-ER/N,
    F1_score = 2*TP/(2*TP+FP+FN) + 2*TN/(2*TN+FP+FN),
    MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),
    HDR_plot = HDR_plot))
}

# Dirichlet copula function ---> OK 
# the density of the copula is given by 
d_cop_u = function(u, alpha_vec){
  # Density function of the copula of Dirichlet of dimension n (as a function of u)
  # pseudo-observations u \in [0,1]^n
  # parameters alpha_vec \in R+^{n+1} 
  
  if( any(u > 1) | any(u < 0) ){
    res = 0
  } else {
    n = length(u)
    cons_num = 1
    sum_num = c()
    den = 1
    for(i in 1:n){
      cons_num = cons_num*gamma(sum(alpha_vec)-alpha_vec[i])
      den = den*(1-qbeta(u[i], alpha_vec[i], sum(alpha_vec)-alpha_vec[i]))^(sum(alpha_vec)-alpha_vec[i]-1)
      sum_num = c(sum_num, qbeta(u[i], alpha_vec[i], sum(alpha_vec)-alpha_vec[i]))
    }
    
    cons_den = gamma(sum(alpha_vec))*gamma(alpha_vec[n+1])
    K = cons_num/cons_den
    
    num = (1-sum(sum_num))^(alpha_vec[n+1]-1)
    
    res = K*num/den
    
    if(sum(sum_num)>1) res = 0 # vincolo in due variabili 
    
  }
  return(res)
}

# Parametric estimation of margins and copula ---> OK 
Pfit_margins_cop = function(data, margin_family){
  
  margin_params = list()
  U = matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])
  for(m in 1:length(margin_family)){
    if(margin_family[m]=="norm"){
      fit_m = fitdistr(data[,m],"normal")
      margin_params[[m]] = list(mean = fit_m$estimate[1], sd = fit_m$estimate[2])
      U[,m] = pnorm(data[,m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd)
    } else if (margin_family[m]=="t"){
      fit_m = fitdistr(data[,m],"t", start = list(m=mean(data[,m]),s=sd(data[,m]), df=3), lower=c(-1, 0.001,1))
      margin_params[[m]] = list(df = fit_m$estimate[3])
      U[,m] = pt(data[,m], df = margin_params[[m]]$df)
    } else if (margin_family[m]=="mixnorm"){
      # mixR::mixfit(data[,m], ncomp = 2) -----> da errore per qualche motivo incomprensibile [sembra essere dovuto al loop su m]
      fit_m = mclust::Mclust(data[,m], modelNames = mclustModelNames("V"), G = 2)$param
      margin_params[[m]] = list(mean = fit_m$mean, sd = sqrt(fit_m$variance$sigmasq), pro = fit_m$pro)
      U[,m] = KScorrect::pmixnorm(data[,m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd, pro = margin_params[[m]]$pro)
    } else if (margin_family[m]=="beta"){
      # adjust this
      fit_m = fitdistr(data[,m], "beta", start=list(shape1=1/2, shape2=1/2)) # in library(fitdistplus) start can be omitted
      margin_params[[m]] = list(shape1 = fit_m$estimate[1], shape2 = fit_m$estimate[2])
      U[,m] = pbeta(data[,m], shape1 = margin_params[[m]]$shape1, shape2 = margin_params[[m]]$shape2)
    } else {
      stop("Execution stopped. Only 'norm', 't', 'beta', and 'mixnorm' specifications are allowed for the marginals")
    }
  }
  
  estCopula = VineCopula::BiCopSelect(U[,1], U[,2])
  
  return(list(margin_family = margin_family, margin_params = margin_params, copula_model = estCopula))
}
#Pfit_margins_cop(data = data_sim, margin_family = margin_family)


# Proposed metrics -------------------

# |- M2. NP_cop (KDE) ---> OK ---------------
NP_cop = function(data){
  U = pobs(data)
  fit_cop = kdecop(U, method = "TLL2nn") # from kdecopula package: https://cran.r-project.org/web/packages/kdecopula/vignettes/kdecopula.pdf
  # here I use the TLL2nn transformation method (see 2.3. of the above and Fig 7 for a comparison among different methods)
  
  fit_x <- kde1d(data[,1]) # estimate marginal 1
  fit_y <- kde1d(data[,2]) # estimate marginal 2
  
  # combine copula and marginals
  npC_est = dkdecop(U, obj = fit_cop)*dkde1d(data[,1], fit_x)*dkde1d(data[,2], fit_y)
  
  return(list(U = U, fit_cop = fit_cop, dist_mv = npC_est))
}


# |- M3. P_Cop ---> OK ---------------
# library(kdevine)
# fit = kdevine(data_S1)
# np_est = dkdevine(as.matrix(data_S1), obj = fit)

P_cop = function(data, margin_family, margin_params, copula_model){
  est_copP = VC2copula::BiCop2copula(family = copula_model$family, par = copula_model$par, par2 = copula_model$par2)
  mvd_est = suppressMessages(mvdc(copula = est_copP,
                 margins = margin_family,
                 paramMargins = margin_params))

  pC_est = dMvdc(data, mvdc = mvd_est)
  
  return(P_est = pC_est)
}


# |- M4. knn_dist ---> OK ---------------
knn_dist = function(ax, ay, data, k){
  all_dist = mapply(function(bx, by) sqrt((ax-bx)^2+(ay-by)^2), data[,1], data[,2])
  k_dist = sum(sort(all_dist)[2:(k+1)])
  return(k_dist)
}


# |- M5. knn_CDFdist ---> OK ---------------

# this is a weighted sum, with weights given by the inverse of the Eucl distance: sum ( num / den )
# the weighted sum works better than the sum (num) / sum (den)
# Note: sum ( den / num ) gives a different (better result). This is what we use
# Apparently taking the distance and weighting it by the cdf is better than the opposite
knn_CDFdist = function(ax, ay, data, k, cdf_x, cdf_y){
  
  # The function can be modified so as to let it estimate the cdf (empirical cdf) when this is not given, i.e. cdf_x = NULL, cdf_y = NULL
  # if(is.null(cdf_x)){
  #   cdf_x = ecdf(data[,1])
  # }
  # 
  # if(is.null(cdf_y)){
  #   cdf_y = ecdf(data[,2])
  # }
  
  all_dist = mapply(function(bx, by) sqrt((ax-bx)^2+(ay-by)^2), data[,1], data[,2])
  #k_dist = sum(sort(all_dist)[1:(k+1)])
  idx = order(all_dist)[1:(k+1)]
  
  k_CDFdist = mapply(function(bx, by) sqrt((cdf_x(ax)-cdf_x(bx))^2+(cdf_y(ay)-cdf_y(by))^2), data[idx,1], data[idx,2])
  
  num = sort(all_dist)[2:(k+1)]
  den = k_CDFdist[2:(k+1)]
  
  return(sum(num/den))
}


# |- M6. d_CDFdist_uv ---> OK ---------------

# Method 6 (in forse)
# CDF distance in a d-neighborhood region (~PDF in a d-region) with UNI ECDFs --> UV reasoning: Does not work very well
d_CDFdist_uv = function(ax, ay, d, data, cdf_x, cdf_y){
  
  # The function can be modified so as to let it estimate the cdf (empirical cdf) when this is not given, i.e. cdf_x = NULL, cdf_y = NULL
  # if(is.null(cdf_x)){
  #   cdf_x = ecdf(data[,1])
  # }
  # 
  # if(is.null(cdf_y)){
  #   cdf_y = ecdf(data[,2])
  # }
  
  if(ax-d[1] < min(data[,1])){
    mydist_X = (cdf_x(ax+d[1])-cdf_x(ax))/d[1]
  } else if (ax+d[1] > max(data[,1])){
    mydist_X = (1-cdf_x(ax-d[1]))/d[1]
  }else{
    mydist_X = (cdf_x(ax+d[1])-cdf_x(ax-d[1]))/(2*d[1])
  }
  
  #mydist_X = (cdf_x(ax+d1)-cdf_x(ax-d1))^2
  #mydist_X = (cdf_x(ax+d1)-cdf_x(ax-d1)/(2*d1))^2
  
  if(ay-d[2] < min(data[,2])){
    mydist_Y = (cdf_y(ay+d[2])-cdf_y(ay))/d[2]
  } else if (ay+d[2] > max(data[,2])){
    mydist_Y = (1-cdf_y(ay-d[2]))/d[2]
  }else{
    mydist_Y = (cdf_y(ay+d[2])-cdf_y(ay-d[2]))/(2*d[2])
  }
  
  #mydist_Y = (cdf_y(ay+d2)-cdf_y(ay-d2))^2
  #mydist_Y = (cdf_y(ay+d2)-cdf_y(ay-d2)/(2*d2))^2
  
  #return(sqrt(mydist_X+mydist_Y)/(2*d1+2*d2))
  return(sqrt(mydist_X^2+mydist_Y^2))
}


# |- M7. d_CDFdist ---> OK ---------------

d_CDFdist = function(ax, ay, d, data){
  colnames(data) = c("V1", "V2")
  a_CDF_mv = empirical_cdf(x=data.table(data), ubounds = data.table(V1 = ax+d[1], V2 = ay+d[2]))$CDF
  b_CDF_mv = empirical_cdf(x=data.table(data), ubounds = data.table(V1 = ax+d[1], V2 = ay-d[2]))$CDF
  c_CDF_mv = empirical_cdf(x=data.table(data), ubounds = data.table(V1 = ax-d[1], V2 = ay+d[2]))$CDF
  d_CDF_mv = empirical_cdf(x=data.table(data), ubounds = data.table(V1 = ax-d[1], V2 = ay-d[2]))$CDF
  
  # area = 4*d1*d2 # we don't need this; it is constant wrt to the distance; see Eq. (4)
  # k_CDFdist_mv = (a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv)/area
  k_CDFdist_mv = a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv
  
  return(k_CDFdist_mv)
}


# |- M8. d_NPCopCDFdist ---> OK -----------
d_NPCopCDFdist = function(ax, ay, d, data, cdf_x, cdf_y, fit_cop){
  
  Ux_r = cdf_x(ax+d[1])
  Ux_l = cdf_x(ax-d[1])
  Uy_r = cdf_y(ay+d[2])
  Uy_l = cdf_y(ay-d[2])
  
  a_CDF_mv = pkdecop(c(Ux_r, Uy_r), obj = fit_cop)
  b_CDF_mv = pkdecop(c(Ux_r, Uy_l), obj = fit_cop)
  c_CDF_mv = pkdecop(c(Ux_l, Uy_r), obj = fit_cop)
  d_CDF_mv = pkdecop(c(Ux_l, Uy_l), obj = fit_cop)
  
  # area = 4*d1*d2 # we don't need this; it is constant wrt to the distance; see Eq. (4)
  # k_CDFdist_mv = (a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv)/area
  k_CDFdist_mv = a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv
  
  return(k_CDFdist_mv)
}


# |- M9. d_PCopCDFdist ---> OK ---------------

d_PCopCDFdist = function(ax, ay, d, margin_family, margin_params, copula_model){
  # x = c(x1, x2) is a point in 2 dimensions
  # d = c(d1, d2) is a distance in 2 dimensions
  # data is a df with 2 columns
  x = c(ax, ay)
  U_r = U_l = c()
  for(m in 1:length(margin_family)){
    if(margin_family[m]=="norm"){
      U_r[m] = pnorm(x[m]+d[m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd)
      U_l[m] = pnorm(x[m]-d[m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd)
    } else if (margin_family[m]=="t"){
      U_r[m] = pt(x[m]+d[m], df = margin_params[[m]]$df)
      U_l[m] = pt(x[m]-d[m], df = margin_params[[m]]$df)
    } else if (margin_family[m]=="mixnorm"){
      U_r[m] = KScorrect::pmixnorm(x[m]+d[m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd, pro = margin_params[[m]]$pro)
      U_l[m] = KScorrect::pmixnorm(x[m]-d[m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd, pro = margin_params[[m]]$pro)
    } else if (margin_family[m]=="beta"){
      U_r[m] = pbeta(x[m]+d[m], shape1 = margin_params[[m]]$shape1, shape2 = margin_params[[m]]$shape2)
      U_l[m] = pbeta(x[m]-d[m], shape1 = margin_params[[m]]$shape1, shape2 = margin_params[[m]]$shape2)
    } else {
      stop("Execution stopped. Only 'norm', 't', 'beta', and 'mixnorm' specifications are allowed for the marginals")
    }
  }
  
  # F(x,y) = C(F(x), F(y))
  a_CDF_mv = BiCopCDF(U_r[1], U_r[2], copula_model)
  b_CDF_mv = BiCopCDF(U_r[1], U_l[2], copula_model)
  c_CDF_mv = BiCopCDF(U_l[1], U_r[2], copula_model)
  d_CDF_mv = BiCopCDF(U_l[1], U_l[2], copula_model)
  
  # area = 4*d1*d2 # we don't need this; it is constant wrt to the distance; see Eq. (4)
  # k_CDFdist_mv = (a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv)/area
  k_CDFdist_mv = a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv
  
  return(k_CDFdist_mv)
}


# Main function evaluating all metrics with final plot --------------------

get_ALL_err = function(mvdc = mvd_def_sD3, true_level = true_q_sD3, Dirichlet_alpha = NA, 
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
  
  # This function evaluates the following 9 measures: M0-M8
  measures = c("M0: Known density","M1: KDE", "M2: KDE-NPCop", "M3: KDE-PCop", 
               "M4: kNN-Euclidean", "M5: kNN-CDFuv", #"M6:eps-CDF-uv", 
               "M6: Ɛ-CDFmv", "M7: Ɛ-CDF-NPCop","M8: Ɛ-CDF-PCop")
  
  if(any(is.na(Dirichlet_alpha)==T)){
    # Simulate the data from an mvdc object and get the true density
    data_sim = rMvdc(n=Nsample, mvdc = mvdc)
    fxy_true = dMvdc(x = data_sim, mvdc = mvdc)
    margin_family = mvdc@margins
    main_title = paste("Scenario:", paste(mvdc@margins,collapse=", "), "margins &", class(mvdc@copula)[1])
  }else{
    data_sim = rdirichlet(Nsample, alpha = Dirichlet_alpha)
    fxy_true = ddirichlet(x = data_sim, alpha = Dirichlet_alpha)
    data_sim = data_sim[,1:(length(Dirichlet_alpha)-1)]
    margin_family = c("beta", "beta")
    main_title = paste0("Scenario: Dirichlet (", paste(Dirichlet_alpha,collapse=", "), ")")
  }
  
  if(plot == T){
    op <- par(mfrow = c(3,3),
              oma = c(2,2,2,1) + 0.1,
              mar = c(2,2,2,1) + 0.1,
              mai = c(1, 0.1, 0.1, 0.1),
              bty = 'n')
    
    if(any(is.na(Dirichlet_alpha)==T)){
      mat = outer(seq(min(data_sim[,1])-1, max(data_sim[,1])+1, length.out = 100), 
                  seq(min(data_sim[,2])-1, max(data_sim[,2])+1, length.out = 100), 
                  FUN = Vectorize(function(x, y){dMvdc(c(x,y), mvdc = mvdc)}))
    } else {
      mat = outer(seq(min(data_sim[,1])-1, max(data_sim[,1])+1, length.out = 100), 
                  seq(min(data_sim[,2])-1, max(data_sim[,2])+1, length.out = 100), 
                  FUN = Vectorize(function(x, y){ddirichlet(c(x,y, 1-x-y), alpha = Dirichlet_alpha)}))
    }
  }
  
  # M0: true density ---> OK
  res_M0 = HDR_class(data = data_sim, true_density = fxy_true, est_density = fxy_true, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[1], ")"))
  
  # M1: KDE ---> OK
  # Default: myh = ks::Hpi.diag(data_S2,binned=TRUE) 
  # OPTIMAL bandwidth: Chacón, J.E., Duong, T., Wand, M.P. (2011), Asymptotics for General Multivariate Kernel Density Derivative Estimators, Statistica Sinica, 21, 807–840.
  myh_opt = (4/(dim(data_sim)[2]+4))^(1/(dim(data_sim)[2]+6))*dim(data_sim)[1]^(-1/(dim(data_sim)[2]+6))*abs(cov(data_sim))^(1/2)
  fit_kde = kde(data_sim, eval.points = data_sim, H = myh_opt) # from ks package
  res_M1 = HDR_class(data = data_sim, true_density = fxy_true, est_density = fit_kde$estimate, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[2], ")"))
  
  # M2: nonparametric copula ---> OK
  npC_est = NP_cop(data_sim)
  res_M2 = HDR_class(data = data_sim, true_density = fxy_true, est_density = npC_est$dist_mv, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[3], ")"))
  
  # M3: parametric copula ---> OK
  parametric_fit = Pfit_margins_cop(data = data_sim, margin_family = margin_family)
  pC_est = P_cop(data = data_sim, margin_family = parametric_fit$margin_family, 
                 margin_params = parametric_fit$margin_params, copula_model = parametric_fit$copula_model) 
  res_M3 = HDR_class(data = data_sim, true_density = fxy_true, est_density = pC_est, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[4], ")"))
  
  # M4: k-neighborhood distances ---> OK
  # Optimal k: sqrt(myn/2). Rule of thumb based on simulations: sqrt(myn*tau) / sqrt(myn*rho) / sqrt(myn/2)
  myk = round(sqrt(Nsample*.5)) # Rule of thumb based on simulations: sqrt(Nsample*tau) / sqrt(Nsample*rho) 
  dist_M4 = pbmapply(knn_dist, data_sim[,1], data_sim[,2], MoreArgs = list(k = myk, data = data_sim))
  # pay attention to the measure being a concentration/sparsity measure 
  res_M4 = HDR_class(data = data_sim, true_density = fxy_true, est_density = 1/dist_M4, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[5], ")"))
  
  # M5: Weighted sum of CDF distances from the k nearest neighbors (UV CDF estimation) ---> OK
  # Variation of https://core.ac.uk/download/pdf/30276753.pdf
  # Not a Metric (triangular inequality satisfied only in specif ordering between three points)
  myk = 30 # quite robust to k 
  est_cdfx = ecdf(data_sim[,1])
  est_cdfy = ecdf(data_sim[,2])
  dist_M5 = pbmapply(knn_CDFdist, data_sim[,1], data_sim[,2], MoreArgs = list(k = myk, data = data_sim, cdf_x=est_cdfx, cdf_y=est_cdfy))
  res_M5 = HDR_class(data = data_sim, true_density = fxy_true, est_density = 1/dist_M5, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[6], ")"))
  
  # # M: CDF distance in a d-neighborhood region (~PDF in a d-region) with UNI ECDFs 
  # myd = c(sum(abs(range(data_sim[,1])))/sqrt(Nsample), sum(abs(range(data_sim[,2])))/sqrt(Nsample))
  # dist_M = pbmapply(d_CDFdist_uv, data_sim[,1], data_sim[,2], MoreArgs = list(d = myd, data = data_sim, cdf_x = est_cdfx, cdf_y = est_cdfy))
  # res_M = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M, true_level = true_level, 
  #                  coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = paste0(coverage_prob*100,"% HDR (", measures[6], ")"))
  
  # M6: CDF distance in a d-neighborhood region (~PDF in a d-region) with Mv ECDF ---> OK
  # Optimal d?
  myd = ifelse(any(is.na(Dirichlet_alpha)==T), exp(2.13-0.3*log(Nsample)), 0.10)
  dist_M6 = pbmapply(d_CDFdist, data_sim[,1], data_sim[,2], MoreArgs = list(data = data_sim, d = rep(myd, 2)))
  res_M6 = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M6, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = bquote(bold(.(coverage_prob*100)*"% HDR (M6: "*epsilon*"-CDFmv)")))
  
  # M8: CDF distance in a d-neighborhood region (~PDF in a d-region) with NP Copula ---> OK
  # from kdecopula package: https://cran.r-project.org/web/packages/kdecopula/vignettes/kdecopula.pdf
  # here I use the TLL2nn transformation method (see 2.3. of the above and Fig 7 for a comparison among different methods)
  # Optimal d? ---> fitted nonlinear regression with exponential decay (otherwise myd = 1.5)
  myd = ifelse(any(is.na(Dirichlet_alpha)==T), exp(1.74-0.26*log(Nsample)), exp(-1.22-0.23*log(Nsample)))
  dist_M7 = pbmapply(d_NPCopCDFdist, data_sim[,1], data_sim[,2], 
                     MoreArgs = list(data = data_sim, d = rep(myd, 2), cdf_x = est_cdfx, cdf_y = est_cdfy, fit_cop = npC_est$fit_cop))
  res_M7 = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M7, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = bquote(bold(.(coverage_prob*100)*"% HDR (M7: "*epsilon*"-CDF-NPCop)")))
  
  # M9: CDF distance in a d-neighborhood region with Parametric Copula ---> OK
  # Optimal d? ---> fitted nonlinear regression with exponential decay (otherwise myd = 0.7). For Dirichlet, take small values < 0.05
  myd = ifelse(any(is.na(Dirichlet_alpha)==T), exp(1.60 - 0.41*log(Nsample)), 0.02)
  dist_M8 = pbmapply(d_PCopCDFdist, data_sim[,1], data_sim[,2], 
                     MoreArgs = list(d = rep(myd, 2), margin_family = parametric_fit$margin_family, 
                                     margin_params = parametric_fit$margin_params, copula_model = parametric_fit$copula_model))
  res_M8 = HDR_class(data = data_sim, true_density = fxy_true, est_density = dist_M8, true_level = true_level, 
                     coverage_prob = coverage_prob, build_plot = plot, mat_plot = mat, main_title = bquote(bold(.(coverage_prob*100)*"% HDR (M6: "*epsilon*"-CDF-PCop)")))
  
  if(plot == T){title(main_title, line = 0.5, cex.main = 2, col.main = "blue", outer = TRUE)}
  
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
  rownames(Error_mat) = c("Total Error Rate", "FPR", "FNR", "Accuracy", "F1 Score", "MCC")
  colnames(Error_mat) = measures
  
  return(Error_mat)
}



