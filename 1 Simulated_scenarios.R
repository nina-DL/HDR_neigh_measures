# Load Packages ---------------------------------------------------------------
library(rgl) 
library(MASS) 
library(VC2copula) 
library(ks)
library(kdecopula)
library(kde1d)
library(VineCopula)
library(pbapply)
library(KScorrect)
library(data.table)
library(gsl) # for Debye function
library(gtools) #for Dirichlet

# hypersetup
myN = 1e7 # this is a big N in order to provide a good approx of the contour and the level that ensures a given coverage prob
coverage_prob = 0.95  # our contour of interest is defined on this coverage probability

# see here for an overview of copula models and relationship with tau:
# https://www.sciencedirect.com/science/article/pii/S0377042710006138#f000010

# Define copulas and marginals  -------------

# Gaussian parameters (All scenarios, except xy"t")
mu_x = 0; mu_y = 1
sd_x = 2; sd_y = 2

# Student-t parameters (Only scenarios with xy"t")
my_df1 = my_df2 = 2

# Dirichlet parameters
myalpha = c(2, 2, 2) # if we use c(1,1,1) we have a unif distribution. how would an HDR be in this case

# Define copula parameters: all copulas are based on a Kendall tau = 0.5
my_tau = 0.5

# Gaussian copula parameters (scenarios A)
par_gaus = sin(pi/2*my_tau) # see e.g., http://www.scielo.org.co/pdf/rce/v43n1/0120-1751-rce-43-01-3.pdf
cop_gaus= BiCop2copula(family = 1, par = par_gaus) # Gauss cop
tau(cop_gaus)

# Student-t copula parameters (scenarios B)
par_t = sin(pi/2*my_tau)
cop_t = BiCop2copula(family = 2, par = par_t, par2 = 6) 
tau(cop_t)

# Frank copula parameters (scenarios C)
debye_1(10, give=FALSE, strict=TRUE) # library(gsl)
par_frank = 5.75
1-4/par_frank*(1-debye_1(par_frank, give=FALSE, strict=TRUE))
cop_frank = BiCop2copula(family = 5, par = par_frank) 
tau(cop_frank)

# Clayton copula parameters (scenarios D)
par_clay = 2*my_tau/(1-my_tau)
cop_clay = BiCop2copula(family = 3, par = par_clay) 
tau(cop_clay)

# graphical parameters
op <- par(mfrow = c(4,4),
          oma = c(1,2,2,1) + 0.1,
          mar = c(1,2,2,1) + 0.1,
          bty = 'n')

cols <- hcl.colors(10, "YlOrRd")

# Scenario sA1t (unimodal): Student-t margins - Gaussian copula ----------------------
mvd_def_sA1t = mvdc(copula = cop_gaus,
                    margins = c("t", "t"),
                    paramMargins = list(list(df = my_df1),
                                        list(df = my_df2)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sA1t)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sA1t)
true_q_sA1t = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sA1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), main = "Unimodal: \nStudent-t marginals", 
        col = cols)
mtext(side = 2, bquote("Gaussian copula ("*theta*" = 0.7)"), line = 2, cex = .9)
contour(mvd_def_sA1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sA1t, add = T)

# use this to build the contour after copula
# DEPRECATED: using now the one provided by contour() with copula which is more accurate
# mat = outer(seq(min(mvt_draws[,1])-1, max(mvt_draws[,1])+1, length.out = 100), 
#             seq(min(mvt_draws[,2])-1, max(mvt_draws[,2])+1, length.out = 100), 
#             FUN = Vectorize(function(x, y){dMvdc(c(x,y), mvdc = mvd_def_sA1t)}))
# 
# contour(seq(min(mvt_draws[,1])-1, max(mvt_draws[,1])+1, length.out = 100), 
#         seq(min(mvt_draws[,2])-1, max(mvt_draws[,2])+1, length.out = 100), 
#         mat, levels = true_q_sA1t, add = T, col = "lightblue", lwd = 2)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sA1 (unimodal): Gaussian margins - Gaussian copula ----------------------

mvd_def_sA1 = mvdc(copula = cop_gaus,
                   margins = c("norm", "norm"),
                   paramMargins = list(list(mean = mu_x, sd = sd_x),
                                       list(mean = mu_y, sd = sd_y)))
set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sA1)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sA1)
true_q_sA1 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sA1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), main = "Unimodal: \nGaussian marginals",
        col = cols)
contour(mvd_def_sA1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sA1, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sA2 (bimodal): Gaussian-MixGaussian margins - Gaussian copula ----------------------

mean_mat <- rbind(c(mu_x,mu_y), c(mu_x+5,mu_y+5))
myw = c(0.5, 0.5)

mvd_def_sA2 = mvdc(copula = cop_gaus,
                   margins = c("norm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[1,1], sd = c(sd_x)),
                                       list(mean = mean_mat[,2], sd = c(sd_x,sd_y), pro = myw)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sA2)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sA2)
true_q_sA2 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sA2, dMvdc, xlim=c(-6, 6), ylim=c(-6, 13), main = "Bimodal: \nGaussian & GaussianMixture marginals", col = cols)
contour(mvd_def_sA2, dMvdc, xlim=c(-7, 7), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sA2, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sA3 (quadrimodal): MixGaussian margins - Gaussian copula ----------------------

mean_mat <- rbind(c(mu_x,mu_y), c(mu_x+9,mu_y+7))
myw = c(0.5, 0.5)

mvd_def_sA3 = mvdc(copula = cop_gaus,
                   margins = c("mixnorm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[,1], sd = c(sd_x, sd_y), pro = myw),
                                       list(mean = mean_mat[,2], sd = c(sd_x, sd_y), pro = myw)))
set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sA3)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sA3)
true_q_sA3 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sA3, dMvdc, xlim=c(-5, 15), ylim=c(-4, 15), main = "Quadrimodal: \nGaussianMixture marginals", col = cols)
contour(mvd_def_sA3, dMvdc, xlim=c(-6, 16), ylim=c(-6, 16), col = "blue", lwd = 2,
        levels = true_q_sA3, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sB1t (unimodal): Student-t margins - t copula ----------------------

mvd_def_sB1t = mvdc(copula = cop_t,
                    margins = c("t", "t"),
                    paramMargins = list(list(df = my_df1),
                                        list(df = my_df2)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sB1t)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sB1t)
true_q_sB1t = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sB1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
mtext(side = 2, bquote("Student-t copula ("*theta*" = 0.7; "*nu*" = 6)"), line = 2, cex = .9)
contour(mvd_def_sB1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sB1t, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sB1 (unimodal): Gaussian margins - t copula ----------------------

mvd_def_sB1 = mvdc(copula = cop_t,
                   margins = c("norm", "norm"),
                   paramMargins = list(list(mean = mu_x, sd = sd_x),
                                       list(mean = mu_y, sd = sd_y)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sB1)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sB1)
true_q_sB1 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sB1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
contour(mvd_def_sB1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sB1, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sB2 (bimodal): Gaussian-MixGaussian margins - t copula ----------------------

mvd_def_sB2 = mvdc(copula = cop_t,
                   margins = c("norm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[1,1], sd = c(sd_x)),
                                       list(mean = mean_mat[,2], sd = c(sd_x,sd_y), pro = myw)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sB2)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sB2)
true_q_sB2 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sB2, dMvdc, xlim=c(-6, 6), ylim=c(-6, 13), col = cols)
contour(mvd_def_sB2, dMvdc, xlim=c(-7, 7), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sB2, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sB3 (quadrimodal): MixGaussian margins - t copula ----------------------

mvd_def_sB3 = mvdc(copula = cop_t,
                   margins = c("mixnorm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[,1], sd = c(sd_x, sd_y), pro = myw),
                                       list(mean = mean_mat[,2], sd = c(sd_x, sd_y), pro = myw)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sB3)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sB3)
true_q_sB3 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sB3, dMvdc, xlim=c(-5, 15), ylim=c(-4, 15), col = cols)
contour(mvd_def_sB3, dMvdc, xlim=c(-6, 16), ylim=c(-5, 16), col = "blue", lwd = 2,
        levels = true_q_sB3, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sC1t (unimodal): Student-t margins - Frank copula ----------------------

mvd_def_sC1t = mvdc(copula = cop_frank,
                    margins = c("t", "t"),
                    paramMargins = list(list(df = my_df1),
                                        list(df = my_df2)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sC1t)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sC1t)
true_q_sC1t = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sC1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
mtext(side = 2, bquote("Frank copula ("*theta*" = 5.75)"), line = 2, cex = .9)
contour(mvd_def_sC1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sC1t, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sC1 (unimodal): Gaussian margins - Frank copula ----------------------

mvd_def_sC1 = mvdc(copula = cop_frank,
                   margins = c("norm", "norm"),
                   paramMargins = list(list(mean = mu_x, sd = sd_x),
                                       list(mean = mu_y, sd = sd_y)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sC1)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sC1)
true_q_sC1 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sC1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
contour(mvd_def_sC1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sC1, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sC2 (bimodal): Gaussian-MixGaussian margins - Frank copula ----------------------

mvd_def_sC2 = mvdc(copula = cop_frank,
                   margins = c("norm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[1,1], sd = c(sd_x)),
                                       list(mean = mean_mat[,2], sd = c(sd_x,sd_y), pro = myw)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sC2)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sC2)
true_q_sC2 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sC2, dMvdc, xlim=c(-6, 6), ylim=c(-6, 13), col = cols)
contour(mvd_def_sC2, dMvdc, xlim=c(-7, 7), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sC2, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sC3 (quadrimodal): MixGaussian margins - Frank copula ----------------------

mvd_def_sC3 = mvdc(copula = cop_frank,
                   margins = c("mixnorm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[,1], sd = c(sd_x, sd_y), pro = myw),
                                       list(mean = mean_mat[,2], sd = c(sd_x, sd_y), pro = myw)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sC3)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sC3)
true_q_sC3 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sC3, dMvdc, xlim=c(-6, 15), ylim=c(-6, 13), col = cols)
contour(mvd_def_sC3, dMvdc, xlim=c(-7, 16), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sC3, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sD1t (unimodal): Student-t margins - Clayton copula ----------------------

mvd_def_sD1t = mvdc(copula = cop_clay,
                    margins = c("t", "t"),
                    paramMargins = list(list(df = my_df1),
                                        list(df = my_df2)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sD1t)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sD1t)
true_q_sD1t = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sD1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
mtext(side = 2, bquote("Clayton copula ("*theta*" = 2)"), line = 2, cex = .9)
contour(mvd_def_sD1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sD1t, add = T)
# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sD1 (unimodal): Gaussian margins - Clayton copula ----------------------

mvd_def_sD1 = mvdc(copula = cop_clay,
                   margins = c("norm", "norm"),
                   paramMargins = list(list(mean = mu_x, sd = sd_x),
                                       list(mean = mu_y, sd = sd_y)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sD1)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sD1)
true_q_sD1 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sD1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
contour(mvd_def_sD1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sD1, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sD2 (unimodal): Gaussian-MixGaussian margins - Clayton copula ----------------------

mvd_def_sD2 = mvdc(copula = cop_clay,
                   margins = c("norm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[1,1], sd = c(sd_x)),
                                       list(mean = mean_mat[,2], sd = c(sd_x,sd_y), pro = myw)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sD2)
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sD2)
true_q_sD2 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sD2, dMvdc, xlim=c(-6, 5), ylim=c(-6, 13), col = cols)
contour(mvd_def_sD2, dMvdc, xlim=c(-7, 6), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sD2, add = T)

# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package

# Scenario sD3 (quadrimodal): MixGaussian margins - Clayton copula ----------------------

mvd_def_sD3 = mvdc(copula = cop_clay,
                   margins = c("mixnorm", "mixnorm"), #mixnorm margin will use qmixnorm from KScorrect package
                   paramMargins = list(list(mean = mean_mat[,1], sd = c(sd_x, sd_y), pro = myw),
                                       list(mean = mean_mat[,2], sd = c(sd_x, sd_y), pro = myw)))

set.seed(12345)
mvt_draws = rMvdc(n=myN, mvdc = mvd_def_sD3)
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sD3)
true_q_sD3 = quantile(fxy_true, probs = (1-coverage_prob))

contour(mvd_def_sD3, dMvdc, xlim=c(-6, 15), ylim=c(-6, 13), col = cols)
contour(mvd_def_sD3, dMvdc, xlim=c(-7, 16), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sD3, add = T)
# 3D plot
# Biv.kde <- kde2d(mvt_draws[,1], mvt_draws[,2], n = 100)   # from MASS package
# col.kde <- heat.colors(length(Biv.kde$z))[rank(Biv.kde$z)]
# persp3d(x=Biv.kde, col = col.kde) # from rgl package


# Scenario sDir (compact): Dirichlet on K=2 simplex ----------------------
myalpha = c(1,1,1)

set.seed(12345)
data_sDir3 = rdirichlet(myN, alpha = myalpha)
fxy_true_sDir = ddirichlet(x = data_sDir3, alpha = myalpha)
true_q_sDir_112 = quantile(fxy_true_sDir, probs = (1-coverage_prob))

# For the contours of the copula density
U1 = pbeta(data_sDir3[,1], myalpha[1], sum(myalpha)-myalpha[1])
U2 = pbeta(data_sDir3[,2], myalpha[2], sum(myalpha)-myalpha[2])
U = cbind(U1, U2)

# check constraint (IGNORE)
# apply(data_sDir3[,1:2], 1, sum) <= 1
# sqrt(1-U1) + sqrt(1-U2) >= 1
# (1-U1)^(1/3) + (1-U2)^(1/3)  >= 1


# after loading d_cop_u() function (the copula density of the Dirichlet)
c_true_sDir = apply(cbind(U1, U2), 1, d_cop_u, alpha_vec = myalpha)
true_qC_sDir_112 = quantile(c_true_sDir, probs = (1-coverage_prob))

# Plot Dir distr in k=2, True HDR, Copula
par(mfrow=c(1,2), bty = 'n')
plot(data_sDir3[1:1000,1:2], cex = 0.2, col = "lightgray", 
     main = bquote("S17: Dirichlet ("*.(myalpha[1])*","*.(myalpha[2])*","*.(myalpha[3])*")"),
     xlab = expression(x[1]), ylab = expression(x[2]), xlim = c(0,1), ylim = c(0,1))

# use this to build the contour 
mat = outer(seq(min(data_sDir3[,1])-1, max(data_sDir3[,1])+1, length.out = 100), 
            seq(min(data_sDir3[,2])-1, max(data_sDir3[,2])+1, length.out = 100), 
            FUN = Vectorize(function(x, y){ddirichlet(c(x,y, 1-x-y), alpha = myalpha)}))

contour(seq(min(data_sDir3[,1])-1, max(data_sDir3[,1])+1, length.out = 100), 
        seq(min(data_sDir3[,2])-1, max(data_sDir3[,2])+1, length.out = 100), 
        mat, levels = true_q_sDir_112, add = T, col = "blue", lwd = 2)

p_cont = quantile(fxy_true_sDir, probs = seq(0.1, 1, 0.1))
for(i in 1:length(cols)){
        contour(seq(min(data_sDir3[,1])-1, max(data_sDir3[,1])+1, length.out = 100), 
                seq(min(data_sDir3[,2])-1, max(data_sDir3[,2])+1, length.out = 100), 
                mat, levels = p_cont[i], add = T, col = cols[i], lwd = 1)
}
plot(U[1:1000,], cex = 0.5, col = "lightgray", main = bquote("Dirichlet ("*.(myalpha[1])*","*.(myalpha[2])*","*.(myalpha[3])*") Copula"),
     xlab = expression(u[1]), ylab = expression(u[2]), xlim = c(0,1), ylim = c(0,1))

# use this to build the contour 
lim = 0.2
mat_cop = outer(seq(min(U1)-lim, max(U2)+lim, length.out = 100), 
                seq(min(U1)-lim, max(U2)+lim, length.out = 100), 
                FUN = Vectorize(function(x, y){d_cop_u(c(x,y), alpha = myalpha)}))

contour(seq(min(U1)-lim, max(U2)+lim, length.out = 100), 
        seq(min(U1)-lim, max(U2)+lim, length.out = 100), 
        mat_cop, levels = true_qC_sDir_112, add = T, col = "blue", lwd = 2)

p_cont = quantile(c_true_sDir, probs = seq(0.1, 1, 0.1))
for(i in 1:length(cols)){
        contour(seq(min(U1)-lim, max(U2)+lim, length.out = 100), 
                seq(min(U1)-lim, max(U2)+lim, length.out = 100), 
                mat_cop, levels = p_cont[i], add = T, col = cols[i], lwd = 1)
}


# Save simulated scenarios ------------------------------------------------

typeline <- function(msg="Enter text: ") {
        if (interactive() ) {
                txt <- readline(msg)
        } else {
                cat(msg);
                txt <- readLines("stdin",n=1);
        }
        return(txt)
}

wanna_save=typeline("Do you want to save the obtained results? (Type either Yes or No) ")

if(wanna_save=="Yes"){
        save(mvd_def_sA1t, mvd_def_sA1, mvd_def_sA2, mvd_def_sA3, 
             mvd_def_sB1t, mvd_def_sB1, mvd_def_sB2, mvd_def_sB3, 
             mvd_def_sC1t, mvd_def_sC1, mvd_def_sC2, mvd_def_sC3, 
             mvd_def_sD1t, mvd_def_sD1, mvd_def_sD2, mvd_def_sD3, 
             true_q_sA1t, true_q_sA1, true_q_sA2, true_q_sA3,
             true_q_sB1t, true_q_sB1, true_q_sB2, true_q_sB3,
             true_q_sC1t, true_q_sC1, true_q_sC2, true_q_sC3,
             true_q_sD1t, true_q_sD1, true_q_sD2, true_q_sD3, 
             true_q_sDir_112, true_qC_sDir_112, true_q_sDir_222, true_qC_sDir_222, myN,
             file = "Def_scenarios.RData")
}

