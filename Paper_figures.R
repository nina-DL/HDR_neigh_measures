library(ggplot2)
library(gridExtra)
library(ggpubr)

# Figure 1: HDR vs Symmetric Quantile interval ----------------------------

# Figure HDR and quantile interval
set.seed(10)
#generate data
stuff<-c(rnorm(10000), rnorm(4000,mean=7))
# kernel density estimate
out.d <- density(stuff)
# Determine W that cuts off p = 0.95
out.diff <- sapply(seq(0.35, 0.5,.001), function(x) sum(out.d$y[out.d$y>=quantile(out.d$y, x)])/sum(out.d$y)-0.95)
which.one <- which.min(out.diff[out.diff>0])
cut.off <- seq(0.35,0.5,.001)[which.one]

quantile(stuff, probs = c((1-0.95)/2, 1-(1-0.95)/2))

# credible intervals
my.x<-out.d$x[out.d$y>=quantile(out.d$y, cut.off)]
my.y<-out.d$y[out.d$y>=quantile(out.d$y, cut.off)]
my.x2 <- out.d$x[out.d$x<=quantile(stuff,.975) & out.d$x>=quantile(stuff,.025)]

# identify regions under curve
grp1.x <- my.x[my.x >=range(my.x[my.x<3])[1] & my.x <=range(my.x[my.x<3])[2]]
grp1.y <- my.y[my.x >=range(my.x[my.x<3])[1] & my.x <=range(my.x[my.x<3])[2]]
grp2.x <- my.x[my.x >=range(my.x[my.x>3])[1] & my.x <=range(my.x[my.x>3])[2]]
grp2.y <- my.y[my.x >=range(my.x[my.x>3])[1] & my.x <=range(my.x[my.x>3])[2]]

plot(density(stuff), xlab="X", main='',cex.lab=1.2, col='grey70', bty="n", ylab = "")
rug(my.x,side=1,col='dodgerblue1',line=-.5)
rug(my.x2,side=1,col=2)
#shade area under curve
polygon(c(grp1.x, grp1.x[length(grp1.x)],grp1.x[1]),c(grp1.y,0,0), col='lightblue', border=NA, density=50)
polygon(c(grp2.x, grp2.x[length(grp2.x)],grp2.x[1]),c(grp2.y,0,0), col='lightblue', border=NA, density=50)
segments(grp1.x[1],quantile(out.d$y, cut.off),grp2.x[length(grp2.x)],quantile(out.d$y, cut.off),lty=2)

arrows(3,.15,0.5,.1,code=2,angle=45,length=.1)
arrows(3.3,.15,5.5,.07,code=2,angle=45,length=.1)
text( 4.85,.155, expression(paste("Coverage = 100(1 - ", alpha, ")% = 95%")),pos=3,offset=0.1)


legend('topright',c(expression(paste(100(1-alpha),"% HDR")),expression(paste("[", q[alpha/2], ", ", q[1-alpha/2], "]"))),col=c('dodgerblue1',2), 
       pch=15, pt.cex=1.2, cex=.9, bty='n', title=expression(bold('Statistical interval')))
lines(out.d$x, out.d$y)
text(grp2.x[length(grp2.x)],quantile(out.d$y, cut.off),expression(f[alpha]),pos=4)



# Figure 2: HDR for different bandwidths h --------------------------------

library(kdecopula)
library(kde1d)
library(VineCopula)
library(VC2copula) 
library(hdrcde)

# load simulated scenarios and the true "estimated" HDR on a big ss
load("Sim Results/Simulated Scenarios/Def_scenarios.RData")

#library(hdrcde) # https://github.com/robjhyndman/hdrcde/blob/master/R/hdr.boxplot.2d.R
coverage_prob = 0.95
mvt_q = rMvdc(n=1e7, mvdc = mvd_def_sA1)
fxy_q = dMvdc(x = mvt_q, mvdc = mvd_def_sA1)
true_q_sA1 = quantile(fxy_q, probs = (1-coverage_prob))

set.seed(1300)
mvt_draws = rMvdc(n=200, mvdc = mvd_def_sA1)
#cor.test(mvt_draws[1:100,1], mvt_draws[1:100,2], method="kendall")
fxy_true = dMvdc(x = mvt_draws, mvdc = mvd_def_sA1)

op <- par(mfrow = c(1,3),
          oma = c(3,1,1,1) + 0.1,
          mar = c(3,3,2,1) + 0.1)

myh = ks::Hpi.diag(mvt_draws ,binned=TRUE) 
hdr_mv = hdr.2d(mvt_draws[,1], mvt_draws[,2], prob = (1-coverage_prob), kde.package = "ks", h=diag(myh))
plot(hdr_mv, pch = 20, show.points = T, outside.points = F, shadecols = "lightblue", pointcol = "cadetblue")
contour(mvd_def_sA1, dMvdc, col = "blue", xlim=c(-7, 7), ylim=c(-7, 7), lwd = 2, levels = true_q_sA1, add = T)
# outliers according to kde 
index <- (hdr_mv$fxy <= 0.99999*min(hdr_mv$falpha))
index1 = index
points(hdr_mv$x[index], hdr_mv$y[index], col="purple", pch = 4, cex = 1.5)
points(hdr_mv$x[117], hdr_mv$y[117], col="red", pch = 4, cex = 2)
title(xlab=expression(X[1]), line=2, cex.lab=.8)
title(ylab=expression(X[2]), line=2, cex.lab=.8)
title(xlab="h = (0.42, 0.42); Wand & Jones, 1994", line=-20, cex.lab=1.2)
title("95% HDR with KDE", line = 0, cex.main = 1.5, col.main = "black", outer = TRUE)

library(etasFLP)
myh_opt_Silv = ((4/(dim(mvt_draws )[2]+2))^(1/(dim(mvt_draws )[2]+4))*dim(mvt_draws)[1]^(-1/(dim(mvt_draws )[2]+4))*sd(mvt_draws))^2
myh_opt_Silv1 = bwd.nrd(mvt_draws[,1], w=replicate(length(mvt_draws[,1]),1), d = 2)
myh_opt_Silv2 = bwd.nrd(mvt_draws[,2], w=replicate(length(mvt_draws[,2]),1), d = 2)
hdr_mv = hdr.2d(mvt_draws[,1], mvt_draws[,2], prob = (1-coverage_prob), kde.package = "ks", h=c(myh_opt_Silv1, myh_opt_Silv2))
plot(hdr_mv, pch = 20, show.points = T, outside.points = F, shadecols = "lightblue", pointcol = "cadetblue")
contour(mvd_def_sA1, dMvdc, col = "blue", xlim=c(-7, 7), ylim=c(-7, 7), lwd = 2, levels = true_q_sA1, add = T)
# outliers according to kde 
index <- (hdr_mv$fxy <= 0.99999*min(hdr_mv$falpha))
index2 = index
mean(index)
points(hdr_mv$x[index], hdr_mv$y[index], col="purple", pch = 4, cex = 1.5)
points(hdr_mv$x[97], hdr_mv$y[97], col="red", pch = 4, cex = 2)
title(xlab=expression(X[1]), line=2, cex.lab=.8)
title(ylab=expression(X[2]), line=2, cex.lab=.8)
title(xlab="h = (0.96, 0.96); Silverman, 1986", line=-20, cex.lab=1.2)

myh_opt_Chacon = (4/(dim(mvt_draws )[2]+4))^(1/(dim(mvt_draws)[2]+6))*dim(mvt_draws )[1]^(-1/(dim(mvt_draws )[2]+6))*abs(cov(mvt_draws ))^(1/2)
hdr_mv = hdr.2d(mvt_draws[,1], mvt_draws[,2], prob = (1-coverage_prob), kde.package = "ks", h=diag(myh_opt_Chacon))
plot(hdr_mv, pch = 20, show.points = T, outside.points = F, shadecols = "lightblue", pointcol = "cadetblue")
contour(mvd_def_sA1, dMvdc, col = "blue", xlim=c(-7, 7), ylim=c(-7, 7), lwd = 2, levels = true_q_sA1, add = T)
# outliers according to kde 
index <- (hdr_mv$fxy <= 0.99999*min(hdr_mv$falpha))
index3 = index
points(hdr_mv$x[index], hdr_mv$y[index], col="purple", pch = 4, cex = 1.5)
points(hdr_mv$x[189], hdr_mv$y[189], col="red", pch = 4, cex = 2)
title(xlab=expression(X[1]), line=2, cex.lab=.8)
title(ylab=expression(X[2]), line=2, cex.lab=.8)
title(xlab="h = (1.12, 1.12); Chacon et al., 2011", line=-20, cex.lab=1.2)

idx = unique(c(which(index1 != index2), which(index1 != index3), which(index2 != index3)))

# Figure 3: CDF distance 3 points -----------------------------------------
library(EnvStats)

ggplot(data.frame(x = c(-4, 4)), aes(x = x)) +
  stat_function(fun = function(x){dnormMix(x, mean1 = -2, sd1 = .5, mean2 = 1, sd2 = .9, p.mix = 0.5)}) + 
  geom_vline(xintercept = -2, linetype="dashed", color = "black") +
  geom_vline(xintercept = -1.6, linetype="dashed", color = "black")+
  #geom_vline(xintercept = -1.14, linetype="dashed", color = "black") +
  geom_vline(xintercept = -0.72, linetype="dashed", color = "black") +
  geom_vline(xintercept = 0.5, linetype="dashed", color = "black") +
  stat_function(fun = function(x){dnormMix(x, mean1 = -2, sd1 = .5, mean2 = 1, sd2 = .9, p.mix = 0.5)}, 
                xlim = c(-2,-1.6),
                geom = "area",
                fill = "violet", 
                alpha = .5) +
  stat_function(fun = function(x){dnormMix(x, mean1 = -2, sd1 = .5, mean2 = 1, sd2 = .9, p.mix = 0.5)}, 
                xlim = c(-1.6,-1.14),
                geom = "area",
                fill = "blue", 
                alpha = .4) +
  stat_function(fun = function(x){dnormMix(x, mean1 = -2, sd1 = .5, mean2 = 1, sd2 = .9, p.mix = 0.5)}, 
                xlim = c(-1.14,-0.72),
                geom = "area",
                fill = "blue", 
                alpha = .4) +
  stat_function(fun = function(x){dnormMix(x, mean1 = -2, sd1 = .5, mean2 = 1, sd2 = .9, p.mix = 0.5)}, 
                xlim = c(-0.72, 0.5),
                geom = "area",
                fill = "cyan", 
                alpha = .5) +
  xlab(" ") + ylab("PDF") + 
  scale_x_continuous(breaks=c(-4,-3,-2,-1.6, -0.72, 0.5, 1,2, 3, 4),
                     labels=c(-4,-3,expression(x[1]*" = -2"),expression(x[2]*" = -1.6"), expression(x[3]*" = -0.72"), expression(x[4]*" = 0.5"), 1,2, 3, 4)) + theme_bw()



# Figure 4: Plot all scenarios ------------------------------------------------------
load("Sim Results/Simulated Scenarios/Def_scenarios.RData")
library(KScorrect)

# graphical parameters
op <- par(mfrow = c(4,4),
          oma = c(1,2,2,1) + 0.1,
          mar = c(1,2,2,1) + 0.1,
          bty = 'n')

cols <- hcl.colors(10, "YlOrRd")

contour(mvd_def_sA1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), main = "Unimodal: \nStudent-t marginals", 
        col = cols)
mtext(side = 2, bquote("Gaussian copula ("*theta*" = 0.7)"), line = 2, cex = .9)
contour(mvd_def_sA1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sA1t, add = T)
legend("bottomright", legend = "S1", bty='n', cex = 1.5)

contour(mvd_def_sA1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), main = "Unimodal: \nGaussian marginals",
        col = cols)
contour(mvd_def_sA1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sA1, add = T)
legend("bottomright", legend = "S2", bty='n', cex = 1.5)

contour(mvd_def_sA2, dMvdc, xlim=c(-6, 6), ylim=c(-6, 13), main = "Bimodal: \nGaussian & GaussianMixture marginals", col = cols)
contour(mvd_def_sA2, dMvdc, xlim=c(-7, 7), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sA2, add = T)
legend("bottomright", legend = "S3", bty='n', cex = 1.5)

contour(mvd_def_sA3, dMvdc, xlim=c(-5, 15), ylim=c(-4, 15), main = "Quadrimodal: \nGaussianMixture marginals", col = cols)
contour(mvd_def_sA3, dMvdc, xlim=c(-6, 16), ylim=c(-6, 16), col = "blue", lwd = 2,
        levels = true_q_sA3, add = T)
legend("bottomright", legend = "S4", bty='n', cex = 1.5)


contour(mvd_def_sB1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
mtext(side = 2, bquote("Student-t copula ("*theta*" = 0.7; "*nu*" = 6)"), line = 2, cex = .9)
contour(mvd_def_sB1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sB1t, add = T)
legend("bottomright", legend = "S5", bty='n', cex = 1.5)

contour(mvd_def_sB1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
contour(mvd_def_sB1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sB1, add = T)
legend("bottomright", legend = "S6", bty='n', cex = 1.5)

contour(mvd_def_sB2, dMvdc, xlim=c(-6, 6), ylim=c(-6, 13), col = cols)
contour(mvd_def_sB2, dMvdc, xlim=c(-7, 7), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sB2, add = T)
legend("bottomright", legend = "S7", bty='n', cex = 1.5)

contour(mvd_def_sB3, dMvdc, xlim=c(-5, 15), ylim=c(-4, 15), col = cols)
contour(mvd_def_sB3, dMvdc, xlim=c(-6, 16), ylim=c(-5, 16), col = "blue", lwd = 2,
        levels = true_q_sB3, add = T)
legend("bottomright", legend = "S8", bty='n', cex = 1.5)


contour(mvd_def_sC1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
mtext(side = 2, bquote("Frank copula ("*theta*" = 5.75)"), line = 2, cex = .9)
contour(mvd_def_sC1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sC1t, add = T)
legend("bottomright", legend = "S9", bty='n', cex = 1.5)

contour(mvd_def_sC1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
contour(mvd_def_sC1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sC1, add = T)
legend("bottomright", legend = "S10", bty='n', cex = 1.5)

contour(mvd_def_sC2, dMvdc, xlim=c(-6, 6), ylim=c(-6, 13), col = cols)
contour(mvd_def_sC2, dMvdc, xlim=c(-7, 7), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sC2, add = T)
legend("bottomright", legend = "S11", bty='n', cex = 1.5)

contour(mvd_def_sC3, dMvdc, xlim=c(-6, 15), ylim=c(-6, 13), col = cols)
contour(mvd_def_sC3, dMvdc, xlim=c(-7, 16), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sC3, add = T)
legend("bottomright", legend = "S12", bty='n', cex = 1.5)


contour(mvd_def_sD1t, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
mtext(side = 2, bquote("Clayton copula ("*theta*" = 2)"), line = 2, cex = .9)
contour(mvd_def_sD1t, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sD1t, add = T)
legend("bottomright", legend = "S13", bty='n', cex = 1.5)

contour(mvd_def_sD1, dMvdc, xlim=c(-6, 6), ylim=c(-6, 6), col = cols)
contour(mvd_def_sD1, dMvdc, xlim=c(-7, 7), ylim=c(-7, 7), col = "blue", lwd = 2,
        levels = true_q_sD1, add = T)
legend("bottomright", legend = "S14", bty='n', cex = 1.5)

contour(mvd_def_sD2, dMvdc, xlim=c(-6, 5), ylim=c(-6, 13), col = cols)
contour(mvd_def_sD2, dMvdc, xlim=c(-7, 6), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sD2, add = T)
legend("bottomright", legend = "S15", bty='n', cex = 1.5)

contour(mvd_def_sD3, dMvdc, xlim=c(-6, 15), ylim=c(-6, 13), col = cols)
contour(mvd_def_sD3, dMvdc, xlim=c(-7, 16), ylim=c(-7, 14), col = "blue", lwd = 2,
        levels = true_q_sD3, add = T)
legend("bottomright", legend = "S16", bty='n', cex = 1.5)

# Figure 5: Dirichlet scenario ------------------------------------------------------
# load("Sim Results/Simulated Scenarios/Def_scenarios.RData")
source("2 HDR_functions.R")
cols <- hcl.colors(10, "YlOrRd")
library(gtools) #for Dirichlet

myN = 1e7 # this is a big N in order to provide a good approx of the contour and the level that ensures a given coverage prob
coverage_prob = 0.95  # our contour of interest is defined on this coverage probability
myalpha = c(1,1,2)

set.seed(12345)
data_sDir3 = rdirichlet(myN, alpha = myalpha)
fxy_true_sDir = ddirichlet(x = data_sDir3, alpha = myalpha)
true_q_sDir_112 = quantile(fxy_true_sDir, probs = (1-coverage_prob))

# For the contours of the copula density
U1 = pbeta(data_sDir3[,1], myalpha[1], sum(myalpha)-myalpha[1])
U2 = pbeta(data_sDir3[,2], myalpha[2], sum(myalpha)-myalpha[2])
U = cbind(U1, U2)

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


load("/Users/blackmamba/Library/CloudStorage/GoogleDrive-nina.deliu@uniroma1.it/Il mio Drive/Nina - MEMOTEF/Multivariate Statistical Intervals & HDRs/R codes – HDRs/Sim Results/Error Evaluation/Res_sDir.RData")


# Total Error Plot (S17: Dirichlet) ---------------------------------------------------------------

Nsamples = c(50, 100, 250, 500, 1000)

measures = c(expression(M[0]*":KDE"), expression(M[1]*":kNN-Eucl"),
             expression(M[2]*":kNN-CDF"), expression(M[3]*":"*epsilon*"-CDF"),
             expression(M[0]^NPCop*":DE"), expression(M[0]^PCop*":DE"), 
             expression(M[3]^NPCop*":"*epsilon*"-CDF"), expression(M[3]^PCop*":"*epsilon*"-CDF"))

measures1 = c("M0:KDE", "M1:kNN-Eucl", "M2:kNN-CDF", "M3:Ɛ-CDF", 
              "MCop0:DE-NP","MCop0:DE-P", "MCop3:CDF-NP","MCop3:CDF-P")

metrics = c("Total Error Rate", "FPR", "FNR", "Accuracy", "F1 Score", "MCC")

n_metrics = length(metrics); n_methods = length(measures); n_ss = length(Nsamples)

# |-- Total Error -------
mydf_TotErr = data.frame("meanEr" = c(ErrMean$'50'[4,c(2,5, 6,7,3,4,8,9)],ErrMean$'100'[4,c(2,5, 6,7,3,4,8,9)], ErrMean$'250'[4,c(2,5, 6,7,3,4,8,9)],
                                      ErrMean$'500'[4,c(2,5, 6,7,3,4,8,9)], ErrMean$'1000'[4,c(2,5, 6,7,3,4,8,9)]), 
                         "SDEr" = c(ErrSD$'50'[4,c(2,5, 6,7,3,4,8,9)],ErrSD$'100'[4,c(2,5, 6,7,3,4,8,9)], ErrSD$'250'[4,c(2,5, 6,7,3,4,8,9)],
                                    ErrSD$'500'[4,c(2,5, 6,7,3,4,8,9)], ErrSD$'1000'[4,c(2,5, 6,7,3,4,8,9)]),
                         "method" = rep(measures1, n_ss),
                         "ss" = c(rep(Nsamples[1],n_methods), rep(Nsamples[2],n_methods), rep(Nsamples[3],n_methods),
                                  rep(Nsamples[4],n_methods), rep(Nsamples[5],n_methods)))

mydf_TotErr$ss <- as.factor(mydf_TotErr$ss)
mydf_TotErr$method <- as.factor(mydf_TotErr$method)
summary(mydf_TotErr)


plot_ERR = ggplot(mydf_TotErr, aes(ss, meanEr)) +
  geom_errorbar(
    aes(ymin = meanEr-SDEr, ymax = meanEr+SDEr, color = method),
    position = position_dodge(0.3), width = 0.2
  )+
  geom_point(aes(color = method, shape = method), position = position_dodge(0.3)) +
  geom_line(aes(linetype = method, group = method, color = method), position = position_dodge(0.3))+
  scale_linetype_manual(" ", values = 1:8, breaks=measures1, labels=measures, name = " ")+
  scale_shape_manual(values = 1:8, breaks = measures1, labels = measures, name = " ") +
  scale_colour_hue(breaks = measures1, labels = measures, name = " ") + 
  xlab('Sample Size (n)') + #ylim(0.2, 0.5) + 
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(color = 'grey85'),
        panel.grid.minor = element_line(color = 'grey91')) +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(size=9)) +
  ggtitle("Accuracy") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
  theme(legend.position="bottom",
        legend.title=element_blank()) #+
#scale_x_continuous(sec.axis=sec_axis(trans=~ . * 1, name="Batch (T)"))

plot_ERR 

# |-- MCC -------
mydf_MCC = data.frame("meanEr" = c(ErrMean$'50'[6,c(2,5, 6,7,3,4,8,9)],ErrMean$'100'[6,c(2,5, 6,7,3,4,8,9)], ErrMean$'250'[6,c(2,5, 6,7,3,4,8,9)],
                                   ErrMean$'500'[6,c(2,5, 6,7,3,4,8,9)], ErrMean$'1000'[6,c(2,5, 6,7,3,4,8,9)]), 
                      "SDEr" = c(ErrSD$'50'[6,c(2,5, 6,7,3,4,8,9)],ErrSD$'100'[6,c(2,5, 6,7,3,4,8,9)], ErrSD$'250'[6,c(2,5, 6,7,3,4,8,9)],
                                 ErrSD$'500'[6,c(2,5, 6,7,3,4,8,9)], ErrSD$'1000'[6,c(2,5, 6,7,3,4,8,9)]),
                      "method" = rep(measures1, n_ss),
                      "ss" = c(rep(Nsamples[1],n_methods), rep(Nsamples[2],n_methods), rep(Nsamples[3],n_methods),
                               rep(Nsamples[4],n_methods), rep(Nsamples[5],n_methods)))

mydf_MCC$ss <- as.factor(mydf_MCC$ss)
mydf_MCC$method <- as.factor(mydf_MCC$method)
summary(mydf_MCC)


plot_MCC = ggplot(mydf_MCC, aes(ss, meanEr)) +
  geom_errorbar(
    aes(ymin = meanEr-SDEr, ymax = meanEr+SDEr, color = method),
    position = position_dodge(0.3), width = 0.2
  )+
  geom_point(aes(color = method, shape = method), position = position_dodge(0.3)) +
  geom_line(aes(linetype = method, group = method, color = method), position = position_dodge(0.3))+
  scale_linetype_manual(" ", values = 1:8, breaks=measures1, labels=measures, name = " ")+
  scale_shape_manual(values = 1:8, breaks = measures1, labels = measures, name = " ") +
  scale_colour_hue(breaks = measures1, labels = measures, name = " ") + 
  xlab('Sample Size (n)') + #ylim(0.2, 0.5) + 
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(color = 'grey85'),
        panel.grid.minor = element_line(color = 'grey91')) +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(size=9)) +
  ggtitle("Matthews Correlation Coefficient (MCC)") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
  theme(legend.position="bottom",
        legend.title=element_blank()) #+
#scale_x_continuous(sec.axis=sec_axis(trans=~ . * 1, name="Batch (T)"))

plot_MCC

# |-- False Positive Rate -------
mydf_FPR = data.frame("meanEr" = c(ErrMean$'50'[2,c(2,5, 6,7,3,4,8,9)],ErrMean$'100'[2,c(2,5, 6,7,3,4,8,9)], ErrMean$'250'[2,c(2,5, 6,7,3,4,8,9)],
                                   ErrMean$'500'[2,c(2,5, 6,7,3,4,8,9)], ErrMean$'1000'[2,c(2,5, 6,7,3,4,8,9)]), 
                      "SDEr" = c(ErrSD$'50'[2,c(2,5, 6,7,3,4,8,9)],ErrSD$'100'[2,c(2,5, 6,7,3,4,8,9)], ErrSD$'250'[2,c(2,5, 6,7,3,4,8,9)],
                                 ErrSD$'500'[2,c(2,5, 6,7,3,4,8,9)], ErrSD$'1000'[2,c(2,5, 6,7,3,4,8,9)]),
                      "method" = rep(measures1, n_ss),
                      "ss" = c(rep(Nsamples[1],n_methods), rep(Nsamples[2],n_methods), rep(Nsamples[3],n_methods),
                               rep(Nsamples[4],n_methods), rep(Nsamples[5],n_methods)))

mydf_FPR$ss <- as.factor(mydf_FPR$ss)
mydf_FPR$method <- as.factor(mydf_FPR$method)
summary(mydf_FPR)


plot_FPR = ggplot(mydf_FPR, aes(ss, meanEr)) +
  geom_errorbar(
    aes(ymin = meanEr-SDEr, ymax = meanEr+SDEr, color = method),
    position = position_dodge(0.3), width = 0.2
  )+
  geom_point(aes(color = method, shape = method), position = position_dodge(0.3)) +
  geom_line(aes(linetype = method, group = method, color = method), position = position_dodge(0.3))+
  scale_linetype_manual(" ", values = 1:8, breaks=measures1, labels=measures, name = " ")+
  scale_shape_manual(values = 1:8, breaks = measures1, labels = measures, name = " ") +
  scale_colour_hue(breaks = measures1, labels = measures, name = " ") + 
  xlab('Sample Size (n)') + #ylim(0.2, 0.5) + 
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(color = 'grey85'),
        panel.grid.minor = element_line(color = 'grey91')) +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(size=9)) +
  ggtitle("False Positive Rate (FPR)") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
  theme(legend.position="bottom",
        legend.title=element_blank()) #+
#scale_x_continuous(sec.axis=sec_axis(trans=~ . * 1, name="Batch (T)"))

plot_FPR

# |-- False Negative Rate -------
mydf_FNR = data.frame("meanEr" = c(ErrMean$'50'[3,c(2,5, 6,7,3,4,8,9)],ErrMean$'100'[3,c(2,5, 6,7,3,4,8,9)], ErrMean$'250'[3,c(2,5, 6,7,3,4,8,9)],
                                   ErrMean$'500'[3,c(2,5, 6,7,3,4,8,9)], ErrMean$'1000'[3,c(2,5, 6,7,3,4,8,9)]), 
                      "SDEr" = c(ErrSD$'50'[3,c(2,5, 6,7,3,4,8,9)],ErrSD$'100'[3,c(2,5, 6,7,3,4,8,9)], ErrSD$'250'[3,c(2,5, 6,7,3,4,8,9)],
                                 ErrSD$'500'[3,c(2,5, 6,7,3,4,8,9)], ErrSD$'1000'[3,c(2,5, 6,7,3,4,8,9)]),
                      "method" = rep(measures1, n_ss),
                      "ss" = c(rep(Nsamples[1],n_methods), rep(Nsamples[2],n_methods), rep(Nsamples[3],n_methods),
                               rep(Nsamples[4],n_methods), rep(Nsamples[5],n_methods)))

mydf_FNR$ss <- as.factor(mydf_FNR$ss)
mydf_FNR$method <- as.factor(mydf_FNR$method)
summary(mydf_FNR)


plot_FNR = ggplot(mydf_FNR, aes(ss, meanEr)) +
  geom_errorbar(
    aes(ymin = meanEr-SDEr, ymax = meanEr+SDEr, color = method),
    position = position_dodge(0.3), width = 0.2
  )+
  geom_point(aes(color = method, shape = method), position = position_dodge(0.3)) +
  geom_line(aes(linetype = method, group = method, color = method), position = position_dodge(0.3))+
  scale_linetype_manual(" ", values = 1:8, breaks=measures1, labels=measures, name = " ")+
  scale_shape_manual(values = 1:8, breaks = measures1, labels = measures, name = " ") +
  scale_colour_hue(breaks = measures1, labels = measures, name = " ") + 
  xlab('Sample Size (n)') + #ylim(0.2, 0.5) + 
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(color = 'grey85'),
        panel.grid.minor = element_line(color = 'grey91')) +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(size=9)) +
  ggtitle("False Negative Rate (FNR)") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
  theme(legend.position="bottom",
        legend.title=element_blank()) #+
#scale_x_continuous(sec.axis=sec_axis(trans=~ . * 1, name="Batch (T)"))

plot_FNR

Combined_Dir = ggarrange(plot_ERR, plot_MCC, plot_FPR, plot_FNR, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")
annotate_figure(Combined_Dir, top = text_grob("S17: Dirichlet(1,1,2)", face = "bold", size = 14))
