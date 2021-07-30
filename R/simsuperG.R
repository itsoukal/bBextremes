#' Fast simulation of stationary univariate processes with any marginal distribution and autocorrelation structure
#'
#' @param n number of realizations
#' @param ACF A vector with the target autocorrelation structure (including lag-0, i.e., 1). The length of the ACF vector should be equal to that of the time series to be generated.
#' @param dist A string indicating the quantile function of the target marginal distribution (i.e., the ICDF).
#' @param params 	A named list with the parameters of the target distribution.
#' @param NatafIntMethod A string ("GH", "Int", or "MC"), indicating the intergation method, to resolve the Nataf integral. Currently only "GH" is implemented.
#' @param NoEval 	A scalar indicating (default: 9) the number of evaluation points for the integration methods.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then another curve is fitted.
#' @param stand Boolean indicator (TRUE or default: FALSE) whether or not standardization is performed.
#'
#' @return A list with 5 elements. Two (2) that concern the Nataf transofrmation (i.e., nataf and acfz), and three (3) that regard the generated time series (in vector format): X: The final time series at the actual domain with the target marginal distribution and correlation structure; Z: The auxiliary Gaussian time series at the Gaussian domain and; U: The auxiliary uniform time series at the Copula domain (i.e., in [0,1]).
#' @export
#'
#' @examples
simsuperG=function(n=1, ACF, dist, params, NatafIntMethod = "GH",
                   NoEval = 9, polydeg = 0, stand=0) {

  nataf=anySim::NatafInvD(targetrho = ACF,
                          fx = dist, fy = dist,
                          paramlistfx = params, paramlistfy = params,
                          NatafIntMethod = NatafIntMethod, NoEval = NoEval, polydeg = polydeg)

  acfz=c(1, nataf$rzEq[-1])
  nataf$rzEq=acfz
  Z=SuperGauss::rSnorm(n = n, acf = acfz[1:(length(acfz)-1)], fft = FALSE)

  if (stand==0) {
    U=pnorm(q = Z, mean = mean(Z), sd=sd(Z))
  } else {
    U=pnorm(q = Z)
  }

  X=qFX(p = U, FX = dist, params)
  sim=list()

  sim$nataf=nataf
  sim$acfz=acfz

  sim$Z=Z
  sim$U=U
  sim$X=X

  return(sim)
}

# N=6*24*365*100
#
# dist='qzi'
# params=list(p0=0.95, Distr=qgamma, shape=0.1, scale=2)
#
# ACF <- anySim::acsCAS(param = c(2, 0.1), N)
#
# sim=simsuperG(n = 1, ACF = ACF, dist = dist, params = params)
# X=sim$X
#
# source("~/Dropbox/R_Home/Extremes/blockmaxxer.R")
#
# AM10=blockmaxxer(x = sim$X, k = 6*24*365, n = 100, fun = max)
# plot(1/(1-ppoints(AM10, a=0)), sort(AM10), log='xy')
#
# min=6
# d=c(1,2,3,5,10,30,60, 90)
# scales=c(min, min*24*d)
#
# dfagg=list()
# for (i in 1:length(scales) ) {
#   dfagg[[i]]=SMATAPKG::aggregateTS(x = X, k = scales[i], Div = 0)
# }
# names(dfagg)=c('1h', '1d', '2d', '3d', '5d,', '10d', '30d', '60d', '90d')
# unlist(lapply(dfagg, SMATAPKG::pdry))
#
#
# source("~/Dropbox/R_Home/Extremes/pbbextremes.R")
# source("~/Dropbox/R_Home/Extremes/acfPOT.R")
# source("~/Dropbox/R_Home/Extremes/pbivnorm.R")
#
#
# ACS=anySim::acsCAS(param = c(3.1, 0.01), lag =  (max(scales)-1))
# PD=0.95
# pkbb=pkcop=pkAR1=rep(NA, length(scales))
# for (i in 1:length(scales)) {
#   pkbb[i]=pbbextremes(u = PD, ACFn = ACS[1:scales[i]])$Pnew
#
# }
#
# plot(c(1, scales), c(0.95, unlist(lapply(dfagg, SMATAPKG::pdry))), t='o', lwd=2, col='red', log='xy', ylim=c(0.001,1))
# points(scales, pkbb, lwd=2, col='black')
#
# ACS=anySim::acsCAS(param = c(3.1, 0.02), lag =  (max(scales)-1))
# plot(0:10000, ACS[1:10001])
# lines(0:10000, sim$acfz[1:10001], col='red')
#
# z=seq(0, 5, by = 0.1)
# u=pnorm(z)
# x=qFX(p = u, FX = dist, PFX = list(p0=0.95, Distr=qgamma, shape=0.1, scale=2))
# ubb=pbbextremes(u = u, ACFn = anySim::acsCAS(param = c(4.8, 0.0165), lag = 6*24*365), fast = 1)$Pnew
#
#
# plot(1/(1-ubb), x, log='xy',t='l', xlim=c(1,100))
# points(1/(1-ppoints(AM10, a=0)), sort(AM10), pch=19, col='orange')
