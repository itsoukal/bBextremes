#' Cumulative distribution funtcion of the distribution of k-length block maxima under the assumption of temporal independence
#'
#' @param q  a scalar denoting the return levels
#' @param FX the cumulative distribution function of distribution of the parent process
#' @param k  a scalar denoting the length of block (e.g., 365 for 1 year, in the case of a daily parent process)
#' @param asRP a boolean indicating whether or not the function will return probabilities or return periods
#' @param ... additional arguments (parameter) that regards the distribution of the parent process
#'
#' @return vector of probabilities (or return periods, if asRP=1)
#' @export
#'
#' @examples
pmaxiid=function(q, FX, k, asRP=0, ...) {
  # u=pburrDK2(q = q, scale=params$scale, shape1=params$shape1, shape2 = params$shape2, PW=params$PW)
  u=FX(q=q, ...)
  p=u^k

  if (asRP==1) {
    p=1/(1-u^k)
  }
  return(p)
}


#' Quantile function (Inverse cumulative distribution funtcion) of the distribution of k-length block maxima under the assumption of temporal independence
#'
#' @param Rp a scalar denoting the return period (in years)
#' @param QFX the quantile function of distribution of the parent process
#' @param k a scalar denoting the length of block (e.g., 365 for 1 year, in the case of a daily parent process)
#' @param ... additional arguments (parameter) that regards the distribution of the parent process
#'
#' @return vector of return levels
#' @export
#'
#' @examples
qmaxiid=function(Rp, QFX, k, ...) {
  # Nr=number of Wet events
  p=(1-(1/Rp))^(1/k)
  # q=qburrDK2(p = p, scale=params$scale, shape1=params$shape1, shape2 = params$shape2, PW=params$PW)
  q=QFX(p=p, ...)
  return(q)
}

# source('C:/Users/jtsou/Dropbox/R_Home/Extremes/bBextremes/R/ExtremeCompleteFunctions.R')
#
# k=360*24
# n=10^3
# params=list(scale=1, shape1=0.5, shape2 = 0.02, PW=0.3)
# x=rburrDK2(n = k*n, scale=params$scale, shape1=params$shape1, shape2 = params$shape2, PW=params$PW)
#
# scales=c(1,2,3,6,12,24, 24*2, 24*3, 24*5, 24*10, 24*30)
# kk=k/scales
#
# AM1=blockmaxxer(x = x, k = kk[1], n = length(x)/kk[1], fun = max)
# Rp=1/(1-ppoints(AM1, a=0))
# plot(Rp, sort(AM1), log='xy', col='orange', pch=19, ylim=c(1, 300))
# q=seq(1, 300)
# RpTheor=piidmaxDK2(q = q, params = params, k = k, asRP = T)
# lines(RpTheor, q, col='red')
# q=qiidmaxDK2(Rp = RpTheor, params = params, k = k);q
# lines(RpTheor, q, col='blue')
#
#
# RpCompEmp=round(Tcomplete_emp(x = x, Nr = k, PW = params$PW), 3);
# points(RpComp, q, col='black', pch=19)
# RpComp=round(Tcomplete_theor(q = q, params = params, Nr = k), 1)
# lines(RpComp, q, col='green', lwd=2)
#
# AMk=list()
# aggx=list()
# for (i in 2:length(kk)) {
#   temp=SMATAPKG::aggregateTS(x = x, k = scales[i], Div = 1)
#   aggx[[i]]=temp
#   AMk[[i]]=blockmaxxer(x = temp, k = kk[i], n = length(temp)/kk[i], fun = max)
#   Rp=1/(1-ppoints(AMk[[i]], a=0))
#   lines(Rp, sort(AMk[[i]]), log='xy', col='red', pch=19)
# }
#
#
#
#
#
#
