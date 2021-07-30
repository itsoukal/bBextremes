#' Distribution of extremes (maxima) of a stochastic process with AR(1) autocorrelation structure and any marginal distribution
#'
#' @param u vector of probabilities
#' @param dist quantile function (character, e.g., qgamma)
#' @param params a list with distribution's CDF (as function) and parameters
#' @param rhoz A scalar, with the lag-1 equivelant (i.e., Gaussian) correlation coefficient
#' @param k duration (e.g., 365) or block length
#'
#' @return vector of probabilities
#' @export
#'
#' @examples
pAR1=function(q, dist, params, rhoz, k){
  #rhoz: lag-1 correlation
  #k = duration (e.g., 365)
  u=pFX(q = q, FX = dist, PFX = params)

  norm.cop <- normalCopula(dim = 2, rhoz, dispstr = 'toep')
  p=rep(NA, length(u))
  for (i in 1:length(u)) {
    H2=pCopula(rep(u[i], 2), norm.cop)
    FX=u[i]
    p[i]=FX*(H2/FX)^(k-1)
    # p[i]=(FX^2/H2)*(H2/FX)^k # identical
  }

  return(p)
}
