#' Distribution of extremes (maxima) of a stochastic process with AR(p) autocorrelation structure and any marginal distribution
#'
#' @param u vector of probabilities
#' @param ACSn A vector, with the equivelant (i.e., Gaussian) autocorrelation structure.
#' @param k duration (e.g., 365) or block length
#' @param p AR(p) model order
#' @param fast Boolean, TRUE or FALSE, indicating whether to use or not an approximation of the bivariate Gaussian copula (default valueis FALSE)
#'
#' @return vector of probabilities
#' @export
#'
#' @examples
pARp=function(q, dist, params, ACSn, k, p, fast=0) {
  # k>=p>=2

  u=pFX(q = q, FX = dist, PFX = params)

  probs=rep(NA, length(u))
  for (i in 1:length(u)) {
    A=pbbextremes(q = q[i], dist = dist, params = params, ACSn = ACSn[1:(p+1)], fast)$Pnew
    B=pbbextremes(q = q[i], dist = dist, params = params, ACSn = ACSn[1:(p)], fast)$Pnew
    probs[i]=B*(A/B)^(k-p)
  }
  return(probs)
}
