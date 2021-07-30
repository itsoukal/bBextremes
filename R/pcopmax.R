#' Distribution of extremes (maxima) of a stochastic process with any autocorrelation structure (ACS) using Copulas
#'
#' @param u vector of probabilities
#' @param ACSn A vector, with the equivelant (i.e., Gaussian) autocorrelation structure.
#'
#' @return vector of probabilities
#' @export
#'
#' @examples
pcopmax=function(q, dist, params, ACSn= ARTApar$ACSn) {
  library(copula)
  k=length(ACSn)

  u=pFX(q = q, FX = dist, PFX = params)
  Pcop=rep(NA, length(u))

  for (i in 1:length(u) ) {
    norm.cop <- normalCopula(dim = k, param = ACSn[-1], dispstr = 'toep')
    Pcop[i]=pCopula(u = rep(u[i], k), copula = norm.cop)
    # Pcop[i]=mvtnorm::pmvnorm(upper = rep(z[i], k), mean = rep(0, k), sigma = toeplitz(ARTApar$ACFn))[1]
    print(paste0('Probability ',i,'/',length(u)))
  }
  return(Pcop)
}
