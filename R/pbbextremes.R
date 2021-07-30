#' Distribution of extremes (maxima) of a stochastic process with any autocorrelation structure (ACS)
#'
#' @param u vector of probabilities
#' @param ACSn A vector, with the equivelant (i.e., Gaussian) autocorrelation structure.
#' @param fast Boolean, TRUE or FALSE, indicating whether to use or not an approximation of the bivariate Gaussian copula (default valueis FALSE)
#' @description Approximation based on the beta-Binomial distribution. For more details see, Serinaldi F, Lombardo F, Kilsby CG (2020) All in order: Distribution of serially correlated order statistics with applications to hydrological extremes. Adv Water Resour 144:103686. doi: 10.1016/j.advwatres.2020.103686.
#' @return vector of probabilities
#' @export
#'
#' @examples
pbbextremes=function(q, dist, params, ACSn= ARTApar$ACSn, fast=0){

  u=pFX(q = q, FX = dist, PFX = params)
  k=length(ACSn)
  rbbk=Pnew=rep(NA, length(u))

  for (i in 1:length(u) ) {
    acf01=c(1, acfOT(PD = u[i], ACSn = ACSn, fast = fast))
    # Cmat=toeplitz(acf01)
    # rbb=(sum(Cmat)-k)/(k*(k-1))
    # rbbk[i]=rbb

    rbbk[i]=sum(SuperGauss::toep.mult(acf = acf01, rep(1,k)) -1)/(k*(k-1))

    Pnew[i]=VGAM::pbetabinom(0, size = k, prob = (1-u[i]), rho =rbbk[i])
  }
  return(list(Pnew=Pnew, rbbk=rbbk))
}
