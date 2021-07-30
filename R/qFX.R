#' Quantile function (Inverse cumulative distribution function)
#'
#' @param p vector of probabilities
#' @param FX quantile function (character, e.g., qgamma)
#' @param PFX a list with distribution's ICDF (as function) and parameters
#'
#' @return vector of quantiles
#' @export
#'
#' @examples
qFX=function(p=ppoints(1000, a=0), FX, PFX) {
  # FX='qgamma'
  # PFX=list(shape=0.5, scale=1)

  PPFX=PFX
  PPFX$p=p
  X=as.vector(do.call(FX, PPFX))
  return(X)
}
