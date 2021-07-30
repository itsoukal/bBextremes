#' cummulative distribution funcion (CDF)
#'
#' @param q vector of quantiles
#' @param FX CDF function (character, e.g., pgamma)
#' @param PFX a list with distribution's CDF (as function) and parameters
#'
#' @return vector of probabilities
#' @export
#'
#' @examples
pFX=function(q, FX, PFX) {
  # FX='pgamma'
  # PFX=list(shape=0.5, scale=1)
  substr(x = FX, start = 1, stop = 1)='p'
  PPFX=PFX
  PPFX$q=q
  U=as.vector(do.call(FX, PPFX))
  return(U)
}

