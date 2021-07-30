#' Fast estimation of equivelant correlation coefficients for Bernoulli processes
#'
#' @param P parameter of Bernoulli process
#' @param ACSn A vector, with the equivelant (i.e., Gaussian) autocorrelation structure.
#'
#' @return vector of equivelant correlation coefficients
#' @export
#'
#' @examples
Gaussian2binary=function(P, ACSn) {
  # min(PD, 1-PD)
  a=0.2567+1.3821*P^(0.3035+0.1287*P)*exp(-0.4837*P)
  b=1.9718+84.2456*P^(1.5883+8.6981*P)*exp(-21.1004*P)
  ACSb=(1 - (1 - ACSz)^(1/b) )^(1/a)
  return(ACSb)
}
