#' Estimate the Intra-cluster (IC) correlation coefficient
#'
#' @param ACFn A vector with the autocorrelation structure (length kmax<=m)
#' @param kmax Maximum scale (kmax<=m)
#' @param SV Boolean (TRUE or default: FALSE) indicating if IC will be calculated in the intermidiate time scales
#'
#' @return The Intra-cluster (IC) correlation coefficient
#' @export
#'
#' @examples
#'
#' k=365
#' ACF=anySim::acsCAS(param = c(3, 0.1), lag = k-1)
#' sum(SuperGauss::toep.mult(acf = ACF, rep(1,k)) -1)/(k*(k-1))
#' (sum(toeplitz(ACF))-k)/(k*(k-1))
#' pracma::tic();ICcorr(ACFn = ACF, kmax = k, SV = 1);pracma::toc()
#'
ICcorr=function(ACFn=acsHurst(H = 0.6, lag=99), kmax=100, SV=0) {
  m=length(ACFn)
  if (kmax>m){
    print("Warning: kmax larger than maximum possible estimation from ACS")
    break()
  }

  rbbK=rep(NA, kmax)
  if (SV==0) {
  temp=toeplitz(ACFn[1:kmax])
  }

  if (SV==0) {
    for (k in 2:kmax) {
      Cmat=temp[1:k, 1:k]
      rbb=(sum(Cmat)-k)/(k*(k-1))
      rbbK[k]=rbb
    }
  }else {
    # Cmat=temp[1:kmax, 1:kmax]
    # rbb=(sum(Cmat)-kmax)/(kmax*(kmax-1))
    # rbbK=rbb
    rbbK=sum(SuperGauss::toep.mult(acf = ACFn, rep(1,m)) -1)/(m*(m-1))
  }
  return(rbbK)
}
