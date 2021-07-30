#' ACF of a dichotomized process
#'
#' @param PD A scalar (0<=PD<1), denoting probability level of the threshold at which dichotomization is applied.
#' @param ACSn A vector, with the equivelant (i.e., Gaussian) autocorrelation structure.
#' @param fast Boolean, TRUE or FALSE, indicating whether to use or not an approximation of the bivariate Gaussian copula (default valueis FALSE)
#'
#' @return A vector (with length equal to ACSn) with the ACF of the dichotomized process.
#' @export
#'
#' @examples
#'  # Define an AR(1) ACF, with rho=0.6, up to lag 100.
#'  ACS=0.6^(0:100)
#'
#'  # Estimate the over-threshold ACF using the "fast", and the "slow" method.
#'  acfOTfast=acfOT(PD = 0.7, ACSn = ACS, fast = 1)
#'  acfOT=acfOT(PD = 0.7, ACSn = ACS, fast = 0)
#'
#'  # Compare the fast", and the "slow" method.
#'  plot(acfOTfast, acfOT)
#'  abline(0,1,col='red')
#'
acfOT=function(PD, ACSn, fast=0) {

  PW=1-PD
  maxlag=length(ACSn)-1
  ACSn=ACSn
  ACFbin=rep(NA, maxlag)

  for (i in 1:maxlag) {

    if (fast==0) {
      pp=1- 2*pnorm(qnorm(PD)) + mvtnorm::pmvnorm(upper=rep(qnorm(PD), 2), lower = rep(-Inf, 2),
                                                  sigma = matrix(c(1, ACSn[(i+1)], ACSn[(i+1)], 1), ncol=2, byrow = TRUE) )
      ACFbin[i]=(pp[1]-PW^2)/(PW*PD)
    } else if (fast==1){
      pp=1- 2*pnorm(qnorm(PD)) + pbivnorm(q = qnorm(PD), rho = ACSn[(i+1)])
      # mvtnorm::pmvnorm(upper=rep(qnorm(1-PD),2), lower = rep(-Inf,2),
      #                                           sigma = matrix(c(1,ACSn[(i+1)],ACSn[(i+1)],1),ncol=2) )
      ACFbin[i]=(pp[1]-PW^2)/(PW*PD)
    } else {
      # notice the min(PD, 1-PD)
      ACFbin=Gaussian2binary(P = min(PD, 1-PD), ACSn)[-1]
    }
    # # This is the acf of the binary sequence
    # ACFbin[i]=(pp[1]-PW^2)/(PW-PW^2)


  }
  return(ACFbin)
}


# library(anySim)
#
# ACS=acsHurst(H=0.7, lag=500, var=1)
# PD=0.8
# PW=1-PD
# fx='qbinom'
# pfx=list(size=1, prob=PW)
#
# ARTApar=EstARTAp(ACF=ACS, maxlag=0, dist=fx, params=pfx,
#                  NatafIntMethod ='Int', NoEval=9, polydeg=0)
#
# Sim=SimARTAp(ARTApar = ARTApar, burn = 1000, steps = 10^5, stand = 0)
# acf(Sim$X)
# lines(0:(length(ACS)-1), ACS)
# plot(Sim$X[1:1000], type='l', col='red')
#
# ACSz=ARTApar$ACFn
#
#
# PD=0.9999
# # Estimate the over-threshold ACF (extremogram) using the "fast", and the "slow" method.
# acfOTfast1=acfOT(PD = PD, ACSn = ACSz, fast = 1)
# acfOTfast2=acfOT(PD =  PD, ACSn = ACSz, fast = 2)
# acfOTv=acfOT(PD = PD, ACSn = ACSz, fast = 0)
#
# # Compare the fast", and the "slow" method.
# plot(acfOTfast1, acfOTv)
# abline(0,1,col='red')
#
# plot(acfOTfast2[-1], acfOTv)
# abline(0,1,col='red')
#
#
#
# acf(Sim$X)
# lines(0:(length(ACS)-1), ACS)
# lines(0:(length(ACS)-1), c(1, acfOTfast), col='red')
# lines(0:(length(ACS)-1), c(1, acfOTv), col='blue')
#
# Nataf=NatafGH(rho = ACS, fx = 'qbinom', fy = 'qbinom',
#               paramlistfx=list(size=1, prob=1-PD),
#               paramlistfy=list(size=1, prob=1-PD), nodes = 151)
# plot(ACS, Nataf, t='p')
# abline(0,1,col='red')
#
# plot(ACS, c(0, acfOTv), t='p')
# abline(0,1,col='red')
# natafbin= NatafBin2(P = min(PD, 1-PD), ACSz = ACS) # See Nataf_Bernoulli2.R
# # plot(Nataf, ACS, t='p')
# lines(ACS,natafbin, col='blue')
#
# lines(ACS[-1], acfOTv, col='green')
#
# plot(Nataf[-1], natafbin[-1])
# abline(0,1,col='red')
#
#
# # compare alternative ways to estimate RhoBB
# acs=acs[1:100]
# bBextremes::ICcorr(acs, SV = 1 )
#
# CC=toeplitz(acs)
# SS=0
# for (i in 1:100) {
#   for (j in 1:100) {
#     SS=SS+CC[i,j]
#   }
# }
# (SS-100)/(100*(100-1))
#
# ONES=(rep(1,100))
# (t(ONES)%*%CC%*%ONES - 100)/(100*(100-1))
#
# X=rgamma(1000000, shape=0.5, scale=10)
# xu=1; p=1-pgamma(q = xu, shape=0.5, scale=10); p
#
# Bin=ifelse(X>xp, 1, 0)
# mean(Bin)

