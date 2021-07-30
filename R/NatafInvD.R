#' @title Direct estimation of equivalent correlation coefficients.
#'
#' @description Direct estimation of equivalent correlation coefficients (i.e., in the Gaussian domain).
#'
#' @param targetrho A scalar or vector of target correlation coefficients.
#' @param fx A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param fy A string indicating the quantile function of the distribution (i.e., the ICDF).
#' @param paramlistfx A named list with the parameters of the distribution.
#' @param paramlistfy A named list with parameters of the distribution.
#' @param NoEval A scalar indicating (default: 9) the number of evaluation points for the integration methods.
#' @param NatafIntMethod A string ("GH", "Int", or "MC"), indicating the integration method, to resolve the Nataf integral.
#' @param polydeg A scalar indicating the order of the fitted polynomial. If polydeg=0, then an alternative two-parameter curve is fitted (see, Papalexiou, 2018). Alternative parametric curves are used if polydeg={-1, -2, -3, -4}. Note that polydeg={-2, -4} fit three-parameter curves.
#' @param ... Additional named arguments for the selected "NatafIntMethod" method.
#'
#' @note Avoid the use of the "GH" method (i.e., NatafIntMethod='GH'), when the marginals are discrete.
#' @return A named list with two elements:
#' dfnataf: A dataframe that contains the pairs of Gaussian and resulting correlation coefficients, upon which the curve (polynomial or other) was fitted.
#' rzEq: A vector with the equivalent correlation coefficients, that result into the target ones (i.e., targetrho).
#' rxmax: A scalar, indicating the maximum possible correlation coefficient admissable by the Nataf model.
#' fitpars: (Optional) In case of polydeg<=0, it returns a vector containing the optimum parameters.
#'
#' @export
#'
#' @examples
#' ## The case of two identrical zero-inflated (i.e., mixed) distributions,
#' ## with p0=0.9 a Gamma distribution
#' ## for the continuous part with shape=0.1 and scale=1.
#'\dontrun{
#' fx=fy='qzi'
#' pfx=pfy=list(Distr=qgamma, p0=0.9, shape=0.1, scale=1)
#' rhoz=seq(from=0, to=1 , length.out = 21)
#' rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhox,rhoz, col='red', pch=19); abline(0,1)
#' rhotarget=seq(from=0.0001, to=0.9999 , length.out = 210)
#' req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'GH', polydeg=8, NoEval = 9)
#' points(rhotarget, req, col='blue', pch=17);
#'
#' req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'GH', polydeg=8, NoEval = 9)
#' points(rhotarget, req, col='green')
#'}
#' ## The case with identical Bernoulli distributions, with size=1 and prob=0.2.
#'
#'\dontrun{
#' fx=fy='qbinom'
#' pfx=pfy=list(size=1, prob=0.2)
#' rhoz=seq(from=0, to=1 , length.out = 21)
#' rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
#' plot(rhox,rhoz, col='red', pch=19); abline(0,1)
#' rhotarget=seq(from=0.0001, to=0.9999 , length.out = 210)
#' req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'Int', polydeg=8, NoEval = 9)$rzEq
#' points(rhotarget, req, col='blue', pch=17);
#'
#' req2=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#'               NatafIntMethod = 'Int', polydeg=8, NoEval = 9)$rzEq
#' points(rhotarget, req2, col='green')
#'}
#'
NatafInvD=function(targetrho, fx, fy, paramlistfx, paramlistfy, NatafIntMethod='GH', NoEval=19, polydeg=8, ...){

  rmax=ifelse(max(targetrho)>0,1,0)
  rmin=ifelse(min(targetrho)<0,-1,0)

  rz=seq(rmin, rmax, length.out = NoEval)
  Index0=which(rz==0)

  if (NatafIntMethod=='GH') {
    rx=NatafGH(rho = rz, fx = fx, fy = fy, paramlistfx = paramlistfx, paramlistfy = paramlistfy, ...)
  } else if (NatafIntMethod=='MC') {
    rx=NatafMC(rho = rz, fx = fx, fy = fy, paramlistfx = paramlistfx, paramlistfy = paramlistfy, ...)
  } else if (NatafIntMethod=='Int') {
    rx=NatafInt(rho = rz, fx = fx, fy = fy, paramlistfx = paramlistfx, paramlistfy = paramlistfy)
  } else {
    print('Select an appropriate method to resolve Nataf integral. That is: "GH, "MC" or "Int"')
  }
  rx[Index0]=0
  dfnataf=data.frame(rz=rz, rx=rx)


  if (polydeg>0){
    lb=min(rz);
    ub=max(rz)
    x=seq(lb, ub, 0.0002)
    fitZ2Y=pracma::polyfit(x = rz, y = rx, n = polydeg)
    y=pracma::polyval(fitZ2Y,x)
    X=(cbind(y,x))
    rzEq=approx(X[,1], X[,2], targetrho)$y
  }else if (polydeg==0) {
    rxmax=rx[length(rx)]
    rxmax=ifelse(rxmax>1,1,rxmax)

    DE=DEoptim::DEoptim(fn = TwoParfunObj, lower = c(0.0001,0.0001),upper = c(1000, 100),
                        rx=rx, rz=rz, rmax=rxmax,
                        control = DEoptim.control(trace = 0,NP=30,itermax =1000))

    b=DE$optim$bestmem[1]
    c=DE$optim$bestmem[2]
    targetrho=ifelse(targetrho>rxmax,NA,targetrho)
    rzEq=TwoParfun(b, c, rx=targetrho, rmax = rxmax)

  }else if (polydeg==-1) {

    rxmax=rx[length(rx)]
    rxmax=ifelse(rxmax>1,1,rxmax)

    DE=DEoptim::DEoptim(fn = TwoParfun2Obj, lower = c(0.000001, 0.000001),upper = c(10000, 100),
                        rx=rx, rz=rz, rmax=rxmax,
                        control = DEoptim.control(trace = 0, NP=30, itermax = 1000))

    b=DE$optim$bestmem[1]
    k=DE$optim$bestmem[2]

    targetrho=ifelse(targetrho>rxmax,NA,targetrho)
    rzEq=TwoParfun2(b, k, rx=targetrho, rmax = rxmax)

  }else if (polydeg==-2) {

    rxmax=rx[length(rx)]
    rxmax=ifelse(rxmax>1,1,rxmax)

    DE=DEoptim::DEoptim(fn = ThreeParfunObj, lower = c(0.0001, 0.0001, 0.0001),upper = c(10000, 100, 100),
                        rx=rx, rz=rz, rmax=rxmax,
                        control = DEoptim.control(trace = 0, NP=60, itermax = 1000))

    b=DE$optim$bestmem[1]
    c=DE$optim$bestmem[2]
    k=DE$optim$bestmem[3]

    targetrho=ifelse(targetrho>rxmax,NA,targetrho)
    rzEq=ThreeParfun(b, c, k, rx=targetrho, rmax = rxmax)

  }else if (polydeg==-3) {

    rxmax=rx[length(rx)]
    rxmax=ifelse(rxmax>1,1,rxmax)

    DE=DEoptim::DEoptim(fn = TwoParfun2ObjInv, lower = c(0.0001, 0.0001),upper = c(1000, 100),
                        rx=rx, rz=rz, rmax=rxmax,
                        control = DEoptim.control(trace = 0, NP=60, itermax = 1000))

    b=DE$optim$bestmem[1]
    c=DE$optim$bestmem[2]


    targetrho=ifelse(targetrho>rxmax,NA,targetrho)
    rzEq=TwoParfun2(b, c, rx=targetrho, rmax = rxmax)

  } else if (polydeg==-4) {

    rxmax=rx[length(rx)]
    rxmax=ifelse(rxmax>1,1,rxmax)

    DE=DEoptim::DEoptim(fn = ThreeParfunObjInv, lower = c(0.0001, 0.0001, 0.0001),upper = c(1000, 100, 100),
                        rx=rx, rz=rz, rmax=rxmax,
                        control = DEoptim.control(trace = 0, NP=60, itermax = 1000))

    b=DE$optim$bestmem[1]
    c=DE$optim$bestmem[2]
    k=DE$optim$bestmem[3]

    targetrho=ifelse(targetrho>rxmax,NA,targetrho)
    rzEq=ThreeParfun(b, c, k, rx=targetrho, rmax = rxmax)
  }

  if (polydeg<=0) {
    temp=list('dfnataf'=dfnataf, 'rzEq'=rzEq, 'fitpars'=DE$optim$bestmem, 'rxmax'=rxmax)
  } else {
    temp=list('dfnataf'=dfnataf, 'rzEq'=rzEq, 'rxmax'=rmax)
  }
  return(temp)
}

TwoParfun<-function(b,c,rx,rmax=1){
  A=(-1+((1+b*rx)^(1-c)))
  B=(-1+(1+b*rmax)^(1-c))
  rzhat=A/B
  return(rzhat)
}

TwoParfunInv<-function(b,c,rz,rmax=1){
  rx=(-((1 - ((1 - (1 + b* rmax)^(1 - c))* (1/(1 - (1 + b* rmax)^(1 - c)) - rz))^(1/(1 - c)))/b))
  # TO-DO
  return(rx)
}

TwoParfunObj<-function(par,rz,rx,rmax=1){
  b=par[1];  c=par[2];
  rzhat=TwoParfun(b = b, c = c,rx = rx,rmax = rmax)

  SSE=sum((rz-rzhat)^2)
  SSE=ifelse(is.nan(SSE),10^6,SSE)
  return(SSE)
}

TwoParfun2<-function(b,c,rx,rmax=1){

  # rzhat=(-1 + (1 + b*rx)^(c))/(-1 + (1 + b*rmax)^(c))
  rzhat=(1 - (1 + b*rx)^(c))/(1 - (1 + b*rmax)^(c)) #lomax
  return(rzhat)
}

TwoParfun2Inv<-function(b,c,rz,rmax=1){

  # rx=((-1 + (1 + (-1 + (1 + b* rmax)^c) *rz)^(1/c))/b)
  rx=(-1 + (1 + (-1 + (1 + b* rmax)^c)* rz)^(1/c))/b #lomax
  return(rx)
}

TwoParfun2Obj<-function(par,rz,rx,rmax=1){
  b=par[1]; c=par[2]
  rzhat=TwoParfun2(b = b, c = c, rx = rx, rmax = rmax)

  SSE=sum((rz-rzhat)^2)
  SSE=ifelse(is.nan(SSE),10^6,SSE)
  return(SSE)
}

TwoParfun2ObjInv<-function(par,rz,rx,rmax=1){
  b=par[1]; c=par[2]
  rxhat=TwoParfun2Inv(b = b, c = c, rz = rz, rmax = rmax)

  SSE=sum((rx-rxhat)^2)
  SSE=ifelse(is.nan(SSE),10^6,SSE)
  return(SSE)
}


ThreeParfun<-function(b,c,k,rx,rmax=1){

  rzhat=(-1 + (1 + b*rx^k)^(c))/(-1 + (1 + b*rmax^k)^(c))
  return(rzhat)
}

ThreeParfunInv<-function(b,c,k,rz,rmax=1){

  rx=((-1 + (1 + (-1 + (1 + b * rmax^k)^c) * rz)^(1/c))/b)^(1/k)
  return(rx)
}

ThreeParfunObj<-function(par,rz,rx,rmax=1){
  b=par[1];  c=par[2]; k=par[3]
  rzhat=ThreeParfun(b = b, c = c, k = k, rx = rx, rmax = rmax)

  SSE=sum((rz-rzhat)^2)
  SSE=ifelse(is.nan(SSE),10^6,SSE)
  return(SSE)
}

ThreeParfunObjInv<-function(par,rz,rx,rmax=1){
  b=par[1];  c=par[2]; k=par[3]
  rxhat=ThreeParfunInv(b = b, c = c, k = k, rz = rz, rmax = rmax)

  SSE=sum((rx-rxhat)^2)
  SSE=ifelse(is.nan(SSE),10^6,SSE)
  return(SSE)
}

# fx=fy='qzi'
# pfx=list(Distr=qgamma, p0=0.97, shape=0.01, scale=1)
# pfy=list(Distr=qgamma, p0=0.97, shape=0.01, scale=1)
# rhoz=seq(from=-1, to=1 , length.out = 21)
# rhox=NatafInt(rho = rhoz, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy)
# plot(rhox,rhoz, col='red', pch=19); abline(0,1)
#
# rhotarget=seq(from=-1, to=1 , by=0.01)
# req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#               NatafIntMethod = 'GH', polydeg=8, NoEval = 51)
# lines(rhotarget, req$rzEq, col='blue', lwd=2);
#
# req=NatafInvD(targetrho = rhotarget, fx = fx, fy = fy, paramlistfx = pfx, paramlistfy = pfy,
#              NatafIntMethod = 'GH', polydeg=-3, NoEval = 9)
# plot(rhoz,rhox, col='red', pch=19); abline(0,1)
# points(req$rzEq,rhotarget, col='green')




