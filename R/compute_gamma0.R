##' Compute the value of gamma_0 from gamma1,...,gamma_m
##'
##' Note: assume we must compute the distibution of W in the controls
##' @title Compute gamma0, given haplotype frequencies
##' @param N0 number of samples with Y=0
##' @param N1 number of samples with Y=1
##' @param W The causal SNPs
##' @param gamma.CV The odds ratios of effect for each CV
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @return The value of gamma0
##' @author Mary Fortune
compute_gamma0<-function(N0,N1,W,gamma.W,freq){
  eta.CV<-rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma.W,`*`))
  e.eta.CV<-exp(eta.CV)
  GenoProbW<-make_GenoProbList(W[1],W,freq)[[1]]
  PW<-GenoProbW[[1]]+GenoProbW[[2]]+GenoProbW[[3]]
  denom<-e.eta.CV*PW
  gamma0<-log(N1)-log(N0)-log(sum(denom))
  return(gamma0)
}


#for when we know PW

##' Compute the value of gamma_0 from gamma1,...,gamma_m
##'
##' Note: assume we know the distibution of W in the controls
##' @title Compute gamma0 given distribution of W
##' @param N0 number of samples with Y=0
##' @param N1 number of samples with Y=1
##' @param W The causal SNPs
##' @param gamma.CV The odds ratios of effect for each CV
##' @param PWgY0 the distribution of causal snps in controls
##' @return The value of gamma0
##' @author Mary Fortune
compute_gamma0_PW<-function(N0,N1,W,gamma.CV,PWgY0){
  e.eta.CV<-exp(rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma.CV,`*`)))
  gamma0<-log(N1)-log(N0)-log(sum(e.eta.CV*PWgY0))
  return(gamma0)
}
