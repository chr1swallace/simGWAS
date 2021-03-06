% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vbeta.R
\name{simulated_vbeta}
\alias{simulated_vbeta}
\title{Compute a simulated var(beta)}
\usage{
simulated_vbeta(
  N0,
  N1,
  snps,
  W,
  gamma.W,
  freq,
  GenoProbList = make_GenoProbList(snps = snps, W = W, freq = freq),
  nrep = 1
)
}
\arguments{
\item{N0}{The number of Y=0}

\item{N1}{The number of Y=1}

\item{snps}{The snps at which we wish to compute the expected Z Score}

\item{W}{The true causal SNPs (these need not be in "snps")}

\item{gamma.W}{The log odds ratios of effect of the true causal SNPs (not including gamma0, the intercept term)}

\item{freq}{Haplotype frequencies as a data.frame, with column Probability indicating relative frequency in controls.}

\item{GenoProbList}{An list of objects giving the probability of seeing each {X,W} genotype vector.  This can be calculated within the function if no value supplied, or you can pass a pre-calculated version}

\item{nrep}{Number of replicates (simulated vectors of Z scores)
under this scenario.  Default=1}
}
\value{
A simulated variance of beta for each SNP X, assuming the causal SNPs are W
}
\description{
Simulate var(beta)
}
\details{
Assumes we have a list, GenoProbList, giving the GenoProb values for each X.
}
\examples{
freq=fake_freq(nhaps=100,nsnps=5) # fake haplotype frequency data
VB=simulated_vbeta(N0=1000,N1=2000,snps=paste0("s",1:5),
                    W="s1",gamma.W=log(1.5),freq=freq)
VB # causal variant is SNP 1, with OR 1.5
}
\author{
Mary Fortune and Chris Wallace
}
