% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GE_translate_inputs_old.R
\name{GE_translate_inputs_old}
\alias{GE_translate_inputs_old}
\title{GE_translate_inputs_old.R}
\usage{
GE_translate_inputs_old(beta_list, rho_list, prob_G, cov_Z = NULL,
  cov_W = NULL)
}
\arguments{
\item{rho_list}{A list of the 6 pairwise covariances between the
covariates.  These should be in the order (1) cov_GE (2) cov_GZ (3) cov_EZ
(4) cov_GW (5) cov_EW (6) cov_ZW. If Z or M are vectors then terms like cov_GZ should be vectors 
(in the appropriate order).
If Z or M are vectors, then cov_ZW should be a vector in the order (cov(Z_1,W_1),...,cov(Z_1,W_q),
cov(Z_2,W_1),........,cov(Z_p,W_q) where Z is a vector of length p and W is a vector of length q.}

\item{prob_G}{Probability that each allele is equal to 1.  Since each SNP has
two alleles, the expectation of G is 2*prob_G.}

\item{cov_Z}{Only used if Z is a vector, gives the covariance matrix of Z (remember by assumption
Z has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
of the vector.}

\item{cov_W}{Only used if W is a vector, gives the covariance matrix of W (remember by assumption
W has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
of the vector.}
}
\value{
A list with the elements:
\item{sig_mat_total}{The sigma parameter for rmvnorm call to generate our data.}
\item{sig_mat_ZZ}{The covariance matrix of Z, i.e. E[ZZ^T]}
\item{sig_mat_WW}{The covariance matrix of W, i.e. E[WW^T]}
}
\description{
Mostly for internal use, function called by GE_bias_normal_old() and GE_scoreeq_sim_old()
to translate the rho_list inputs and return a total covariance matrix for simulation/
checking validity of covariance structure.  If invalid covariance structure, will stop
and return an error message.
}
\examples{
GE_translate_inputs_old( beta_list=as.list(runif(n=6, min=0, max=1)), 
						rho_list=as.list(rep(0.3,6)), prob_G=0.3)
}
\keyword{internal}
