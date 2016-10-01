#' GE_translate_inputs.R
#'
#' Mostly for internal use, function called by GE_bias_normal() and GE_scoreeq_sim()
#' to translate the rho_list inputs and return a total covariance matrix for simulation/
#' checking validity of covariance structure.  If invalid covariance structure, will stop
#' and return an error message.
#' 
#' @param rho_list A list of the 6 pairwise covariances between the
#' covariates.  These should be in the order (1) cov_GE (2) cov_GZ (3) cov_EZ
#' (4) cov_GW (5) cov_EW (6) cov_ZW. If Z or M are vectors then terms like cov_GZ should be vectors 
#' (in the appropriate order).
#' If Z or M are vectors, then cov_ZW should be a vector in the order (cov(Z_1,W_1),...,cov(Z_1,W_q),
#' cov(Z_2,W_1),........,cov(Z_p,W_q) where Z is a vector of length p and W is a vector of length q.
#' @param cov_Z Only used if Z is a vector, gives the covariance matrix of Z (remember by assumption
#' Z has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
#' of the vector.
#' @param cov_W Only used if W is a vector, gives the covariance matrix of W (remember by assumption
#' W has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
#' of the vector.
#' @param prob_G Probability that each allele is equal to 1.  Since each SNP has
#' two alleles, the expectation of G is 2*prob_G.
#'
#' @return A list with the elements:
#' \item{sig_mat_total}{The sigma parameter for rmvnorm call to generate our data.}
#' \item{sig_mat_ZZ}{The covariance matrix of Z, i.e. E[ZZ^T]}
#' \item{sig_mat_WW}{The covariance matrix of W, i.e. E[WW^T]}
#'
#' @keywords internal
#' @export
#' @examples 
#' GE_translate_inputs( beta_list=as.list(runif(n=6, min=0, max=1)), 
#'							rho_list=as.list(rep(0.3,6)), prob_G=0.3)

GE_translate_inputs <- function(beta_list, rho_list, prob_G, cov_Z=NULL, cov_W=NULL)
{
	# First, make sure we got good inputs
  	if (length(beta_list) != 6 | length(rho_list) != 6 | class(beta_list) != 'list' | class(rho_list) != 'list')
  	{
  	  stop('Input vectors not the right size!')
  	}

	# How long are these vectors?  Remember that W is the same length as M by assumption.
  	num_Z <- length(beta_list[[5]])
  	num_W <- length(beta_list[[6]])
  	num_rho <- 2*(num_Z+num_W) + num_Z*num_W + 1
  	
  	# Make sure we have compatible lengths for rho_list
    if (length(rho_list[[2]]) != num_Z | length(rho_list[[3]]) != num_Z | length(rho_list[[4]]) != num_W
   			| length(rho_list[[5]]) != num_W | length(rho_list[[6]]) != num_Z*num_W) {
   		stop('Incompatible number of elements in beta/rho_list')
   	}
    	
   	# Fill in our covariances.
   	rho_GE <- rho_list[[1]]; rho_GZ <- rho_list[[2]]; rho_EZ <- rho_list[[3]]
   	rho_GW <- rho_list[[4]]; rho_EW <- rho_list[[5]]; rho_ZW <- rho_list[[6]]
   	
   	
   	
   	################################################################
    # Build our covariance matrix in steps.
    # The 3x3 in the top left is always the same to build, vectors or not.
    w <- qnorm(1-prob_G)					# Threshold for generating G
    r_GE <- rho_GE / (2*dnorm(w))	
    sig_mat_GE <- matrix(data=c(1, 0, r_GE, 0, 1, r_GE, r_GE, r_GE, 1), nrow=3) 
    
    # Build the p*3 matrix that describes Z with G1,G2,E
    sig_mat_Z_column <- matrix(data=NA, nrow=num_Z, ncol=3)
    r_GZ <- rho_GZ / (2*dnorm(w))
    for (i in 1:num_Z) {
    	sig_mat_Z_column[i,] <- c(r_GZ[i], r_GZ[i], rho_EZ[i])
    }
    
    # Build the q*3 matrix that describes W with G1,G2,E
    sig_mat_W_column <- matrix(data=NA, nrow=num_W, ncol=3)
    r_GW <- rho_GW / (2*dnorm(w))
    for (i in 1:num_W) {
    	sig_mat_W_column[i,] <- c(r_GW[i], r_GW[i], rho_EW[i])
    }
    
    # Build the p*q matrix that describes Z with W
    sig_mat_Z_W <- matrix(data=NA, nrow=num_Z, ncol=num_W)
    for (i in 1:num_Z) {
    	start_ind <- (i-1)*num_W+1
    	end_ind <- i*num_W
    	sig_mat_Z_W[i,] <- rho_ZW[start_ind:end_ind]
    }
    
    # If Z or W vectorized, build the ZZ and WW covariance matrices too
    if (num_Z > 1) {
    	sig_mat_ZZ <- matrix(data=0, nrow=num_Z, ncol=num_Z)
    	sig_mat_ZZ[upper.tri(sig_mat_ZZ)] <- cov_Z
    	sig_mat_ZZ <- sig_mat_ZZ + t(sig_mat_ZZ)
    	diag(sig_mat_ZZ) <- 1
    } else {
    	sig_mat_ZZ <- matrix(data=1, nrow=1, ncol=1)
    }
    
    if (num_W > 1) {
    	sig_mat_WW <- matrix(data=0, nrow=num_W, ncol=num_W)
    	sig_mat_WW[upper.tri(sig_mat_WW)] <- cov_W
    	sig_mat_WW <- sig_mat_WW + t(sig_mat_WW)
    	diag(sig_mat_WW) <- 1
    } else {
    	sig_mat_WW <- matrix(data=1, nrow=1, ncol=1)
    }
   
    # Now put it all together
    sig_mat_total <- matrix(data=NA, nrow=(3+num_Z+num_W), ncol=(3+num_Z+num_W))
    sig_mat_total[1:3, 1:3] <- sig_mat_GE
    sig_mat_total[4:(3+num_Z), 1:3] <- sig_mat_Z_column
    sig_mat_total[1:3, 4:(3+num_Z)] <- t(sig_mat_Z_column)
    sig_mat_total[(4+num_Z):(3+num_Z+num_W), 1:3] <- sig_mat_W_column
    sig_mat_total[1:3, (4+num_Z):(3+num_Z+num_W)] <- t(sig_mat_W_column)
    sig_mat_total[4:(3+num_Z), 4:(3+num_Z)] <- sig_mat_ZZ
    sig_mat_total[(4+num_Z):(3+num_Z+num_W), (4+num_Z):(3+num_Z+num_W)] <- sig_mat_WW
    sig_mat_total[4:(3+num_Z), (4+num_Z):(3+num_Z+num_W)] <- sig_mat_Z_W
    sig_mat_total[(4+num_Z):(3+num_Z+num_W), 4:(3+num_Z)] <- t(sig_mat_Z_W)
    if (!isSymmetric(sig_mat_total)) {stop("Problem building covariance matrix!")}
    
    # Now make sure we can actually generate data with this structure
    test_data <- tryCatch(mvtnorm::rmvnorm(n=1, sigma=sig_mat_total), 
    				warning=function(w) w, error=function(e) e)
    if (class(test_data)[1] != 'matrix') {stop('You specified an impossible covariance matrix!')}
    
    return(list(sig_mat_total=sig_mat_total, sig_mat_ZZ=sig_mat_ZZ, sig_mat_WW=sig_mat_WW))
}




