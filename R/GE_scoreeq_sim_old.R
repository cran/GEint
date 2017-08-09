#' GE_scoreeq_sim_old.R
#'
#' Here we perform simulation to verify that we have solved for
#' the correct alpha values in GE_bias_norm_squaredmis_old().
#' Make the same assumptions as in GE_bias_norm_squaredmis_old().
#'
#' @param num_sims The number of simulations to run, we suggest 5000.
#' @param num_sub The number of subjects to generate in every simulation, we suggest 2000.
#' @param beta_list A list of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' If Z or M is a vector, then beta_Z and beta_M should be vectors.
#' @param rho_list A list of the 6 pairwise covariances between the
#' covariates.  These should be in the order (1) cov_GE (2) cov_GZ (3) cov_EZ
#' (4) cov_GW (5) cov_EW (6) cov_ZW.
#' Again if Z or W are vectors then terms like cov_GZ should be vectors (in the order
#' cov(G,Z_1),...,cov(G,Z_p)) where Z is of dimension p, and similarly for W.
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
#' @return A list of the fitted values alpha
#'
#' @export
#' @examples 
#' GE_scoreeq_sim_old( num_sims=10, beta_list=as.list(runif(n=6, min=0, max=1)), 
#'							rho_list=as.list(rep(0.3,6)), prob_G=0.3)

GE_scoreeq_sim_old <- function(num_sims=5000, num_sub=2000, beta_list, prob_G, rho_list, cov_Z=NULL, cov_W=NULL)
{
	# Need survival function
	surv <- function(x) {1-pnorm(x)}
  	
  # Fill in our covariances.
  rho_GE <- rho_list[[1]]; rho_GZ <- rho_list[[2]]; rho_EZ <- rho_list[[3]]
  rho_GW <- rho_list[[4]]; rho_EW <- rho_list[[5]]; rho_ZW <- rho_list[[6]]
  	
  # Quantities necessary for calculating higher order moments
  w <- qnorm(1-prob_G)					
  r_GE <- rho_GE / (2*dnorm(w))	
  r_GZ <- rho_GZ / (2*dnorm(w))
  r_GW <- rho_GW / (2*dnorm(w))
   
   # Get the total covariance matrix (also some basic validity checks)
	translated_inputs <- GE_translate_inputs_old(beta_list=beta_list, rho_list=rho_list, 
									prob_G=prob_G, cov_Z=cov_Z, cov_W=cov_W)
	sig_mat <- translated_inputs$sig_mat_total
	sig_mat_ZZ <- translated_inputs$sig_mat_ZZ
	sig_mat_WW <- translated_inputs$sig_mat_WW
	
	# Different from GE_bias_normal_squaredmis in that we need num_Z/num_W=0 if they don't exist.
	if (is.null(sig_mat_ZZ)) {
	  num_build_Z <- 0
	} else {
	  num_build_Z <- length(beta_list[[5]])
	}
  if (is.null(sig_mat_WW)) {
    num_build_W <- 0
  }	else {
    num_build_W <- length(beta_list[[6]])
  }

	results <- matrix(data=NA, nrow=num_sims, ncol=(4+num_build_Z+num_build_W))
	for(i in 1:num_sims)
	{
		# Sim covariates
		sim_data <- mvtnorm::rmvnorm(n=num_sub, mean=c(rep(0,3+num_build_Z+num_build_W)), 
		                    sigma=sig_mat)		
		snp1 <- as.numeric(sim_data[,1]>w)
		snp2 <- as.numeric(sim_data[,2]>w)
		G <- snp1 + snp2 - 2*prob_G
		E <- sim_data[,3]
		
		# Build the design matrix depending on if Z and M/W exist.
		if (num_build_Z > 0) { 
		  Z <- sim_data[,4:(3+num_build_Z)] 
		  if (num_build_W > 0)  {         # Normal
		    W <- sim_data[,(4+num_build_Z):(3+num_build_Z+num_build_W)]
		    d_right <- cbind(1, G, E^2, G*E^2, Z, W^2)
		    d_wrong <- cbind(1, G, E, G*E, Z, W)
		  } else if (num_build_W == 0) {
		    d_right <- cbind(1, G, E^2, G*E^2, Z)
		    d_wrong <- cbind(1, G, E, G*E, Z)
		  }
		} else if (num_build_Z == 0) {
		  if (num_build_W > 0)  {         
		    W <- sim_data[,(4+num_build_Z):(3+num_build_Z+num_build_W)]
		    d_right <- cbind(1, G, E^2, G*E^2, W^2)
		    d_wrong <- cbind(1, G, E, G*E, W)
		  } else if (num_build_W == 0) {          # No other covariates in model
		    d_right <- cbind(1, G, E^2, G*E^2)
		    d_wrong <- cbind(1, G, E, G*E)
		  }
		}
		
		# Remove 0s from beta vector
		true_beta <- beta_list
		if (num_build_Z == 0 & num_build_W == 0) {
		  true_beta <- true_beta[-c(5,6)]
		} else if (num_build_Z == 0 & num_build_W != 0) {
		  true_beta <- true_beta[-5]
		} else if (num_build_Z != 0 & num_build_W == 0) {
		  true_beta <- true_beta[-6]
		}

		Y <- d_right %*% unlist(true_beta) + rnorm(num_sub)
		
		# Solve for beta_hat
		b_hat <- solve(t(d_wrong) %*% d_wrong) %*% t(d_wrong) %*% Y
		results[i, 1:length(b_hat)] <- b_hat
		
		# Checkpoint
		if (i%%1000 == 0) {cat(i, "done\n")}
	}
	
	# Get the column names correct
	if (num_build_Z == 0 & num_build_W == 0) {
	  colnames(results) <- c('alpha_0', 'alpha_G', 'alpha_E', 'alpha_I')
	} else if (num_build_Z == 0 & num_build_W != 0) {
	  colnames(results) <- c('alpha_0', 'alpha_G', 'alpha_E', 'alpha_I', rep('ALPHA_W',num_build_W))
	} else if (num_build_Z != 0 & num_build_W == 0) {
	  colnames(results) <- c('alpha_0', 'alpha_G', 'alpha_E', 'alpha_I', rep('ALPHA_Z',num_build_Z))
	} else {
	  colnames(results) <- c('alpha_0', 'alpha_G', 'alpha_E', 'alpha_I', rep('ALPHA_Z',num_build_Z), rep('ALPHA_W',num_build_W))
	}
	
	sim_alpha <- apply(results, 2, mean)
	
	return(sim_alpha)
}

