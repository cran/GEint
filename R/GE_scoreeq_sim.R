#' GE_scoreeq_sim.R
#'
#' Here we perform simulation to verify that we have solved for
#' the correct alpha values in GE_bias_norm_squaredmis().
#' Make the same assumptions as in GE_bias_norm_squaredmis().
#'
#' @param num_sims The number of simulations to run, we suggest 5000.
#' @param num_sub The number of subjects to generate in every simulation, we suggest 2000.
#' @param beta_list A list of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' If G or Z or M is a vector, then beta_G/beta_Z/beta_M should be vectors.
#' If Z and/or M/W do not exist in your model, then set beta_Z and/or beta_M = 0.
#' @param rho_list A list of expectations (which happen to be covariances if all covariates
#' are centered at 0) in the order specified by GE_enumerate_inputs().
#' If Z and/or M/W do not exist in your model, then treat them as constants 0. For example,
#' if Z doesn't exist and W includes 2 covariates, then set cov(EZ) = 0 and cov(ZW) = (0,0).
#' If describing expectations relating two vectors, i.e. Z includes two covariates and W
#' includes three covariates, sort by the first term and then the second. Thus in the 
#' example, the first three terms of cov(ZW) are cov(Z_1,W_1),cov(Z_1,W_2), cov(Z_1,W_3), 
#' and the last three terms are cov(Z_3,W_1), cov(Z_3,W_2), cov(Z_3,W_3).
#' @param prob_G Probability that each allele is equal to 1.  Since each SNP has
#' two alleles, the expectation of G is 2*prob_G. Should be a d*1 vector.
#' @param cov_Z Should be a matrix equal to cov(Z) or NULL if no Z. 
#' @param cov_W Should be a matrix equal to cov(W) or NULL if no W. 
#' @param corr_G Should be a matrix giving the *pairwise correlations* between each SNP
#' in the set, or NULL. Must be specified if G is a vector.  For example, the [2,3] element
#' of the matrix would be the pairwise correlation between SNP2 and SNP3.
#'
#' @return A list of the fitted values alpha
#'
#' @export
#' @examples 
#' GE_scoreeq_sim( num_sims=10, num_sub=1000, beta_list=as.list(runif(n=6, min=0, max=1)), 
#' rho_list=as.list(rep(0.3,6)), prob_G=0.3, cov_Z=1, cov_W=1)

GE_scoreeq_sim <- function(num_sims=5000, num_sub=2000, beta_list, rho_list, prob_G, cov_Z=NULL, cov_W=NULL, corr_G=NULL)
{
  # Need survival function.
  surv <- function(x) {1-pnorm(x)}
  
  # For thresholding
  w_vec <- qnorm(1-prob_G)
  
  # Record some initial quantities
  rho_GE <- rho_list[[1]]; rho_GZ <- rho_list[[2]]; rho_EZ <- rho_list[[3]]
  rho_GW <- rho_list[[4]]; rho_EW <- rho_list[[5]]; rho_ZW <- rho_list[[6]]

  beta_0 <- beta_list[[1]]; beta_G <- beta_list[[2]]; beta_E <- beta_list[[3]]
  beta_I <- beta_list[[4]]; beta_Z <- beta_list[[5]]; beta_M <- beta_list[[6]]
  
  # How long are they
  num_G <- length(beta_list[[2]])
  num_build_W <- length(beta_list[[6]])
  if (num_build_W == 1 & beta_list[[6]][1] == 0) {
    num_build_W = 0
  } 
  num_build_Z <- length(beta_list[[5]])
  if (num_build_Z == 1 & beta_list[[5]][1] == 0) {
    num_build_Z = 0
  }
  
  # Some error checking, make sure the covariance matrix is ok
  translated_inputs <- GE_translate_inputs(beta_list=beta_list, rho_list=rho_list, 
                                               prob_G=prob_G, cov_Z=cov_Z, cov_W=cov_W, corr_G=corr_G)
  sig_mat <- translated_inputs$sig_mat_total

  # Loop through simulations
	results <- matrix(data=NA, nrow=num_sims, ncol=(2+2*num_G+num_build_Z+num_build_W))
	for(i in 1:num_sims)
	{
		# Sim covariates
		sim_data <- mvtnorm::rmvnorm(n=num_sub, mean=c(rep(0,1+2*num_G+num_build_Z+num_build_W)), 
		                    sigma=sig_mat)	
		
		G <- matrix(data=NA, nrow=num_sub, ncol=num_G)
		for (j in 1:num_G) {
		  snp1 <- as.numeric(sim_data[, (j*2-1)] > w_vec[j])
		  snp2 <- as.numeric(sim_data[, (j*2)] > w_vec[j])
		  G[, j] <- snp1 + snp2 - 2*prob_G[j]
		}
		E <- sim_data[, (2*num_G + 1)]
		
		# Build the design matrix depending on if Z and M/W exist.
		if (num_build_Z > 0) { 
		  Z <- sim_data[,(2*num_G+2):(2*num_G+1+num_build_Z)] 
		  if (num_build_W > 0)  {         # Normal
		    W <- sim_data[, (2*num_G+2+num_build_Z):(2*num_G+1+num_build_Z+num_build_W)]
		    d_right <- cbind(1, G, E^2, G*E^2, Z, W^2)
		    d_wrong <- cbind(1, G, E, G*E, Z, W)
		  } else if (num_build_W == 0) {
		    d_right <- cbind(1, G, E^2, G*E^2, Z)
		    d_wrong <- cbind(1, G, E, G*E, Z)
		  }
		} else if (num_build_Z == 0) {
		  if (num_build_W > 0)  {         
		    W <- sim_data[,(2*num_G+2):(2*num_G+1+num_build_W)]
		    d_right <- cbind(1, G, E^2, G*E^2, W^2)
		    d_wrong <- cbind(1, G, E, G*E, W)
		  } else if (num_build_W == 0) {          # No other covariates in model
		    d_right <- cbind(1, G, E^2, G*E^2)
		    d_wrong <- cbind(1, G, E, G*E)
		  }
		}
		
		# Remove 0s from beta_list (it's a list)
		true_beta <- beta_list
		if (num_build_Z == 0 & num_build_W == 0) {
		  true_beta <- true_beta[-c(5,6)]
		} else if (num_build_Z == 0 & num_build_W != 0) {
		  true_beta <- true_beta[-5]
		} else if (num_build_Z != 0 & num_build_W == 0) {
		  true_beta <- true_beta[-6]
		}

		# Simulate outcome
		Y <- d_right %*% unlist(true_beta) + rnorm(num_sub)
		
		# Solve for beta_hat
		b_hat <- solve(t(d_wrong) %*% d_wrong) %*% t(d_wrong) %*% Y
		results[i, 1:length(b_hat)] <- b_hat
		
		# Checkpoint
		if (i%%1000 == 0) {cat(i, "done\n")}
	}
	
	# Get the column names correct
	colnames(results) <- c('alpha_0', rep('alpha_G', num_G), 'alpha_E', 
	                       rep('alpha_I', num_G), rep('alpha_Z', num_build_Z), 
	                       rep('alpha_W', num_build_W)) 
	
	# Get the simulated alphas
	sim_alpha <- apply(results, 2, mean)
	
	return(sim_alpha)
}

