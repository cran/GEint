#' GE_nleqslv.R
#'#'
#' Uses package nleqslv to get a numerical solution to the score equations, which
#' we can use to check our direct solution from GE_bias().
#' 
#' @param beta_list A list of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' If G or Z or M is a vector, then beta_G/beta_Z/beta_M should be vectors.
#' If Z and/or M/W do not exist in your model, then set beta_Z and/or beta_M = 0.
#' @param cov_list A list of expectations (which happen to be covariances if all covariates
#' are centered at 0) in the order specified by GE_enumerate_inputs().
#' If Z and/or M/W do not exist in your model, then treat them as constants 0. For example,
#' if Z doesn't exist and W includes 2 covariates, then set cov(EZ) = 0 and cov(ZW) = (0,0).
#' If describing expectations relating two vectors, i.e. Z includes two covariates and W
#' includes three covariates, sort by the first term and then the second. Thus in the 
#' example, the first three terms of cov(ZW) are cov(Z_1,W_1),cov(Z_1,W_2), cov(Z_1,W_3), 
#' and the last three terms are cov(Z_3,W_1), cov(Z_3,W_2), cov(Z_3,W_3).
#' @param cov_mat_list  A list of matrices of expectations as specified by GE_enumerate_inputs().
#' @param mu_list A list of means as specified by GE_enumerate_inputs().
#' @param HOM_list A list of higher order moments as specified by GE_enumerate_inputs().
#' 
#' @return A list of the fitted coefficients alpha
#'
#' @export
#' @examples 
#' solutions <- GE_bias_normal_squaredmis( beta_list=as.list(runif(n=6, min=0, max=1)),
#' rho_list=as.list(rep(0.3,6)), prob_G=0.3, cov_Z=1, cov_W=1)
#' GE_nleqslv(beta_list=solutions$beta_list, solutions$cov_list, solutions$cov_mat_list,
#' solutions$mu_list, solutions$HOM_list)

GE_nleqslv <- function(beta_list, cov_list, cov_mat_list, mu_list, HOM_list)
{

  # Record some initial quantities
  beta_0 <- beta_list[[1]]; beta_G <- beta_list[[2]]; beta_E <- beta_list[[3]]
  beta_I <- beta_list[[4]]; beta_Z <- beta_list[[5]]; beta_M <- beta_list[[6]]
  
  # Some means
  mu_f <- mu_list[[1]]; mu_h <- mu_list[[2]]; mu_Z <- as.matrix(mu_list[[3]], ncol=1)
  mu_M <- as.matrix(mu_list[[4]], ncol=1); mu_W <- as.matrix(mu_list[[5]], ncol=1)
  
  # Some covariances
  mu_GE <- as.matrix(cov_list[[1]], ncol=1)
  mu_Gf <- as.matrix(cov_list[[2]], ncol=1)
  mu_Gh <- as.matrix(cov_list[[3]], ncol=1)
  mu_EE <- as.matrix(cov_list[[4]], ncol=1)
  mu_Ef <- as.matrix(cov_list[[5]], ncol=1)
  mu_EZ <- as.matrix(cov_list[[6]], ncol=1)
  mu_EM <- as.matrix(cov_list[[7]], ncol=1)
  mu_EW <- as.matrix(cov_list[[8]], ncol=1)
  mu_fZ <- as.matrix(cov_list[[9]], ncol=1)
  mu_fW <- as.matrix(cov_list[[10]], ncol=1)
  
  # Matrix covariances
  mu_GG <- cov_mat_list[[1]]
  mu_GZ <- cov_mat_list[[2]]; mu_ZG <- t(mu_GZ)
  mu_GM <- cov_mat_list[[3]]
  mu_GW <- cov_mat_list[[4]]; mu_WG <- t(mu_GW)
  mu_ZZ <- cov_mat_list[[5]]
  mu_ZW <- cov_mat_list[[6]]; mu_WZ <- t(mu_ZW)
  mu_ZM <- cov_mat_list[[7]]
  mu_WW <- cov_mat_list[[8]]
  mu_WM <- cov_mat_list[[9]]
  
  # Some higher order moments
  mu_GEG <- HOM_list[[1]]
  mu_GhG <- HOM_list[[2]]
  mu_GEE <- HOM_list[[3]]
  mu_GEf <- HOM_list[[4]]
  mu_GEh <- HOM_list[[5]]
  mu_GEZ <- HOM_list[[6]]; mu_ZEG <- t(mu_GEZ)
  mu_GEM <- HOM_list[[7]]
  mu_GEW <- HOM_list[[8]]; mu_WEG <- t(mu_GEW)
  mu_GhW <- HOM_list[[9]]; mu_WhG <- t(mu_GhW)
  mu_GhZ <- HOM_list[[10]]; mu_ZhG <- t(mu_GhZ)
  mu_GEEG <- HOM_list[[11]]
  mu_GEfG <- HOM_list[[12]]
  mu_GEhG <- HOM_list[[13]]
	
	#################################################################			
	# Define the set of score equations we will be solving
	# Different score equations for no Z and no M/W
	num_Z <- length(beta_Z); if(num_Z == 1 & beta_Z[1] == 0) {num_Z <- 0}
	num_W <- length(beta_M); if(num_W == 1 & beta_M[1] == 0) {num_W <- 0}
	num_G <- length(beta_G)
	
	score_eqs_all <- function(x)
	{
	 	alpha_0 <- x[1]
	 	alpha_G <- x[2:(2+num_G-1)]
	 	alpha_E <- x[(2+num_G)]
	 	alpha_I <- x[(3+num_G):(3+2*num_G-1)]
	 	alpha_Z <- x[(3+2*num_G):(3+2*num_G+num_Z-1)]
	 	alpha_W <- x[(3+2*num_G+num_Z):(3+2*num_G+num_Z+num_W-1)]
	 	y <- numeric(2+num_Z+num_W+2*num_G)
	 	
	 	# For alpha_0
	 	y[1] = alpha_0 + t(mu_GE) %*% alpha_I + t(mu_Z) %*% alpha_Z + t(mu_W) %*% alpha_W - beta_0 - 
	  	  beta_E*mu_f - t(mu_Gh) %*% beta_I - t(mu_Z) %*% beta_Z - t(mu_M) %*% beta_M
	
	 	# For alpha_G
	  y[2:(2+num_G-1)] = mu_GG %*% alpha_G + mu_GE %*% alpha_E + mu_GEG %*% alpha_I +
	     mu_GZ %*% alpha_Z + mu_GW %*% alpha_W - mu_GG %*% beta_G - 
	     beta_E*mu_Gf - mu_GhG %*% beta_I - mu_GZ %*% beta_Z - mu_GM %*% beta_M
	  
	  # For alpha_E
		y[(2+num_G)] = t(mu_GE) %*% alpha_G + alpha_E*mu_EE +  t(mu_GEE) %*% alpha_I +
		    t(mu_EZ) %*% alpha_Z + t(mu_EW) %*% alpha_W - t(mu_GE) %*% beta_G -
		    beta_E*mu_Ef - t(mu_GEh) %*% beta_I - t(mu_EZ) %*% beta_Z - t(mu_EM) %*% beta_M
	
		# For alpha_I
		y[(3+num_G):(3+2*num_G-1)] = alpha_0*mu_GE + mu_GEG %*% alpha_G + alpha_E*mu_GEE + 
		    mu_GEEG %*% alpha_I + mu_GEZ %*% alpha_Z + mu_GEW %*% alpha_W - beta_0*mu_GE -
		    mu_GEG %*% beta_G - beta_E*mu_GEf - mu_GEhG %*% beta_I - mu_GEZ %*% beta_Z - 
		    mu_GEM %*% beta_M
	
		# For alpha_Z
	  y[(3+2*num_G):(3+2*num_G+num_Z-1)] = alpha_0*mu_Z + mu_ZG %*% alpha_G + 
	      alpha_E*mu_EZ + mu_ZEG %*% alpha_I + mu_ZZ %*% alpha_Z + 
		    mu_ZW %*% alpha_W - beta_0*mu_Z - mu_ZG %*% beta_G - beta_E*mu_fZ - 
	      mu_ZhG %*% beta_I - mu_ZZ %*% beta_Z - mu_ZM %*% beta_M
	
	 # For alpha_W
	  y[(3+2*num_G+num_Z):(3+2*num_G+num_Z+num_W-1)] = alpha_0*mu_W + mu_WG %*% alpha_G + 
	    alpha_E*mu_EW + mu_WEG %*% alpha_I + mu_WZ %*% alpha_Z + mu_WW %*% alpha_W - 
	    beta_0*mu_W - mu_WG %*% beta_G - beta_E*mu_fW - mu_WhG %*% beta_I - 
	    mu_WZ %*% beta_Z - mu_WM %*% beta_M
	
	 # Output y
		y
	}

	# Score equations no W
	score_eqs_no_W <- function(x)
	{
	  alpha_0 <- x[1]
	  alpha_G <- x[2:(2+num_G-1)]
	  alpha_E <- x[(2+num_G)]
	  alpha_I <- x[(3+num_G):(3+2*num_G-1)]
	  alpha_Z <- x[(3+2*num_G):(3+2*num_G+num_Z-1)]
	  alpha_W <- 0; beta_W <- 0; beta_M <- 0
	  y <- numeric(2+2*num_G+num_Z)
	  
	  # For alpha_0
	  y[1] = alpha_0 + t(mu_GE) %*% alpha_I + t(mu_Z) %*% alpha_Z + t(mu_W) %*% alpha_W - beta_0 - 
	    beta_E*mu_f - t(mu_Gh) %*% beta_I - t(mu_Z) %*% beta_Z - t(mu_M) %*% beta_M
	  
	  # For alpha_G
	  y[2:(2+num_G-1)] = mu_GG %*% alpha_G + mu_GE %*% alpha_E + mu_GEG %*% alpha_I +
	    mu_GZ %*% alpha_Z + mu_GW %*% alpha_W - mu_GG %*% beta_G - 
	    beta_E*mu_Gf - mu_GhG %*% beta_I - mu_GZ %*% beta_Z - mu_GM %*% beta_M
	  
	  # For alpha_E
	  y[(2+num_G)] = t(mu_GE) %*% alpha_G + alpha_E*mu_EE +  t(mu_GEE) %*% alpha_I +
	    t(mu_EZ) %*% alpha_Z + t(mu_EW) %*% alpha_W - t(mu_GE) %*% beta_G -
	    beta_E*mu_Ef - t(mu_GEh) %*% beta_I - t(mu_EZ) %*% beta_Z - t(mu_EM) %*% beta_M
	  
	  # For alpha_I
	  y[(3+num_G):(3+2*num_G-1)] = alpha_0*mu_GE + mu_GEG %*% alpha_G + alpha_E*mu_GEE + 
	    mu_GEEG %*% alpha_I + mu_GEZ %*% alpha_Z + mu_GEW %*% alpha_W - beta_0*mu_GE -
	    mu_GEG %*% beta_G - beta_E*mu_GEf - mu_GEhG %*% beta_I - mu_GEZ %*% beta_Z - 
	    mu_GEM %*% beta_M
	  
	  # For alpha_Z
	  y[(3+2*num_G):(3+2*num_G+num_Z-1)] = alpha_0*mu_Z + mu_ZG %*% alpha_G + 
	    alpha_E*mu_EZ + mu_ZEG %*% alpha_I + mu_ZZ %*% alpha_Z + 
	    mu_ZW %*% alpha_W - beta_0*mu_Z - mu_ZG %*% beta_G - beta_E*mu_fZ - 
	    mu_ZhG %*% beta_I - mu_ZZ %*% beta_Z - mu_ZM %*% beta_M
	    
	  # Output y
	  y
	}

	# Score equations no Z
	score_eqs_no_Z <- function(x)
	{
	  alpha_0 <- x[1]
	  alpha_G <- x[2:(2+num_G-1)]
	  alpha_E <- x[2+num_G]
	  alpha_I <- x[(3+num_G):(3+2*num_G-1)]
	  alpha_Z <- 0; beta_Z <- 0
	  alpha_W <- x[(3+2*num_G):(3+2*num_G+num_W-1)]
	  y <- numeric(2+2*num_G+num_W)
	  
	  # For alpha_0
	  y[1] = alpha_0 + t(mu_GE) %*% alpha_I + t(mu_Z) %*% alpha_Z + t(mu_W) %*% alpha_W - beta_0 - 
	    beta_E*mu_f - t(mu_Gh) %*% beta_I - t(mu_Z) %*% beta_Z - t(mu_M) %*% beta_M
	  
	  # For alpha_G
	  y[2:(2+num_G-1)] = mu_GG %*% alpha_G + mu_GE %*% alpha_E + mu_GEG %*% alpha_I +
	    mu_GZ %*% alpha_Z + mu_GW %*% alpha_W - mu_GG %*% beta_G - 
	    beta_E*mu_Gf - mu_GhG %*% beta_I - mu_GZ %*% beta_Z - mu_GM %*% beta_M
	 
	  # For alpha_E
	  y[(2+num_G)] = t(mu_GE) %*% alpha_G + alpha_E*mu_EE +  t(mu_GEE) %*% alpha_I +
	    t(mu_EZ) %*% alpha_Z + t(mu_EW) %*% alpha_W - t(mu_GE) %*% beta_G -
	    beta_E*mu_Ef - t(mu_GEh) %*% beta_I - t(mu_EZ) %*% beta_Z - t(mu_EM) %*% beta_M
	  
	  # For alpha_I
	  y[(3+num_G):(3+2*num_G-1)] = alpha_0*mu_GE + mu_GEG %*% alpha_G + alpha_E*mu_GEE + 
	    mu_GEEG %*% alpha_I + mu_GEZ %*% alpha_Z + mu_GEW %*% alpha_W - beta_0*mu_GE -
	    mu_GEG %*% beta_G - beta_E*mu_GEf - mu_GEhG %*% beta_I - mu_GEZ %*% beta_Z - 
	    mu_GEM %*% beta_M
	 
	  # For alpha_W
	  y[(3+2*num_G):(3+2*num_G+num_W-1)] = alpha_0*mu_W + mu_WG %*% alpha_G + 
	    alpha_E*mu_EW + mu_WEG %*% alpha_I + mu_WZ %*% alpha_Z + mu_WW %*% alpha_W - 
	    beta_0*mu_W - mu_WG %*% beta_G - beta_E*mu_fW - mu_WhG %*% beta_I - 
	    mu_WZ %*% beta_Z - mu_WM %*% beta_M
	  
	  # Output y
	  y
	}

	# Score equations no Z and no W
	score_eqs_no_Z_no_W <- function(x)
	{
    alpha_0 <- x[1]
    alpha_G <- x[2:(2+num_G-1)]
    alpha_E <- x[2+num_G]
	  alpha_I <- x[(3+num_G):(3+2*num_G-1)]
	  alpha_Z <- 0; beta_Z <- 0
	  alpha_W <- 0; alpha_M <- 0; beta_M <- 0
    y <- numeric(2+2*num_G)
    
    # For alpha_0
    y[1] = alpha_0 + t(mu_GE) %*% alpha_I + t(mu_Z) %*% alpha_Z + t(mu_W) %*% alpha_W - beta_0 - 
      beta_E*mu_f - t(mu_Gh) %*% beta_I - t(mu_Z) %*% beta_Z - t(mu_M) %*% beta_M
    
    # For alpha_G
    y[2:(2+num_G-1)] = mu_GG %*% alpha_G + mu_GE %*% alpha_E + mu_GEG %*% alpha_I +
      mu_GZ %*% alpha_Z + mu_GW %*% alpha_W - mu_GG %*% beta_G - 
      beta_E*mu_Gf - mu_GhG %*% beta_I - mu_GZ %*% beta_Z - mu_GM %*% beta_M
    
    # For alpha_E
    y[(2+num_G)] = t(mu_GE) %*% alpha_G + alpha_E*mu_EE +  t(mu_GEE) %*% alpha_I +
      t(mu_EZ) %*% alpha_Z + t(mu_EW) %*% alpha_W - t(mu_GE) %*% beta_G -
      beta_E*mu_Ef - t(mu_GEh) %*% beta_I - t(mu_EZ) %*% beta_Z - t(mu_EM) %*% beta_M
    
    # For alpha_I
    y[(3+num_G):(3+2*num_G-1)] = alpha_0*mu_GE + mu_GEG %*% alpha_G + alpha_E*mu_GEE + 
      mu_GEEG %*% alpha_I + mu_GEZ %*% alpha_Z + mu_GEW %*% alpha_W - beta_0*mu_GE -
      mu_GEG %*% beta_G - beta_E*mu_GEf - mu_GEhG %*% beta_I - mu_GEZ %*% beta_Z - 
      mu_GEM %*% beta_M
      
    # Output y
	  y
	}
	
	# Decide which score equation to use
	if ( num_Z > 0 & num_W > 0 ) {
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,rep(0, num_G),0,rep(0, num_G), rep(0, (num_Z+num_W))), fn=score_eqs_all)
	} else if ( num_Z > 0 & num_W == 0 ) {
	  mu_WW <- 0
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,rep(0, num_G),0,rep(0, num_G), rep(0, num_Z)), fn=score_eqs_no_W)
	} else if ( num_Z == 0 & num_W > 0 ) {
	  mu_ZZ <- 0
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,rep(0, num_G),0,rep(0, num_G), rep(0, num_W)), fn=score_eqs_no_Z)
	} else if ( num_Z == 0 & num_W == 0) {
	  mu_ZZ <- 0; mu_WW <- 0;
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,rep(0, num_G),0,rep(0, num_G)), fn=score_eqs_no_Z_no_W)
	}
	
	return(solved_scoreeqs)
}


