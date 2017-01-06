#' GE_nleqslv
#'#'
#' Uses package nleqslv to get a numerical solution to the score equations, which
#' we can use to check our direct solution from GE_bias().
#' 
#' @param beta_list A list of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' If Z or M is a vector, then beta_Z and beta_M should be vectors.
#' If Z and/or M/W do not exist in your model, then set beta_Z and/or beta_M = 0.
#' @param cov_list A list of expectations (which happen to be covariances if all covariates
#' are centered at 0) in the order specified by GE_enumerate_inputs().
#' If Z and/or M/W do not exist in your model, then treat them as constants 0. For example,
#' set cov(EZ) = 0 and cov(ZW) = 0.
#' @param cov_mat_list  A list of matrices of expectations as specified by GE_enumerate_inputs().
#' @param mu_list A list of means as specified by GE_enumerate_inputs().
#' @param HOM_list A list of higher order moments as specified by GE_enumerate_inputs().
#' 
#' @return A list of the fitted coefficients alpha
#'
#' @export
#' @examples 
#' solutions <- GE_bias_normal_squaredmis( beta_list=as.list(runif(n=6, min=0, max=1)), 
#'							rho_list=as.list(rep(0.3,6)), prob_G=0.3)
#' GE_nleqslv(beta_list=solutions$beta_list, solutions$cov_list, solutions$cov_mat_list, 
#'						solutions$mu_list, solutions$HOM_list)

GE_nleqslv <- function(beta_list, cov_list, cov_mat_list, mu_list, HOM_list)
{

	# Here we extract the relevant parameters from the inputs
	beta_0 <- beta_list[[1]]
	beta_G <- beta_list[[2]]
	beta_E <- beta_list[[3]]
	beta_I <- beta_list[[4]]
	BETA_Z <- beta_list[[5]]
	BETA_M <- beta_list[[6]]

	mu_f <- mu_list[[1]]
	mu_h <- mu_list[[2]]
	MU_Z <- mu_list[[3]]
	MU_M <- mu_list[[4]]
	MU_W <- mu_list[[5]]
	
	num_Z <- length(MU_Z)
	num_W <- length(MU_W)
	
	mu_GG <- cov_list[[1]]
	mu_GE <- cov_list[[2]]
	mu_Gf <- cov_list[[3]]
	mu_Gh <- cov_list[[4]]
	MU_GZ <- cov_list[[5]]
	MU_GM <- cov_list[[6]]
	MU_GW <- cov_list[[7]]
	mu_EE <- cov_list[[8]]
	mu_Ef <- cov_list[[9]]
	MU_EZ <- cov_list[[10]]
	MU_EM <- cov_list[[11]]
	MU_EW <- cov_list[[12]]
	MU_fZ <- cov_list[[13]]
	MU_fW <- cov_list[[14]]
	
	MU_ZZ <- cov_mat_list[[1]]
 	MU_WW <- cov_mat_list[[2]]
	MU_ZW <- cov_mat_list[[3]]
  MU_WZ <- cov_mat_list[[4]]
  MU_ZM <- cov_mat_list[[5]]	
  MU_WM <- cov_mat_list[[6]]		
	
	mu_GGE <- HOM_list[[1]]
	mu_GGh <- HOM_list[[2]]
	mu_GEE <- HOM_list[[3]]
	mu_GEf <- HOM_list[[4]]
	mu_GEh <- HOM_list[[5]]
	MU_GEZ <- HOM_list[[6]]
	MU_GEM <- HOM_list[[7]]
	MU_GEW <- HOM_list[[8]]
	MU_GhW <- HOM_list[[9]]
	MU_GhZ <- HOM_list[[10]]
	mu_GGEE <- HOM_list[[11]]
	mu_GGEf <- HOM_list[[12]]
	mu_GGEh <- HOM_list[[13]]
	
	#################################################################			
	# Define the set of score equations we will be solving
	# Different score equations for no Z and no M/W
	score_eqs_all <- function(x)
	{
	 	alpha_0 <- x[1]
	 	alpha_G <- x[2]
	 	alpha_E <- x[3]
	 	alpha_I <- x[4]
	 	ALPHA_Z <- x[5:(4+num_Z)]
	 	ALPHA_W <- x[(5+num_Z):(4+num_Z+num_W)]
	 	y <- numeric(4+num_Z+num_W)
	 	y[1] = alpha_0 + alpha_I*mu_GE + t(MU_Z) %*% ALPHA_Z + t(MU_W) %*% ALPHA_W - beta_0 - 
	  	beta_E*mu_f - beta_I*mu_Gh - t(MU_Z) %*% BETA_Z - t(MU_M) %*% BETA_M
	
	  y[2] = alpha_G*mu_GG + alpha_E*mu_GE + alpha_I*mu_GGE + t(MU_GZ) %*% ALPHA_Z + t(MU_GW) %*% ALPHA_W - 
	  	beta_G*mu_GG - beta_E*mu_Gf - beta_I*mu_GGh - t(MU_GZ) %*% BETA_Z - t(MU_GM) %*% BETA_M
	
		y[3] = alpha_G*mu_GE + alpha_E*mu_EE + alpha_I*mu_GEE + t(MU_EZ) %*% ALPHA_Z + t(MU_EW) %*% ALPHA_W - 
			beta_G*mu_GE - beta_E*mu_Ef - beta_I*mu_GEh - t(MU_EZ) %*% BETA_Z - t(MU_EM) %*% BETA_M
	
		y[4] = alpha_0*mu_GE + alpha_G*mu_GGE + alpha_E*mu_GEE + alpha_I*mu_GGEE + t(MU_GEZ) %*% ALPHA_Z + 
			t(MU_GEW) %*% ALPHA_W - beta_0*mu_GE - beta_G*mu_GGE - beta_E*mu_GEf - beta_I*mu_GGEh -
	  	t(MU_GEZ) %*% BETA_Z - t(MU_GEM) %*% BETA_M
	
	  y[5:(4+num_Z)] = alpha_0*MU_Z + alpha_G*MU_GZ + alpha_E*MU_EZ + alpha_I*MU_GEZ + MU_ZZ %*% ALPHA_Z + 
		  MU_ZW %*% ALPHA_W - beta_0*MU_Z - beta_G*MU_GZ - beta_E*MU_fZ - beta_I*MU_GhZ - 
		  MU_ZZ %*% BETA_Z - MU_ZM %*% BETA_M
	
	  y[(5+num_Z):(4+num_Z+num_W)] = alpha_0*MU_W + alpha_G*MU_GW + alpha_E*MU_EW + alpha_I*MU_GEW + MU_WZ %*% ALPHA_Z + 
		  MU_WW %*% ALPHA_W - beta_0*MU_W - beta_G*MU_GW - beta_E*MU_fW - beta_I*MU_GhW -
		  MU_WZ %*% BETA_Z - MU_WM %*% BETA_M
	
		y
	}

	# Score equations no W
	score_eqs_no_W <- function(x)
	{
	  alpha_0 <- x[1]
	  alpha_G <- x[2]
	  alpha_E <- x[3]
	  alpha_I <- x[4]
	  ALPHA_Z <- x[5:(4+num_Z)]
	  ALPHA_W <- 0
	  y <- numeric(4+num_Z)
	  y[1] = alpha_0 + alpha_I*mu_GE + t(MU_Z) %*% ALPHA_Z + t(MU_W) %*% ALPHA_W - beta_0 - 
	      beta_E*mu_f - beta_I*mu_Gh - t(MU_Z) %*% BETA_Z - t(MU_M) %*% BETA_M
	    
	  y[2] = alpha_G*mu_GG + alpha_E*mu_GE + alpha_I*mu_GGE + t(MU_GZ) %*% ALPHA_Z + t(MU_GW) %*% ALPHA_W - 
	      beta_G*mu_GG - beta_E*mu_Gf - beta_I*mu_GGh - t(MU_GZ) %*% BETA_Z - t(MU_GM) %*% BETA_M
	    
	  y[3] = alpha_G*mu_GE + alpha_E*mu_EE + alpha_I*mu_GEE + t(MU_EZ) %*% ALPHA_Z + t(MU_EW) %*% ALPHA_W - 
	      beta_G*mu_GE - beta_E*mu_Ef - beta_I*mu_GEh - t(MU_EZ) %*% BETA_Z - t(MU_EM) %*% BETA_M
	    
	  y[4] = alpha_0*mu_GE + alpha_G*mu_GGE + alpha_E*mu_GEE + alpha_I*mu_GGEE + t(MU_GEZ) %*% ALPHA_Z + 
	      t(MU_GEW) %*% ALPHA_W - beta_0*mu_GE - beta_G*mu_GGE - beta_E*mu_GEf - beta_I*mu_GGEh -
	      t(MU_GEZ) %*% BETA_Z - t(MU_GEM) %*% BETA_M
	    
	  y[5:(4+num_Z)] = alpha_0*MU_Z + alpha_G*MU_GZ + alpha_E*MU_EZ + alpha_I*MU_GEZ + MU_ZZ %*% ALPHA_Z + 
	      MU_ZW %*% ALPHA_W - beta_0*MU_Z - beta_G*MU_GZ - beta_E*MU_fZ - beta_I*MU_GhZ - 
	      MU_ZZ %*% BETA_Z - MU_ZM %*% BETA_M
	    
	  y
	}

	# Score equations no Z
	score_eqs_no_Z <- function(x)
	{
	  alpha_0 <- x[1]
	  alpha_G <- x[2]
	  alpha_E <- x[3]
	  alpha_I <- x[4]
	  ALPHA_Z <- 0
	  ALPHA_W <- x[(5):(4+num_W)]
	  y <- numeric(4+num_W)
	  y[1] = alpha_0 + alpha_I*mu_GE + t(MU_Z) %*% ALPHA_Z + t(MU_W) %*% ALPHA_W - beta_0 - 
	      beta_E*mu_f - beta_I*mu_Gh - t(MU_Z) %*% BETA_Z - t(MU_M) %*% BETA_M
	    
	  y[2] = alpha_G*mu_GG + alpha_E*mu_GE + alpha_I*mu_GGE + t(MU_GZ) %*% ALPHA_Z + t(MU_GW) %*% ALPHA_W - 
	      beta_G*mu_GG - beta_E*mu_Gf - beta_I*mu_GGh - t(MU_GZ) %*% BETA_Z - t(MU_GM) %*% BETA_M
	    
	  y[3] = alpha_G*mu_GE + alpha_E*mu_EE + alpha_I*mu_GEE + t(MU_EZ) %*% ALPHA_Z + t(MU_EW) %*% ALPHA_W - 
	      beta_G*mu_GE - beta_E*mu_Ef - beta_I*mu_GEh - t(MU_EZ) %*% BETA_Z - t(MU_EM) %*% BETA_M
	    
	  y[4] = alpha_0*mu_GE + alpha_G*mu_GGE + alpha_E*mu_GEE + alpha_I*mu_GGEE + t(MU_GEZ) %*% ALPHA_Z + 
	      t(MU_GEW) %*% ALPHA_W - beta_0*mu_GE - beta_G*mu_GGE - beta_E*mu_GEf - beta_I*mu_GGEh -
	      t(MU_GEZ) %*% BETA_Z - t(MU_GEM) %*% BETA_M
	    
	  y[(5):(4+num_W)] = alpha_0*MU_W + alpha_G*MU_GW + alpha_E*MU_EW + alpha_I*MU_GEW + MU_WZ %*% ALPHA_Z + 
	      MU_WW %*% ALPHA_W - beta_0*MU_W - beta_G*MU_GW - beta_E*MU_fW - beta_I*MU_GhW -
	      MU_WZ %*% BETA_Z - MU_WM %*% BETA_M
	    
    y
	}

	# Score equations no Z and no W
	score_eqs_no_Z_no_W <- function(x)
	{
    alpha_0 <- x[1]
    alpha_G <- x[2]
    alpha_E <- x[3]
	  alpha_I <- x[4]
	  ALPHA_Z <- 0
	  ALPHA_W <- 0
    y <- numeric(4)
	  y[1] = alpha_0 + alpha_I*mu_GE + t(MU_Z) %*% ALPHA_Z + t(MU_W) %*% ALPHA_W - beta_0 - 
	      beta_E*mu_f - beta_I*mu_Gh - t(MU_Z) %*% BETA_Z - t(MU_M) %*% BETA_M
	    
	  y[2] = alpha_G*mu_GG + alpha_E*mu_GE + alpha_I*mu_GGE + t(MU_GZ) %*% ALPHA_Z + t(MU_GW) %*% ALPHA_W - 
	      beta_G*mu_GG - beta_E*mu_Gf - beta_I*mu_GGh - t(MU_GZ) %*% BETA_Z - t(MU_GM) %*% BETA_M
	    
	  y[3] = alpha_G*mu_GE + alpha_E*mu_EE + alpha_I*mu_GEE + t(MU_EZ) %*% ALPHA_Z + t(MU_EW) %*% ALPHA_W - 
	      beta_G*mu_GE - beta_E*mu_Ef - beta_I*mu_GEh - t(MU_EZ) %*% BETA_Z - t(MU_EM) %*% BETA_M
	    
    y[4] = alpha_0*mu_GE + alpha_G*mu_GGE + alpha_E*mu_GEE + alpha_I*mu_GGEE + t(MU_GEZ) %*% ALPHA_Z + 
	      t(MU_GEW) %*% ALPHA_W - beta_0*mu_GE - beta_G*mu_GGE - beta_E*mu_GEf - beta_I*mu_GGEh -
	      t(MU_GEZ) %*% BETA_Z - t(MU_GEM) %*% BETA_M

	    y
	}
	
	# Decide which score equation to use
	if ( !is.null(MU_ZZ) & !is.null(MU_WW) ) {
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,0,0,0, rep(0, (num_Z+num_W))), fn=score_eqs_all)
	} else if ( !is.null(MU_ZZ) & is.null(MU_WW) ) {
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,0,0,0, rep(0, num_Z)), fn=score_eqs_no_W)
	} else if ( is.null(MU_ZZ) & !is.null(MU_WW) ) {
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,0,0,0, rep(0, num_W)), fn=score_eqs_no_Z)
	} else if ( is.null(MU_ZZ) & is.null(MU_WW) ) {
	  solved_scoreeqs = nleqslv::nleqslv(x=c(0,0,0,0), fn=score_eqs_no_Z_no_W)
	}
	
	return(solved_scoreeqs)
}


