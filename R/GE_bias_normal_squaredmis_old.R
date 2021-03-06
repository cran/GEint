#' GE_bias_normal_squaredmis_old.R
#'
#' A function to calculate the bias in testing for GxE interaction, making many more
#' assumptions than GE_bias_old().  The additional assumptions are added to simplify the process
#' of calculating/estimating many higher order moments which the user may not be familiar with. \cr
#' The following assumptions are made: \cr
#' (1) All fitted covariates besides G (that is, E, all Z, and all W) have a marginal standard 
#' normal distribution with mean 0 and variance 1.  This corresponds to the case of the researcher
#' standardizing all of their fitted covariates. \cr
#' (2) G is generated by means of thresholding two independent normal RVs and is centered to have mean 0.
#' (3) The joint distributions of E, Z, W, and the thresholded variables underlying G can be described
#' by a multivariate normal distribution. \cr
#' (4) The misspecification is of the form f(E)=h(E)=E^2, and M_j=W_j^2 for all j. In particular,
#' W always has the same length as M here. \cr
#' 
#' @param beta_list A list of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' If Z or M is a vector, then beta_Z and beta_M should be vectors.
#' If Z or M is not in the model (i.e. all covariates other than G+E 
#' have been specified incorrectly, or all covariates other than G+E have been 
#' specified correctly, or the only covariates are G+E), then set beta_Z=0 and/or beta_M=0.
#' @param rho_list A list of the 6 pairwise covariances between the
#' covariates.  These should be in the order (1) cov_GE (2) cov_GZ (3) cov_EZ
#' (4) cov_GW (5) cov_EW (6) cov_ZW.
#' Again if Z or W are vectors then terms like cov_GZ should be vectors (in the order
#' cov(G,Z_1),...,cov(G,Z_p)) where Z is of dimension p, and similarly for W.
#' If Z or M are vectors, then cov_ZW should be a vector in the order (cov(Z_1,W_1),...,cov(Z_1,W_q),
#' cov(Z_2,W_1),........,cov(Z_p,W_q) where Z is a vector of length p and W is a vector of length q.
#' If Z or M are not in the model then treat them as the constant 0.
#' So for example if Z is not in the model and M (and therefore W) is a vector of length 2, we would
#' have cov_EZ=0 and cov(ZW) = (0,0).
#' @param cov_Z Only specify this if Z is a vector, gives the covariance matrix of Z (remember by assumption
#' Z has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
#' of the vector.
#' @param cov_W Only specify this if W is a vector, gives the covariance matrix of W (remember by assumption
#' W has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
#' of the vector.
#' @param prob_G Probability that each allele is equal to 1.  Since each SNP has
#' two alleles, the expectation of G is 2*prob_G.
#' 
#' @return A list with the elements:
#' \item{alpha_list}{The asymptotic values of the fitted coefficients alpha.}
#' \item{beta_list}{The same beta_list that was given as input.}
#' \item{cov_list}{The list of all covariances (both input and calculated) for use with GE_nleqslv() 
#'	and GE_bias().}
#' \item{cov_mat_list}{List of additionally calculated covariance matrices for use with GE_nleqslv()
#' and GE_bias().}
#' \item{mu_list}{List of calculated means for f(E), h(E), Z, M, and W for use with GE_nleqslv() 
#' and GE_bias().}
#' \item{HOM_list}{List of calculated Higher Order Moments for use with GE_nleqslv() and GE_bias().}
#'
#' @export
#' @examples 
#' GE_bias_normal_squaredmis_old( beta_list=as.list(runif(n=6, min=0, max=1)), 
#'							rho_list=as.list(rep(0.3,6)), prob_G=0.3)

GE_bias_normal_squaredmis_old <- function(beta_list, rho_list, prob_G, cov_Z=NULL, cov_W=NULL)
{
  # Need survival function.
  surv <- function(x) {1-pnorm(x)}
  
  # Record some initial quantities
  rho_GE <- rho_list[[1]]; rho_GZ <- rho_list[[2]]; rho_EZ <- rho_list[[3]]
  rho_GW <- rho_list[[4]]; rho_EW <- rho_list[[5]]; rho_ZW <- rho_list[[6]]
  w <- qnorm(1-prob_G)					
  r_GE <- rho_GE / (2*dnorm(w))	
  r_GZ <- rho_GZ / (2*dnorm(w))
  r_GW <- rho_GW / (2*dnorm(w))
  	
  beta_0 <- beta_list[[1]]; beta_G <- beta_list[[2]]; beta_E <- beta_list[[3]]
  beta_I <- beta_list[[4]]; BETA_Z <- beta_list[[5]]; BETA_M <- beta_list[[6]]
  
  # Even if Z or M/W is 0 keep the num at 1.
  num_W <- length(beta_list[[6]])
  num_Z <- length(beta_list[[5]])
  
  
  # Some error checking, make sure the covariance matrix is ok
  translated_inputs <- GE_translate_inputs_old(beta_list, rho_list, prob_G, cov_Z, cov_W)
  sig_mat <- translated_inputs$sig_mat_total
  sig_mat_ZZ <- translated_inputs$sig_mat_ZZ
  sig_mat_WW <- translated_inputs$sig_mat_WW
  
  # Check for Z and W.
  # Some obvious and by assumption means.
  if (is.null(sig_mat_ZZ)) {
    MU_Z <- 0
  } else {
    MU_Z <- rep(0, num_Z)
  }
  if (is.null(sig_mat_WW)) {
    MU_W <- 0
    MU_M <- 0
  } else {
    MU_M <- rep(1, num_W)		
    MU_W <- rep(0, num_W)	
  }
	
  # More obvious and by assumption.
  mu_f <- 1
  mu_h <- 1
  
  # Now calculate other harder, necessary terms that have 
  # been determined by our assumptions + inputs
  ########################
  # Covariances
  mu_GE <- rho_GE
  mu_Gf <- 2*r_GE^2*w*dnorm(w) + 2*surv(w) - 2*prob_G
  mu_Gh <- mu_Gf
  mu_GG <- 2*prob_G*(1-prob_G)
  mu_EE <- 1
  mu_Ef <- 0
  MU_GZ <- rho_GZ  	# Vector
  MU_GW <- rho_GW		# Vector
  MU_EM <- 	rep(0, num_W)				# Vector, in particular because third moment of W is 0
  MU_fW <- 	rep(0, num_W)		# Vector
  MU_EW <- rho_EW			# Vector
  MU_EZ <- rho_EZ	
  MU_fZ <- 	rep(0, num_Z)	
  
  # Depends on if M exists.
  if (is.null(sig_mat_WW)) {
    MU_GM <- 0
  } else {
    MU_GM <- 	2*r_GW^2*w*dnorm(w) + 2*surv(w) - 2*prob_G	# Vector, see gen_cor_bin_normal for explanation
  } 

  ########################
  # Matrix covariances
  # MU_ZW is not the same as MU_WZ because the dimensions of the matrix are not the same!
  # Remember the covariances in rho_ZW are in the order cov(Z_1,W_1), cov(Z_1,W_2), ..., cov(Z_2,W_1),...
  MU_ZW <- matrix(data=rho_ZW, nrow=num_Z, ncol=num_W, byrow=TRUE)	# Matrix	 
  MU_WZ <- t(MU_ZW)
  MU_ZM <- matrix(data=0, nrow=num_Z, ncol=num_W) 		# Matrix
  MU_WM <- matrix(data=0, nrow=num_W, ncol=num_W) 			# Matrix
  MU_ZZ <- sig_mat_ZZ 		# Matrix
  MU_WW <- sig_mat_WW		# Matrix
  
  
  ########################
  # Higher order moments with G+E
  # We need as intermediate quantities E[G_1E], E[G_1E^2], E[G_1E^3], E[G_1G_2E], E[G_1G_2E^2], E[G1EZ], 
  # E[G1EW], E[G1EW^2], E[G1WE^2], E[G1ZE^2], E[G1G2E^2]
  mu_G1_E <- r_GE*dnorm(w)
  mu_G1_EE <- r_GE^2*w*dnorm(w) + surv(w)
  mu_G1_EEE <- r_GE^3*w^2*dnorm(w) - r_GE^3*dnorm(w) + 3*r_GE*dnorm(w)
  
  # E[G1G2E] requires numerical integration
  temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  f_G1_G2_E <- function(x,w,r_GE) {
  	x*dnorm(x)*mvtnorm::pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  }
  mu_G1_G2_E <- pracma::quadinf(f=f_G1_G2_E, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  # E[G1G2E^2] requires numerical integration
  temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  f_G1_G2_EE <- function(x,w,r_GE) {
  	x^2*dnorm(x)*mvtnorm::pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  }
  mu_G1_G2_EE <- pracma::quadinf(f=f_G1_G2_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  # E[G1G2E^3] requires numerical integration
  temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  f_G1_G2_EEE <- function(x,w,r_GE) {
  	x^3*dnorm(x)*mvtnorm::pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  }
  mu_G1_G2_EEE <- pracma::quadinf(f=f_G1_G2_EEE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
    
  # See gen_cor_bin_normal to see how to do these
  mu_GGE <- 2*mu_G1_E + 2*mu_G1_G2_E - 8*prob_G*mu_G1_E
  mu_GGh <- 2*mu_G1_EE + 2*mu_G1_G2_EE + 4*prob_G^2*1 - 8*prob_G*mu_G1_EE
  mu_GEE <- mu_Gf
  mu_GEf <- 2*(r_GE^3*w^2*dnorm(w) - r_GE^3*dnorm(w) + 3*r_GE*dnorm(w))
  mu_GEh <- mu_GEf
  
  mu_GGEE <- 2*mu_G1_EE + 2*mu_G1_G2_EE + 4*prob_G^2*1 - 8*prob_G*mu_G1_EE
  mu_GGEf <- 2*mu_G1_EEE + 2*mu_G1_G2_EEE + 4*prob_G^2*0 - 8*prob_G*mu_G1_EEE
  mu_GGEh <- mu_GGEf
  
  ##########################################
  # Harder ones involving Z and W
  
  # E[G1EZ] requires numerical integration
  f_G1_E_Z <- function(x, w, r_EZ, r_GE, r_GZ) {
  	( r_EZ * x * surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / sqrt(1-r_GE^2) ) * 
  		(r_GZ-r_GE*r_EZ) / sqrt(1-r_GE^2) ) * x* dnorm(x)
  }
  if (is.null(sig_mat_ZZ)) {
    mu_G1_E_Z <- 0
  } else {
    mu_G1_E_Z <- rep(NA, num_Z)
    for (i in 1:num_Z) {
  	  mu_G1_E_Z[i] <- pracma::quadinf(f= f_G1_E_Z, xa=-Inf, xb=Inf, w=w, r_EZ=rho_EZ[i], r_GE=r_GE, r_GZ=r_GZ[i])$Q
    }
  }
  
  # E[G1EW] requires numerical integration
  f_G1_E_W <- function(x, w, r_EW, r_GE, r_GW) {
  	( r_EW * x * surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / sqrt(1-r_GE^2) ) * 
  		(r_GW-r_GE*r_EW) / sqrt(1-r_GE^2) ) * x* dnorm(x)
  }
  if (is.null(sig_mat_WW)) {
    mu_G1_E_W <- 0
  } else {
    mu_G1_E_W <- rep(NA, num_W)
    for (i in 1:num_W) {
  	  mu_G1_E_W[i] <- pracma::quadinf(f= f_G1_E_W, xa=-Inf, xb=Inf, w=w, r_EW=rho_EW[i], r_GE=r_GE, r_GW=r_GW[i])$Q
    }
  }
  
  # E[G1EW^2] requires numerical integration
  f_G1_E_WW <- function(x, w, r_GE, r_GW, r_EW) {
  	( r_EW * x* surv( (w-x*r_GW) / sqrt(1-r_GW^2) ) + dnorm( (w-r_GW*x) / 
  		sqrt(1-r_GW^2) ) * (r_GE-r_GW*r_EW) / sqrt(1-r_GW^2) ) * x^2 * dnorm(x)
  }
  if (is.null(sig_mat_WW)) {
    mu_G1_E_WW <- 0
  } else {
    mu_G1_E_WW <- rep(NA, num_W)
    for (i in 1:num_W) {
  	  mu_G1_E_WW[i] <- pracma::quadinf(f=f_G1_E_WW, xa=-Inf, xb=Inf, w=w , r_GE=r_GE, r_GW=r_GW[i], r_EW=rho_EW[i])$Q
    }
  }
  
  # E[G1WE^2] requires numerical integration
  f_G1_W_EE <- function(x, w, r_GE, r_GW, r_EW) {
  	( r_EW * x* surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / 
  		sqrt(1-r_GE^2) ) * (r_GW-r_GE*r_EW) / sqrt(1-r_GE^2) ) * x^2 * dnorm(x)
  }
  if (is.null(sig_mat_WW)) {
    mu_G1_W_EE <- 0
  } else {
    mu_G1_W_EE <- rep(NA, num_W)
    for (i in 1:num_W) {
  	  mu_G1_W_EE[i] <- pracma::quadinf(f=f_G1_W_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE, r_GW=r_GW[i], r_EW=rho_EW[i])$Q
    }
  }
  
  # E[G1ZE^2] requires numerical integration
  f_G1_Z_EE <- function(x, w, r_GE, r_GZ, r_EZ) {
  	( r_EZ * x* surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / 
  		sqrt(1-r_GE^2) ) * (r_GZ-r_GE*r_EZ) / sqrt(1-r_GE^2) ) * x^2 * dnorm(x)
  }
  if (is.null(sig_mat_ZZ)) {
    mu_G1_Z_EE <- 0
  } else {
    mu_G1_Z_EE <- rep(NA, num_Z)
    for (i in 1:num_Z) {
  	  mu_G1_Z_EE[i] <- pracma::quadinf(f=f_G1_Z_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE, r_GZ=r_GZ[i], r_EZ=rho_EZ[i])$Q
    }
  }

  # See gen_cor_bin_normal to see how to do these (vectors)
  MU_GEZ <- 2*mu_G1_E_Z - 2*prob_G*rho_EZ			# Vector
  MU_GEW <- 2*mu_G1_E_W	- 2*prob_G*rho_EW			# Vector
  MU_GEM <-	2*mu_G1_E_WW							# Vector
  MU_GhW <- 2*mu_G1_W_EE
  MU_GhZ <- 2*mu_G1_Z_EE
 
  
  ########################
  # Some shortcut quantities
  A <- (mu_GE * MU_GZ / mu_GG - MU_EZ) / (mu_EE - mu_GE^2/mu_GG)
  B <- (mu_GE * MU_GW / mu_GG - MU_EW) / (mu_EE - mu_GE^2/mu_GG)
  
  # O will be set to 0 if no Z
  if (is.null(MU_ZZ)) {
    O <- 0
    solve_O <- 0
  } else {
    O <- MU_Z%*%t(MU_Z) + MU_GZ%*%t(MU_GZ)/mu_GG - MU_ZZ - A %*% t(MU_EZ - MU_GZ*mu_GE/mu_GG)
    solve_O <- solve(O)
  }
  
  C <- (B %*% t(MU_EZ - MU_GZ*mu_GE/mu_GG) - MU_W%*%t(MU_Z) - MU_GW%*%t(MU_GZ)/mu_GG + MU_WZ) %*% solve_O
  
  # Q will be 0 if no W
  if ( is.null(MU_WW) ) {
    Q <- 0
    solve_Q <- 0
  } else {
    Q <- MU_W%*%t(MU_W) + MU_GW%*%t(MU_GW)/mu_GG - MU_WW + B %*% t(MU_GW*mu_GE/mu_GG - MU_EW) + 
      C %*% ( MU_Z%*%t(MU_W) + MU_GZ%*%t(MU_GW)/mu_GG - MU_ZW + A %*% t(MU_GW*mu_GE/mu_GG - MU_EW) )
    solve_Q <- solve(Q)
  }
  
  D <- (mu_GE * mu_GGE / mu_GG - mu_GEE) / (mu_EE - mu_GE^2 / mu_GG)
  E <- t(MU_GEZ - MU_Z*mu_GE - MU_GZ*mu_GGE/mu_GG + D*(MU_EZ - MU_GZ*mu_GE/mu_GG)) %*% solve_O
  EFF <- ( t(MU_W*mu_GE + MU_GW*mu_GGE/mu_GG - MU_GEW + D*(MU_GW * mu_GE / mu_GG - MU_EW)) + 
             E %*% (A %*% t(MU_GW*mu_GE/mu_GG - MU_EW) + MU_Z%*%t(MU_W) + MU_GZ%*%t(MU_GW)/mu_GG - MU_ZW) ) %*% solve_Q
  
  
  # Solve for \alpha_I
  alpha_I_num <- beta_E * (-mu_f*mu_GE - mu_Gf*mu_GGE/mu_GG + mu_GEf + D * (mu_Ef - mu_Gf*mu_GE/mu_GG)) +
    beta_E * E %*% (-mu_f*MU_Z - MU_GZ*mu_Gf/mu_GG + MU_fZ + A * (mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_I * (-mu_Gh*mu_GE - mu_GGh*mu_GGE/mu_GG + mu_GGEh + D * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    beta_I * E %*% (MU_GhZ -mu_Gh*MU_Z - MU_GZ*mu_GGh/mu_GG + A * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    t(MU_GEM - MU_M*mu_GE - MU_GM*mu_GGE/mu_GG + D * (MU_EM - MU_GM*mu_GE/mu_GG)) %*% BETA_M + 
    E %*% (A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M - 
    beta_E * EFF %*% (-mu_f*MU_W - MU_GW*mu_Gf/mu_GG + MU_fW + B %*% as.matrix(mu_Ef - mu_Gf*mu_GE/mu_GG)) - 
    beta_E * EFF %*% C %*% (-mu_f*MU_Z - mu_Gf*MU_GZ/mu_GG + MU_fZ + A %*% as.matrix(mu_Ef - mu_Gf*mu_GE/mu_GG)) - 
    beta_I * EFF %*% (-mu_Gh*MU_W - MU_GW*mu_GGh/mu_GG + MU_GhW + B %*% as.matrix(mu_GEh - mu_GGh*mu_GE/mu_GG)) -
    beta_I * EFF %*% C %*% (MU_GhZ - MU_Z*mu_Gh - MU_GZ*mu_GGh/mu_GG + A %*% as.matrix(mu_GEh - mu_GGh*mu_GE/mu_GG)) - 
    EFF %*% ( -MU_W%*%t(MU_M) - MU_GW%*%t(MU_GM)/mu_GG + MU_WM + B %*% t(MU_EM - MU_GM*mu_GE/mu_GG) ) %*% BETA_M - 
    EFF %*% C %*% ( A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M
  
  alpha_I_denom <- EFF %*% ( MU_W*mu_GE + MU_GW*mu_GGE/mu_GG - MU_GEW + B * (mu_GGE*mu_GE/mu_GG - mu_GEE) ) +
    EFF %*% C %*% ( MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG - MU_GEZ + A * (mu_GGE*mu_GE/mu_GG - mu_GEE) ) - 
    ( mu_GE^2 + mu_GGE^2/mu_GG - mu_GGEE + D * (mu_GGE*mu_GE/mu_GG - mu_GEE) ) - 
    E %*% ( MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG - MU_GEZ + A * (mu_GGE*mu_GE/mu_GG - mu_GEE) )
  
  alpha_I <- alpha_I_num / alpha_I_denom
  
  R <- beta_E * (-MU_W*mu_f - MU_GW*mu_Gf/mu_GG + MU_fW + B * (mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_E * C %*% (-mu_f*MU_Z - MU_GZ*mu_Gf/mu_GG + MU_fZ + A %*% as.matrix(mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_I * (-MU_W*mu_Gh - MU_GW*mu_GGh/mu_GG + MU_GhW + B * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    beta_I * C %*% (MU_GhZ - MU_Z*mu_Gh - MU_GZ*mu_GGh/mu_GG + A * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    ( B %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_W%*%t(MU_M) - MU_GW%*%t(MU_GM)/mu_GG + MU_WM) %*% BETA_M + 
    C %*% ( A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M + 
    alpha_I * (MU_W*mu_GE + MU_GW*mu_GGE/mu_GG - MU_GEW + B * (mu_GGE*mu_GE/mu_GG - mu_GEE)) + 
    as.numeric(alpha_I) * C %*% (MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG - MU_GEZ + A*(mu_GGE*mu_GE/mu_GG - mu_GEE))
  
  ALPHA_W <- - solve_Q %*% R
  
  P <- beta_E * (-MU_Z*mu_f - MU_GZ*mu_Gf/mu_GG + MU_fZ + A * (mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_I * (MU_GhZ - MU_Z*mu_Gh - MU_GZ*mu_GGh/mu_GG + A * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    alpha_I * (MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG - MU_GEZ + A * (mu_GGE*mu_GE/mu_GG - mu_GEE)) + 
    ( A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M + 
    ( A %*% t(MU_GW*mu_GE/mu_GG - MU_EW) + MU_Z%*%t(MU_W) + MU_GZ%*%t(MU_GW)/mu_GG - MU_ZW) %*% ALPHA_W
  
  Bz_Az <- solve_O %*% P
  ALPHA_Z <- BETA_Z - solve_O %*% P
  
  alpha_E <- ( beta_E * (mu_Ef - mu_Gf*mu_GE/mu_GG) + beta_I * (mu_GEh - mu_GGh*mu_GE/mu_GG) +  
                 alpha_I * (mu_GGE*mu_GE/mu_GG - mu_GEE) + t(MU_EZ - MU_GZ*mu_GE/mu_GG) %*% Bz_Az + 
                 t(MU_EM - MU_GM*mu_GE/mu_GG) %*% BETA_M + t(MU_GW*mu_GE/mu_GG - MU_EW) %*% ALPHA_W ) /
    (mu_EE - mu_GE^2/mu_GG)
  
  Bg_Ag <- ( alpha_E*mu_GE - beta_E*mu_Gf + alpha_I*mu_GGE - beta_I*mu_GGh - t(MU_GZ) %*% Bz_Az + 
               t(MU_GW) %*% ALPHA_W - t(MU_GM) %*% BETA_M ) / mu_GG
  alpha_G <- beta_G - Bg_Ag
  
  
  alpha_0 <- beta_0 + beta_E*mu_f + beta_I*mu_Gh - alpha_I*mu_GE + t(MU_Z) %*% Bz_Az +
    t(MU_M) %*% BETA_M - t(MU_W) %*% ALPHA_W
  
  # Return 
  return(list(alpha_list=list(alpha_0, alpha_G, alpha_E, alpha_I, ALPHA_Z, ALPHA_W),
              beta_list = list(beta_0, beta_G, beta_E, beta_I, BETA_Z, BETA_M),
              cov_list = list(mu_GG, mu_GE, mu_Gf, mu_Gh, MU_GZ, MU_GM, MU_GW, mu_EE,
              mu_Ef, MU_EZ, MU_EM, MU_EW, MU_fZ, MU_fW),
              cov_mat_list = list(MU_ZZ, MU_WW, MU_ZW, MU_WZ, MU_ZM, MU_WM),
              mu_list = list(mu_f, mu_h, rep(0,num_Z), MU_M, rep(0,num_W)),
              HOM_list = list(mu_GGE, mu_GGh, mu_GEE, mu_GEf, mu_GEh, MU_GEZ, MU_GEM, MU_GEW,
              MU_GhW, MU_GhZ, mu_GGEE, mu_GGEf, mu_GGEh)))
}