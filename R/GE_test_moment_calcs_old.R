#' GE_test_moment_calcs_old.R
#'
#' Test function mostly for internal use to ensure the higher order moments (covariances)
#' calculated in GE_bias_normal_squaredmis_old() are correct. Will give warning messages if
#' some calculations appear to be incorrect.  If receive warning messages, run again, and
#' if still receive the same warning messages, something may be wrong.
#'
#' @param rho_list A list of the 6 pairwise covariances between the
#' covariates.  These should be in the order (1) cov_GE (2) cov_GZ (3) cov_EZ
#' (4) cov_GW (5) cov_EW (6) cov_ZW.  If Z and/or W include multiple covariates, then
#' terms like cov_GZ should be a vector.  
#' If Z and/or M/W do not exist in your model, then treat them as constants 0. For example,
#' if Z doesn't exist and W includes 2 covariates, then set cov(EZ) = 0 and cov(ZW) = (0,0).
#' @param cov_Z Only used if Z is a vector, gives the covariance matrix of Z (remember by assumption
#' Z has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
#' of the vector. If Z or M are vectors, then cov_ZW should be a vector in the order
#' (cov(Z_1,W_1),cov(Z_1,W_2),...,cov(Z_1,W_q),cov(Z_2,W_1),........,cov(Z_p,W_q) where Z is 
#' a vector of length p and W is a vector of length q.
#' @param cov_W Only used if W is a vector, gives the covariance matrix of W (remember by assumption
#' W has mean 0 and variance 1).  The (i,j) element of the matrix should be the (i-1)(i-2)/2+j element
#' of the vector.
#' @param prob_G Probability that each allele is equal to 1.  Since each SNP has
#' two alleles, the expectation of G is 2*prob_G.
#' @param num_sum Number of subjects to do the simulation with.
#' @param test_threshold How much margin for error on tests?
#'
#' @return Nothing
#'
#' @export
#' @keywords internal
#' @examples 
#' GE_test_moment_calcs_old(beta_list=as.list(runif(n=6, min=0, max=1)), 
#'			rho_list=as.list(rep(0.3,6)), prob_G=0.3)

GE_test_moment_calcs_old <- function(beta_list, rho_list, prob_G, cov_Z=NULL, cov_W=NULL, num_sub=2000000, test_threshold=0.003)
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
									
	  # Generate test data
	  test_data <- mvtnorm::rmvnorm(n=num_sub, sigma=sig_mat)
	
	  # Get the individual components
  	G1 <- as.numeric(test_data[,1] > w)
  	G2 <- as.numeric(test_data[,2] > w)
  	G <- G1 + G2 - 2*prob_G
	  E <- test_data[,3]
	  
	  # Figure out how many Z and W we have.
	  # Even if Z or M/W = 0 then we keep their number at 1 to facilitate ease of testing.
	  num_Z <- length(beta_list[[5]])
	  num_W <- length(beta_list[[6]])
	  # And how to build their covariate vectors.
	  if (is.null(sig_mat_ZZ)) {
	    num_build_Z <- 0            # Introduce this here just for building the covariate vectors 
	    Z <- matrix(data=0, nrow=length(E), ncol=1)
	  } else {
	    Z <- as.matrix(test_data[,4:(3+num_Z)])
	    num_build_Z <- num_Z
	  }
	  if (is.null(sig_mat_WW)) {
	    num_build_W <- 0
	    W <- matrix(data=0, nrow=length(E), ncol=1)
	  } else {
	    num_build_W <- num_W
	    W <- as.matrix(test_data[,(4+num_build_Z):(3+num_build_Z+num_W)])
	  }
	 
  	# Start testing with calculated G covariances
  	temp_test <- rho_GE - cov(G,E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with rho_GE')}
  	
  	for (i in 1:num_Z)
  	{
  		temp_test <- rho_GZ[i] - cov(G,Z[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with rho_GZ')}
  	}
  	
  	for (i in 1:num_W)
  	{
  		temp_test <- rho_GW[i] - cov(G,W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with rho_GW')}
  	}
  	
  	 ########################
 	  # More covariances
  	mu_GE <- rho_GE
 	  mu_Gf <- 2*r_GE^2*w*dnorm(w) + 2*surv(w) - 2*prob_G
 	  mu_Gh <- mu_Gf
 	  mu_GG <- 2*prob_G*(1-prob_G)
 	  MU_GZ <- rho_GZ  	# Vector
 	  MU_GW <- rho_GW		# Vector
 	  MU_EM <- 	rep(0, num_W)				# Vector, in particular because third moment of W is 0
 	  MU_EZ <- rho_EZ			# Vector
 	  MU_EW <- rho_EW			# Vector
 	  mu_EE <- 1
 	  mu_Ef <- 0
 	  MU_fZ <- 	rep(0, num_Z)	# Vector
 	  MU_fW <- 	rep(0, num_W)		# Vector
 	  
 	  # Depends on if M exists
 	  if (num_build_W != 0) {
 	    MU_GM <- 	2*r_GW^2*w*dnorm(w) + 2*surv(w) - 2*prob_G	# Vector, see gen_cor_bin_normal for explanation
 	  } else {
 	    MU_GM <- 0
 	  }
 	  
 	  # Test the above
  	temp_test <- mu_GE - mean(G*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GE')}
  	
  	temp_test <- mu_Gf - mean(G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_Gf')}
  	
  	temp_test <- mu_Gh - mean(G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_Gh')}
  	
  	temp_test <- mu_GG - mean(G*G)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GG')}
  	
  	for (i in 1:num_Z)
  	{
  		temp_test <- MU_GZ[i] - mean(G*Z[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_GZ')}
  	}
  	
  	for (i in 1:num_W)
  	{
  		temp_test <- MU_GW[i] - mean(G*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_GW')}
  	}
  	
  	for (i in 1:num_W)
  	{
  		temp_test <- MU_GM[i]- mean(G*W[,i]*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_GM')}
  	}
  	
  	for (i in 1:num_W)
  	{
  		temp_test <- MU_EM[i] - mean(E*W[,i]*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_EM')}
  	}
  	
  	for (i in 1:num_Z)
  	{
  		temp_test <- MU_EZ[i] - mean(E*Z[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_EZ')}
  	}
  	
  	for (i in 1:num_W)
  	{
  		temp_test <- MU_EW[i] - mean(E*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_EW')}
  	}
  	
  	for (i in 1:num_Z)
  	{
  		temp_test <- MU_fZ[i] - mean(E*E*Z[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_fZ')}
  	}
  	
  	for (i in 1:num_W)
  	{
  		temp_test <- MU_fW[i] - mean(E*E*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_fW')}
  	}
  	
   	########################
 	  # Matrix covariances
 	  MU_ZW <- matrix(data=rho_ZW, nrow=num_Z, ncol=num_W, byrow=TRUE)	# Matrix	 
  	MU_WZ <- t(MU_ZW)	 
 	  MU_ZM <- matrix(data=0, nrow=num_Z, ncol=num_W) 		# Matrix
 	  MU_WM <- matrix(data=0, nrow=num_W, ncol=num_W) 			# Matrix
  	MU_ZZ <- sig_mat_ZZ 		# Matrix
  	MU_WW <- sig_mat_WW		# Matrix
  	
  	# Check them
  	for (i in 1:num_Z)
  	{
  		for (j in 1:num_W)
  		{
  			temp_test <- MU_ZW[i,j] - mean(Z[,i]*W[,j])
  			if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZW')}
  		}
  	}
  	
  	for (i in 1:num_Z)
  	{
  		for (j in 1:num_W)
  		{
  			temp_test <- MU_ZM[i,j] - mean(Z[,i]*W[,j]*W[,j])
  			if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZM')}
  		}
  	}
  	  
  	for (i in 1:num_W)
  	{
  		for (j in 1:num_W)
  		{
  			temp_test <- MU_WM[i,j] - mean(W[,i]*W[,j]*W[,j])
  			if (abs(temp_test) > test_threshold) {warning('Problem with mu_WM')}
  		}
  	}	
  	
  	# If Z exists.
  	if ( !is.null(MU_ZZ) ) {
  	  for (i in 1:num_Z)
  	  {
  	  	for (j in 1:num_Z)
  	  	{
  	  		temp_test <- MU_ZZ[i,j] - mean(Z[,i]*Z[,j])
  	  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZZ')}
  	  	}
  	  }	
  	}
  	
  	# If W exists.
  	if ( !is.null(MU_WW) ) {
  	  for (i in 1:num_W)
  	  {
  	  	for (j in 1:num_W)
  	  	{
  	  		temp_test <- MU_WW[i,j] - mean(W[,i]*W[,j])
  	  		if (abs(temp_test) > test_threshold) {warning('Problem with mu_WM')}
  	  	}
  	  }
  	}
  	
    ########################
  	# Higher order moments, intermediate quantities
  	mu_G1_E <- r_GE*dnorm(w)
  	mu_G1_EE <- r_GE^2*w*dnorm(w) + surv(w)
  	mu_G1_EEE <- r_GE^3*w^2*dnorm(w) - r_GE^3*dnorm(w) + 3*r_GE*dnorm(w)
  	
  	temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  	f_G1_G2_E <- function(x,w,r_GE) {
  		x*dnorm(x)*mvtnorm::pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  	}
  	mu_G1_G2_E <- pracma::quadinf(f=f_G1_G2_E, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  	temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  	f_G1_G2_EE <- function(x,w,r_GE) {
  		x^2*dnorm(x)*mvtnorm::pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  	}
  	mu_G1_G2_EE <- pracma::quadinf(f=f_G1_G2_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  	temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  	f_G1_G2_EEE <- function(x,w,r_GE) {
  		x^3*dnorm(x)*mvtnorm::pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  	}
  	mu_G1_G2_EEE <- pracma::quadinf(f=f_G1_G2_EEE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  	
  	# Higher order moments, see gen_cor_bin_normal to see how to do these
  	mu_GGE <- 2*mu_G1_E + 2*mu_G1_G2_E - 8*prob_G*mu_G1_E
  	mu_GGh <- 2*mu_G1_EE + 2*mu_G1_G2_EE + 4*prob_G^2*1 - 8*prob_G*mu_G1_EE
  	mu_GEE <- mu_Gf
  	mu_GEf <- 2*(r_GE^3*w^2*dnorm(w) - r_GE^3*dnorm(w) + 3*r_GE*dnorm(w))
  	mu_GEh <- mu_GEf
  
  	mu_GGEE <- 2*mu_G1_EE + 2*mu_G1_G2_EE + 4*prob_G^2*1 - 8*prob_G*mu_G1_EE
  	mu_GGEf <- 2*mu_G1_EEE + 2*mu_G1_G2_EEE + 4*prob_G^2*0 - 8*prob_G*mu_G1_EEE
  	mu_GGEh <- mu_GGEf
  
 	  ##################################
  	# Check	
  	temp_test <- mu_G1_E - mean(G1*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E')}
  	
  	temp_test <- mu_G1_EE - mean(G1*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_EE')}
  	
  	temp_test <- mu_G1_EEE - mean(G1*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_EEE')}
  	
  	temp_test <- mu_G1_G2_E - mean(G1*G2*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_G2_E')}
  	
  	temp_test <- mu_G1_G2_EE - mean(G1*G2*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_G2_EE')}
  	
  	temp_test <- mu_G1_G2_EEE - mean(G1*G2*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_G2_EEE')}

  	temp_test <- mu_GGE - mean(G*G*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GGE')}
  	
  	temp_test <- mu_GGh - mean(G*G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GGh')}
  	
  	temp_test <- mu_GEE - mean(G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEE')}
  	
  	temp_test <- mu_GEf - mean(G*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEf')}
  	
  	temp_test <- mu_GEh - mean(G*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEh')}
  
  
  
 	  ##############
  	# Harder intermediate quantities involving Z and W
  	f_G1_E_Z <- function(x, w, r_EZ, r_GE, r_GZ) {
  		( r_EZ * x * surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / sqrt(1-r_GE^2) ) * 
  			(r_GZ-r_GE*r_EZ) / sqrt(1-r_GE^2) ) * x* dnorm(x)
  	}
  	if ( !is.null(sig_mat_ZZ) ) {
  	  mu_G1_E_Z <- rep(NA, num_Z)
  	  for (i in 1:num_Z) {
  		  mu_G1_E_Z[i] <- pracma::quadinf(f= f_G1_E_Z, xa=-Inf, xb=Inf, w=w, r_EZ=rho_EZ[i], r_GE=r_GE, r_GZ=r_GZ[i])$Q
  		  temp_test <- mu_G1_E_Z[i] - mean(G1*E*Z[,i])
  		  if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E_Z')}
  	  }
  	} else {
  	  mu_G1_E_Z <- 0
  	}
  
  	mu_G1_E_W <- rep(NA, num_W)
  	f_G1_E_W <- function(x, w, r_EW, r_GE, r_GW) {
  		( r_EW * x * surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / sqrt(1-r_GE^2) ) * 
  			(r_GW-r_GE*r_EW) / sqrt(1-r_GE^2) ) * x* dnorm(x)
  	}
  	if ( !is.null(sig_mat_WW) ) {
  	  for (i in 1:num_W) {
  		  mu_G1_E_W[i] <- pracma::quadinf(f= f_G1_E_W, xa=-Inf, xb=Inf, w=w, r_EW=rho_EW[i], r_GE=r_GE, r_GW=r_GW[i])$Q
  		  temp_test <- mu_G1_E_W[i] - mean(G1*E*W[,i])
  		  if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E_W')}
  	  }
  	} else {
  	  mu_G1_E_W <- 0
  	}

	  mu_G1_E_WW <- rep(NA, num_W)
  	f_G1_E_WW <- function(x, w, r_GE, r_GW, r_EW) {
  		( r_EW * x* surv( (w-x*r_GW) / sqrt(1-r_GW^2) ) + dnorm( (w-r_GW*x) / 
  			sqrt(1-r_GW^2) ) * (r_GE-r_GW*r_EW) / sqrt(1-r_GW^2) ) * x^2 * dnorm(x)
  	}
  	if ( !is.null(sig_mat_WW) ) {
  	  for (i in 1:num_W) {
  		  mu_G1_E_WW[i] <- pracma::quadinf(f=f_G1_E_WW, xa=-Inf, xb=Inf, w=w , r_GE=r_GE, r_GW=r_GW[i], r_EW=rho_EW[i])$Q
  		  temp_test <- mu_G1_E_WW[i] - mean(G1*E*W[,i]*W[,i])
  		  if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E_WW')}
  	  }
  	} else {
  	  mu_G1_E_WW <- 0
  	}
  
    mu_G1_W_EE <- rep(NA, num_W)
  	f_G1_W_EE <- function(x, w, r_GE, r_GW, r_EW) {
  		( r_EW * x* surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / 
  			sqrt(1-r_GE^2) ) * (r_GW-r_GE*r_EW) / sqrt(1-r_GE^2) ) * x^2 * dnorm(x)
  	}
  	if ( !is.null(sig_mat_WW) ) {
  	  for (i in 1:num_W) {
  		  mu_G1_W_EE[i] <- pracma::quadinf(f=f_G1_W_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE, r_GW=r_GW[i], r_EW=rho_EW[i])$Q
  		  temp_test <- mu_G1_W_EE[i] - mean(G1*W[,i]*E*E)
  		  if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_W_EE')}
  	  }
  	} else {
  	  mu_G1_W_EE <- 0
  	}

    mu_G1_Z_EE  <- rep(NA, num_Z)
  	f_G1_Z_EE <- function(x, w, r_GE, r_GZ, r_EZ) {
  		( r_EZ * x* surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / 
  			sqrt(1-r_GE^2) ) * (r_GZ-r_GE*r_EZ) / sqrt(1-r_GE^2) ) * x^2 * dnorm(x)
  	}	
  	if ( !is.null(sig_mat_ZZ) ) {
  	  for (i in 1:num_Z) {
  		  mu_G1_Z_EE[i] <- pracma::quadinf(f=f_G1_Z_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE, r_GZ=r_GZ[i], r_EZ=rho_EZ[i])$Q
  		  temp_test <- mu_G1_Z_EE[i] - mean(G1*Z[,i]*E*E)
  		  if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_Z_EE')}
  	  }
  	} else {
  	  mu_G1_Z_EE <- 0
  	}
  	
  	##############
  	# Check higher orders with Z and W
  	MU_GEZ <- 2*mu_G1_E_Z - 2*prob_G*rho_EZ			# Vector
  	MU_GEW <- 2*mu_G1_E_W	- 2*prob_G*rho_EW		# Vector
  	MU_GEM <-	2*mu_G1_E_WW				# Vector
  	MU_GhW <- 2*mu_G1_W_EE
  	MU_GhZ <- 2*mu_G1_Z_EE
  
  	# Test
 	  for (i in 1:num_Z) { 	
  		temp_test <- MU_GEZ[i] - mean(G*E*Z[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with MU_GEZ')}
  	}
  	
  	for (i in 1:num_W) {
  		temp_test <- MU_GEW[i] - mean(G*E*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with MU_GEW')}
  	}
  	
  	for (i in 1:num_W) {
  		temp_test <- MU_GEM[i] - mean(G*E*W[,i]*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with MU_GEM')}
  	}
  	
  	for (i in 1:num_W) {
  		temp_test <- MU_GhW[i] - mean(G*E*E*W[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with MU_GhW')}
  	}
  	
  	for (i in 1:num_Z) {
  		temp_test <- MU_GhZ[i] - mean(G*E*E*Z[,i])
  		if (abs(temp_test) > test_threshold) {warning('Problem with MU_GhZ')}
	}
  	
  	cat('Done')
}




