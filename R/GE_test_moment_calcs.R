#' GE_test_moment_calcs.R
#'
#' Test function mostly for internal use to ensure the higher order moments (covariances)
#' calculated in GE_bias_normal_squaredmis() are correct. Will give warning messages if
#' some calculations appear to be incorrect.  If receive warning messages, run again, and
#' if still receive the same warning messages, something may be wrong.
#'
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
#' @param cov_Z Should be a matrix equal to cov(Z) or NULL. Must be specified if Z is a vector.  
#' @param cov_W Should be a matrix equal to cov(W) or NULL. Must be specified if W is a vector.  
#' @param corr_G Should be a matrix giving the *pairwise correlations* between each SNP
#' in the set, or NULL. Must be specified if G is a vector.  For example, the [2,3] element
#' of the matrix would be the pairwise correlation between SNP2 and SNP3.
#' @param num_sub Number of subjects to simulate/test with.  
#' @param test_threshold How close does our simulated value have to be to the analytic value? 
#'
#' @return Nothing
#'
#' @export
#' @keywords internal
#' @examples 
#' GE_test_moment_calcs(beta_list=as.list(runif(n=6, min=0, max=1)), 
#' rho_list=as.list(rep(0.3,6)), prob_G=0.3, cov_Z=1, cov_W=1)

GE_test_moment_calcs <- function(beta_list, rho_list, prob_G, cov_Z=NULL, cov_W=NULL, num_sub=2000000, test_threshold=0.003, corr_G=NULL)
{
    num_G <- length(beta_list[[2]])
	  num_Z <- length(beta_list[[5]])
	  num_W <- length(beta_list[[6]])
    normal_squaredmis_results <- GE_bias_normal_squaredmis(beta_list=beta_list, 
                                                               rho_list=rho_list, 
                                                               prob_G=prob_G, 
                                                               cov_Z=cov_Z, 
                                                               cov_W=cov_W, 
                                                               corr_G=corr_G)
		# What we think the moments are
    mu_list <- normal_squaredmis_results$mu_list
    cov_list <- normal_squaredmis_results$cov_list
    cov_mat_list <- normal_squaredmis_results$cov_mat_list
    HOM_list <- normal_squaredmis_results$HOM_list
    
    # Fill in our covariances.
    rho_GE <- rho_list[[1]]; rho_GZ <- rho_list[[2]]; rho_EZ <- rho_list[[3]]
    rho_GW <- rho_list[[4]]; rho_EW <- rho_list[[5]]; rho_ZW <- rho_list[[6]]
    rho_ZZ <- cov_Z; rho_WW <- cov_W
    
    # The means
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
	  mu_GG <- as.matrix(cov_mat_list[[1]])
	  mu_GZ <- cov_mat_list[[2]]; mu_ZG <- t(mu_GZ)
	  mu_GM <- cov_mat_list[[3]]
	  mu_GW <- cov_mat_list[[4]]; mu_WG <- t(mu_GW)
	  mu_ZZ <- as.matrix(cov_mat_list[[5]])
	  mu_ZW <- cov_mat_list[[6]]; mu_WZ <- t(mu_ZW)
	  mu_ZM <- cov_mat_list[[7]]
	  mu_WW <- as.matrix(cov_mat_list[[8]])
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

	  ############################################################
	  # Generate test data
	  sig_mat <- normal_squaredmis_results$sig_mat
	  test_data <- mvtnorm::rmvnorm(n=num_sub, sigma=sig_mat)
	  
	  # For thresholding
	  w_vec <- qnorm(1-prob_G)
	  
	  # Get E, test it
	  E <- test_data[, 2*num_G+1]
	  
	  temp_test <- mu_EE - mean(E*E)
	  if (abs(temp_test) > test_threshold) {warning('Problem with mu_EE')}
	  
	  temp_test <- mu_Ef - mean(E*E*E)
	  if (abs(temp_test) > test_threshold) {warning('Problem with mu_Ef')}
	  
	  # Get G, test G and GE
	  G_mat <- matrix(data=NA, nrow=num_sub, ncol=num_G)
	  for (i in 1:num_G) {
	    G_mat[, i] <- as.numeric(test_data[, i*2-1] > w_vec[i]) + 
	      as.numeric(test_data[, i*2] > w_vec[i]) - 2*prob_G[i]
	    
	    # Test GE
	    temp_test <- rho_GE[i] - cov(G_mat[, i], E)
	    if (abs(temp_test) > test_threshold) {warning('Problem with rho_GE')}
	    
	    temp_test <- mu_GE[i] - mean(G_mat[, i]*E)
	    if (abs(temp_test) > test_threshold) {warning('Problem with mu_GE')}
	    
	    temp_test <- mu_Gf[i] - mean(G_mat[, i]*E*E)
	    if (abs(temp_test) > test_threshold) {warning('Problem with mu_Gf')}
	    
	    temp_test <- mu_Gh[i] - mean(G_mat[, i]*E*E)
	    if (abs(temp_test) > test_threshold) {warning('Problem with mu_Gh')}
	    
	    temp_test <- mu_GEE[i] - mean(G_mat[, i]*E*E)
	    if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEE')}
	    
	    temp_test <- mu_GEf[i] - mean(G_mat[, i]*E*E*E)
	    if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEf')}
	    
	    temp_test <- mu_GEh[i] - mean(G_mat[, i]*E*E*E)
	    if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEh')}
	  }
	  
	  # Higher order GE moments
	  for (i in 1:num_G) {
	    for (j in 1:num_G) {
	      temp_G1 <- G_mat[, i]
	      temp_G2 <- G_mat[, j]
	      
	      temp_test <- mu_GG[i, j] - mean(temp_G1*temp_G2)
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_GG')}
	      
	      temp_test <- mu_GEG[i, j] - mean(temp_G1*temp_G2*E)
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEG')}
	      
	      temp_test <- mu_GEEG[i, j] - mean(temp_G1*temp_G2*E*E)
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEEG')}
	      
	      temp_test <- mu_GhG[i, j] - mean(temp_G1*temp_G2*E*E)
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_GhG')}
	      
	      temp_test <- mu_GEfG[i, j] - mean(temp_G1*temp_G2*E*E*E)
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEfG')}
	      
	      temp_test <- mu_GEhG[i, j] - mean(temp_G1*temp_G2*E*E*E)
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEhG')}
	    }
	  }
	  
	  # Now do G and E with Z
	  if (num_Z > 0) {
	    Z_mat <- as.matrix(test_data[, (2+2*num_G):(1+2*num_G+num_Z)])
	    for (i in 1:num_Z) {
	      
	      temp_test <- rho_EZ[i] - cov(E,Z_mat[,i])
	      if (abs(temp_test) > test_threshold) {warning('Problem with rho_EZ')}
	      
	      temp_test <- mu_EZ[i] - mean(E*Z_mat[,i])
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_EZ')}
	      
	      temp_test <- mu_fZ[i] - mean(E*E*Z_mat[,i])
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_fZ')}
	      
	      for (j in 1:num_Z) {
	        temp_test <- mu_ZZ[i, j] - mean(Z_mat[, i]*Z_mat[, j])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZZ')}
	      }
	      
	      for (j in 1:num_G) {
	        temp_G <- G_mat[, j]
	        
	        temp_test <- rho_GZ[((j-1)*num_Z+i)] - cov(temp_G, Z_mat[,i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with rho_GZ')}
	        
	        temp_test <- mu_GZ[j, i] - mean(temp_G*Z_mat[,i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GZ')}
	        
	        temp_test <- mu_GEZ[j, i] - mean(temp_G*E*Z_mat[,i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEZ')}
	        
	        temp_test <- mu_ZEG[i, j] - mean(temp_G*E*Z_mat[,i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZEG')}
	        
	        temp_test <- mu_GhZ[j, i] - mean(temp_G*E*E*Z_mat[,i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GhZ')}
	      }
	    }
	  }
	  
	  
	  if (num_W > 0) {
	    W_mat <- as.matrix(test_data[, (2+2*num_G+num_Z):(1+2*num_G+num_Z+num_W)])
	    
	    for (i in 1:num_W) {
	      
	      temp_test <- rho_EW[i] - cov(E, W_mat[, i])
	      if (abs(temp_test) > test_threshold) {warning('Problem with rho_EW')}
	      
	      temp_test <- mu_EW[i] - mean(E*W_mat[, i])
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_EW')}
	      
	      temp_test <- mu_fW[i] - mean(E*E*W_mat[, i])
	      if (abs(temp_test) > test_threshold) {warning('Problem with mu_fW')}
	      
	      for (j in 1:num_W) {
	        temp_test <- mu_WW[i, j] - mean(W_mat[, i]*W_mat[, j])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_WW')}
	        
	        temp_test <- mu_WM[i, j] - mean(W_mat[, i]*W_mat[, j]^2)
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_WM')}
	      }
	      
	      for (j in 1:num_G) {
	        temp_G <- G_mat[, j]
	        
	        temp_test <- rho_GW[((j-1)*num_W+i)] - cov(temp_G, W_mat[, i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with rho_GW')}
	        
	        temp_test <- mu_GW[j, i] - mean(temp_G*W_mat[, i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GW')}
	        
	        temp_test <- mu_GM[j, i] - mean(temp_G*W_mat[, i]^2)
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GM')}
	        
	        temp_test <- mu_GEW[j, i] - mean(temp_G*E*W_mat[, i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEW')}
	        
	        temp_test <- mu_WEG[i, j] - mean(temp_G*E*W_mat[, i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_WEG')}
	        
	        temp_test <- mu_GEM[j, i] - mean(temp_G*E*W_mat[, i]^2)
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEM')}
	        
	        temp_test <- mu_GhW[j, i] - mean(temp_G*E*E*W_mat[, i])
	        if (abs(temp_test) > test_threshold) {warning('Problem with mu_GhW')}
	      }
	      
	      if (num_Z > 0) {
	        for (j in 1:num_Z) {
	          temp_test <- rho_ZW[((j-1)*num_W+i)] - mean(Z_mat[, j]*W_mat[, i])
	          if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZW')}
	          
	          temp_test <- mu_ZW[j, i] - mean(Z_mat[, j]*W_mat[, i])
	          if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZW')}
	          
	          temp_test <- mu_ZM[j, i] - mean(Z_mat[, j]*W_mat[, i]^2)
	          if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZM')}
	        }
	      }
	    }
	  }

  	cat('Done')
}




