#' GE_translate_inputs.R
#'
#' Mostly for internal use, function called by GE_bias_normal() and GE_scoreeq_sim()
#' to translate the rho_list inputs and return a total covariance matrix for simulation/
#' checking validity of covariance structure.  If invalid covariance structure, will stop
#' and return an error message.
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
#' @param prob_G Probability that each allele is equal to 1.  Since each SNP has
#' two alleles, the expectation of G is 2*prob_G.
#' @param cov_Z Should be a matrix equal to cov(Z) or NULL if no Z.  
#' @param cov_W Should be a matrix equal to cov(W) or NULL if no W.
#' @param corr_G Should be a matrix giving the *pairwise correlations* between each SNP
#' in the set, or NULL. Must be specified if G is a vector.  For example, the [2,3] element
#' of the matrix would be the pairwise correlation between SNP2 and SNP3. Diagonal
#' should be 1.

#'
#' @return A list with the elements:
#' \item{sig_mat_total}{The sigma parameter for rmvnorm call to generate our data.}
#'
#' @keywords internal
#' @export
#' @examples 
#' GE_translate_inputs( beta_list=as.list(runif(n=6, min=0, max=1)),
#' rho_list=as.list(rep(0.3,6)), prob_G=0.3, cov_Z=1, cov_W=1)

GE_translate_inputs <- function(beta_list, rho_list, prob_G, cov_Z=NULL, cov_W=NULL, corr_G=NULL)
{
	  # First, make sure we got good inputs
  	if (length(beta_list) != 6 | length(rho_list) != 6 | !inherits(beta_list, "list") | !inherits(rho_list, "list"))
  	{
  	  stop('Input vectors not the right size!')
  	}

	  # How long are these vectors?  Remember that W is the same length as M by assumption.
    # The number is 0 if an element of beta_list is NULL
    num_G <- length(beta_list[[2]])
    num_I <- length(beta_list[[4]])
    num_Z <- length(beta_list[[5]]); if(num_Z == 1 & beta_list[[5]][1] == 0) {num_Z <- 0}
    num_W <- length(beta_list[[6]]); if(num_W == 1 & beta_list[[6]][1] == 0) {num_W <- 0}
    
    # Make sure we have the same number of effect sizes for G and G*E in true model.
    # Also same as the length of the MAF for G.
    if (num_G != num_I | num_G != length(prob_G)) {
      stop('Discordance between effect sizes for G and probabilities of G')
    }

    # Fill in our covariances.
    rho_GE <- rho_list[[1]]; rho_GZ <- rho_list[[2]]; rho_EZ <- rho_list[[3]]
    rho_GW <- rho_list[[4]]; rho_EW <- rho_list[[5]]; rho_ZW <- rho_list[[6]]

  	# Make sure we have compatible lengths for rho_list
    if (length(rho_GE) != num_G) {
      stop('Incompatible number of elements in beta/rho_list')
    }
    if (num_Z > 0) {
      if (length(rho_GZ) != num_Z*num_G | length(rho_EZ) != num_Z)
      {
        stop('Incompatible number of elements in beta/rho_list')
      }
      if (num_W > 0) {
        if (length(rho_ZW) != num_Z*num_W) {
          stop('Incompatible number of elements in beta/rho_list')
        }
      }
    }
    if (num_W > 0) {
      if (length(rho_GW) != num_G*num_W | length(rho_EW) != num_W) {
        stop('Incompatible number of elements in beta/rho_list')
      }
    }
    
    if (num_G > 1) {
      if (ncol(corr_G) != num_G | nrow(corr_G) != num_G) {
        stop('Incompatible number of elements in beta/rho_list')
      }
    }
    
   	################################################################
    # Build our covariance matrix in steps.
    
    # (1) First step is to build the 2d*2d upper left corner for the G
    # Use bindata to do this
    if (num_G > 1) {
      #G_bin_struct <- matrix(data=0, nrow=num_G, ncol=num_G)
      #G_bin_struct[upper.tri(G_bin_struct)] <- rho_GG
      #G_bin_struct <- G_bin_struct + t(G_bin_struct)
      G_bin_struct <- corr_G
      
      cprob <- tryCatch(bindata::bincorr2commonprob(margprob = prob_G, bincorr=G_bin_struct),
                        error=function(e) e)
      if ('error' %in% class(cprob)) {
        stop ('You specified an invalid corr_G structure')
      }
      sigma_struct <- tryCatch(bindata::commonprob2sigma(commonprob=cprob), 
                               error=function(e) e)
      if ('error' %in% class(sigma_struct)) {
        stop ('You specified an invalid corr_G structure')
      }
      
      # Now that we have the structure, we need to account for two thresholded
      # normals making up every G, so insert an extra d terms.
      G_segment <- matrix(data=0, nrow=2*num_G, ncol=2*num_G)
      for (i in 1:num_G) {
        odd_seq <- seq(from=1, to=(2*num_G-1), by=2)
        even_seq <- odd_seq + 1
        G_segment[i*2-1, odd_seq] <- sigma_struct[i, ]
        G_segment[i*2, even_seq] <- sigma_struct[i, ]
      }
    } else {
      G_segment <- diag(x=1, nrow=2, ncol=2)
    }
    
    # Now generate the 2d+1 row/column, which describes GE
    w_vec <- qnorm(1-prob_G)
    r_GE <- rho_GE / (2*dnorm(w_vec))
    GE_segment <- rep(r_GE, each=2)
    
    # Now make the GZ part
    if (num_Z != 0) {
      GZ_segment <- matrix(data=NA, nrow=2*num_G, ncol=num_Z) 
      for (i in 1:num_G) {
        r_GZ <-  rho_GZ[((i-1)*num_Z+1):(i*num_Z)] / (2*dnorm(w_vec[i]))
        GZ_segment[i*2-1, ] <- r_GZ
        GZ_segment[i*2, ] <- r_GZ
      }
    }
    
    # Build the GW part
    if (num_W != 0) {
      GW_segment <- matrix(data=NA, nrow=2*num_G, ncol=num_W) 
      for (i in 1:num_G) {
        r_GW <- rho_GW[((i-1)*num_W+1):(i*num_W)] / (2*dnorm(w_vec[i]))
        GW_segment[i*2-1, ] <- r_GW
        GW_segment[i*2, ] <- r_GW
      }
    }
   
    # Now put it all together
    # Add the G segment first
    MVN_sig_tot <- matrix(data=NA, nrow=(2*num_G+num_Z+num_W+1), ncol=(2*num_G+num_Z+num_W+1))
    MVN_sig_tot[1:(2*num_G), 1:(2*num_G)] <- G_segment
    
    # Then the GE segment
    MVN_sig_tot[(2*num_G+1), 1:(2*num_G)] <- GE_segment
    MVN_sig_tot[1:(2*num_G), (2*num_G+1)] <- GE_segment
    MVN_sig_tot[(2*num_G+1), (2*num_G+1)] <- 1
    
    # Most complicated is if both Z and W exist
    if (num_Z > 0 & num_W > 0) {
      # Then the GZ segment
      MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), 1:(2*num_G)] <- t(GZ_segment)
      MVN_sig_tot[1:(2*num_G), (2*num_G+2):(2*num_G+num_Z+1)] <- GZ_segment
      
      # Then the GW segment
      MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), 1:(2*num_G)] <- t(GW_segment)
      MVN_sig_tot[1:(2*num_G), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- GW_segment
      
      # Then the EZ segment
      MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), (2*num_G+1)] <- rho_EZ
      MVN_sig_tot[(2*num_G+1), (2*num_G+2):(2*num_G+num_Z+1)] <- rho_EZ
      
      # Then the ZZ segment
      if (num_Z == 1) {
        MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), (2*num_G+2):(2*num_G+num_Z+1)] <- 1
      } else {
        MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), (2*num_G+2):(2*num_G+num_Z+1)] <- cov_Z
      }
      
      # Then the EW segment
      MVN_sig_tot[(2*num_G+1), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- rho_EW
      MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), (2*num_G+1)] <- rho_EW
      
      # Build the ZW part
      ZW_segment <- matrix(data=NA, nrow=num_Z, ncol=num_W) 
      for (i in 1:num_Z) {
        ZW_segment[i, ] <- rho_ZW[((i-1)*num_W+1):(i*num_W)]
      }
      
      # Then the ZW segment
      MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- ZW_segment
      MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), (2*num_G+2):(2*num_G+num_Z+1)] <- t(ZW_segment)
      
      # Then the WW segment
      if (num_W == 1) {
        MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- 1
      } else {
        MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- cov_W
      }
    }
    
    # If only Z and no W
    if (num_Z > 0 & num_W == 0) {
      # Then the GZ segment
      MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), 1:(2*num_G)] <- t(GZ_segment)
      MVN_sig_tot[1:(2*num_G), (2*num_G+2):(2*num_G+num_Z+1)] <- GZ_segment
      
      # Then the EZ segment
      MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), (2*num_G+1)] <- rho_EZ
      MVN_sig_tot[(2*num_G+1), (2*num_G+2):(2*num_G+num_Z+1)] <- rho_EZ
      
      # Then the ZZ segment
      if (num_Z == 1) {
        MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), (2*num_G+2):(2*num_G+num_Z+1)] <- 1
      } else {
        MVN_sig_tot[(2*num_G+2):(2*num_G+num_Z+1), (2*num_G+2):(2*num_G+num_Z+1)] <- cov_Z
      }
      
    }
    
    # If only W and no Z
    if (num_Z == 0 & num_W > 0) {
      # Then the GW segment
      MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), 1:(2*num_G)] <- t(GW_segment)
      MVN_sig_tot[1:(2*num_G), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- GW_segment
    
      # Then the EW segment
      MVN_sig_tot[(2*num_G+1), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- rho_EW
      MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), (2*num_G+1)] <- rho_EW
      
      # Then the WW segment
      if (num_W == 1) {
        MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- 1
      } else {
        MVN_sig_tot[(2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1), (2*num_G+num_Z+2):(2*num_G+num_Z+num_W+1)] <- cov_W
      }
    }
    
  
    # Final sanity check that the building went correctly.
    if (!isSymmetric(MVN_sig_tot)) {stop("Problem building covariance matrix!")}
    
    # Now make sure we can actually generate data with this structure
    test_data <- tryCatch(mvtnorm::rmvnorm(n=1, sigma=MVN_sig_tot), 
    				warning=function(w) w, error=function(e) e)
    if (class(test_data)[1] != 'matrix') {stop('You specified an impossible covariance matrix!')}
    
    return(list(sig_mat_total=MVN_sig_tot))
}




