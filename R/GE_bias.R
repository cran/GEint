#' GE_bias.R
#'
#' A function to calculate the bias in testing for GxE interaction.
#' 
#' @param beta_list A list of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' If Z or M is a vector, then beta_Z and beta_M should be vectors.
#' If Z and/or M/W do not exist in your model, then set beta_Z and/or beta_M = 0.
#' @param cov_list A list of expectations (which happen to be covariances if all covariates
#' are centered at 0) in the order specified by GE_enumerate_inputs().
#' If Z and/or M/W do not exist in your model, then treat them as constants 0. For example,
#' if Z doesn't exist and W includes 2 covariates, then set cov(EZ) = 0 and cov(ZW) = (0,0).
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
#' GE_bias(beta_list=solutions$beta_list, solutions$cov_list, solutions$cov_mat_list, 
#'						solutions$mu_list, solutions$HOM_list)

GE_bias <- function(beta_list, cov_list, cov_mat_list, mu_list, HOM_list)
{
  # Record some initial quantities
  beta_0 <- beta_list[[1]]; beta_G <- beta_list[[2]]; beta_E <- beta_list[[3]]
  beta_I <- beta_list[[4]]; BETA_Z <- beta_list[[5]]; BETA_M <- beta_list[[6]]


  # Some means
  mu_f <- mu_list[[1]]; mu_h <- mu_list[[2]]; MU_Z <- mu_list[[3]]
  MU_M <- mu_list[[4]]; MU_W <- mu_list[[5]]
 
  # Some covariances
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


  ########################
  # Matrix covariances
  # MU_ZW is not the same as MU_WZ because the dimensions of the matrix are not the same!
  MU_ZZ <- cov_mat_list[[1]]
  MU_WW <- cov_mat_list[[2]]
  MU_ZW <- cov_mat_list[[3]]
  MU_WZ <- cov_mat_list[[4]]
  MU_ZM <- cov_mat_list[[5]]	
  MU_WM <- cov_mat_list[[6]]			
  
  # Some higher order moments
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
  return(list(alpha_0, alpha_G, alpha_E, alpha_I, ALPHA_Z, ALPHA_W))
}