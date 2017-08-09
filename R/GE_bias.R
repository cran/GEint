#' GE_bias.R
#'
#' A function to calculate the bias in testing for GxE interaction.
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
#' GE_bias(beta_list=solutions$beta_list, solutions$cov_list, solutions$cov_mat_list, 
#'						solutions$mu_list, solutions$HOM_list)

GE_bias <- function(beta_list, cov_list, cov_mat_list, mu_list, HOM_list)
{
  # Record some initial quantities
  beta_0 <- beta_list[[1]]; beta_G <- beta_list[[2]]; beta_E <- beta_list[[3]]
  beta_I <- beta_list[[4]]; beta_Z <- beta_list[[5]]; beta_M <- beta_list[[6]]
  
  # Numbers
  num_G <- length(beta_G)
  num_Z <- length(beta_Z); if(num_Z == 1 & beta_Z[1] == 0) {num_Z <- 0}
  num_W <- length(beta_M); if(num_W == 1 & beta_M[1] == 0) {num_W <- 0}
  
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

  ########################################################################
  # Start solving
  # Some shortcut quantities
  U <- solve(mu_GG)
  A <- (mu_ZG %*% U %*% mu_GE - mu_EZ) / as.numeric(mu_EE - t(mu_GE) %*% t(U) %*% mu_GE)
  B <- (mu_WG %*% U %*% mu_GE - mu_EW) / as.numeric(mu_EE - t(mu_GE) %*% t(U) %*% mu_GE)
 
  # O will be set to 0 if no Z
  if (num_Z == 0) {
    O <- as.matrix(0)
    solve_O <- as.matrix(0)
  } else {
    O <- mu_Z %*% t(mu_Z) + mu_ZG %*% U %*% mu_GZ - mu_ZZ - 
            A %*% (t(mu_EZ) - t(mu_GE) %*% U %*% mu_GZ)
    solve_O <- solve(O)
  }
  
  C <- (-mu_W %*% t(mu_Z) - mu_WG %*% U %*% mu_GZ + mu_WZ + 
          B %*% (t(mu_EZ) - t(mu_GE) %*% U %*% mu_GZ)) %*% solve_O
 
  # Q will be 0 if no W
  if ( num_W == 0 ) {
    Q <- as.matrix(0)
    solve_Q <- as.matrix(0)
    # Need it once for later
    mu_WW <- as.matrix(0)
  } else {
    Q <- mu_W %*% t(mu_W) + mu_WG %*% U %*% mu_GW - mu_WW + 
      B %*% (t(mu_GE) %*% U %*% mu_GW - t(mu_EW)) + 
      C %*% (A %*% (t(mu_GE) %*% U %*% mu_GW - t(mu_EW)) + 
      mu_Z %*% t(mu_W) + mu_ZG %*% U %*% mu_GW - mu_ZW)
    solve_Q <- solve(Q)
  }
  
  D <- (mu_GEG %*% U %*% mu_GE - mu_GEE) / as.numeric(mu_EE - t(mu_GE) %*% t(U) %*% mu_GE)
  E <- (mu_GEZ - mu_GE %*% t(mu_Z) - mu_GEG %*% U %*% mu_GZ + 
            D %*% (t(mu_EZ) - t(mu_GE) %*% U %*% mu_GZ)) %*% solve_O
  EFF <- (mu_GE %*% t(mu_W) + mu_GEG %*% U %*% mu_GW - mu_GEW + 
            D %*% (t(mu_GE) %*% U %*% mu_GW - t(mu_EW)) + 
            E %*% (A %*% (t(mu_GE) %*% U %*% mu_GW - t(mu_EW)) +
            mu_Z %*% t(mu_W) + mu_ZG %*% U %*% mu_GW - mu_ZW)) %*% solve_Q
  
  # Solve for \alpha_I
  S <- EFF %*% (mu_W %*% t(mu_GE) + mu_WG %*% U %*% mu_GEG - mu_WEG + 
            B %*% (t(mu_GE) %*% U %*% mu_GEG - t(mu_GEE))) + 
        EFF %*% C %*% (mu_Z %*% t(mu_GE) + mu_ZG %*% U %*% mu_GEG - mu_ZEG + 
            A %*% (t(mu_GE) %*% U %*% mu_GEG - t(mu_GEE))) -
        (mu_GE %*% t(mu_GE) + mu_GEG %*% U %*% mu_GEG - mu_GEEG + 
            D %*% (t(mu_GE) %*% U %*% mu_GEG - t(mu_GEE))) - 
        E %*% (mu_Z %*% t(mu_GE) + mu_ZG %*% U %*% mu_GEG - mu_ZEG + 
            A %*% (t(mu_GE) %*% U %*% mu_GEG - t(mu_GEE)))
  TEE <- beta_E * (-mu_f * mu_GE - mu_GEG %*% U %*% mu_Gf + mu_GEf + 
            D %*% (mu_Ef - t(mu_Gf) %*% t(U) %*% mu_GE)) + 
          beta_E * E %*% (-mu_f * mu_Z - mu_ZG %*% U %*% mu_Gf + mu_fZ + 
            A %*% (mu_Ef - t(mu_Gf) %*% t(U) %*% mu_GE)) + 
        (-mu_GE %*% t(mu_Gh) - mu_GEG %*% U %*% mu_GhG + mu_GEhG + 
            D %*% (t(mu_GEh) - t(mu_GE) %*% U %*% mu_GhG)) %*% beta_I + 
          E %*% (mu_ZhG - mu_Z %*% t(mu_Gh) - mu_ZG %*% U %*% mu_GhG + 
            A %*% (t(mu_GEh) - t(mu_GE) %*% U %*% mu_GhG)) %*% beta_I + 
        (mu_GEM - mu_GE %*% t(mu_M) - mu_GEG %*% U %*% mu_GM + 
            D %*% (t(mu_EM) - t(mu_GE) %*% U %*% mu_GM)) %*% beta_M + 
        E %*% (A %*% (t(mu_EM) - t(mu_GE) %*% U %*% mu_GM) - 
            mu_Z %*% t(mu_M) - mu_ZG %*% U %*% mu_GM + mu_ZM) %*% beta_M - 
        beta_E * EFF %*% (-mu_f * mu_W - mu_WG %*% U %*% mu_Gf + mu_fW + 
            B %*% (mu_Ef - t(mu_Gf) %*% t(U) %*% mu_GE)) - 
        beta_E * EFF %*% C %*% (-mu_f * mu_Z - mu_ZG %*% U %*% mu_Gf + mu_fZ + 
            A %*% (mu_Ef - t(mu_Gf) %*% t(U) %*% mu_GE)) - 
        EFF %*% (-mu_W %*% t(mu_Gh) - mu_WG %*% U %*% mu_GhG + mu_WhG + 
            B %*% (t(mu_GEh) - t(mu_GE) %*% U %*% mu_GhG)) %*% beta_I - 
        EFF %*% C %*% (mu_ZhG - mu_Z %*% t(mu_Gh) - mu_ZG %*% U %*% mu_GhG +
            A %*% (t(mu_GEh) - t(mu_GE) %*% U %*% mu_GhG)) %*% beta_I - 
        EFF %*% (-mu_W %*% t(mu_M) - mu_WG %*% U %*% mu_GM + mu_WM + 
            B %*% (t(mu_EM) - t(mu_GE) %*% U %*% mu_GM)) %*% beta_M -
        EFF %*% C %*% (A %*% (t(mu_EM) - t(mu_GE) %*% U %*% mu_GM) - 
            mu_Z %*% t(mu_M) - mu_ZG %*% U %*% mu_GM + mu_ZM) %*% beta_M
  

  alpha_I <- solve(S) %*% TEE
  
  R <- beta_E * (-mu_f * mu_W - mu_WG %*% U %*% mu_Gf + mu_fW + 
          B %*% (mu_Ef - t(mu_Gf) %*% t(U) %*% mu_GE)) + 
        beta_E * C %*% (-mu_f * mu_Z - mu_ZG %*% U %*% mu_Gf + mu_fZ + 
          A %*% ((mu_Ef) - t(mu_Gf) %*% t(U) %*% mu_GE)) + 
        (-mu_W %*% t(mu_Gh) - mu_WG %*% U %*% mu_GhG + mu_WhG + 
          B %*% (t(mu_GEh) - t(mu_GE) %*% U %*% mu_GhG)) %*% beta_I + 
        C %*% (mu_ZhG - mu_Z %*% t(mu_Gh) - mu_ZG %*% U %*% mu_GhG + 
          A %*% (t(mu_GEh) - t(mu_GE) %*% U %*% mu_GhG)) %*% beta_I + 
        (-mu_W %*% t(mu_M) - mu_WG %*% U %*% mu_GM + mu_WM + 
          B %*% (t(mu_EM) - t(mu_GE) %*% U %*% mu_GM)) %*% beta_M + 
        C %*% (A %*% (t(mu_EM) - t(mu_GE) %*% U %*% mu_GM) - 
          mu_Z %*% t(mu_M) - mu_ZG %*% U %*% mu_GM + mu_ZM) %*% beta_M + 
        (mu_W %*% t(mu_GE) + mu_WG %*% U %*% mu_GEG - mu_WEG + 
          B %*% (t(mu_GE) %*% U %*% mu_GEG - t(mu_GEE))) %*% alpha_I + 
        C %*% (mu_Z %*% t(mu_GE) + mu_ZG %*% U %*% mu_GEG - mu_ZEG + 
          A %*% (t(mu_GE) %*% U %*% mu_GEG - t(mu_GEE))) %*% alpha_I
  Q <- mu_W %*% t(mu_W) + mu_WG %*% U %*% mu_GW - mu_WW + 
          B %*% (t(mu_GE) %*% U %*% mu_GW - t(mu_EW)) + 
        C %*% (A %*% (t(mu_GE) %*% U %*% mu_GW - t(mu_EW)) + 
          mu_Z %*% t(mu_W) + mu_ZG %*% U %*% mu_GW - mu_ZW)
  
  alpha_W <- - solve_Q %*% R
  
  P <- beta_E * (-mu_f * mu_Z - mu_ZG %*% U %*% mu_Gf + mu_fZ + 
          A %*% (mu_Ef - t(mu_Gf) %*% t(U) %*% mu_GE)) + 
        (mu_ZhG - mu_Z %*% t(mu_Gh) - mu_ZG %*% U %*% mu_GhG + 
          A %*% (t(mu_GEh) - t(mu_GE) %*% U %*% mu_GhG)) %*% beta_I + 
        (mu_Z %*% t(mu_GE) + mu_ZG %*% U %*% mu_GEG - mu_ZEG + 
          A %*% (t(mu_GE) %*% U %*% mu_GEG - t(mu_GEE))) %*% alpha_I + 
        (A %*% (t(mu_EM) - t(mu_GE) %*% U %*% mu_GM) - 
           mu_Z %*% t(mu_M) - mu_ZG %*% U %*% mu_GM + mu_ZM) %*% beta_M + 
        (A %*% (t(mu_GE) %*% U %*% mu_GW - t(mu_EW)) + 
            mu_Z %*% t(mu_W) + mu_ZG %*% U %*% mu_GW - mu_ZW) %*% alpha_W
  
  Bz_Az <- solve_O %*% P
  Az_Bz <- -Bz_Az
  alpha_Z <- beta_Z - solve_O %*% P
  
  alpha_E <- as.numeric(( beta_E * (mu_Ef - t(mu_Gf) %*% t(U) %*% mu_GE) + 
              t(beta_I) %*% (mu_GEh - t(mu_GhG) %*% t(U) %*% mu_GE) + 
              t(alpha_I) %*% (t(mu_GEG) %*% t(U) %*% mu_GE - mu_GEE) + 
              t(Bz_Az) %*% (mu_EZ - t(mu_GZ) %*% t(U) %*% mu_GE) + 
              t(beta_M) %*% (mu_EM - t(mu_GM) %*% t(U) %*% mu_GE) + 
              t(alpha_W) %*% (t(mu_GW) %*% t(U) %*% mu_GE - mu_EW) ) / 
                  (mu_EE - t(mu_GE) %*% t(U) %*% mu_GE))
  
  Bg_Ag <- U %*% (alpha_E * mu_GE - beta_E * mu_Gf + mu_GEG %*% alpha_I - 
            mu_GhG %*% beta_I + mu_GZ %*% (Az_Bz) + mu_GW %*% alpha_W - mu_GM %*% beta_M)
  alpha_G <- beta_G - Bg_Ag
  
  alpha_0 <- as.numeric(beta_0 + beta_E * mu_f + t(beta_I) %*% mu_Gh - t(alpha_I) %*% mu_GE +  
                          t(Bz_Az) %*% mu_Z + t(beta_M) %*% mu_M) - t(alpha_W) %*% mu_W
  
  # Return 
  return(list(alpha_0=alpha_0, alpha_G=alpha_G, alpha_E=alpha_E, 
              alpha_I=alpha_I, alpha_Z=alpha_Z, alpha_W=alpha_W))
}