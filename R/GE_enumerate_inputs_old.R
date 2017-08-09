#' GE_enumerate_inputs_old.R
#'
#' Call this function to display the necessary inputs for GE_bias.  If you see a term like
#' mu_Gf, that means it is the scalar E[G*f(E)].  If you see a term like MU_GM, that means it
#' is the vector c(E[G*M_1], E[G*M_2], ..., E[G*M_q]) where M is of dimension q.  For terms in 
#' cov_mat_list like E[ZZ], these should be matrices where the (i,j) element is E[Z_i*Z_j].
#'
#' @return Nothing
#'
#' @export
#' @examples 
#' GE_enumerate_inputs_old()

GE_enumerate_inputs_old <- function()
{
	beta_help <- paste("beta_list should contain the following 6 inputs (in order):",
					"beta_0",
					"beta_G",
					"beta_E",
					"beta_I",
					"BETA_Z",
					"BETA_M", sep='\n')

	mu_help <- paste("mu_list should contain the following 5 inputs (in order):",
					"mu_f",
					"mu_h",
					"MU_Z",
					"MU_M",
					"MU_W", sep='\n')
					
	cov_help <- paste("cov_list should contain the following 14 inputs (in order):",
					"mu_GG",
					"mu_GE",
					"mu_Gf",
					"mu_Gh",
					"MU_GZ",
					"MU_GM",
					"MU_GW",
					"mu_EE",
					"mu_Ef",
					"MU_EZ",
					"MU_EM",
					"MU_EW",
					"MU_fZ",
					"MU_fW", sep='\n')
					
	mat_help <- paste("cov_mat_list should be a list of the matrices:",
					"MU_ZZ = Z%*%t(Z)",
					"MU_WW",
					"MU_ZW",
					"MU_WZ",
					"MU_ZM",
					"MU_WM", sep='\n')
					
	HOM_help <- paste("HOM_list should be a list of the following higher order moments:",
					"mu_GGE",
					"mu_GGh",
					"mu_GEE",
					"mu_GEf",
					"mu_GEh",
					"MU_GEZ",
					"MU_GEM",
					"MU_GEW",
					"MU_GhW",
					"MU_GhZ",
					"mu_GGEE",
					"mu_GGEf",
					"mu_GGEh", sep='\n')
					
	cat(beta_help, '\n\n')
	cat(mu_help, '\n\n')
	cat(cov_help, '\n\n')
	cat(mat_help, '\n\n')
	cat(HOM_help, '\n\n')
}

