#' GE_BICS.R
#'
#' A function to perform inference on the GxE interaction regression coefficient.
#' Shows better small sample performance than comparable methods.
#' 
#' @param outcome The outcome vector
#' @param design_mat The design matrix of covariates
#' @param num_boots The number of bootstrap resamples to perform - we suggest 1000
#' @param desired_coef The column in the design matrix holding the interaction covariate
#' @param outcome_type Either 'D' for dichotomous outcome or 'C' for continuous outcome
#' 
#' @return The p-value for the interaction effect
#'
#' @import stats
#'
#' @export
#' @examples 
#' E <- rnorm(n=500)
#' G <- rbinom(n=500, size=2, prob=0.3)
#' design_mat <- cbind(1, G, E, G*E)
#' outcome <- rnorm(500)
#' GE_BICS(outcome=outcome, design_mat=design_mat, desired_coef=4, outcome_type='C')

GE_BICS <- function(outcome, design_mat, num_boots=1000, desired_coef, outcome_type)
{
	colnames(design_mat) <- 1:ncol(design_mat)
	n <- length(outcome)
	
	# Fit the initial mod
	if (outcome_type == 'C') {
		init_mod <- geepack::geeglm(outcome~design_mat - 1, id=1:n, std.err='san.se')
	} else if (outcome_type == 'D') {
		init_mod <- geepack::geeglm(outcome~design_mat - 1, family=binomial(link='logit'), 	
				id=1:n, std.err='san.se')
	} else {
		stop('Invalid outcome type!')
	}
	beta_init <- summary(init_mod)$coefficients[desired_coef,1]
	se_init <- summary(init_mod)$coefficients[desired_coef,2]

	# Bootstrapping
	b_vec <- rep(NA, num_boots)
	z_vec <- rep(NA, num_boots)
	for (i in 1:num_boots)
	{
		# Get bootstrap samples
		samp_index <- sample(1:n, size=n, replace=TRUE)
		temp_Y <- outcome[samp_index]
		temp_X <- design_mat[samp_index,]
		
		# If n small, need to check and make sure we can do the fitting/d_mat nonsingular.
		# This runs fast when n small so ok to use a little more computing power.
		if (n <= 500)
		{
			svd_min <- svd(temp_X)$d[ncol(temp_X)]
			if (svd_min < 10^(-4)) {next}
		} 
		
		# Get the sandwich estimator fast
		if (outcome_type == 'C') {
			bread <- solve(crossprod(temp_X))
			b_hat <- ( bread %*% crossprod(temp_X,temp_Y) )
			fit_y <- temp_X %*% b_hat
			resid_sq <- as.numeric( (temp_Y - fit_y)^2 )
			meat <- crossprod(temp_X * resid_sq, temp_X)
			b_k <- b_hat[desired_coef]
			s_k <- sqrt( (bread %*% meat %*% bread)[desired_coef,desired_coef] )	
		} else {
			boot_mod <- tryCatch(speedglm::speedglm.wfit(y=temp_Y, X=temp_X, family=binomial()), 
				warning=function(w) w, error=function(e) e)
			if (length(class(boot_mod)) > 1) {next}
			
			fitted <- as.numeric( rje::expit( temp_X %*% boot_mod$coefficients ) )
			bread <- solve( crossprod(temp_X * fitted * (1-fitted), temp_X) )
			meat <- crossprod(temp_X*(temp_Y-fitted)^2, temp_X)
			s_k <- sqrt( (bread %*% meat %*% bread)[desired_coef,desired_coef] )
			b_k <- boot_mod$coefficients[desired_coef]
		}
		
		b_vec[i] <- b_k
		z_vec[i] <- ((b_k-beta_init)/s_k)^2
	}
	
	# Match moments
	mean_Z <- mean(z_vec, na.rm=TRUE)
	var_Z <- var(z_vec, na.rm=TRUE)
	cee <- var_Z / (2*mean_Z)
	a <- mean_Z / cee
	
	# Calculate p-value
	test_stat <- (beta_init / se_init)^2 / cee
	p_value <- 1-pgamma(test_stat, shape=a/2, scale=2)
	
	
	return( p_value )
}

