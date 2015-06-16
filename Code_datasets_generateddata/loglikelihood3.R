# function to caclulate the loglikelyhood of a set of visits.
# used by the optim function to estimate the score function
# loglikelihood2() uses count data while loglikelihood() uses binary data

loglikelihood3 <- function(Phi, U, A, lambda, nMW, nPlants){

	Prob = array(0, dim = c(nPlants, nMW))
	LL = 0
	col_sum = 0

	remove_watches = vector()

	# To calculate probability
	for (mw in 1:nMW) {
		#Calculating the denominator of probability P(v_t(j) = i)
		col_sum = t(A[,mw])%*%exp(Phi)

		#calculating probability
		Prob[,mw] = (A[,mw]*exp(Phi))/col_sum[1][1]

		if(any(is.infinite(Prob[,mw]))){
			Prob[,mw] = 0
		}

	} # mw

	# Calculating Log Likelihood
	temp = 0
	for (mw in 1:nMW) {
		k = 0
		term1_denom = 0
		term2 = 0
		for (pl in 1:nPlants) {
			use = U[pl,mw]
			# if use > 0 so that we don't have issues with log(Prob)
			# being zero
			if (use > 0){
				k = k + use
				term1_denom = term1_denom + lfactorial(use)
				term_t = term2
				term2 = term2 + use * log(Prob[pl,mw])
			}
		} # pl
		temp = temp + lfactorial(k)- term1_denom + term2
	} # mw

	f = sqrt(t(Phi) %*% Phi)
	LL = temp - (lambda * f)

	return(-LL) 

} # function end
