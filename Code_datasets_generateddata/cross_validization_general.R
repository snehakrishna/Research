# cross_validization function

source("loglikelihood3.R")
library(bipartite)

cross_validization <- function(U, A, nMW, num_division, lambda_options, nPlants){

	Zero = array(0, dim = c(nPlants, 1))

	#random_MW used to test first time
#	random_MW = c(5, 6, 1, 5, 7, 3, 3, 3, 2, 4, 6, 5, 1, 6, 5, 3, 1, 1, 6, 7, 2,
#		2, 1, 7, 7, 6, 1, 5, 6, 4, 2, 5, 1, 7, 3, 4, 2, 4, 7, 3, 1, 2, 6, 2, 1, 
#		7, 6, 4, 2, 4, 6, 3, 4, 3, 5, 4, 4, 2, 5, 5, 7, 3, 1, 2, 3, 7, 3, 5, 4, 
#		6, 6, 4, 4, 3, 5, 7, 2, 1, 1, 2, 7, 6, 5)
	random_MW = array(0, dim = c(nMW))
	nremove = ceiling(nMW/num_division)
	count = 1
	for (i in 1:num_division){
		for (j in 1:nremove){
			if (count <= nMW){
				random_MW[count] = i
				count = count + 1
			}
		} #j
	} #i
	random_MW = sample(random_MW, replace = FALSE)
	
	num_lambda = length(lambda_options)
	votes = array(0, dim = c(num_lambda, 1))
	
	for (div in 1:num_division){
		print(div)
		A_new = A
		U_new = U
	
		watches = vector()
		for (i in 1:nMW){
			if (random_MW[i] == div){
				watches = append(watches, i)
			}
		} #i

		A_new = A_new[,-watches]
		U_new = U_new[,-watches]
		meadow_watch = meadow_watch[-watches]
		#print(meadow_watch[32])
	
		A_left = A[,watches]
		U_left = U[,watches]

		min_LL = Inf
		vote_lambda = 0
		for (lambda in 1:num_lambda){
			print(c("\t", lambda))
			Phi = optim(Zero, loglikelihood3, gr = NULL, U_new, A_new, lambda_options[lambda], nMW-length(watches), nPlants, method = "BFGS")
			LL = loglikelihood3(Phi$par, U_left, A_left, lambda_options[lambda], length(watches), nPlants)
			if(LL < min_LL){
				min_LL = LL
				vote_lambda = lambda
			}
		} #lambda
		votes[vote_lambda] = votes[vote_lambda] + 1
	} #div

	max_votes = 0
	max_lambda = 0
	for (n in 1:length(votes)){
		if (votes[n] >= max_votes){
			max_votes = votes[n]
			max_lambda = n
		}
	}

	chosen_lambda = lambda_options[max_lambda]

	return(chosen_lambda)
}

# .25 for normal and random_MW = 
# [1] 5 6 1 5 7 3 3 3 2 4 6 5 1 6 5 3 1 1 6 7 2 2 1 7 7 6 1 5 6 4 2 5 1 7 3 4 2 4
#[39] 7 3 1 2 6 2 1 7 6 4 2 4 6 3 4 3 5 4 4 2 5 5 7 3 1 2 3 7 3 5 4 6 6 4 4 3 5 7
#[77] 2 1 1 2 7 6 5
