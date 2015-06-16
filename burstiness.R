# function to caclulate the loglikelihood of a set of visits.
# based on the burstiness model

# Phi_d will have Phi scores and d where d is the last element
# length(Phi_d) = nPlants + 1
burstiness <- function(Phi_d, U, A, nMW, nPlants){

	Phi = Phi_d[-length(Phi_d)]
	d = Phi_d[length(Phi_d)]

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

		if(col_sum[1][1] == 0){
			Prob[,mw] = 0
			remove_watches = append(mw, remove_watches)
		}
	} # mw

	for (mw in 1:nMW) {
		V = vector()
		list = vector(mode = "list")
		lh = vector() #likelihoods

		# change U to a list of visits v
		for(i in 1:length(U[,mw])){
			V = append(rep(i, U[i,mw]), V)
			if(U[i,mw] > 0){
				list[list(length)+1] = i
			}
		} #i

		total_use = sum(U[,mw])

		#likelihood calculations
		for(e in list){
			lh = append(log(Prob[e,mw]),lh)
		} #e

		for(int in 2:total_use){
			list_temp = vector()
			lh_temp = vector()
			for(e in 1:length(list)){
				utemp = U[,mw]
				for(pl in 1:nPlants){
					utemp[pl] = utemp[pl] - length(list[e][which(list[e]==pl)])
					if(utemp[pl] > 0){
						if(list[e][1] == pl){
							trans = log((1-d) + d*Prob[pl, mw])
						}
						else{
							trans = log(d*Prob[pl,mw])
						}
						temp = append(pl, list[e])
						list_temp = append(temp, list_temp)
						lh_temp = append(lh[e]+trans, lh_temp)
					}
				} #pl
			} #e
			list = list_temp
			lh = lh_temp
		} #int
	} #mw
	

	return(-LL) 

} # function end
