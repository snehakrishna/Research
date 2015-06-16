# goodness of fit tests

setwd(path.expand("~//Research//thesis_stuff"))
source("loglikelihood3.R")
source("cross_validization_general.R")
library(bipartite)

#allData = read.csv("Final_Interactions.csv")
#anthesis = read.csv("Final_Anthesis.csv")

allData_org = read.csv("Interactions_Final.csv")
anthesis_org = read.csv("Anthesis_Final.csv")

allData_org$MeadowW = paste(allData_org$MEADOW, allData_org$WATCH, sep = '')
anthesis_org$MeadowW = paste(anthesis_org$MEADOW, anthesis_org$WATCH, sep = '')

allData_org$MeadowWY = paste(allData_org$MeadowW, allData_org$year, sep = '')
anthesis_org$MeadowWY = paste(anthesis_org$MeadowW, anthesis_org$YEAR%%100, sep = '')

years = unique(allData_org$year)

for(year in years){
	allData = allData_org[which(allData_org$year == year),]
	anthesis = anthesis_org[which(anthesis_org$YEAR%%100 == year),]

	plantNames = levels(allData$PLTSP_NAME)
	plantNames = plantNames[-c(1)]
	pollNames = levels(allData$VISSP_NAME) # pollinators
	pollNames = pollNames[-c(1)]
	nPlants = length(plantNames)
	nPolls = length(pollNames)

	meadow_watch = unique(allData$MeadowW)
	meadow_watch = meadow_watch[-c(1)]
	nMW = length(meadow_watch)

	meadow_watch_year = unique(allData$MeadowWY)
	meadow_watch_year = meadow_watch[-c(1)]
	nMWY = length(meadow_watch_year)

	A = matrix(0,nrow=nPlants,ncol=nMW,dimnames=list(plantNames, meadow_watch)) #Colums represent how many meadows
	U = array(0, dim = c(nPlants, nPolls, nMW))
	Zero = array(0, dim = c(nPlants, 1))
	Prob = array(0, dim = c(nPlants, nMW))
	Prob_zero = array(0, dim = c(nPlants, nMW))

	# Filling 2D array with availability and use data
	for (mw in 1:nMW){
		mdata = allData[which(allData$MeadowW==meadow_watch[mw]),] 
		plantNamesM = levels(mdata$PLTSP_NAME)
		plantNamesM = plantNamesM[-c(1)]
		anthesis_m = anthesis[which(anthesis$MeadowW==meadow_watch[mw]),]

		for (pl in 1:nPlants) {
			# get subset where plant exists.
			m_plants = mdata[which(mdata$PLTSP_NAME==plantNames[pl]),]
			anthesisT = anthesis_m[which(anthesis_m$SPP_NAME==plantNames[pl]),]

			# Availability
			if(dim(anthesisT)[1]>0){
				fl = as.numeric(as.character(anthesisT$NO_FLW))
				st = as.numeric(as.character(anthesisT$NO_STALK))
				temp = fl*st
				if(any(is.na(temp))){
					temp = na.omit(temp)
				}
		      	A[pl, mw] = sum(temp)
			}
			if(dim(m_plants)[1]>0){ 
				for (po in 1:nPolls){
					int = m_plants[which(m_plants$VISSP_NAME==pollNames[po]),]
					ints = as.numeric(as.character(int$NO_INT))
					if(any(is.na(ints))){
						ints = na.omit(ints)
					}
					if(A[pl,mw] == 0){
						U[pl,po,mw] = 0
					}
					else{
						U[pl,po,mw] = sum(ints)
					}
				} #po
			}
		} # pl
	} # mw

	# starting with po=50 (Apis mellifera)
	# working with U[,50,]
	po = 50
	U_true = U[,po,]
	count_true = vector()
	meadow_use = colSums(U_true)

	# find regularization + phi function
	num_division = 7
	lambda_options = c(0, .25, .5, 1, 2)
	#lambda = cross_validization(U_true, A, nMW, num_division, lambda_options, nPlants)
	lambda = .25
	Phi = optim(Zero, loglikelihood3, gr = NULL, U_true, A, lambda, nMW, nPlants, method = "BFGS")

	# Expected use data
	for (mw in 1:nMW) {
		#Calculating the denominator of probability P(v_t(j) = i)
		col_sum = t(A[,mw])%*%exp(Phi$par)
		#calculating probability
		Prob[,mw] = (A[,mw]*exp(Phi$par))/col_sum[1][1]
	} # mw

	E = array(0, dim = c(nPlants, nMW))
	for (pl in 1:nPlants){
		for (mw in 1:nMW){
			E[pl, mw] = meadow_use[mw]*Prob[pl, mw]/meadow_use[mw]
		}
	}

	for (mw in 1:nMW) {
		#Calculating the denominator of probability P(v_t(j) = i)
		col_sum = t(A[,mw])%*%exp(Zero)
		#calculating probability
		Prob_zero[,mw] = (A[,mw]*exp(Zero))/col_sum[1][1]
	} # mw

	E_zero = array(0, dim = c(nPlants, nMW))
	for (pl in 1:nPlants){
		for (mw in 1:nMW){
			E_zero[pl, mw] = meadow_use[mw]*Prob_zero[pl, mw]/meadow_use[mw]
		}
	}

	flips = 10
	U_syn = array(0, dim = c(nPlants, nMW))
	for(mw in 1:nMW){
		U_syn[,mw] = rmultinom(1, flips, Prob[,mw])
	}

	U_zero = array(0, dim = c(nPlants, nMW))
	for(mw in 1:nMW){
		U_zero[,mw] = rmultinom(1, flips, Prob_zero[,mw])
	}
print("here2")

	count_syn = vector()
	count_true = vector()
	count_zero = vector()
	E_vect = vector()
	E_zero_vect = vector()
	Prob_vect = vector()
	Prob_zero_vect = vector()
	for (pl in 1:nPlants) {
		for (mw in 1:nMW){
			if (A[pl, mw] > 0 && meadow_use[mw] > 0){ # && meadow_use[mw] > 0
				count_syn = append(count_syn, c(U_syn[pl, mw]))
				count_true = append(count_true, c(U_true[pl,mw]))
				E_vect = append(E_vect, c(E[pl,mw]))
				Prob_vect = append(Prob_vect, c(Prob[pl,mw]))
				count_zero = append(count_zero, c(U_zero[pl,mw]))
				E_zero_vect = append(E_zero_vect, c(E_zero[pl,mw]))
				Prob_zero_vect = append(Prob_zero_vect, c(Prob_zero[pl,mw]))
			}
		}
	}
print("here1")

	#Calculate probability for "true" data
	#A_temp = array(1, dim = c(nPlants, 1))
	#col_sum = t(A_temp)%*%exp(Phi$par)
	#prob_true = (A_temp * exp(Phi$par))/col_sum[1][1]
	col_true = sum(count_true)
	prob_true = count_true/col_true[1][1]

	#Calculate probability for synthetic data
	#count_plant_syn = rowSums(U_syn)
	#count_total_syn = sum(count_plant)
	#prob_syn = count_plant/count_total
	col_syn = sum(count_syn)
	prob_syn = count_syn/col_syn[1][1]
	prob_e = E_vect/sum(E_vect)

print("here")

	#chi^2 with Zero Phi
	obs = count_zero
	exp = E_zero_vect
	df = length(obs) - length(which(meadow_use[] !=0)) - nPlants - 1
	x2stat_zero <- sum(((obs - exp)^2)/exp)
	goodness_zero <- 1 - pchisq(x2stat_zero,df = df)
	print("zero test")
	print(x2stat_zero)

	# want to use count_true because we want to use the field data as our
	# observed data. The prob_e is our expected probability
	# perform goodness of fit tests
	chisq.test(count_true, p = prob_e)
	obs = count_true
	exp = E_vect
	df = length(obs) - length(which(meadow_use[] !=0)) - nPlants - 1
	x2stat <- sum(((obs - exp)^2)/exp)
	result <- 1 - pchisq(x2stat,df,lower.tail=FALSE)
	print("regular test")
	print(x2stat)
	print(result)
}