#evals with cross validation

# Plots to evaluate the multinomial model with count data

setwd(path.expand("~//Research//thesis_stuff"))
source("likelihood_traits.R")
source("cross_validization_traits.R")
library(bipartite)

allData_ori = read.csv("Interactions_Final.csv")
anthesis_ori = read.csv("Anthesis_Final.csv")
traits = read.csv("PlantCovars.csv")
	traits$generalShapeExclusion = as.factor(traits$generalShapeExclusion)
	traits$ExclusionClosed = as.factor(traits$ExclusionClosed)
	traits$EXVisible_Sneha = as.factor(traits$EXVisible_Sneha)
	traits$ExclusionPendant = as.factor(traits$ExclusionPendant)
	traits$LifeForm_Sneha = as.factor(traits$LifeForm_Sneha)
	traits$FlowerForm_Sneha = as.factor(traits$FlowerForm_Sneha)
	traits$EXDiel_Sneha = as.factor(traits$EXDiel_Sneha)
	traits$EXHandle_Sneha = as.factor(traits$EXHandle_Sneha)
	traits$ExclusionFeeble = as.factor(traits$ExclusionFeeble)
	traits$ExclusionPlatform = as.factor(traits$ExclusionPlatform)
traits = traits[c(2,3,5,7,9,12,15,17,18,20,22,23)]

years = unique(allData_ori$year)

for(year in years){
	allData = allData_ori[which(allData_ori$year == years[1]),]
	anthesis = anthesis_ori[which(anthesis_ori$YEAR == 2011),]

	allData$MeadowW = paste(allData$MEADOW, allData$WATCH, sep = '')
	anthesis$MeadowW = paste(anthesis$MEADOW, anthesis$WATCH, sep = '')

	plantNames = levels(allData$PLTSP_NAME)
	plantNames = plantNames[-c(1)]
	pollNames = levels(allData$VISSP_NAME) # pollinators
	pollNames = pollNames[-c(1)]
	nPlants = length(plantNames)
	nPolls = length(pollNames)

	meadow_watch = unique(allData$MeadowW)
	meadow_watch = meadow_watch[-c(1)]
	nMW = length(meadow_watch)

	A = array(0, dim = c(nPlants, nMW))#nrow=nPlants,ncol=nMWY,dimnames=list(plantNames, meadow_watch_year)) #Colums represent how many meadows
	U = array(0, dim = c(nPlants, nMW))
	Zero = array(0, dim = c(nPlants, 1))
	Prob = array(0, dim = c(nPlants, nMW))

	normal = runif(nPlants,-2,2) 
	random = runif(nPlants,-2,2)
	uniform = matrix(2, nPlants, 1)
	true_special = matrix(-2, nPlants, 1)
	true_special[1] = 2
	half_special = matrix(-2, nPlants, 1)
	half_special[1] = 2
	half_special[6] = 2
	half_special[77] = 2

	# Filling 2D array with availability data
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
		} #pl
	} # mw

	traits2 = data.frame(Biom.flwr = numeric(),
				PlantName = character(),
				LifeForm_Sneha = character(),
				FlowerForm_Sneha = character(),
				ExclusionFeeble = character(),
				generalShapeExclusion = character(),
				ExclusionClosed = character(),
				EXVisible_Sneha = character(),
				ExclusionPendant = character(),
				EXDiel_Sneha = character(),
				EXHandle_Sneha = character(),
				ExclusionPlatform = character())

	traits = traits[c(2,3,5,7,9,12,15,17,18,20,22,23)]

	#get only traits of flowers in study
	for(pl in 1:nPlants){
		traits_temp = traits[which(traits$PlantName == plantNames[pl]),]
		traits2 = rbind(traits2, traits_temp)
	}
	traits = traits2[order(traits2$PlantName),]

	flip_options = c(50)
	num_flip_options = length(flip_options)
	num_starts = 1
	num_true = 3
	trials = 10
	#c(options for flips, options for Phi_test, options for Phi)
	cor_data_kendall = array(0, dim = c(num_flip_options, num_true, trials))
	cor_data_spearman = array(0, dim = c(num_flip_options, num_true, trials))
	cor_data_pearsons = array(0, dim = c(num_flip_options, num_true, trials))
	data_squared_error = array(0, dim = c(num_flip_options, num_true, trials))
	data_chosen_lambda = array(0, dim = c(num_flip_options, num_true, trials))
	
	# in order to create a data frame:
	score_distribution = vector()
	flip_number = vector()
	optim_start = vector()
	cor_kendall = vector()
	cor_spearman = vector()
	cor_pearsons = vector()
	squared_error = vector()
	chosen_lambda = vector()

	for (i in 1:3){
		#true phi function
		if (i == 1){
			Phi_test = true_special
			Phi_name = "true_special"
		} else if (i == 2){
			Phi_test = half_special
			Phi_name = "half_special"
		} else if (i == 3){
			Phi_test = normal
			Phi_name = "normal"
		}

		col_sum = 0
		# To calculate probability
		for (mw in 1:nMW) {
			#Calculating the denominator of probability P(v_t(j) = i)
			col_sum = t(A[,mw])%*%exp(Phi_test)
	
			#calculating probability
			Prob[,mw] = (A[,mw]*exp(Phi_test))/col_sum[1][1]
			if(any(is.na(Prob[,mw]))){
				print(col_sum)
				print(mw)
			}
		} # mw

		for (j in 1:num_flip_options){
			flips = flip_options[j]
	
			#print(c("Phi ", i))
			Phi = Zero
			for (l in 1:trials){
				print(c("trial ", l))

				# Generating Use data
				for(mw in 1:nMW){
					U[,mw] = rmultinom(1, flips, Prob[,mw])	#rmultinom(#trials, #'flips" in each trial, probability vector)
				}

				num_division = 7
				lambda_options = c(0, .25, .5, 1, 2)
				lambda = cross_validization(U, A, nMW, num_division, lambda_options, nPlants)

				result = optim(Phi, likelihood_traits, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")
				cor_data_kendall[j,i,l] = cor(Phi_test, result$par, method= "kendall")
				cor_data_spearman[j,i,l] = cor(Phi_test, result$par, method= "spearman")
				cor_data_pearsons[j,i,l] = cor(Phi_test, result$par, method= "pearson")
				data_chosen_lambda[j,i,l] = lambda

				score_distribution = append(score_distribution, c(Phi_name))
				flip_number = append(flip_number, c(flips))
				cor_kendall = append(cor_kendall, c(cor_data_kendall[j,i,l]))
				cor_spearman = append(cor_spearman, c(cor_data_spearman[j,i,l]))
				cor_pearsons = append(cor_pearsons, c(cor_data_pearsons[j,i,l]))
				chosen_lambda = append(chosen_lambda, c(data_chosen_lambda[j,i,l]))
				print("Trial Done.")
			}
		}
	} #all data gathered for one year

	trial_data = data.frame(score_distribution, 
				flip_number, cor_kendall, cor_spearman,
				cor_pearsons, chosen_lambda)

	filename = paste("MWtraitsCV_", as.character(year), ".RData", sep = '')
	save(trial_data, file=filename)

} # all data gathered for all years
