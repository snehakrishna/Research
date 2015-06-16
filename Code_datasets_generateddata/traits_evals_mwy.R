#evals with no cross validation

# Plots to evaluate the multinomial model with count data

setwd(path.expand("~//Research//thesis_stuff"))
source("likelihood_traits.R")
library(bipartite)

allData = read.csv("Interactions_Final.csv")
anthesis = read.csv("Anthesis_Final.csv")
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

allData$MeadowW = paste(allData$MEADOW, allData$WATCH, sep = '')
anthesis$MeadowW = paste(anthesis$MEADOW, anthesis$WATCH, sep = '')

allData$MeadowWY = paste(allData$MeadowW, allData$year, sep = '')
anthesis$MeadowWY = paste(anthesis$MeadowW, anthesis$YEAR%%100, sep = '')

plantNames = levels(allData$PLTSP_NAME)
plantNames = plantNames[-c(1)]
pollNames = levels(allData$VISSP_NAME) # pollinators
pollNames = pollNames[-c(1)]
nPlants = length(plantNames)
nPolls = length(pollNames)

meadow_watch = unique(allData$MeadowWY)
meadow_watch = meadow_watch[-c(1)]
nMW = length(meadow_watch)

A = array(0, dim = c(nPlants, nMW))

# Filling 2D array with availability data
remove_watch = vector()
for (mw in 1:nMW){
	mdata = allData[which(allData$MeadowWY==meadow_watch[mw]),] 
	plantNamesM = levels(mdata$PLTSP_NAME)
	plantNamesM = plantNamesM[-c(1)]
	anthesis_m = anthesis[which(anthesis$MeadowWY==meadow_watch[mw]),]

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
	if(sum(A[,mw]) == 0){
		remove_watch = append(mw, remove_watch)
	}
} # mw

A = A[,-remove_watch]
nMW = dim(A)[2]

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

#get only traits of flowers in study
for(pl in 1:nPlants){
	traits_temp = traits[which(traits$PlantName == plantNames[pl]),]
	traits2 = rbind(traits2, traits_temp)
}
traits = traits2[order(traits2$PlantName),]

testmatrix = model.matrix(~ traits$ExclusionClosed
	+ traits$ExclusionPendant + traits$EXVisible_Sneha
	+ traits$generalShapeExclusion + traits$LifeForm_Sneha
	+ traits$FlowerForm_Sneha + traits$EXDiel_Sneha
	+ traits$EXHandle_Sneha + traits$ExclusionFeeble
	+ traits$ExclusionPlatform)[,-1]

traits$intercept = 1

tmatrix = as.matrix(data.frame(traits$intercept, traits$Biom.flwr, testmatrix))

U = array(0, dim = c(nPlants, nMW))
Zero = array(0, dim = c(dim(tmatrix)[2], 1))
Prob = array(0, dim = c(nPlants, nMW))

#normal = runif(dim(tmatrix)[2]+1, -1, 1) # assigning random weights
normal = runif(nPlants,-2,2) 
random = runif(nPlants,-2,2)
uniform = matrix(2, nPlants, 1)
true_special = matrix(-2, nPlants, 1)
true_special[1] = 2
half_special = matrix(-2, nPlants, 1)
half_special[1] = 2
half_special[6] = 2
half_special[77] = 2

flip_options = c(50)
num_flip_options = 1
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

Weight1 = Zero
Weight1[3] = 2
Weight1[13] = 2
Weight2 = Zero
Weight2[2] = .5
Weight2[5] = 2
Weight3 = runif(dim(tmatrix)[2],-1.5,1.5)
Weight3[2] = .5

for (i in 1:3){
	#true phi function
	if (i == 1){
		Weight_test = Weight1
		Phi_test = tmatrix %*% Weight_test
		Phi_name = "true_special"
		#Weight_test = solve(t(tmatrix) %*% tmatrix) %*% t(tmatrix) %*% Phi_test
	} else if (i == 2){
		Weight_test = Weight2
		Phi_test = tmatrix %*% Weight_test
		Phi_name = "half_special"
		#Weight_test = solve(t(tmatrix) %*% tmatrix) %*% t(tmatrix) %*% Phi_test
	} else if (i == 3){
		Weight_test = Weight3
		Phi_test = tmatrix %*% Weight_test
		Phi_name = "Normal Distribution"
		#Weight_test = solve(t(tmatrix) %*% tmatrix) %*% t(tmatrix) %*% Phi_test
	}

	col_sum = 0
	# To calculate probability

	for (mw in 1:nMW) {
		#Calculating the denominator of probability P(v_t(j) = i)
		col_sum = t(A[,mw])%*% exp(Phi_test)

		#calculating probability
		Prob[,mw] = (A[,mw]*exp(Phi_test))/col_sum[1][1]
		if(any(is.infinite(Prob[,mw]))){
			Prob[,mw] = 0
		}
	} # mw

	for (j in 1:num_flip_options){
		flips = flip_options[j]

		#print(c("Weight ", i))
		for (l in 1:trials){
			print(c("trial ", l))

			# Generating Use data
			for(mw in 1:nMW){
				U[,mw] = rmultinom(1, flips, Prob[,mw])	#rmultinom(#trials, #'flips" in each trial, probability vector)
			}

			lambda = 0
			result = optim(Zero, likelihood_traits, gr = NULL, U, A, tmatrix, lambda, nMW, nPlants, method = "BFGS")

			Phi_result = tmatrix %*% result$par

			cor_data_kendall[j,i,l] = cor(Phi_test, Phi_result, method= "kendall")
			cor_data_spearman[j,i,l] = cor(Phi_test, Phi_result, method= "spearman")
			cor_data_pearsons[j,i,l] = cor(Phi_test, Phi_result, method= "pearson")

			score_distribution = append(score_distribution, c(Phi_name))
			flip_number = append(flip_number, c(flips))
			cor_kendall = append(cor_kendall, c(cor_data_kendall[j,i,l]))
			cor_spearman = append(cor_spearman, c(cor_data_spearman[j,i,l]))
			cor_pearsons = append(cor_pearsons, c(cor_data_pearsons[j,i,l]))
			print("Trial Done.")
		}
	}
} #all data gathered for one year

trial_data = data.frame(score_distribution, 
				flip_number, cor_kendall, cor_spearman,
				cor_pearsons)

filename = "MWYtraitsCVno.RData"
save(trial_data, file=filename)

names_vect = c("traits.intercept","traits.Biom.flwr", 
	"traits.ExclusionClosedopen.access",
	"traits.ExclusionPendantsuspended", 
	"traits.EXVisible_Snehanotbright", 
	"traits.generalShapeExclusionmoderate.exclusion",
	"traits.generalShapeExclusionpoor.exclusion",
	"traits.generalShapeExclusionsevere.exclusion",
	"traits.LifeForm_SnehaP",
	"traits.FlowerForm_Snehaexclusion",
	"traits.EXDiel_Snehayesnight",
	"traits.EXHandle_Snehaokay",
	"traits.ExclusionFeeblestrong",
	"traits.ExclusionPlatformplatform",
	"traits.ExclusionPlatformweak"   )

ntraits = dim(tmatrix)[2]

x = seq(-2,2,by=0.01)
y = x

#result1 -- Weight1
#result2 -- Weight2
#result3 -- Weight3

plot(range(Weight_test), range(result$par), 
	main="Comparison between True and Fitted Weight Vectors",
  	xlab="True Weights", ylab="Fitted Weights", pch=19, col = "white")
for(i in 1:ntraits){
	points(Weight1[i], result1[i], col = "blue", pch = 19)
	points(Weight2[i], result2[i], col = "green", pch = 19)
	points(Weight3[i], result3[i], col = "black", pch = 19)
}
lines(x, y)

plot(range(Phi1, Phi2, Phi3), range(rPhi1, rPhi2, rPhi3), 
	main="Comparison between True and Fitted Weight Vectors",
  	xlab="True Weights", ylab="Fitted Weights", pch=19, col = "white")
for(i in 1:ntraits){
	points(Weight1[i], result1[i], col = "blue", pch = 19)
	points(Weight2[i], result2[i], col = "green", pch = 19)
	points(Weight3[i], result3[i], col = "black", pch = 19)
}
lines(x, y)