# goodness of fit tests

setwd(path.expand("~//Research//thesis_stuff"))
source("likelihood_traits.R")
source("cross_validization_traits.R")
library(bipartite)

#allData = read.csv("Final_Interactions.csv")
#anthesis = read.csv("Final_Anthesis.csv")

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

years = unique(allData$year)

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

meadow_watch = unique(allData$MeadowW)
meadow_watch = meadow_watch[-c(1)]
nMW = length(meadow_watch)

meadow_watch_year = unique(allData$MeadowWY)
meadow_watch_year = meadow_watch_year[-c(1)]
nMWY = length(meadow_watch_year)

A = matrix(0,nrow=nPlants,ncol=nMWY,dimnames=list(plantNames, meadow_watch_year)) #Colums represent how many meadows
U = array(0, dim = c(nPlants, nPolls, nMWY))
Zero = array(0, dim = c(nPlants, 1))
Prob = array(0, dim = c(nPlants, nMWY))

# Filling 2D array with availability and use data
remove_watch = vector()
for (mw in 1:nMWY){
	mdata = allData[which(allData$MeadowWY==meadow_watch_year[mw]),] 
	plantNamesM = levels(mdata$PLTSP_NAME)
	plantNamesM = plantNamesM[-c(1)]
	anthesis_m = anthesis[which(anthesis$MeadowWY==meadow_watch_year[mw]),]

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
			for (po in 50:50){
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

#A = A[,-remove_watch]
#U = U[,,-remove_watch]
#nMW = dim(A)[2]

# starting with po=50 (Apis mellifera)
# working with U[,50,]
po = 50
U_true = U[,po,]
count_true = vector()
meadow_use = colSums(U_true)

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

Zero = array(0, dim = c(dim(tmatrix)[2], 1))

# find regularization + phi function
num_division = 7
lambda = 0
Weights = optim(Zero, likelihood_traits, gr = NULL, U_true, A, tmatrix, lambda, nMW, nPlants, method = "BFGS")

Phi = tmatrix %*% Weights$par

# Expected use data
for (mw in 1:nMW) {
	#Calculating the denominator of probability P(v_t(j) = i)
	col_sum = t(A[,mw])%*%exp(Phi)
	#calculating probability
	Prob[,mw] = (A[,mw]*exp(Phi))/col_sum[1][1]
} # mw

E = array(0, dim = c(nPlants, nMW))
for (pl in 1:nPlants){
	for (mw in 1:nMW){
		E[pl, mw] = meadow_use[mw]*Prob[pl, mw]/meadow_use[mw]
	}
}

flips = 10
U_syn = array(0, dim = c(nPlants, nMW))
for(mw in 1:nMW){
	U_syn[,mw] = rmultinom(1, flips, Prob[,mw])
}

count_syn = vector()
count_true = vector()
E_vect = vector()
Prob_vect = vector()
for (pl in 1:nPlants) {
	for (mw in 1:nMW){
		if (A[pl, mw] > 0 && meadow_use[mw] > 0){ # && meadow_use[mw] > 0
			count_syn = append(count_syn, c(U_syn[pl, mw]))
			count_true = append(count_true, c(U_true[pl,mw]))
			E_vect = append(E_vect, c(E[pl,mw]))
			Prob_vect = append(Prob_vect, c(Prob[pl,mw]))
		}
	}
}

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

# want to use count_true because we want to use the field data as our
# observed data. The prob_e is our expected probability
# perform goodness of fit tests
chisq.test(count_true, p = prob_e)
obs = count_true
exp = E_vect
df = length(obs) - length(which(meadow_use[] !=0)) - nPlants - 1
x2stat <- sum(((obs - exp)^2)/exp)
result <- pchisq(x2stat,df,lower.tail=FALSE)
print(result)