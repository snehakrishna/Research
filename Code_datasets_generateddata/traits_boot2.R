sneha_boot = function(U, nMW, nPlants, Zero, A, tmatrix){
	data = array(0, dim = dim(U))
	print("here")
	for (mw in 1:dim(U)[2]){
		if(sum(U[,mw]) == 0){
			next
		}
		S = vector()
		for(pl in 1:dim(U)[1]){
			S = append(S, rep(pl, U[pl,mw]))
		}
		if(length(S)==1){
			V = S
		}
		else{
			V = sample(S, replace=T)
		}

		for(pl in 1:dim(U)[1]){
			data[pl,mw] = length(V[which(V==pl)])
		}
		
	}
print("here1")
	lambda = 0
	result = optim(Zero, likelihood_traits, gr = NULL, data, A, 
				tmatrix, lambda, nMW, nPlants, method = "BFGS")
print("here2")
	print(result)
	return(result$par)
}

sneha_boot_results = function(U, nMW, nPlants, Zero, A, n, tmatrix){
	ntraits = dim(tmatrix)[2]
	res = array(0, dim = c(ntraits, n))
	for(i in 1:n){
		res[,i] = sneha_boot(U, nMW, nPlants, Zero, A, tmatrix)
	}
	upper = vector()
	lower = vector()
	for(t in 1:ntraits){
		s = sd(res[t,])
		a = mean(res[t,])
		error = qnorm(0.975)*s/sqrt(n)
		lower = append(lower, a-error)
		upper = append(upper, a + error)
	}
	return(data.frame(upper, lower))
}

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
U = array(0, dim = c(nPlants, nMW))

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
		if(dim(m_plants)[1]>0){ 
			int = m_plants[which(m_plants$VISSP_NAME==pollNames[50]),]
			ints = as.numeric(as.character(int$NO_INT))
			if(any(is.na(ints))){
				ints = na.omit(ints)
			}
			if(A[pl,mw] == 0){
				U[pl,mw] = 0
			}
			else { 
				U[pl,mw] = sum(ints)
			}
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

n = 200
sboot = sneha_boot_results(U, nMW, nPlants, Zero, A, n, tmatrix)

save(sboot, file = "tboot_200.RData")
