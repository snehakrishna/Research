#traits_trial

# goodness of fit tests

setwd(path.expand("~//Research//thesis_stuff"))
source("loglikelihood3.R")
source("cross_validization.R")
library(bipartite)
library(glmnet)

allData = read.csv("Interactions_Final.csv")
anthesis = read.csv("Anthesis_Final.csv")
traits = read.csv("PlantCovars.csv")
#traits = as.factor(traits) #may want to factor just some cols
#traits$Biom.flwr
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

years = unique(allData$year)
year = 2011
allData = allData[which(allData$year == years[1]),]
anthesis = anthesis[which(anthesis$YEAR == 2011),]

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

A = matrix(0,nrow=nPlants,ncol=nMW,dimnames=list(plantNames, meadow_watch)) #Colums represent how many meadows
U = array(0, dim = c(nPlants, nPolls, nMW))
Zero = array(0, dim = c(nPlants, 1))
Prob = array(0, dim = c(nPlants, nMW))

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
		if(dim(m_plants)[1]>0){ 
			#for (po in 1:nPolls){
			#ignoring loop for just apis mellifera for now
			po = 50
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
			#} #po
		}
	} # pl
	if(sum(A[,mw]) == 0){
		# find a mw w/o any plants available
		print(mw)
		zero_use = append(zero_use, mw)
	}
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

# starting with po=50 (Apis mellifera)
# working with U[,50,]
po = 50
U_true = U[,po,]
count_true = vector()
meadow_use = colSums(U_true)

# find regularization + phi function
num_division = 7
lambda_options = c(0, .25, .5, 1, 2)
#lambda = cross_validization(U_true, A, nMW, num_division, 
#					lambda_options, nPlants)

#if no time to do cross validation
lambda = .25

Phi = optim(Zero, loglikelihood3, gr = NULL, U_true, 
		A, lambda, nMW, nPlants, method = "BFGS")

result = lm(Phi$par ~ traits$Biom.flwr + traits$ExclusionClosed
	+ traits$ExclusionPendant + traits$EXVisible_Sneha
	+ traits$generalShapeExclusion + traits$LifeForm_Sneha
	+ traits$FlowerForm_Sneha + traits$EXDiel_Sneha
	+ traits$EXHandle_Sneha + traits$ExclusionFeeble
	+ traits$ExclusionPlatform)

testmatrix = model.matrix(Phi$par ~ traits$ExclusionClosed
	+ traits$ExclusionPendant + traits$EXVisible_Sneha
	+ traits$generalShapeExclusion + traits$LifeForm_Sneha
	+ traits$FlowerForm_Sneha + traits$EXDiel_Sneha
	+ traits$EXHandle_Sneha + traits$ExclusionFeeble
	+ traits$ExclusionPlatform)[,-1]

tmatrix = as.matrix(data.frame(traits$Biom.flwr, testmatrix))

test = glmnet(x = tmatrix, y = Phi$par, family = "gaussian", alpha = 1)
testcv = cv.glmnet(x = tmatrix, y = Phi$par, family = "gaussian", alpha = 1)
coef(testcv, s = "lambda.min")

reg = lm(Phi$par ~ traits$Biom.flwr)
plot(x = log(traits$Biom.flwr), y = Phi$par, main = "Biomass/Flower effect on Phi", 
	xlab = "log(Biomass)", ylab = "Phi")
abline(reg)

plot(x = traits$ExclusionClosed, y = Phi$par), main = "Closed effect on Phi", 
	xlab = "Closed", ylab = "Phi")
plot(x = traits$ExclusionPendant, y = Phi$par, main = "Pendant effect on Phi", 
	xlab = "Pendant", ylab = "Phi")
plot(x = traits$EXVisible_Sneha, y = Phi$par, main = "Visibility effect on Phi", 
	xlab = "Visibility", ylab = "Phi")

reg = lm(Phi$par ~ traits$generalShapeExclusion)
plot(x = traits$generalShapeExclusion, y = Phi$par, main = "General Tube Shape effect on Phi", 
	xlab = "General Tube Shape", ylab = "Phi")
plot(x = traits$LifeForm_Sneha, y = Phi$par, main = "Life Form effect on Phi", 
	xlab = "Life Form", ylab = "Phi")
plot(x = traits$FlowerForm_Sneha, y = Phi$par, main = "Flower Form effect on Phi", 
	xlab = "Flower Form", ylab = "Phi")
plot(x = traits$EXDiel_Sneha, y = Phi$par, main = "Diel effect on Phi", 
	xlab = "Diel", ylab = "Phi")
plot(x = traits$EXHandle_Sneha, y = Phi$par, main = "Pollen Size effect on Phi", 
	xlab = "Pollen Size", ylab = "Phi")

reg = lm(Phi$par ~ traits$ExclusionFeeble)
plot(x = traits$ExclusionFeeble, y = Phi$par, main = "Feebleness effect on Phi", 
	xlab = "Feeble", ylab = "Phi")

#lm(Phi$par ~ traits$Biom.flwr + traits$EXVisible_Sneha)

#L1 regularization using glmnet
#alpha = 1 for lasso
#family = gaussian
# R squared value
# traits do seem to contain some information relating to the preference