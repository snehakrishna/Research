#evals with cross validation

# Plots to evaluate the multinomial model with count data

setwd(path.expand("~//Research//thesis_stuff"))
source("loglikelihood3.R")
source("cross_validization_general.R")
library(bipartite)

allData_ori = read.csv("Interactions_Final.csv")
anthesis_ori = read.csv("Anthesis_Final.csv")

allData_ori$MeadowW = paste(allData_ori$MEADOW, allData_ori$WATCH, sep = '')
anthesis_ori$MeadowW = paste(anthesis_ori$MEADOW, anthesis_ori$WATCH, sep = '')

allData_ori$MeadowWY = paste(allData_ori$MeadowW, allData_ori$year, sep = '')
anthesis_ori$MeadowWY = paste(anthesis_ori$MeadowW, anthesis_ori$YEAR %% 100, sep = '')

years = unique(allData_ori$year)

Phi_vect = vector()
year_vect = vector()
plant_vect = vector()

#for all meadow watch years
plantNames = levels(allData_ori$PLTSP_NAME)
plantNames = plantNames[-c(1)]
pollNames = levels(allData_ori$VISSP_NAME) # pollinators
pollNames = pollNames[-c(1)]
nPlants = length(plantNames)
nPolls = length(pollNames)

meadow_watch= unique(allData_ori$MeadowWY)
meadow_watch= meadow_watch[-c(1)]
nMW = length(meadow_watch)

A = array(0, dim = c(nPlants, nMW))#nrow=nPlants,ncol=nMWY,dimnames=list(plantNames, meadow_watch_year)) #Colums represent how many meadows
#U = array(0, dim = c(nPlants, nPolls, nMW))
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

# Filling 2D array with availability and use data
for (mw in 1:nMW){
	mdata = allData_ori[which(allData_ori$MeadowWY==meadow_watch[mw]),] 
	plantNamesM = levels(mdata$PLTSP_NAME)
	plantNamesM = plantNamesM[-c(1)]
	anthesis_m = anthesis_ori[which(anthesis_ori$MeadowWY==meadow_watch[mw]),]

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
			int = m_plants[which(m_plants$VISSP_NAME==pollNames[170]),]
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

num_division = 7
#lambda_options = c(0, .25, .5, 1, 2)
#lambda = cross_validization(U, A, nMW, num_division, lambda_options, nPlants)

lambda = .25

result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")

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
#U_true = U[,po,]
count_true = vector()
meadow_use = colSums(U)

testmatrix = model.matrix(~ traits$ExclusionClosed
	+ traits$ExclusionPendant + traits$EXVisible_Sneha
	+ traits$generalShapeExclusion + traits$LifeForm_Sneha
	+ traits$FlowerForm_Sneha + traits$EXDiel_Sneha
	+ traits$EXHandle_Sneha + traits$ExclusionFeeble
	+ traits$ExclusionPlatform)[,-1]

traits$intercept = 1

tmatrix = as.matrix(data.frame(traits$intercept, traits$Biom.flwr, testmatrix))

# traits do seem to contain some information relating to the preference

Zero = array(0, dim = c(dim(tmatrix)[2], 1))

weights = optim(Zero, likelihood_traits, gr = NULL, U, A, 
		tmatrix, lambda, nMW, nPlants, method = "BFGS")
print(weights$par)

Phi_result = tmatrix %*% weights$par

print(cor(result$par, Phi_result, method= "pearson"))
print(cor(result$par, Phi_result, method= "kendall"))
print(cor(result$par, Phi_result, method= "spearman"))

L1 = loglikelihood3(Phi_result, U, A, lambda, nMW, nPlants)
L0 = loglikelihood3(result$par, U, A, lambda, nMW, nPlants)
g2stat = 2*(L0-L1)
lrtest = 1-pchisq(g2stat,df = df)
print(lrtest)
