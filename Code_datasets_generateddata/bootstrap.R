setwd(path.expand("~//Research//thesis_stuff"))
source("loglikelihood3.R")
source("cross_validization_general.R")
library(bipartite)
library(boot)

likelihood_boot <- function(data, indices, nMW, nPlants, Zero, A){
	U <- data[indices,] # allows boot to select sample 
	A <- A[indices,]
	lambda = .25
	result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")
	return(result$par)
}

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
po = 50

meadow_watch= unique(allData_ori$MeadowWY)
meadow_watch= meadow_watch[-c(1)]
nMW = length(meadow_watch)

A = array(0, dim = c(nPlants, nMW))#nrow=nPlants,ncol=nMWY,dimnames=list(plantNames, meadow_watch_year)) 
	#Colums represent how many meadows
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
			int = m_plants[which(m_plants$VISSP_NAME==pollNames[po]),]
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

likelihood_cens <- function(data, nMW, nPlants, Zero, A){
	U = data
	lambda = .25
	result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")
	print("done")
	return(result)
}

bootObj = boot(data = U, statistic = likelihood_boot, R = 1, nMW = nMW, nPlants = nPlants, Zero = Zero, A = A)

censbootObj = censboot(data = U, statistic = likelihood_cens, R = 1, nMW = nMW, nPlants = nPlants, Zero = Zero, A = A)