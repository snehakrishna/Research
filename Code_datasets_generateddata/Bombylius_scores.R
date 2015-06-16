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
			int = m_plants[which(m_plants$VISSP_NAME==pollNames[74]),]
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
lambda_options = c(0, .25, .5, 1, 2)
lambda = cross_validization(U, A, nMW, num_division, lambda_options, nPlants)

result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")

for(pl in 1:length(result$par)){
	Phi_vect = append(result$par[pl], Phi_vect)
	year_vect = append(999, year_vect)
	plant_vect = append(pl, plant_vect)
}

#getting vectors for each year
for(year in years){
	allData = allData_ori[which(allData_ori$year == year),]
	anthesis = anthesis_ori[which(anthesis_ori$YEAR%%100 == year),]

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
				int = m_plants[which(m_plants$VISSP_NAME==pollNames[74]),]
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
	lambda_options = c(0, .25, .5, 1, 2)
	lambda = cross_validization(U, A, nMW, num_division, lambda_options, nPlants)

	result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")

	for(pl in 1:length(result$par)){
		Phi_vect = append(result$par[pl], Phi_vect)
		year_vect = append(year, year_vect)
		plant_vect = append(pl, plant_vect)
	}

} # all data gathered for all years

trial_data = data.frame(Phi_vect, year_vect, plant_vect)

filename = "Bombylius_scores.RData"
save(trial_data, file=filename)

load("Bombylius_scores.RData")

test1 = trial_data[which(trial_data$year_vect == 11),]
test2 = trial_data[which(trial_data$year_vect == 12),]
test3 = trial_data[which(trial_data$year_vect == 13),]
test4 = trial_data[which(trial_data$year_vect == 14),]
testall = trial_data[which(trial_data$year_vect == 999),]

test1 = test1[order(test1$Phi_vect),]
test2 = test2[order(test2$Phi_vect),]
test3 = test3[order(test3$Phi_vect),]
test4 = test4[order(test4$Phi_vect),]
testall = testall[order(testall$Phi_vect),]