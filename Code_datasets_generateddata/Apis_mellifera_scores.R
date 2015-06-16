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

a_48 = 0
a_100 = 0
c_48 = 0
c_100 = 0

# Filling 2D array with availability and use data
for (mw in 1:nMW){
	mdata = allData_ori[which(allData_ori$MeadowWY==meadow_watch[mw]),] 
	plantNamesM = levels(mdata$PLTSP_NAME)
	plantNamesM = plantNamesM[-c(1)]
	anthesis_m = anthesis_ori[which(anthesis_ori$MeadowWY==meadow_watch[mw]),]

	if(dim(mdata)[1] < 10){
		print(mw)
		print(dim(mdata))
	}

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
#			if(pl == 48){
#				print(plantNames[pl])
#				print(meadow_watch[mw])
#				print(A[pl,mw])
#				a_48 = a_48 + A[pl,mw]
#			}
#			if(pl == 100){
#				print(plantNames[pl])
#				print(meadow_watch[mw])
#				print(A[pl,mw])
#				a_100 = a_100 + A[pl,mw]
#			}
			if(pl == 88){
				print(meadow_watch[mw])
				print(A[pl,mw])
			}
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
			if(pl == 48){
				c_48 = c_48 + U[pl,mw]
			}
			if(pl == 100){
				c_100 = c_100 + U[pl,mw]
			}
			if(pl == 88){
print("hello")
				print(U[pl,mw])
			}
		}
	} #pl
} # mw

num_division = 7
#lambda_options = c(0, .25, .5, 1, 2)
#lambda = cross_validization(U, A, nMW, num_division, lambda_options, nPlants)

lambda = .25

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

filename = "Apis_mellifera_scores.RData"
save(trial_data, file=filename)

load("Apis_mellifera_scores.RData")

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

range_max = max(test1$Phi_vect, test2$Phi_vect, test3$Phi_vect, test4$Phi_vect, testall$Phi_vect)
range_min = min(test1$Phi_vect, test2$Phi_vect, test3$Phi_vect, test4$Phi_vect, testall$Phi_vect)

par(mfrow=c(5, 1))

plot(x = range(test1$plant_vect), 
	y = range(test1$Phi_vect, test2$Phi_vect, test3$Phi_vect, 
	test4$Phi_vect, testall$Phi_vect),
	main = "Apis Mellifera Scores over 2011",
	xlab = "Plant Species Number",
	ylab = "Plant Score", col = "white")
for (pl in 1:dim(test1)[1]){
	points(x = test1$plant_vect[pl], y = test1$Phi_vect[pl])
}
x = test1$plant_vect
y = array(0, dim = c(length(test1$plant_vect), 1))
lines(x,y)
text(locator(), labels = c("x = 0"))

plot(x = range(test2$plant_vect), 
	y = range(test1$Phi_vect, test2$Phi_vect, test3$Phi_vect, 
	test4$Phi_vect, testall$Phi_vect),
	main = "Apis Mellifera Scores over 2012",
	xlab = "Plant Species Number",
	ylab = "Plant Score", col = "white")
for (pl in 1:dim(test2)[1]){
	points(x = test2$plant_vect[pl], y = test2$Phi_vect[pl])
}
lines(x,y)
text(locator(), labels = c("x = 0"))

plot(x = range(test3$plant_vect), 
	y = range(test1$Phi_vect, test2$Phi_vect, test3$Phi_vect, 
	test4$Phi_vect, testall$Phi_vect),
	main = "Apis Mellifera Scores Over 2013",
	xlab = "Plant Species Number",
	ylab = "Plant Score", col = "white")
for (pl in 1:dim(test3)[1]){
	points(x = test3$plant_vect[pl], y = test3$Phi_vect[pl])
}
lines(x,y)
text(locator(), labels = c("x = 0"))

plot(x = range(test4$plant_vect), 
	y = range(test1$Phi_vect, test2$Phi_vect, test3$Phi_vect, 
	test4$Phi_vect, testall$Phi_vect),
	main = "Apis Mellifera Scores Over 2014",
	xlab = "Plant Species Number",
	ylab = "Plant Score", col = "white")
for (pl in 1:dim(test4)[1]){
	points(x = test4$plant_vect[pl], y = test4$Phi_vect[pl])
}
lines(x,y)
text(locator(), labels = c("x = 0"))

plot(x = range(testall$plant_vect), 
	y = range(test1$Phi_vect, test2$Phi_vect, test3$Phi_vect, 
	test4$Phi_vect, testall$Phi_vect),
	main = "Apis Mellifera Scores Over All Years",
	xlab = "Plant Species Number",
	ylab = "Plant Score", col = "white")
for (pl in 1:dim(testall)[1]){
	points(x = testall$plant_vect[pl], y = testall$Phi_vect[pl])
}
lines(x,y)
text(locator(), labels = c("x = 0"))
