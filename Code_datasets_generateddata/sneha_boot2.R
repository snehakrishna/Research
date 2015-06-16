sneha_boot = function(U, nMW, nPlants, Zero, A){
	data = array(0, dim = dim(U))
	for (mw in 1:dim(U)[2]){
		if(sum(U[,mw]) == 0){
			next
		}
		S = vector()
		for(pl in 1:dim(U)[1]){
			S = append(S, rep(pl, U[pl,mw]))
		}
		print(U[,mw])
		print(S)
		if(length(S)==1){
			V = S
		}
		else{
			V = sample(S, replace=T)
		}

		for(pl in 1:dim(U)[1]){
			data[pl,mw] = length(V[which(V==pl)])
			if(data[pl,mw] >0 && A[pl,mw] == 0){
				#print('inhere')
				#print(pl)
				#print(mw)
			}
		}
		
	}
	lambda = .25
	result = optim(Zero, loglikelihood3, gr = NULL, data, A, lambda, nMW, nPlants, method = "BFGS")
	return(result$par)
}

sneha_boot_results = function(U, nMW, nPlants, Zero, A, n){
	res = array(0, dim = c(nPlants, n))
	for(i in 1:n){
		res[,i] = sneha_boot(U, nMW, nPlants, Zero, A)
	}
	upper = vector()
	lower = vector()
	for(pl in 1:nPlants){
		s = sd(res[pl,])
		a = mean(res[pl,])
		error = qnorm(0.975)*s/sqrt(n)
		lower = append(lower, a-error)
		upper = append(upper, a + error)
	}
	return(data.frame(upper, lower))
}

setwd(path.expand("~//Research//thesis_stuff"))
source("loglikelihood3.R")
source("cross_validization_general.R")

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

n = 200
sboot = sneha_boot_results(U, nMW, nPlants, Zero, A, n)

save(sboot, file = "sboot_200.RData")
