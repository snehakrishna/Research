sneha_boot = function(U, nMW, nPlants, Zero, A){
	print("here")
	data = array(0, dim = dim(U))
print("here3")
	for (mw in 1:dim(U)[2]){
print("here4")
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
		}
		
	}
	lambda = .25
	print("here1")
	result = optim(Zero, loglikelihood3, gr = NULL, data, A, lambda, nMW, nPlants, method = "BFGS")
	print("here2")
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

years = unique(allData_ori$year)

for(year in years){
	allData = allData_ori[which(allData_ori$year == year),]
	anthesis = anthesis_ori[which(anthesis_ori$YEAR%%100 == year),]

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

	n = 100
	sboot = sneha_boot_results(U, nMW, nPlants, Zero, A, n)

	filename = paste("sboot_100_", as.character(year), ".RData", sep = '')

	save(sboot, file = filename)
}