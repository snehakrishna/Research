# activity levels of pollinators

setwd("C:/Users/Sneha Krishna/Documents/Research")
library(bipartite)
allData = read.csv("11to13Interactions_v5.csv")
#allData = as.data.frame.matrix(read.csv("11to13Interactions_v5.csv"))
allData$MeadowW = paste(allData$MEADOW, allData$WATCH, sep = '')

years = unique(allData$year)
nyears = length(years)

high_activity_poll = vector()
high_activity_meadow = vector()
high_activity_year = vector()

for(y in 1:nyears){
	year = years[y]
	yearData = allData[which(allData$year == year),]

	plantNames = levels(yearData$PLTSP_NAME)
	plantNames = plantNames[-c(1)]
	pollNames = levels(yearData$VISSP_NAME) # pollinators
	pollNames = pollNames[-c(1)]
	nPlants = length(plantNames)
	nPolls = length(pollNames)

	meadow_watch = unique(allData$MeadowW)
	meadow_watch = meadow_watch[-c(1)]
	nMW = length(meadow_watch)
	#nMW = 1

	use = matrix(0,nrow=nPolls,ncol=nMW,dimnames=list(pollNames, meadow_watch))
	use_v = vector()

	for (po in 1:nPolls) {
		all_ints = yearData[which(allData$VISSP_NAME==pollNames[po]),]
		if (dim(all_ints)[1]>0) { 
			for(m in 1:nMW){
				mdata = yearData[which(allData$MeadowW==meadow_watch[m]),]
				pollNamesM = levels(mdata$VISSP_NAME)
				pollNamesM = pollNamesM[-c(1)]

				int_meadow = mdata[which(mdata$VISSP_NAME==pollNamesM[po]),]
		
				#m_plants = mdata[which(mdata$Plant==plantNames[pl]),]

				use_v = append(use_v, c(dim(int_meadow)[1]))

				if(dim(int_meadow)[1]>0) {
					use[po, m] = dim(int_meadow)[1]
				}

				if(dim(int_meadow)[1]>100){
					high_activity_poll = append(high_activity_poll, c(m))
					high_activity_meadow = append(high_activity_meadow, c(po))
					high_activity_year = append(high_activity_year, c(year))
				}
			}
		}
	} # po

	print(year)
	if(!is.na(year) && year == 2011){
		print("entered 2011")
		use_11 = use
		use_v_11 = use_v
	}
	else if(!is.na(year) && year == 2012){
		print("entered 2012")
		use_12 = use
		use_v_12 = use_v
	}
	else if(!is.na(year) && year == 2013){
		print("entered 2013")
		use_13 = use
		use_v_13 = use_v
	}
	else{
		use_NA = use
		use_v_NA = use_v
	}
}
hist(use_v, breaks=50)

plot(use_v_11, main="Activity of Pollinators by meadow watch--2011", 
  	xlab="", ylab="Number of interactions per Meadow Watch", pch=19)

plot(use_v_12, main="Activity of Pollinators by meadow watch--2012", 
  	xlab="", ylab="Number of interactions per Meadow Watch", pch=19)

plot(use_v_13, main="Activity of Pollinators by meadow watch--2013", 
  	xlab="", ylab="Number of interactions per Meadow Watch", pch=19)

total_11 = rowSums(use_11)
total_12 = rowSums(use_12)
total_13 = rowSums(use_13)

plot(total_11, main="Activity of Pollinators for 2011", 
  	xlab="", ylab="Number of interactions per Meadow Watch", pch=19)
plot(total_12, main="Activity of Pollinators for 2012", 
  	xlab="", ylab="Number of interactions per Meadow Watch", pch=19)
plot(total_13, main="Activity of Pollinators for 2013", 
  	xlab="", ylab="Number of interactions per Meadow Watch", pch=19)

hist(total_11, 1000)

visweb(use, type = "none")