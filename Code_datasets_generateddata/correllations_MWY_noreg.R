#Graphs 121114_<x>.eps

setwd("C:/Users/Sneha Krishna/Documents/Research/thesis_stuff")
source("loglikelihood3.r")
library(bipartite)
library(Metrics)
library("reshape2")

year = 11
year_long = year + 2000

load("MWYtraitsCVno_weights.RData")
#load("MWYevalsCVno_14.RData")
#load("MWYevalsCVset.RData")
lambda = 0

new_data = trial_data
new_data$flip_number = NULL

new_data = trial_data[which(trial_data$score_distribution == "true_special"),]
mean(new_data$cor_kendall)
sd(new_data$cor_kendall)
range(new_data$cor_kendall)
mean(new_data$cor_spearman)
sd(new_data$cor_spearman)
range(new_data$cor_spearman)
mean(new_data$cor_pearsons)
sd(new_data$cor_pearsons)
range(new_data$cor_pearsons)

#Figure 1
data <- melt(new_data, id.var = "score_distribution")
par(mfrow=c(2, 2))
boxplot(value ~ variable, data = data, main="Comparing 3 types of correlations", 
  	xlab="Type of correlation", ylab="Value of correlation")
boxplot(cor_kendall ~ score_distribution, data = new_data, main="Kendall's tau", 
  	xlab="Type of score function", ylab="Value of correlation")
boxplot(cor_spearman ~ score_distribution, data = new_data, main="Spearman's rho", 
  	xlab="Type of score function", ylab="Value of correlation")
boxplot(cor_pearsons ~ score_distribution, data = new_data, main="Pearson's product-moment", 
  	xlab="Type of score function", ylab="Value of correlation")

#Figure 3, 4, 5
allData = read.csv("Interactions_Final.csv")
anthesis = read.csv("Anthesis_Final.csv")

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

A = array(0, dim = c(nPlants, nMW))#nrow=nPlants,ncol=nMWY,dimnames=list(plantNames, meadow_watch_year)) 
	#Colums represent how many meadows

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
	} #pl
	if(sum(A[,mw]) == 0){
		remove_watch = append(mw, remove_watch)
	}
} # mw

A = A[,-remove_watch]
nMW = dim(A)[2]

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

avail = rowSums(A)

flips = 50

x = seq(-2,2,by=0.01)
y = x

# Figure 3
Phi = normal
for (mw in 1:nMW) {
	#Calculating the denominator of probability P(v_t(j) = i)
	col_sum = t(A[,mw])%*%exp(Phi)
	#calculating probability
	Prob[,mw] = (A[,mw]*exp(Phi))/col_sum[1][1]
} # mw
# Synthetic use data
for(mw in 1:nMW){
	U[,mw] = rmultinom(1, flips, Prob[,mw])
}
use = rowSums(U)
lambda = .25
result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")
par(cex.lab = 1.25)
plot(range(Phi), range(result$par), 
	main=paste("Normal Distribution\nwith Regularization", sep = ''), 
  	xlab="True Phi", ylab="Calculated Phi", pch=19, col = "white")
for(u in 1:length(use)){
	if(use[u] == 0){
		if(avail[u] == 0){
			points(Phi[u], result$par[u], col = "blue", pch = 19)
		}
		else{
			points(Phi[u], result$par[u], col = "green", pch = 19)
		}
	}
	else if(use[u] < 11){
		points(Phi[u], result$par[u], col = "red", pch = 19)
	}
	else{
		points(Phi[u], result$par[u], pch = 19)
	}
}
lines(x, y)
legend("bottomright", inset=.05, title="Legend",
  	c("Use = 0, Availability = 0","Use = 0, Availability != 0",
	"0 < Use < 11", "11 <= Use"), fill=c("blue", "green", "red", "black"), 
	horiz=FALSE)
#plot(normal, result$par, main="Normal Distribution", 
#  	xlab="True Phi", ylab="Calculated Phi", pch=19)


# Figure 4
Phi = half_special
for (mw in 1:nMW) {
	#Calculating the denominator of probability P(v_t(j) = i)
	col_sum = t(A[,mw])%*%exp(Phi)
	#calculating probability
	Prob[,mw] = (A[,mw]*exp(Phi))/col_sum[1][1]
} # mw
# Synthetic use data
for(mw in 1:nMW){
	U[,mw] = rmultinom(1, flips, Prob[,mw])
}
use = rowSums(U)
result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")
plot(range(Phi), range(result$par),
	main=paste("Half Special Distribution\nwithout Regularization", sep = ''), 
  	xlab="True Phi", ylab="Calculated Phi", pch=19)
for(u in 1:length(use)){
	if(use[u] == 0){
		if(avail[u] == 0){
			points(Phi[u], result$par[u], col = "blue", pch = 19)
		}
		else{
			points(Phi[u], result$par[u], col = "green", pch = 19)
		}
	}
	else if(use[u] < 11){
		points(Phi[u], result$par[u], col = "red", pch = 19)
	}
	else{
		points(Phi[u], result$par[u], pch = 19)
	}
}
lines(x, y)
legend("bottomright", inset=.05, title="Legend",
  	c("Use = 0, Availability = 0","Use = 0, Availability != 0",
	"0 < Use < 11", "11 <= Use"), fill=c("blue", "green", "red", "black"), 
	horiz=FALSE)
#plot(Phi, result$par, main="Half Special Distribution", 
#  	xlab="True Phi", ylab="Calculated Phi", pch=19)


#Figure 5
Phi = true_special
for (mw in 1:nMW) {
	#Calculating the denominator of probability P(v_t(j) = i)
	col_sum = t(A[,mw])%*%exp(Phi)
	#calculating probability
	Prob[,mw] = (A[,mw]*exp(Phi))/col_sum[1][1]
} # mw
# Synthetic use data
for(mw in 1:nMW){
	U[,mw] = rmultinom(1, flips, Prob[,mw])
}
use = rowSums(U)
result = optim(Zero, loglikelihood3, gr = NULL, U, A, lambda, nMW, nPlants, method = "BFGS")
plot(range(Phi), range(result$par),
	main=paste("True Special Distribution\nwithout Regularization", sep = ''), 
  	xlab="True Phi", ylab="Calculated Phi", pch=19)
for(u in 1:length(use)){
	if(use[u] == 0){
		if(avail[u] == 0){
			points(Phi[u], result$par[u], col = "blue", pch = 19)
		}
		else{
			points(Phi[u], result$par[u], col = "green", pch = 19)
		}
	}
	else if(use[u] < 11){
		points(Phi[u], result$par[u], col = "red", pch = 19)
	}
	else{
		points(Phi[u], result$par[u], pch = 19)
	}
}
lines(x, y)
legend("bottomright", inset=.05, title="Legend",
  	c("Use = 0, Availability = 0","Use = 0, Availability != 0",
	"0 < Use < 11", "11 <= Use"), fill=c("blue", "green", "red", "black"), 
	horiz=FALSE)
#plot(Phi, result$par, main="True Special Distribution", 
#  	xlab="True Phi", ylab="Calculated Phi", pch=19)

#null.probs = Zero
#chisq.test(result$par, p=null.probs)



#Figure
