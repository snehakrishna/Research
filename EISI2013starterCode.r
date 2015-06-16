# Set your working directory to wherever you put the data:
setwd("C:/Users/Sneha Krishna/Documents/Research")

# make sure the bipartite package is installed, then load it:
library(bipartite)

# need to run set options for bipartite once per session
options(expressions=50000)

# load an RData file containing the interaction survey data from 2011
#  in a data.frame called allData
load("2011datax2.RData")

################################
# Loading and visualizing data #
################################

# explore the data a little bit
names(allData) 
levels(allData$MeadowID) 
unique(allData$plotNum)

# get the names and numbers of the plants, polls, ints 
# (and remove non-names)
plantNames = levels(allData$Plant)
plantNames = plantNames[-c(1)]
pollNames = levels(allData$Andy.Ident) # pollinators
pollNames = pollNames[-c(1)]
intNames = levels(allData$Int)
intNames = intNames[-c(1)]
nPlants = length(plantNames)
nPolls = length(pollNames)
nInts = length(intNames)

# build matrices for bipartite of all the interactions in the dataset.
# this takes awhile since we loop over all plants and all pollinators
allIntMatwithNAs = matrix(0,nrow=nPlants,ncol=nPolls,dimnames=list(plantNames,pollNames))
allIntMatwithoutNAs = matrix(0,nrow=nPlants,ncol=nPolls,dimnames=list(plantNames,pollNames))
allIntMatPA = matrix(0,nrow=nPlants,ncol=nPolls,dimnames=list(plantNames,pollNames))
for (pl in 1:nPlants) {
  for (po in 1:nPolls) { 
    thisPPsubset = allData[intersect(which(allData$Plant==plantNames[pl]),which(allData$Andy.Ident==pollNames[po])),]
    if (length(thisPPsubset)>0) { 
      allIntMatwithoutNAs[pl,po] = sum(thisPPsubset$Interactions,na.rm=T)
      allIntMatwithNAs[pl,po] = sum(thisPPsubset$Interactions)
      allIntMatPA[pl,po] = as.numeric(sum(thisPPsubset$INT.PA.)>0)
    }
  } # po
} # pl

# check out a subset of these matrices: number of interactions
#  observed between plants (on the rows) and pollinators (on the
#  columns)
allIntMatwithNAs[1:10,1:10]
allIntMatwithoutNAs[1:10,1:10]
allIntMatPA[1:10,1:10]

# visualize these interactions (also takes awhile!)
visweb(allIntMatwithNAs)
visweb(allIntMatwithoutNAs)
visweb(allIntMatPA)
# that's hard to look at - try another way
plotweb(allIntMatwithNAs)
plotweb(allIntMatwithoutNAs)
plotweb(allIntMatPA)
# still a lot - narrow to a single meadow? 

# Repeat this process for just one meadow, "M2"
M2data = allData[which(allData$MeadowID=="M2"),]

plantNamesM2 = levels(M2data$Plant)
plantNamesM2 = plantNamesM2[-c(1)]
nPlantsM2 = length(plantNamesM2)

pollNamesM2 = levels(M2data$Andy.Ident)
pollNamesM2 = pollNamesM2[-c(1)]
nPollsM2 = length(pollNamesM2)

M2intMatwithNAs = matrix(0,nrow=nPlantsM2,ncol=nPollsM2,dimnames=list(plantNamesM2,pollNamesM2))
M2intMatwithoutNAs = matrix(0,nrow=nPlantsM2,ncol=nPollsM2,dimnames=list(plantNamesM2,pollNamesM2))
M2intMatPA = matrix(0,nrow=nPlants,ncol=nPolls,dimnames=list(plantNamesM2,pollNamesM2))
total_ints = 0
for (pl in 1:nPlantsM2) {
  for (po in 1:nPollsM2) {
    thisPPsubset = M2Data[intersect(which(M2Data$PLTSP_NAME==plantNamesM2[pl]),which(M2Data$VISSP_NAME==pollNamesM2[po])),]
    if (length(thisPPsubset)>0) { 
      M2intMatwithoutNAs[pl,po] = sum(thisPPsubset$NO_INT,na.rm=T)
      M2intMatwithNAs[pl,po] = sum(thisPPsubset$NO_INT)
      M2intMatPA[pl,po] = as.numeric(sum(thisPPsubset$NO_INT)>0)
	total_ints = total_ints + M2intMatwithoutNAs[pl,po]
    }
  } # po
} # pl

M2intMatwithNAs[1:10,1:10]
M2intMatwithoutNAs[1:10,1:10]
M2intMatPA[1:10,1:10]

# visualize M2 network
visweb(allIntMatwithNAs)
visweb(allIntMatwithoutNAs)
visweb(allIntMatPA) 

plotweb(allIntMatwithNAs)
plotweb(allIntMatwithoutNAs)
plotweb(allIntMatPA)


# examples of how to write visualizations to files (change filenames)
png(file="intMatM2bin.png",height=800,width=1200)
visweb(M2intMatwithoutNAs>0,prednames=FALSE,preynames=FALSE)
dev.off()

jpeg(file="intWebM2abund.png",height=800,width=1200)
plotweb(M2intMatwithoutNAs)
dev.off()

# create images for all subplots for M2
for (pn in 1:10) {
  print(pn)

  M2intMatPlot1 = matrix(0,nrow=nPlantsM2,ncol=nPollsM2,dimnames=list(plantNamesM2,pollNamesM2))

  for (pl in 1:nPlantsM2) {
    for (po in 1:nPollsM2) {
      M2intMatPlot1[pl,po] = sum(M2data[intersect(intersect(which(M2data$Plant==plantNamesM2[pl]),which(M2data$Andy.Ident==pollNamesM2[po])),which(M2data$plotNum==pn)),]$Interactions)
    } # po
  } # pl
  
  pdf(file=paste("intMatPlot",pn,".pdf",sep=""))
  visweb(M2intMatPlot1)
  dev.off()
} # pn

# Potential exercise: build these for another meadow to practice
#  subsetting the data, do binary versions instead of abundance, etc.

# NEXT: check out the network metrics that the bipartite package
#  provides

metricsM2 = networklevel(M2intMatPA)
names(metricsM2)
# check out the documentation to start learning about what these mean:
?networklevel

metricsM2[7] # nestedness

# Exercise idea: Calculate nestedness at the 10 subplots of M2
#  and see how consistent it is.

# NEXT: more plotting

# Here is another matrix plotting function that I grabbed from the
#  web, and modified slightly.  This one does NOT rearrange rows and
#  columns to maximize nestedness, as bipartite does.

# ----- Define a function for plotting a matrix ----- #
myImagePlotAbun <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
#layout(matrix(data=c(1,2), nrow=1, ncol=1), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 #ColorRamp <- rgb( seq(0,1,length=256),  # Red
 #                  seq(0,1,length=256),  # Green
 #                  seq(1,0,length=256))  # Blue
 ColorRamp <- gray(seq(1,0,length=256)) # grayscale version
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 ###########par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 #image(1:length(xLabels), 1:length(yLabels), t(x), col=gray(c(1,0)), xlab="",
 #ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
# ----- END plot function ----- #

# Plot for all of M2
myImagePlotAbun(M2intMat)

# Plot for subplot 10 - the last one we computed above.  Still has all
#  the species in M2, even if they don't appear in the subplot.

myImagePlotAbun(M2intMatPlot1)

# Potential exercise: use this plotting function to look at
#  variability across subplots of M2.

##############
# Clustering #
##############

# load a plot by plant matrix for clustering by veg.
# this data is from the flower abundance surveys in 2011.
# 150 plots: 15 meadows by 10 subplots.
# aggregated plants over the entire summer.

load("flowerSurveyAbundance2011.RData")

flwAbundPlot = t(flwAbundPlot)
myImagePlotAbun(flwAbundPlot)

# remove flowers that were never seen
neverSeenIdxs = which(colSums(flwAbundPlot)==0)
flwAbundPlot = flwAbundPlot[,-neverSeenIdxs]

library(cluster)

dev.new()
# how many clusters?
wss = (nrow(flwAbundPlot)-1)*sum(apply(flwAbundPlot,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(flwAbundPlot,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# Look at solutions with different nClusters
nClusters = 6
fit <- kmeans(flwAbundPlot, nClusters) 
# get cluster means 
aggregate(flwAbundPlot,by=list(fit$cluster),FUN=mean)
# append cluster assignment
data.frame(flwAbundPlot, fit$cluster)
# Cluster Plot against 1st 2 principal components
dev.new()
# vary parameters for most readable graph
clusplot(flwAbundPlot, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

# Next steps/exercises: do this on different subsets of the data, look
#  into clusters (one is consistently CNM, what's going on there?),
#  find other 2-d projections to visualize clusters better, etc.

######################
# Occupancy modeling #
######################

library(unmarked)

load("2011detHists.RData") 
# contains: plantDetHists, pollDetHists, intDetHists

logistic = function(x) {
  y = (1 + exp(-x))^-1
  y
}

# format data for occupancy model for first interaction
occuY = unmarkedFrameOccu(intDetHists[1,,])

# fit constant occupancy model
occuOut = occu(~1 ~1, occuY)
occupancyRate = logistic(occuOut[1]@estimates)
detectionRate = logistic(occuOut[2]@estimates)

# I have not prepared any covariates yet (like weather for detection).

