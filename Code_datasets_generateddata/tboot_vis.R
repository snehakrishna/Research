setwd(path.expand("~//Research//thesis_stuff"))

allData_ori = read.csv("Interactions_Final.csv")
anthesis_ori = read.csv("Anthesis_Final.csv")

#for all meadow watch years
plantNames = levels(allData_ori$PLTSP_NAME)
plantNames = plantNames[-c(1)]
pollNames = levels(allData_ori$VISSP_NAME) # pollinators
pollNames = pollNames[-c(1)]
nPlants = length(plantNames)
nPolls = length(pollNames)

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

testmatrix = model.matrix(~ traits$ExclusionClosed
	+ traits$ExclusionPendant + traits$EXVisible_Sneha
	+ traits$generalShapeExclusion + traits$LifeForm_Sneha
	+ traits$FlowerForm_Sneha + traits$EXDiel_Sneha
	+ traits$EXHandle_Sneha + traits$ExclusionFeeble
	+ traits$ExclusionPlatform)[,-1]

traits$intercept = 1

tmatrix = as.matrix(data.frame(traits$intercept, traits$Biom.flwr, testmatrix))
ntraits = dim(tmatrix)[2]

load("tboot_200.RData")

sboot$ID<-seq.int(nrow(sboot))

sboot = sboot[order(-sboot$upper),]

par(mar = c(13,2,4.1,1))

maximum = max(sboot$upper)
minimum = min(sboot$lower)

newnames = c("Intercept", "Biomass",
	"Closed (open access)",
	"Pendant (suspended)",
	"Visibility (not bright)",
	"Tube Shape (moderate exclusion)",
	"Tube Shape (poor exclusion)",
	"Tube Shape (severe exclusion)",
	"Life Form (Perennial)",
	"Flower Form (exclusion exists)",
	"Diel (yes night)",
	"Pollen Size (average)",
	"Feebleness (strong)",
	"Platform (platform)",
	"Platform (weak)")

plot(c(1, 15), range(maximum, minimum),  type = "n", 
	xaxt = 'n', xlab = "", ylab = "Scores", 
	main = "95% Confindence Intervals using Bootstapping for All Years")
for(pl in 1:ntraits){
	points(pl, sboot$upper[pl], col = "blue")
	points(pl, sboot$lower[pl], col = "red")
	abline(v = pl)
}
#axis(1, 1:ntraits, labels = plantNames[1:109], las = 2)

lablist<- newnames 
axis(1, at=seq(1, ntraits, by=1), labels = FALSE)
text(seq(1, ntraits, by=1), par("usr")[3] - 3, labels = lablist, srt = 90, pos = 1, xpd = TRUE)

legend("bottomleft", inset=.05, title="Legend",
  	c("Upper Bound", "Lower Bound"), fill=c("blue", "red"), 
	horiz=FALSE, bg = "white")

