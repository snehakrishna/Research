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

load("sboot_200.RData")

maximum = max(sboot$upper)
minimum = min(sboot$lower)

sboot$ID<-seq.int(nrow(sboot))

sboot = sboot[order(-sboot$upper),]

par(mar = c(13,2,4.1,0))
x = 1:nPlants

newnames = vector()
plot(c(1, 109), range(maximum, minimum),  type = "n", 
	xaxt = 'n', xlab = "", ylab = "Scores", 
	main = "95% Confindence Intervals using Bootstapping for All Years")
for(pl in 1:nPlants){
	points(pl, sboot$upper[pl], col = "blue")
	points(pl, sboot$lower[pl], col = "red")
	abline(v = pl)
	newnames = append(newnames, plantNames[sboot$ID[pl]])
}
#axis(1, 1:nPlants, labels = plantNames[1:109], las = 2)

lablist<- newnames #as.vector(plantNames[1:nPlants])
axis(1, at=seq(1, nPlants, by=1), labels = FALSE)
text(seq(1, nPlants, by=1) - .5, par("usr")[3] - 9.5, labels = lablist, srt = 90, pos = 1, xpd = TRUE)

legend("bottomleft", inset=.05, title="Legend",
  	c("Upper Bound", "Lower Bound"), fill=c("blue", "red"), 
	horiz=FALSE, bg = "white")

