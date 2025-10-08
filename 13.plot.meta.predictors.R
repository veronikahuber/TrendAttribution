#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

###########################################################
# PLOT SELECTED META-PREDICTORS BY SUBPERIODS
# SUPPLEMENTARY FIG. 3
###########################################################

periodlabel <- sapply(yearlist, function(x) paste(range(x), collapse="-"))

pdf("FigS3.pdf",width=9,height=3)
layout(t(1:3))
par(mar=c(4, 4, 2, 1.5),mgp=c(2.5, 1, 0))

####################################
# PLOT AVERAGE ANNUAL TEMPERATURES
####################################

plot(1:length(yearlist),cityinfo$average.annual.temp[cityinfo$city==cities[1]],type="n",
     ylab= "Average annual temperature (°C)",xlab="",xaxt="n",main="",
     ylim=c(min(cityinfo$average.annual.temp)-1,max(cityinfo$average.annual.temp)+1))

# PLOTXAXIS LABELS
axis(side=1,at=1:length(yearlist),label=NA)
text(1:length(yearlist), rep(min(cityinfo$average.annual.temp)-1.6,length(yearlist)), 
     srt = 35, xpd = TRUE,adj=c(1,0),
     labels = periodlabel, cex = 0.8)

for (i in seq(cities)){
  lines(1:length(yearlist),cityinfo$average.annual.temp[cityinfo$city==cities[i]],col="grey")
}

# PLOT CITY AVERAGE
lines(1:length(yearlist),avgtemppred,col="black",lwd=2)

# PLOT FIGURE LABEL
text(1,max(cityinfo$average.annual.temp)+1,labels="a")

legend("top",c("City-specific","City average"),lty=1,col=c("grey","black"),bty="n")

#############################################
# PLOT TEMPORAL VARIABILITY HEAT ALERT DAYS
#############################################

plot(1:length(yearlist),cityinfo$sum.heat.days[cityinfo$city==cities[1]]/nyper,type="n",
     ylab="Average annual number of heat alerts (days)",xlab="",xaxt="n",
     ylim=c(min(cityinfo$sum.heat.days/nyper)-2,max(cityinfo$sum.heat.days/nyper)+1))

# PLOTXAXIS LABELS
axis(side=1,at=1:length(yearlist),label=NA)
text(1:length(yearlist),rep(min(cityinfo$sum.heat.days/nyper)-4,length(yearlist)), 
     srt = 35, xpd = TRUE,adj=c(1,0),
     labels = periodlabel, cex = 0.8)

for (i in seq(cities)){
  lines(1:length(yearlist),cityinfo$sum.heat.days[cityinfo$city==cities[i]]/nyper,col="grey")
}

# PLOT CITY AVERAGE
lines(1:length(yearlist),alertdaypred/nyper,col="black",lwd=2)

text(1,max(cityinfo$sum.heat.days/nyper)+1,labels="b")

################################
# PLOT AVERAGE POPULATION AGE
################################

plot(1:length(yearlist),cityinfo$average.population.age[cityinfo$city==cities[1]],type="n",
     ylab="Average population age (years)",xlab="",xaxt="n",
     ylim=c(min(cityinfo$average.population.age)-2,max(cityinfo$average.population.age)+1))

# PLOTXAXIS LABELS
axis(side=1,at=1:length(yearlist),label=NA)
text(1:length(yearlist),rep(min(cityinfo$average.population.age)-3,length(yearlist)), 
     srt = 35, xpd = TRUE,adj=c(1,0),
     labels = periodlabel, cex = 0.8)

for (i in seq(cities)){
  lines(1:length(yearlist),cityinfo$average.population.age[cityinfo$city==cities[i]],col="grey")
}

# PLOT CITY AVERAGE
lines(1:length(yearlist),avgagepred,col="black",lwd=2)

text(1,max(cityinfo$average.population.age)+1,labels="c")

dev.off()

####
