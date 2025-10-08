#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

###########################################################################
# PLOT FIG.2
# USING RELATIVE SCALE WITH AVERAGE PERCENTILE TEMPERATURE

# PANEL A: PERIOD-SPECIFIC POOLED EXPOSURE-RESPONSE FUNCTIONS
# PANEL B: LIFE EXPECTANCY BY PERIOD AND CITY
# PANEL C: RR TREND OVER TIME (FIXING META-PREDICTORS AT FIRST PERIOD VALUE)
# PANEL D: EXPOSURE-REPONSE FUNCTION FOR LOW AND HIGH LIFE-EXPECTANCIES
##########################################################################

# NUMBER OF SUBPERIODS
nsub <- length(yearlist)

##################################################
# PREDICT POOLED ASSOCIATIONS
##################################################

# DEFINE THE META-PREDICTOR VALUES (AVERAGES ACROSS CITIES)
yearpred <- sapply(yearlist, mean)
lifeexpectpred <- tapply(cityinfo$life.expectancy,cityinfo$period,mean)
alertdaypred <- tapply(cityinfo$sum.heat.days,cityinfo$period,mean)
avgtemppred <- tapply(cityinfo$average.annual.temp,cityinfo$period,mean)
avgagepred <- tapply(cityinfo$average.population.age,cityinfo$period,mean)

# ASSEMBLE META-PREDICTORS
datapred <- data.frame(year=yearpred,
                       average.annual.temp=avgtemppred,
                       life.expectancy=lifeexpectpred,
                       sum.heat.days=alertdaypred,
                       average.population.age=avgagepred)

# PREDICT COEF/VCOV
pred <- predict(modelfinal, datapred, vcov=T)

##################################################
# PREDICT BASED ON COUNTERFACTUAL META-PREDICTORS
##################################################

# PREDICT ASSOCIATION BY SUBPERIOD WITH LIFE-EXPECTANCY SET TO VALUE OF FIRST PERIOD
lifeexpectpredcf <- rep(lifeexpectpred[1],length(yearlist))
datapredcf1 <- data.frame(year=yearpred,
                          average.annual.temp=avgtemppred,
                          sum.heat.days=alertdaypred,
                          life.expectancy=lifeexpectpredcf,
                          average.population.age=avgagepred)

predcf1 <- predict(modelfinal, datapredcf1, vcov=T)

# PREDICT ASSOCIATION FOR LAST SUBPERIOD FOR LOW AND HIGH LIFE-EXPECTANCY
datapredcf2 <- data.frame(year=c(yearpred[nsub],yearpred[nsub]),
                        average.annual.temp=c(avgtemppred[nsub],avgtemppred[nsub]),
                        sum.heat.days=c(alertdaypred[nsub],alertdaypred[nsub]),
                        life.expectancy=c(lifeexpectpred[1],lifeexpectpred[nsub]),
                        average.population.age=c(avgagepred[nsub],avgagepred[nsub]))

predcf2 <- predict(modelfinal, datapredcf2, vcov=T)
names(predcf2) <- c("low LE","high LE")

# PREDICT ASSOCIATION BY SUBPERIOD WITH AVERAGE ANNUAL TEMPERATURES SET TO VALUE OF FIRST PERIOD
avgtemppredcf <- rep(avgtemppred[1],length(yearlist))
datapredcf3 <- data.frame(year=yearpred,
                          average.annual.temp=avgtemppredcf,
                          sum.heat.days=alertdaypred,
                          life.expectancy=lifeexpectpred,
                          average.population.age=avgagepred)

predcf3 <- predict(modelfinal, datapredcf3, vcov=T)

# PREDICT ASSOCIATION BY SUBPERIOD WITH NUMBER OF HEAT ALERTS SET TO VALUE OF FIRST PERIOD
alertdaypredcf <- rep(alertdaypred[1],length(yearlist))
datapredcf4 <- data.frame(year=yearpred,
                          average.annual.temp=avgtemppred,
                          sum.heat.days=alertdaypredcf,
                          life.expectancy=lifeexpectpred,
                          average.population.age=avgagepred)

predcf4 <- predict(modelfinal, datapredcf4, vcov=T)

# PREDICT ASSOCIATION BY SUBPERIOD WITH AVERAGE POPULATION AGE SET TO VALUE OF FIRST PERIOD
avgagepredcf <- rep(avgagepred[1],length(yearlist))
datapredcf5 <- data.frame(year=yearpred,
                          average.annual.temp=avgtemppred,
                          sum.heat.days=alertdaypred,
                          life.expectancy=lifeexpectpred,
                          average.population.age=avgagepredcf)

predcf5 <- predict(modelfinal, datapredcf5, vcov=T)

#####################################################
# PREDICT FIRST-PERIOD POOLED ASSOCIATION FOR PANEL A 
######################################################

# AVERAGE TMEAN AT 99TH PERCENTILE
tmean99 <- avgtmeansum$tmean[avgtmeansum$perc=="99.0%"]

# DEFINE SPLINE TRANSFORMATION ORIGINALLY USED IN FIRST-STAGE MODELS
setknots <- avgtmeansum$tmean[avgtmeansum$perc %in% paste0(varper, ".0%")]

if(!is.null(vardegree)) {
  bvar <- onebasis(avgtmeansum$tmean, fun=varfun, degree=vardegree, knots=setknots)
}else{
  bvar <- onebasis(avgtmeansum$tmean, fun=varfun, knots=setknots)
}

# DEFINE THE CENTERING POINT 
cen <- avgtmeansum$tmean[avgtmeansum$perc=="50.0%"]

# PLOTTING LABELS
xperc <- c(0,1,5,25,50,75,90,99,100)
xval <- avgtmeansum$tmean[avgtmeansum$perc %in% paste0(xperc, ".0%")]

cp1 <- crosspred(bvar, coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                 model.link="log", at=avgtmeansum$tmean, cen=cen)

####################################################
#  SAVE 99TH RR PER PERIOD/OVER TIME FOR PANEL B
###################################################

# PREDICT ASSOCIATIONS FOR OBSERVED META-PREDICTORS
rr99 <- t(sapply(seq(nrow(datapred)), function(i) {
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=pred[[i]]$fit, vcov=pred[[i]]$vcov,
                  model.link="log", at=tmean99, cen=cen)
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

# PREDICT ASSOCIATIONS FOR FIRST-PERIOD LIFE EXPECTANCY, AVERAGE ANNUAL TEMP, AND ALERT NUMBERS
rr99cf1 <- t(sapply(seq(nrow(datapredcf1)), function(i) {
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf1[[i]]$fit, vcov=predcf1[[i]]$vcov,
                  model.link="log", at=tmean99, cen=cen)
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

rr99cf3 <- t(sapply(seq(nrow(datapredcf3)), function(i) {
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf3[[i]]$fit, vcov=predcf3[[i]]$vcov,
                  model.link="log", at=tmean99, cen=cen)
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

rr99cf4 <- t(sapply(seq(nrow(datapredcf4)), function(i) {

  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf4[[i]]$fit, vcov=predcf4[[i]]$vcov,
                  model.link="log", at=tmean99, cen=cen)
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

rr99cf5 <- t(sapply(seq(nrow(datapredcf5)), function(i) {
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf5[[i]]$fit, vcov=predcf5[[i]]$vcov,
                  model.link="log", at=tmean99, cen=cen)
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

#######################################################
# PREDICT ASSOCIATIONS FOR LOW AND HIGH LIFE EXPECTANCY
########################################################

# DEFINE THE CENTERING POINT (AT POINT OF MINIMUM RISK)
cpcf2low <- crosspred(bvar, coef=predcf2[[1]]$fit, vcov=predcf2[[1]]$vcov,
                 model.link="log", at=avgtmeansum$tmean, cen=cen)
cpcf2high <- crosspred(bvar, coef=predcf2[[2]]$fit, vcov=predcf2[[2]]$vcov,
                      model.link="log", at=avgtmeansum$tmean, cen=cen)

####################
# INITIATE PLOT
####################

# COLORS FOR PANEL A
library(RColorBrewer); library(grDevices)
cl <- hcl.colors(length(yearlist),palette="roma" ,rev = T)

cllight <- vector("character",length(cl))
for (i in seq(cl)){
  cllight[i] <- do.call(rgb,c(as.list(col2rgb(cl[i])),alpha=255/6,max=255))
}

pdf(paste0("Fig2_rel_",fext,".pdf"),width=7,height=6)
layout(matrix(1:4,ncol=2,byrow = T))
par(mar=c(4, 4, 1, 1.5),mgp=c(2.5, 1, 0))

###############
# PLOT PANEL A
###############

plot(cp1, ylim=c(0.9,2.6), xlab="Temperature percentile", ylab="RR",
     lab=c(6,5,7), las=1, lwd=1.5, xaxt="n", mgp=c(2.5,1,0), 
     cex.axis=0.8, col=cl[1],
     ci.arg=list(col=cllight[1]), main="")
axis(1, at=xval, labels=paste0(xperc, "%"), cex.axis=0.9)

# LOOP OVER SECOND TO LAST PERIOD
for (i in seq(yearlist)[-1]){

  cp <- crosspred(bvar, coef=pred[[i]]$fit, vcov=pred[[i]]$vcov,
                  model.link="log", at=avgtmeansum$tmean, cen=cen)
  lines(cp, lwd=1.5, col=cl[i], ci="area", ci.arg=list(col=cllight[i]))
}

abline(v=tmean99, lty=2, col=grey(0.8))

legend("topleft", sapply(yearlist, function(x) paste(range(x), collapse="-")), lwd=1.5, col=cl, bty="n", inset=0.1)

box(lty=1,lwd=1)

# PLOT PANEL LABEL
text(xval[1],2.6,labels="a")

#######################
# PANEL B
#######################

# PLOT LIFE EXPECTANCY
plot(1:length(yearlist),cityinfo$life.expectancy[cityinfo$city==cities[1]],type="n",
     ylab= "Average LE (years)",xlab="",xaxt="n",main="",
     ylim=c(min(cityinfo$life.expectancy)-1,max(cityinfo$life.expectancy)+1))

# PLOTXAXIS LABELS
axis(side=1,at=1:length(yearlist),label=NA)
text(1:length(yearlist), rep(min(cityinfo$life.expectancy)-2.5,length(yearlist)), 
     srt = 35, xpd = TRUE,adj=c(1,0),
     labels = period, cex = 0.9)

for (i in seq(cities)){
  lines(1:length(yearlist),cityinfo$life.expectancy[cityinfo$city==cities[i]],col="grey")
}

# PLOT CITY AVERAGE
lines(1:length(yearlist),lifeexpectpred,col="black",lwd=2)

# PLOT PANEL LABEL
text(1,max(cityinfo$life.expectancy)+1,labels="b")

##################
# PLOT PANEL C
##################

ylab <-expression(RR~at~99^th~percentile)
ylim <- range(c(range(rr99),range(rr99cf1),range(rr99cf3),range(rr99cf4),range(rr99cf5)))*c(0.8,1.0)

# PLOT
plot(yearpred, seq(yearpred), type="n", ylim=ylim, 
     ylab=ylab, xlab="", xaxt="n",bty="l", las=1, mgp=c(2.5,1,0), cex.axis=0.8,
     main="")

# PLOTXAXIS LABELS
axis(1,at=yearpred,labels=NA)

text(yearpred, par("usr")[3]-0.15, 
     srt = 35, adj = 1, xpd = TRUE,
     labels = period, cex = 0.9)

arrows(yearpred, rr99[,2], yearpred, rr99[,3], col="grey",
       code=3, angle=90, length=0.05, lwd=1)
points(yearpred, rr99[,1],  type="o", col="black", pch=19)

# COUNTERFACTUAL LIFE-EXPECTANCY
arrows(yearpred, rr99cf1[,2], yearpred, rr99cf1[,3], col="pink",
       code=3, angle=90, length=0.05, lwd=1)
points(yearpred, rr99cf1[,1],  type="o", col="red", pch=19)

# # COUNTERFACTUAL AVG ANNUAL TEMPERATURES
lines(yearpred, rr99cf3[,1], col="black", lty=2)

# # COUNTERFACTUAL NUMBER OF HEAT ALERTS
lines(yearpred, rr99cf4[,1], col="black", lty=3)

# # COUNTERFACTUAL AVERAGE POPULATION HEALTH
lines(yearpred, rr99cf5[,1], col="black", lty=4)

legend("bottomleft",
       c("Observed","Constant LE",
         "Constant annual mean temperatures",
         "Constant number of heat alerts",
         "Constant average population age"),
       col=c("black","red","black","black","black"),pch=c(19,19,NA,NA,NA), 
       lty=c(1,1,2,3,4),bty="n",cex=0.7)

box(lty=1,lwd=1)

# PLOT PANEL LABEL
text(yearpred[1],ylim[2],labels="c")

#######################
# PANEL D
#######################

plot(cpcf2low, ylim=c(0.8,4.2), xlab="Temperature percentile", ylab="RR",
     lab=c(6,5,7), las=1, lwd=1.5, mgp=c(2.5,1,0), xaxt="n",cex.axis=0.8, col="red",
     ci.arg=list(col=alpha("red",0.3)), main="",cex.main=0.8)
axis(1, at=xval, labels=paste0(xperc, "%"), cex.axis=0.9)

lines(cpcf2high, lwd=1.5, col="black", ci="area", ci.arg=list(col=alpha("black",0.3)))

abline(v=tmean99, lty=2, col=grey(0.8))

legend("top", c(paste0("LE (",period[nsub],") = ",round(lifeexpectpred[nsub],digits=1)," years"),
                    paste0("LE (",period[1],") = ",round(lifeexpectpred[1],digits=1)," years")), 
       lwd=1.5, col=c("black","red"), bty="n",cex=0.9)

box(lty=1,lwd=1)

# PLOT PANEL LABEL
text(xval[1],4.2,labels="d")


dev.off()

####