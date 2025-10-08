#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

###########################################################################
# PLOT SUPPLEMENTARY FIG.9
# USING ABSOLUTE SCALE WITH CITY-AVERAGE OF TEMPERATURES

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

###############################################################
#ABSOLUTE TEMPERATURE SCALE
#USE COUNTRY-AVERAGE TIMESERIES BY PERIOD AS TEMPERATURE BASIS
################################################################

# DETERMINE MINIMUM RISK TEMPERATURE (MRT) AND PERCENTILE (MRP)
# BY PERIOD AND BY LIFE EXPECTANCY LEVEL

mmtp <- mmtpcf1 <- mmtpcf3 <- mmtpcf4 <- mmtpcf5 <-vector("numeric",length=length(yearlist))
mmt <- mmtcf1 <- mmtcf3 <- mmtcf4 <- mmtcf5<- vector("numeric",length=length(yearlist))

mmtpcf2 <-  vector("numeric",length=2)
mmtcf2 <- vector("numeric",length=2)

mmtprange <- 25:99

# FIX KNOTS FOR ALL PERIODS
avgknotsabs <- quantile(tmeancountrysum$tmeancountry, varper / 100, na.rm = T)
setknots <- avgknotsabs

# LOOP OVER SUB-PERIODS

for (i in seq(yearlist)){
  
  # EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
  tdata <- subset(tmeancountrysum,year %in% yearlist[[i]])$tmeancountry
  predvar <- quantile(tdata,mmtprange/100,na.rm=T)

  # DEFINE TEMPERATURE BASIS
  argvar <- list(x = predvar, fun = varfun,knots=setknots)
  if (!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis, argvar)
  
  # DEFINE THE CENTERING POINT (AT POINT OF MINIMUM RISK)
  
  # SENSITIVITY ANALYSIS: USE MODEL 0 COEFFICIENTS TO DETERMINE MRP/MRT
  #mmtp[i] <- (mmtprange)[which.min(bvar%*%coef(model0))]
  #mmt[i] <- mmtcf1[i] <- mmtcf3[i] <- mmtcf4[i] <- quantile(tdata,mmtp[i]/100,na.rm=T)
  #if (i==length(yearlist)){ mmtcf2[1] <- mmtcf2[2] <-quantile(tdata,mmtp[i]/100,na.rm=T)}
  
  #PREDICT: PREDICTED COEFFICIENTS PER PERIOD
  #FACTUAL
  mmtp[i] <- (mmtprange)[which.min(bvar%*%pred[[i]]$fit)]
  mmt[i] <- quantile(tdata,mmtp[i]/100,na.rm=T)

  # COUNTERFACTUAL
  mmtpcf1[i] <- (mmtprange)[which.min(bvar%*%predcf1[[i]]$fit)]
  mmtcf1[i] <- quantile(tdata,mmtpcf1[i]/100,na.rm=T)

  mmtpcf3[i] <- (mmtprange)[which.min(bvar%*%predcf3[[i]]$fit)]
  mmtcf3[i] <- quantile(tdata,mmtpcf3[i]/100,na.rm=T)

  mmtpcf4[i] <- (mmtprange)[which.min(bvar%*%predcf4[[i]]$fit)]
  mmtcf4[i] <- quantile(tdata,mmtpcf4[i]/100,na.rm=T)
  
  mmtpcf5[i] <- (mmtprange)[which.min(bvar%*%predcf5[[i]]$fit)]
  mmtcf5[i] <- quantile(tdata,mmtpcf5[i]/100,na.rm=T)

  # PREDICT MMTs FOR HIGH/LOW LIFE EXPECTANCY
  if (i==length(yearlist)){
  for (j in 1:2){
  mmtpcf2[j] <- (mmtprange)[which.min(bvar%*%predcf2[[j]]$fit)]
  mmtcf2[j] <- quantile(tdata,mmtpcf2[j]/100,na.rm=T)
    }
  }
}

#####################################################
# PREDICT FIRST-PERIOD POOLED ASSOCIATION FOR PANEL A 
######################################################

# EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
tdata <- subset(tmeancountrysum,year %in% yearlist[[1]])$tmeancountry

# DEFINE TEMPERATURE BASIS
argvar <- list(
  x = tdata, fun = varfun,
  knots = setknots
)

if (!is.null(vardegree)) argvar$degree <- vardegree
bvar <- do.call(onebasis, argvar)

# DEFINE THE CENTERING POINT (AT POINT OF MINIMUM RISK)
cen <- mmt[1]

cp1 <- crosspred(bvar, coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                 model.link="log", by=0.1, cen=cen)

# FIXED TEMPERATURE (SECOND-PERIOD 99TH PERCENTILE - IS THE COOLEST PERIOD)
t99fix <- quantile(subset(tmeancountrysum,year %in% yearlist[[2]])$tmeancountry,0.99,na.rm=T)

####################################################
#  SAVE 99TH RR PER PERIOD/OVER TIME FOR PANEL B
###################################################

# PREDICT ASSOCIATIONS FOR OBSERVED META-PREDICTORS
rr99 <- t(sapply(seq(nrow(datapred)), function(i) {
  
  # EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
  tdata <- subset(tmeancountrysum,year %in% yearlist[[i]])$tmeancountry
  
  # DEFINE TEMPERATURE BASIS
  argvar <- list(
    x = tdata, fun = varfun,
    knots = setknots
  )
  if (!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis, argvar)
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=pred[[i]]$fit, vcov=pred[[i]]$vcov,
                  model.link="log", at=t99fix, cen=mmt[i])
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

# PREDICT ASSOCIATIONS FOR FIRST-PERIOD LIFE EXPECTANCY, AVERAGE ANNUAL TEMP, AND ALERT NUMBERS
rr99cf1 <- t(sapply(seq(nrow(datapredcf1)), function(i) {
  
  # EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
  tdata <- subset(tmeancountrysum,year %in% yearlist[[i]])$tmeancountry

  # DEFINE TEMPERATURE BASIS
  argvar <- list(
    x = tdata, fun = varfun,
    knots = setknots
  )
  if (!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis, argvar)
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf1[[i]]$fit, vcov=predcf1[[i]]$vcov,
                  model.link="log", at=t99fix, cen=mmtcf1[i])
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

rr99cf3 <- t(sapply(seq(nrow(datapredcf3)), function(i) {
  
  # EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
  tdata <- subset(tmeancountrysum,year %in% yearlist[[i]])$tmeancountry

  # DEFINE TEMPERATURE BASIS
  argvar <- list(
    x = tdata, fun = varfun,
    knots = setknots
  )
  if (!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis, argvar)
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf3[[i]]$fit, vcov=predcf3[[i]]$vcov,
                  model.link="log", at=t99fix, cen=mmtcf3[i])
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

rr99cf4 <- t(sapply(seq(nrow(datapredcf4)), function(i) {
  
  # EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
  tdata <- subset(tmeancountrysum,year %in% yearlist[[i]])$tmeancountry

  # DEFINE TEMPERATURE BASIS
  argvar <- list(
    x = tdata, fun = varfun,
    knots = setknots
  )
  if (!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis, argvar)
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf4[[i]]$fit, vcov=predcf4[[i]]$vcov,
                  model.link="log", at=t99fix, cen=mmtcf4[i])
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

rr99cf5 <- t(sapply(seq(nrow(datapredcf5)), function(i) {
  
  # EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
  tdata <- subset(tmeancountrysum,year %in% yearlist[[i]])$tmeancountry

  # DEFINE TEMPERATURE BASIS
  argvar <- list(
    x = tdata, fun = varfun,
    knots = setknots
  )
  if (!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis, argvar)
  
  # PREDICT ASSOCIATIONS AT 99TH PERCENTILE
  cp <- crosspred(bvar, coef=predcf5[[i]]$fit, vcov=predcf5[[i]]$vcov,
                  model.link="log", at=t99fix, cen=mmtcf5[i])
  
  # EXTRACT RR ANC CI
  est <- c(with(cp, c(allRRfit, allRRlow, allRRhigh)))
  names(est) <- c("RR","RRlow","RRhigh")
  
  # RETURN
  est
}))

#######################################################
# PREDICT ASSOCIATIONS FOR LOW AND HIGH LIFE EXPECTANCY
########################################################

# EXTRACT COUNTRY-AVERAGE TIME SERIES FOR LAST SUBPERIOD
tdata <- subset(tmeancountrysum,year %in% yearlist[[nsub]])$tmeancountry

# DEFINE TEMPERATURE BASIS
argvar <- list(
  x = tdata, fun = varfun,
  knots = setknots
)

if (!is.null(vardegree)) argvar$degree <- vardegree
bvar <- do.call(onebasis, argvar)

# DEFINE THE CENTERING POINT (AT POINT OF MINIMUM RISK)
cpcf2low <- crosspred(bvar, coef=predcf2[[1]]$fit, vcov=predcf2[[1]]$vcov,
                 model.link="log", by=0.1, cen=mmtcf2[1])
cpcf2high <- crosspred(bvar, coef=predcf2[[2]]$fit, vcov=predcf2[[2]]$vcov,
                      model.link="log", by=0.1, cen=mmtcf2[2])

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

pdf(paste0("FigS9_",fext,".pdf"),width=7,height=6)
layout(matrix(1:4,ncol=2,byrow = T))
par(mar=c(4, 4, 1, 1.5),mgp=c(2.5, 1, 0))

###############
# PLOT PANEL A
###############

plot(cp1, ylim=c(0.9,2.7), xlab="Temperature (°C)", ylab="RR",
     lab=c(6,5,7), las=1, lwd=1.5, mgp=c(2.5,1,0), cex.axis=0.8, col=cl[1],
     ci.arg=list(col=cllight[1]), main="",
     xlim=c(min(tmeancountrysum$tmeancountry),max(tmeancountrysum$tmeancountry)))

# LOOP OVER SECOND TO LAST PERIOD
for (i in seq(yearlist)[-1]){
  
  # EXTRACT PERIOD-SPECIFIC COUNTRY-AVERAGE TIME SERIES
  tdata <- subset(tmeancountrysum,year %in% yearlist[[i]])$tmeancountry

  # DEFINE TEMPERATURE BASIS
  argvar <- list(
    x = tdata, fun = varfun,
    knots = setknots
  )
  if (!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis, argvar)
  
  # DEFINE THE CENTERING POINT (AT POINT OF MINIMUM RISK)
  cen <- mmt[i]
  
  cp <- crosspred(bvar, coef=pred[[i]]$fit, vcov=pred[[i]]$vcov,
                  model.link="log", by=0.1, cen=cen)
  lines(cp, lwd=1.5, col=cl[i], ci="area", ci.arg=list(col=cllight[i]))
}

abline(v=t99fix, lty=2, col=grey(0.8))

legend("topleft", sapply(yearlist, function(x) paste(range(x), collapse="-")), lwd=1.5, col=cl, bty="n", inset=0.1)

box(lty=1,lwd=1)

# PLOT PANEL LABEL
text(min(tmeancountrysum$tmeancountry),2.7,labels="a")

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
     labels = period, cex = 0.8)

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
     labels = period, cex = 0.8)

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

plot(cpcf2low, ylim=c(0.8,4), xlab="Temperature (°C)", ylab="RR",
     lab=c(6,5,7), las=1, lwd=1.5, mgp=c(2.5,1,0), cex.axis=0.8, col="red",
     ci.arg=list(col=alpha("red",0.3)), main="",
     xlim=c(min(tmeancountrysum$tmeancountry),max(tmeancountrysum$tmeancountry)),cex.main=0.8)

lines(cpcf2high, lwd=1.5, col="black", ci="area", ci.arg=list(col=alpha("black",0.3)))

abline(v=t99fix, lty=2, col=grey(0.8))

legend("topleft", c(paste0("LE (",period[nsub],") = ",round(lifeexpectpred[nsub],digits=1)," years"),
                    paste0("LE (",period[1],") = ",round(lifeexpectpred[1],digits=1)," years")), 
       lwd=1.5, col=c("black","red"), bty="n",cex=0.9)

box(lty=1,lwd=1)


# PLOT PANEL LABEL
text(min(tmeancountrysum$tmeancountry),4,labels="d")

dev.off()

####