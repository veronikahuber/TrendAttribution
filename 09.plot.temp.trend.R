#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

##############################################################################
# 1) PLOT ANNUAL MEAN SUMMER TEMPERATURES PER CITY
# AND FIT LINEAR REGRESSION LINES (SUPPLEMENTARY FIG. 1)

# 2) PLOT COUNTRY-AVERAGE TEMPERATURES AND LINEAR TRENDS (FIG. 1a)

# 3) PLOT CITY-SPECIFIC LINEAR REGRESSION SLOPES (FIG. 1b)
##############################################################################

# DEFINE MONTHS TO AVERAGE TEMPERATURES

#DEFAULT
mchoice <- sm

# INITIATE ARRAYS TO STORE RESULTS OF LINEAR REGRESSIONS
# AND COUNTERFACTUAL WARM-SEASON TEMPERATURE VERSIONS
ncities <- length(cities)
temptrend <- temptrendp <- array(NA,dim=c(ncities+1,10),
                                 dimnames=list(c(cities,"Germany"),c("fact","cfact.best",paste0("cfact.ci.",1:8))))
temptrendci <- array(NA,dim=c(ncities+1,10,2),
                                 dimnames=list(c(cities,"Germany"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),c("ci.l","ci.u")))

mst <- array(NA,dim=c(ncities,length(ny),10),dimnames=list(cities,ny,c("fact","cfact.best",paste0("cfact.ci.",1:8))))

# DEFINE PLOTTING PARAMETERS
colors <- c("#1f77b4","#ff7f0e")
colorslight <- vector("character",2)

for (i in seq(colors)){
  colorslight[i] <- do.call(rgb,c(as.list(col2rgb(colors[i])),alpha=255/6,max=255))
}

colorsextralight <- do.call(rgb,c(as.list(col2rgb(colors[2])),alpha=255/15,max=255))

# COMPUTE FACTUAL AND COUNTERFACTUAL MST PER YEAR 
mst[,,1] <- t(sapply(dlist,function(x) tapply(x$tmeanf[x$year %in% ny & x$month %in% mchoice],x$year[x$year %in% ny & x$month %in% mchoice],mean,na.rm=T)))

# LOOP OVER COUNTERFACTUAL DATA VERSIONS
for (i in seq(tmeancfvar)){
  mst[,,i+1] <- t(sapply(dlist,function(x) tapply(x[x$year %in% ny & x$month %in% mchoice,tmeancfvar[i]],x$year[x$year %in% ny & x$month %in% mchoice],mean,na.rm=T)))
}

# MEAN WARMING
mstwarm <- sapply(1:9, function(x) rowMeans(mst[,,1],na.rm=T)-rowMeans(mst[,,x+1],na.rm=T))
dimnames(mstwarm)[[2]] <- c("cfact.best",paste0("cfact.ci.",1:8))
  
# COUNTRY-AVERAGE TEMPERATURE
mstfcountry <- colMeans(mst[,,1])
mstcfcountry <- apply(mst[,,2:10],2:3,mean,na.rm=T)
mstwarmcountry <- mean(mstfcountry)-colMeans(mstcfcountry)

# LABELS FOR RESULTS OF WILCOXON RANK TESTS
wrtlabel <- function(p){
  if (p < 0.001) {
    labels <- "***"
  } else if (p < 0.01){
    labels <- "**"
  } else if (p < 0.05){
    labels <- "*"
  } else {
    labels <- "ns"
  }
}

############################
# PLOT TEMPERATURE BY CITY
#############################

xlab <- "Years"
ylab1 <- "Temperature (°C)"

plotname <- paste0("FigS1_",fext,".pdf")
pdf(plotname,width=9,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)

# LOOP OVER CITIES
for(i in seq(length(dlist))) {
  
  # INITIATE PLOT
  plot(ny,mst[i,,1],main=cities[i],type="n",bty="l",xlab=xlab,ylab=ylab1,
       ylim=c(min(mst[i,,],na.rm=T),max(mst[i,,],na.rm=T)+1))
  
  # PLOT UNCERTAINTY BAND FOR COUNTERFACTUAL TRENDS
  
  for (j in 3:10){
    
    temp <- mst[i,,j]

    #  LINEAR FIT LOWEST AIC AND BIC
    fitlm <- lm(temp ~ ny)
    bestfit <- fitlm
    
    # SAVE LINEAR SLOPE, 95% CI, AND P-VALUE
    temptrend[i,j] <- bestfit$coefficients[2]
    temptrendci[i,j,] <- confint(bestfit,'ny',level=0.95)
    temptrendp[i,j] <- summary(bestfit)$coefficients[2,4]
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(bestfit, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))

    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorsextralight,border=NA)
    
  }
  
  #LOOP OVER TEMPERATURE VERSIONS (FACT, CFACT)
  
  for (j in 1:2){
    
    temp <- mst[i,,j]

    # DERIVE MIN,MAX RANGE FOR COUNTERFACTUAL ESTIMATES
    arrows(ny, apply(mst[i,,2:10],1,min), ny, apply(mst[i,,2:10],1,max), col=colorslight[j],
           code=3, angle=90, length=0.02, lwd=0.7)
    
    points(ny,temp,col=colors[j],pch=17)
    
    # TRY FITTING NON-LINEAR REGRESSION LINE 
    # NB: LINEAR FIT LOWEST AIC AND BIC
    fitlm <- lm(temp ~ ny)
    #fitdf2 <- lm(temp ~ ns(ny,df=2))
    #fitdf3 <- lm(temp ~ ns(ny,df=3))
    #fitdf4 <- lm(temp ~ ns(ny,df=4))
    #AIC(fitlm,fitdf2,fitdf3,fitdf4)
    #BIC(fitlm,fitdf2,fitdf3,fitdf4)
    
    bestfit <- fitlm

    # SAVE LINEAR SLOPE, 95% CI, AND P-VALUE
    temptrend[i,j] <- bestfit$coefficients[2]
    temptrendci[i,j,] <- confint(bestfit,'ny',level=0.95)
    temptrendp[i,j] <- summary(bestfit)$coefficients[2,4]
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(bestfit, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))
    # PLOT SPLINE TREND
    lines(ny,bestfit$fitted.values,col=colors[j],xpd=F,lty=1)
    
    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorslight[j],border=NA)
    
  }
  
  if (i==1) {legend("topleft",c("with CC","w/o CC"),col=colors,lty=1,pch=17,cex=0.8,bty="n")}
}

dev.off()

######################################
# PLOT COUNTRY AVERAGE TEMPERATURE
#######################################

pdf(paste0("Fig1_",fext,".pdf"),width=5,height=2.5)
layout(matrix(1:2,ncol=2,byrow = T))
par(mar=c(3, 4, 2, 1) + 0.1, mgp=c(2.5, 1, 0),cex=0.75)

plot(ny,mstfcountry,type="n",bty="l",xlab="",ylab=ylab1,
     ylim=c(min(mstcfcountry,na.rm=T),max(mstfcountry,na.rm=T)+1))
box("plot",lty=1,col="black")

# PLOT UNCERTAINTY BANDS FOR COUNTERFACTUAL TRENDS

for (j in 2:9){
  
  temp <- mstcfcountry[,j]
  
  #  LINEAR FIT LOWEST AIC AND BIC
  fitlm <- lm(temp ~ ny)
  bestfit <- fitlm
  
  # SAVE LINEAR SLOPE, 95% CI, AND P-VALUE
  temptrend["Germany",j+1] <- bestfit$coefficients[2]
  temptrendci["Germany",j+1,] <- confint(bestfit,'ny',level=0.95)
  temptrendp["Germany",j+1] <- summary(bestfit)$coefficients[2,4]
  
  # COMPUTE ERROR BANDS (2*SE)
  pred = predict(bestfit, newdata = list(ny), se = TRUE)
  se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                               "lower" = fit-2*se.fit))
  
  # PLOT UNCERTAINTY BAND
  polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorsextralight,border=NA)
  
}

# PLOT FACTUAL AND BEST-ESTIMATE COUNTERFACTUAL TEMPERATURES

for (j in c(1,2)){
  
  if (j==1) {temp<-mstfcountry}else{temp<-mstcfcountry[,1]}
  
  # DERIVE MIN,MAX RANGE FOR COUNTERFACTUAL ESTIMATES
  arrows(ny, apply(mstcfcountry,1,min), ny, apply(mstcfcountry,1,max), col=colorslight[j],
         code=3, angle=90, length=0.02, lwd=0.7)
  
  points(ny,temp,col=colors[j],pch=17)
  
  # TRY FITTING NON-LINEAR REGRESSION LINE 
  # NB: LINEAR FIT LOWEST AIC AND BIC
  fitlm <- lm(temp ~ ny)
  #fitdf2 <- lm(temp ~ ns(ny,df=2))
  #fitdf3 <- lm(temp ~ ns(ny,df=3))
  #fitdf4 <- lm(temp ~ ns(ny,df=4))
  #AIC(fitlm,fitdf2,fitdf3,fitdf4)
  #BIC(fitlm,fitdf2,fitdf3,fitdf4)
  
  bestfit <- fitlm
  
  # SAVE LINEAR SLOPE, 95% CI, AND P-VALUE
  temptrend["Germany",j] <- bestfit$coefficients[2]
  temptrendci["Germany",j,] <- confint(bestfit,'ny',level=0.95)
  temptrendp["Germany",j] <- summary(bestfit)$coefficients[2,4]
  
  # COMPUTE ERROR BANDS (2*SE)
  pred = predict(bestfit, newdata = list(ny), se = TRUE)
  se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                               "lower" = fit-2*se.fit))
  # PLOT SPLINE TREND
  lines(ny,bestfit$fitted.values,col=colors[j],xpd=F,lty=1)
  
  # PLOT UNCERTAINTY BAND
  polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorslight[j],border=NA)
  
}

legend("top",c("with CC","w/o CC"),col=colors,lty=1,pch=17, cex=0.8,bty="n")

# PLOT PANEL LABEL
text(ny[1],max(mstfcountry,na.rm=T)+1,labels="a")

###############################################
# PLOT LINEAR REGRESSION SLOPES
###############################################

# COMPUTE WILCOXON-RANK-TEST ON CITY-SPECIFIC TRENDS
# TO DETERMINE WHETHER CLIMATE CHANGE INDUCED A SIGNIFICANT EFFECT

wrt <- wilcox.test(temptrend[1:ncities,1],as.vector(temptrend[1:ncities,2:10]),alternative="two.sided")
labels <- wrtlabel(wrt$p.value)

ylab <- expression(paste(Delta," Temperature ","(",degree,"C/year)"))
ylim <- c(0,0.1)

plot(rep(1,ncities),temptrend[1:ncities,1],type="n",xaxt="n",
     ylab=ylab,xlab="",xlim=c(0,3),ylim=ylim,main="")
axis(side=1,at=c(1,2),labels=c("with CC","w/o CC"),cex.axis=1)

# PLOT FACTUAL TEMPERATURE TRENDS
boxplot(temptrend[1:ncities,1],add=T,at=1,border="black",col=colorslight[1],
        xaxt="n",yaxt="n",ylim=ylim,outpch=18)
points(rep(1,ncities),temptrend[1:ncities,1],col=colors[1],pch=2)
points(1,temptrend[ncities+1,1],pch=4,col="black",cex=1.7)

# PLOT BOXPLOT FOR ALL COUNTERFACTUAL TEMPERATURES
boxplot(as.vector(temptrend[1:ncities,2:10]),add=T,at=2,border="black",col=colorslight[2],
        xaxt="n",yaxt="n",ylim=ylim,outpch= 18)

# ADD BEST-ESTIMATE COUNTERFACTUAL TEMPERATURES
points(rep(2,ncities),temptrend[1:ncities,2],col=colors[2],pch=2)
points(2,temptrend[ncities+1,2],pch=4,col="black",cex=1.7)

# ADD SIGNIFICANCE LEVEL OF DIFFERENCES BETWEEN FACTUAL AND COUNTERFACTUAL TRENDS
yline <- max(temptrend[1:ncities,1])+0.01
lines(c(1,2),c(yline,yline),lty=1,col="black")
lines(c(1,1),c(yline,yline-0.002),lty=1,col="black")
lines(c(2,2),c(yline,yline-0.002),lty=1,col="black")
text(1.5,yline+0.005,labels=labels,cex=1.2)

# PLOT PANEL LABEL
text(0,ylim[2],labels="b")

dev.off()

########################
