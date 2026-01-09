#####################################################################################
# R code for the analysis in 

# Huber, V., Breitner-Busch, S., Feldbusch, H. et al. Improvements in life expectancy 
# mask rising trends in heat-related excess mortality attributable to climate change. 
# Nat Commun 16, 11632 (2025). https://doi.org/10.1038/s41467-025-66681-0
######################################################################################

#####################################################################
# PLOT ATTRIBUTABLE NUMBERS AND DEATH RATES (PER 100 00O POPULATION)
# SUPPLEMENTARY FIG 6
#####################################################################

##################################################################################
# PLOT POOLED AN (DEATH COUNT AND MORTALITY RATES) WITH/WITHOUT LE IMPROVEMENTS
##################################################################################

ylab <- "AN"

# DEFINE POINT TYPES
pch <- c(20,18)

ylim <- c(0,max(c(max(heatantotnoadapt[,basemort,"fact","ci.u"],na.rm=T),max(heatantotvar[,basemort,"fact","ci.u"],na.rm=T))))
  
plotname <- paste0("FigS6_",fext,".pdf")
pdf(plotname,width=5,height=4)
layout(matrix(1:4,ncol=2,byrow = T))
par(mar=c(2, 4.1, 2, 1.5), mgp=c(2.5, 1, 0),cex=0.75)

#########################
#ATTRIBUTABLE NUMBER
#########################

# LOOP OVER SETUP WITH/WITHOUT LE IMPROVEMENTS
for (i in (c(1,2))){
  
  if (i==1){heatantot <- heatantotvar} 
  else {heatantot <- heatantotnoadapt}
  
  # INITIATE PLOT
  plot(ny,heatantot[,basemort,"fact","est"],main=title[i],cex.main=1,
       type="n",bty="l",xlab="",ylab=ylab,
       ylim=ylim)
  box("plot",lty=1,"black")
  
  # LOOP OVER SETUP WITH/WITHOUT CLIMATE CHANGE
  
  for (j in c(1,2)){
    
    arrows(ny, heatantot[,basemort,j,"ci.l"], ny, heatantot[,basemort,j,"ci.u"], col=colorslight[j],
           code=3, angle=90, length=0.02, lwd=1)
    points(ny,heatantot[,basemort,j,"est"],col=colors[j],pch=pch[i])

    # LINEAR REGRESSION LINE 
    # NB: LINEAR FIT LOWEST AIC AND BIC
    fitlm <- lm(heatantot[,basemort,j,"est"] ~ ny)

    bestfit <- fitlm
    
    # SAVE LINEAR SLOPE
    antrend["Germany",j,i] <- bestfit$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    antrendci["Germany",j,i,] <- confint(bestfit,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    antrendp["Germany",j,i] <- summary(bestfit)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(bestfit, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                "lower" = fit-2*se.fit))
    # PLOT TREND LINE
    lines(ny,pred$fit,col=colors[j],xpd=F,lty=1)
    
    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorslight[j],border=NA)

  }
  
  legend("topright",c("with CC","w/o CC"),col=colors,lty=1,pch=pch[i],cex=0.8,bty="n")
  
  # PLOT PANEL LABEL
  text(ny[1],ylim[2],labels=panellabels[i])
}

###############################
# ATTRIBUTABLE MORTALITY RATES
###############################

# DIVIDE AN BY ANNUAL POPULATION

# SUM POPULATION DATA ACROSS CITIES
popsum <- sapply(ny,function(x) sum(metadat$population.corrected[metadat$year==x]))

heatantotvarrate <- apply(heatantotvar,2:4, function(x) x*100000/popsum)
heatantotnoadaptrate <- apply(heatantotnoadapt,2:4, function(x) x*100000/popsum)

ylim <- c(0,max(c(max(heatantotnoadaptrate[,basemort,"fact","ci.u"],na.rm=T),max(heatantotvarrate[,basemort,"fact","ci.u"],na.rm=T))))
ylab <- "AN (per 100 000 population)"

# LOOP OVER SETUP WITH/WITHOUT LE IMPROVEMENTS
for (i in (c(1,2))){
  
  if (i==1){heatantot <- heatantotvarrate} 
  else {heatantot <- heatantotnoadaptrate}
  
  # INITIATE PLOT
  plot(ny,heatantot[,basemort,"fact","est"],main="",
       type="n",bty="l",xlab="",ylab=ylab,
       ylim=ylim)
  box("plot",lty=1,"black")
  
  # LOOP OVER SETUP WITH/WITHOUT CLIMATE CHANGE
  
  for (j in c(1,2)){
    
    arrows(ny, heatantot[,basemort,j,"ci.l"], ny, heatantot[,basemort,j,"ci.u"], col=colorslight[j],
           code=3, angle=90, length=0.02, lwd=1)
    points(ny,heatantot[,basemort,j,"est"],col=colors[j],pch=pch[i])

    # TRY FITTING NON-LINEAR REGRESSION LINE 
    # NB: LINEAR FIT LOWEST AIC AND BIC
    fitlm <- lm(heatantot[,basemort,j,"est"] ~ ny)
    
    bestfit <- fitlm
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(bestfit, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))
    # PLOT TREND LINE
    lines(ny,pred$fit,col=colors[j],xpd=F,lty=1)
    
    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorslight[j],border=NA)
    
  }
  
  # PLOT PANEL LABEL
  text(ny[1],ylim[2],labels=panellabels[2+i])
}

dev.off()


#######
