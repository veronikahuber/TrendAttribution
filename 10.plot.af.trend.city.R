#####################################################################################
# R code for the analysis in 

# Huber, V., Breitner-Busch, S., Feldbusch, H. et al. Improvements in life expectancy 
# mask rising trends in heat-related excess mortality attributable to climate change. 
# Nat Commun 16, 11632 (2025). https://doi.org/10.1038/s41467-025-66681-0
######################################################################################

#####################################################################
# SUPPLEMENTARY FIG.5a,b, FIG. 8a,b
#
# PLOT ATTRIBUTABLE NUMBERS AND FRACTIONS BY CITY
# FIT LINEAR TREND AND SAVE INFORMATION ON REGRESSION COEFFICIENTS
#####################################################################

# INITIATE ARRAYS TO STORE RESULTS OF LINEAR REGRESSIONS
# HEAT-RELATED MORTALITY
morttrend <- morttrendp <- antrend <- antrendp <- array(NA,dim=c(length(cities)+1,10,2),
                                 dimnames=list(c(cities,"Germany"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),c("adapt","noadapt")))
morttrendci <- antrendci <- array(NA,dim=c(length(cities)+1,10,2,2),
                     dimnames=list(c(cities,"Germany"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),c("adapt","noadapt"),c("ci.l","ci.u")))

# HEAT-RELATED MORTALITY ATTRIBUTABLE TO CLIMATE CHANGE
attrmorttrend <- attrmorttrendp <- attrantrend <- attrantrendp <- array(NA,dim=c(length(cities)+1,2,2,9),
                                 dimnames=list(c(cities,"Germany"),c("adapt","noadapt"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8))))
attrmorttrendci <- attrantrendci <- array(NA,dim=c(length(cities)+1,2,2,9,2),
                                         dimnames=list(c(cities,"Germany"),c("adapt","noadapt"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),c("ci.l","ci.u")))

# NAMES OF COUNTERFACTUAL TEMPERATURE VERSIONS WITHOUT BEST ESTIMATE
cfactci <- paste0("cfact.ci.",1:8)

###########################################
# PLOT HEAT-ATTRIBUTABLE FRACTIONS BY CITY
############################################

# DEFINE WHETHER TO USE OBSERVED ("obs") OR AVERAGED ("sim") MORTALITY
basemort <- "obs" 

if (basemort=="sim") {
  fext <- vers[3]
}

#########################
# WITH LE IMPROVEMENTS
#########################

xlab <- "Years"
ylab2 <- "AF (%)"

# DEFINE POINT TYPES
pch <- c(20,18)

plotname <- paste0("FigS5a_",fext,".pdf")
pdf(plotname,width=9,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)

# LOOP OVER CITIES
for(i in seq(length(dlist))) {

  # INITIATE PLOT
  plot(ny,heatafvar[i,,basemort,"fact","est"],main=cities[i],type="n",bty="l",xlab=xlab,ylab=ylab2,
     ylim=c(0,max(heatafvar[i,,basemort,"fact","ci.u"],na.rm=T)))
  
  # LOOP OVER ALTERNATIVE COUNTERFACTUAL TEMPERATURE VERSIONS AND SAVE LINEAR TREND ESTIMATES
  for (c in seq(cfactci)){

    fitlm <- lm(heatafvar[i,,basemort,2+c,"est"] ~ ny)
    
    # SAVE LINEAR SLOPE
    morttrend[i,2+c,"adapt"] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    morttrendci[i,2+c,"adapt",] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    morttrendp[i,2+c,"adapt"] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))

  }

  # PLOT BEST ESTIMATES FOR FACTUAL AND COUNTERFACTUAL TEMPERATURES
  # AND SAVE LINEAR TREND ESTIMATES
  
  for (j in c(1,2)){
  
    arrows(ny, heatafvar[i,,basemort,j,"ci.l"], ny, heatafvar[i,,basemort,j,"ci.u"], col=colorslight[j],
           code=3, angle=90, length=0.02, lwd=0.7)
    points(ny,heatafvar[i,,basemort,j,"est"],col=colors[j],pch=pch[1])

    # TRY FITTING NON-LINEAR REGRESSION LINE 
    # NB: LINEAR FIT LOWEST AIC AND BIC
    fitlm <- lm(heatafvar[i,,basemort,j,"est"] ~ ny, na.action="na.exclude")
    
    # SAVE LINEAR SLOPE
    morttrend[i,j,"adapt"] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    morttrendci[i,j,"adapt",] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    morttrendp[i,j,"adapt"] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))
    # PLOT TREND LINE
    lines(ny,pred$fit,col=colors[j],xpd=F,lty=1)
 
    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorslight[j],border=NA)
    
  }
    
  if (i==1) {legend("topright",c("with CC","w/o CC"),col=colors,lty=1,pch=pch[1],cex=0.8,bty="n")}
}

dev.off()


###########################
# WITHOUT LE IMPROVEMENTS
###########################

plotname <- paste0("FigS5b_",fext,".pdf")
pdf(plotname,width=9,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)

# LOOP OVER CITIES
for(i in seq(length(dlist))) {
  
  # INITIATE PLOT
  plot(ny,heatafnoadapt[i,,basemort,"fact","est"],main=cities[i],type="n",bty="l",xlab=xlab,ylab=ylab2,
       ylim=c(0,max(heatafnoadapt[i,,basemort,"fact","ci.u"],na.rm=T)))
  
  # LOOP OVER ALTERNATIVE COUNTERFACTUAL TEMPERATURE VERSIONS AND SAVE LINEAR TREND ESTIMATES
  for (c in seq(cfactci)){
    
    fitlm <- lm(heatafnoadapt[i,,basemort,2+c,"est"] ~ ny,na.action=na.exclude)
    
    # SAVE LINEAR SLOPE
    morttrend[i,2+c,"noadapt"] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    morttrendci[i,2+c,"noadapt",] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    morttrendp[i,2+c,"noadapt"] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))

  }
  
  # PLOT BEST ESTIMATES FOR FACTUAL AND COUNTERFACTUAL TEMPERATURES
  # AND SAVE LINEAR TREND ESTIMATES
  
  for (j in c(1,2)){
    
    arrows(ny, heatafnoadapt[i,,basemort,j,"ci.l"], ny, heatafnoadapt[i,,basemort,j,"ci.u"], col=colorslight[j],
           code=3, angle=90, length=0.02, lwd=0.7)
    points(ny,heatafnoadapt[i,,basemort,j,"est"],col=colors[j],pch=pch[2])
    
    # TRY FITTING NON-LINEAR REGRESSION LINE 
    # NB: LINEAR FIT LOWEST AIC AND BIC
    fitlm <- lm(heatafnoadapt[i,,basemort,j,"est"] ~ ny, na.action=na.exclude)
    
    # SAVE LINEAR SLOPE
    morttrend[i,j,"noadapt"] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    morttrendci[i,j,"noadapt",] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    morttrendp[i,j,"noadapt"] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))

    # PLOT TREND LINE
    lines(ny,pred$fit,col=colors[j],xpd=F,lty=1)

    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorslight[j],border=NA)
    
  }
  
  if (i==1) {legend("topleft",c("with CC","w/o CC"),col=colors,lty=1,pch=pch[2],cex=0.8,bty="n")}
}

dev.off()

########################################
# ATTRIBUTABLE TO CLIMATE CHANGE 
# ABSOLUTE DIFFERENCES FACT - CFACT
########################################

ylab <- expression('AF'['CC']~'(%)')
col <- c("black","red")
collight <- colextralight <- vector("character",2)

for (i in seq(col)){
  collight[i] <- do.call(rgb,c(as.list(col2rgb(col[i])),alpha=255/6,max=255))
  colextralight[i] <- do.call(rgb,c(as.list(col2rgb(col[i])),alpha=255/18,max=255))
}

plotname <- paste0("FigS8a_",fext,".pdf")
pdf(plotname,width=9,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)

# LOOP OVER CITIES
for(i in seq(length(dlist))) {

  for (a in (c(2,1))){
    
    if (a==1){heataf <- heatafvardif} 
    else {heataf <- heatafnoadaptdif}
    
    if(a==2) {plot(ny,heataf[i,,basemort,"abs","cfact.best","est"],main=cities[i],type="n",bty="l",xlab=xlab,
                   ylab=ylab,
                   ylim=c(0,max(heatafnoadaptdif[i,,basemort,"abs",,"ci.u"],na.rm=T)))}
    
    # LOOP OVER ALTERNATIVE COUNTERFACTUAL TEMPERATURE VERSIONS AND SAVE LINEAR TREND ESTIMATES
    for (c in seq(cfactci)){
      
      fitlm <- lm(heataf[i,,basemort,"abs",1+c,"est"] ~ ny)
      
      # SAVE LINEAR SLOPE
      attrmorttrend[i,a,"abs",1+c] <- fitlm$coefficients[2]
      # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
      attrmorttrendci[i,a,"abs",1+c,] <- confint(fitlm,'ny',level=0.95)
      # AND CORRESPONDING P-VALUE
      attrmorttrendp[i,a,"abs",1+c] <- summary(fitlm)$coefficients[2,4] 
      
      # COMPUTE ERROR BANDS (2*SE)
      pred = predict(fitlm, newdata = list(ny), se = TRUE)
      se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                   "lower" = fit-2*se.fit)) 
    }
    
    arrows(ny, heataf[i,,basemort,"abs","cfact.best","ci.l"], ny, heataf[i,,basemort,"abs","cfact.best","ci.u"], col=collight[a],
           code=3, angle=90, length=0.02, lwd=0.7)
    points(ny,heataf[i,,basemort,"abs","cfact.best","est"],col=col[a],pch=pch[a])

    # INCLUDE LINEAR REGRESSION LINE FOR BEST ESTIMATE
    fitlm <- lm(heataf[i,,basemort,"abs","cfact.best","est"] ~ ny)
    
    # SAVE LINEAR SLOPE
    attrmorttrend[i,a,"abs","cfact.best"] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    attrmorttrendci[i,a,"abs","cfact.best",] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    attrmorttrendp[i,a,"abs","cfact.best"] <- summary(fitlm)$coefficients[2,4] 

    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))
    # PLOT TREND LINE
    lines(ny,pred$fit,col=col[a],xpd=F,lty=1)

    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=collight[a],border=NA)
    
   }
  
  if (i==1) {legend("top",c("with LE improvements","w/o LE improvements"),lty=1,col=col,pch=pch,bty="n",cex=0.7)}
  
}

dev.off()

########################################
# ATTRIBUTABLE TO CLIMATE CHANGE 
# RELATIVE DIFFERENCES FACT - CFACT
########################################

ylab <- expression('P'['CC']~'(%)')
plotname <- paste0("FigS8b_",fext,".pdf")
pdf(plotname,width=9,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)

# LOOP OVER CITIES
for(i in seq(length(dlist))) {
  
    for (a in (c(2,1))){
    
    if (a==1){heataf <- heatafvardif} 
    else {heataf <- heatafnoadaptdif}
    
    ylim <- c(min(heataf[i,,basemort,"rel",,"ci.l"],na.rm=T),max(heatafvardif[i,,basemort,"rel",,"ci.u"],na.rm=T))

    if(a==2) {plot(ny,heataf[i,,basemort,"rel","cfact.best","est"],main=cities[i],type="n",bty="l",xlab=xlab,
                   ylab=ylab,ylim=ylim)}
    
    # LOOP OVER ALTERNATIVE COUNTERFACTUAL TEMPERATURE VERSIONS AND SAVE LINEAR TREND ESTIMATES
    for (c in seq(cfactci)){
      
      fitlm <- lm(heataf[i,,basemort,"rel",1+c,"est"] ~ ny)
      
      # SAVE LINEAR SLOPE
      attrmorttrend[i,a,"rel",1+c] <- fitlm$coefficients[2]
      # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
      attrmorttrendci[i,a,"rel",1+c,] <- confint(fitlm,'ny',level=0.95)
      # AND CORRESPONDING P-VALUE
      attrmorttrendp[i,a,"rel",1+c] <- summary(fitlm)$coefficients[2,4] 
      
      # COMPUTE ERROR BANDS (2*SE)
      pred = predict(fitlm, newdata = list(ny), se = TRUE)
      se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                   "lower" = fit-2*se.fit)) 

    }
    
    arrows(ny, heataf[i,,basemort,"rel","cfact.best","ci.l"], ny, heataf[i,,basemort,"rel","cfact.best","ci.u"], col=collight[a],
           code=3, angle=90, length=0.02, lwd=0.7)
    points(ny,heataf[i,,basemort,"rel","cfact.best","est"],col=col[a],pch=pch[a])
    
    # INCLUDE LINEAR REGRESSION LINE
    fitlm <- lm(heataf[i,,basemort,"rel","cfact.best","est"] ~ ny)
    
    # SAVE LINEAR SLOPE
    attrmorttrend[i,a,"rel","cfact.best"] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    attrmorttrendci[i,a,"rel","cfact.best",] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    attrmorttrendp[i,a,"rel","cfact.best"] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))
    # PLOT TREND LINE
    lines(ny,pred$fit,col=col[a],xpd=F,lty=1)
    
    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=collight[a],border=NA)
    
  }
  
  if (i==1) {legend("bottomleft",c("with LE improvements","w/o LE improvements"),lty=1,col=col,pch=pch,bty="n",cex=0.7)}
  
}

dev.off()


###
