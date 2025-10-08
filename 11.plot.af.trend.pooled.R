#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

########################################################################################
# PLOT POOLED HEAT-RELATED MORTALITY AND LINEAR TRENDS (FIG. 3)
# 
# WITH/WITHOUT CLIMATE CHANGE
# WITH/WITHOUT IMPROVEMENTS IN LIFE EXPECTANCIES
#
# PLOT CORRESPONDING ATTRIBUTABLE MORTALITY MEASURES (FIG. 4, SUPPLEMENTARY FIG. 7)
########################################################################################

#####################################
# PLOT POOLED ATTRIBUTABLE FRACTIONS
#####################################

title <- c("with LE improvements",
           "w/o LE improvements")
ylab <- "AF (%)"
panellabels <- c("a","b","c","d")

ylim <- c(0,max(c(max(heataftotnoadapt[,basemort,"fact","ci.u"],na.rm=T),max(heataftotvar[,basemort,"fact","ci.u"],na.rm=T))))

# DEFINE POINT TYPES
pch <- c(20,18)

plotname <- paste0("Fig3_",fext,".pdf")
pdf(plotname,width=5,height=4.5)
layout(matrix(1:4,ncol=2,byrow = T))
par(mar=c(2, 4.1, 3, 1.5), mgp=c(2.5, 1, 0),cex=0.75)

# LOOP OVER SETUP WITH/WITHOUT LE IMPROVEMENTS
for (i in (c(1,2))){
  
  if (i==1){heataftot <- heataftotvar}  else {heataftot <- heataftotnoadapt}
  
  # INITIATE PLOT
  plot(ny,heataftot[,basemort,"fact","est"],main=title[i],cex.main=1,
       type="n",bty="l",xlab="",ylab=ylab,ylim=ylim)
  box("plot",lty=1,"black")

  # LOOP OVER ALTERNATIVE COUNTERFACTUAL TEMPERATURE VERSIONS AND SAVE LINEAR TREND ESTIMATES
  for (c in seq(cfactci)){
    
    fitlm <- lm(heataftot[,basemort,2+c,"est"] ~ ny)
    
    # SAVE LINEAR SLOPE
    morttrend["Germany",2+c,i] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    morttrendci["Germany",2+c,i,] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    morttrendp["Germany",2+c,i] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit))

  }  
  
  # PLOT BEST ESTIMATES FOR FACTUAL AND COUNTERFACTUAL TEMPERATURES
  # AND SAVE LINEAR TREND ESTIMATES
  
  # LOOP OVER SETUP WITH/WITHOUT CLIMATE CHANGE
  for (j in c(1,2)){
    
    arrows(ny, heataftot[,basemort,j,"ci.l"], ny, heataftot[,basemort,j,"ci.u"], col=colorslight[j],
           code=3, angle=90, length=0.02, lwd=1)
    points(ny,heataftot[,basemort,j,"est"],col=colors[j],pch=pch[i])

    # TRY FITTING NON-LINEAR REGRESSION LINE 
    # NB: LINEAR FIT LOWEST AIC AND BIC
    fitlm <- lm(heataftot[,basemort,j,"est"] ~ ny)
    
    bestfit <- fitlm
    
    # SAVE LINEAR SLOPE
    morttrend["Germany",j,i] <- bestfit$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    morttrendci["Germany",j,i,] <- confint(bestfit,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    morttrendp["Germany",j,i] <- summary(bestfit)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(bestfit, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                "lower" = fit-2*se.fit))
    # PLOT TREND LINE
    lines(ny,pred$fit,col=colors[j],xpd=F,lty=1)
    
    # PLOT UNCERTAINTY BAND
    polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=colorslight[j],border=NA)

  }
  
  legend("topright",c("with CC","w/o CC"),col=colors,lty=1,pch=pch[i], cex=0.8,bty="n")
  
  # PLOT PANEL LABEL
  text(ny[1],ylim[2],labels=panellabels[i])
}

#########################################
# PLOT LINEAR SLOPES BY CITY AND BOXPLOT
##########################################

# TEXT SIZE AND OFFSET FOR ANNOTATION OF SIGNFICANCE TESTS
cex <- c(0.8,1.2)
offset <- c(0.02,0.02)
  
pch=c(1,5)
ylab <- expression(paste(Delta,"AF (%/year)"))
par(mar=c(3, 4.1, 2, 1.5), mgp=c(2.5, 1, 0),cex=0.75)

# DEFINE Y-SCALE
ylim <- range(morttrend)
ylim[2] <- ylim[2]+0.05

# LOOP OVER SETUP WITH/WITHOUT LE IMPROVEMENTS
for (a in 1:2){
  
  # COMPUTE WILCOXON-RANK-TEST ON CITY-SPECIFIC TRENDS
  # TO DETERMINE WHETHER CLIMATE CHANGE INDUCED A SIGNIFICANT EFFECT
  wrt <- wilcox.test(morttrend[1:ncities,"fact",a],as.vector(morttrend[1:ncities,2:10,a]),alternative="two.sided")
  labels <- wrtlabel(wrt$p.value)
  print(wrt$p.value)
  
  plot(rep(1,ncities),morttrend[1:ncities,1,a],type="n",xaxt="n",
       ylab=ylab,xlab="",xlim=c(0,3),ylim=ylim,cex.axis=0.9)
  axis(side=1,at=c(1,2),labels=c("with CC","w/o CC"),cex.axis=0.85)
  
  # PLOT AF TRENDS BASED ON FACTUAL TEMPERATURES
  boxplot(morttrend[1:ncities,"fact",a],add=T,at=1,border="black",col=colorslight[1],
          xaxt="n",yaxt="n",ylim=ylim, outpch = 18)
  points(rep(1,ncities),morttrend[1:ncities,"fact",a],col=colors[1],pch=pch[a])
  points(1,morttrend[ncities+1,"fact",a],pch=4,col="black",cex=1.7)
  
  # PLOT BOXPLOT FOR AF TRENDS BASED ON ALL COUNTERFACTUAL TEMPERATURES
  boxplot(as.vector(morttrend[1:ncities,2:10,a]),add=T,at=2,border="black",col=colorslight[2],
          xaxt="n",yaxt="n",ylim=ylim,outpch= 18)
  
  # ADD BEST-ESTIMATE COUNTERFACTUAL AF TREND
  points(rep(2,ncities),morttrend[1:ncities,"cfact.best",a],col=colors[2],pch=pch[a])
  points(2,morttrend[ncities+1,"cfact.best",a],pch=4,col="black",cex=1.7)
  
  # ADD SIGNIFICANCE LEVEL OF DIFFERENCES 
  yline <- ylim[2]-offset[a]
  lines(c(1,2),c(yline,yline),lty=1,col="black")
  lines(c(1,1),c(yline,yline-offset[a]/4),lty=1,col="black")
  lines(c(2,2),c(yline,yline-offset[a]/4),lty=1,col="black")
  text(1.5,yline+offset[a]/1.5,labels=labels,cex=cex[a])
  
  abline(h=0,lty=2,col="grey")
  
  # PLOT PANEL LABEL
  text(0,ylim[2],labels=panellabels[2+a])
}

dev.off()

########################################
# PLOT ATTRIBUTABLE MORTALITY MEASURES
########################################

# RE-DEFINE POINT TYPES
pch <- c(20,18)
ylab <- expression('AF'['CC']~'(%)')

plotname <- paste0("Fig4_",fext,".pdf")
pdf(plotname,width=2.5,height=4)
layout(matrix(1:2,ncol=1,byrow = T))
par(mar=c(2, 4.1, 1, 1.5), mgp=c(2.5, 1, 0),cex=0.75)

##########################################
# PLOT ABSOLUTE DIFFERENCES FACT - CFACT
##########################################

# LOOP OVER SETUP WITH/WITHOUT LE IMPROVEMENTS
for (i in (c(2,1))){
  
  if (i==1){heataftot <- heataftotvardif} 
  else {heataftot <- heataftotnoadaptdif}
  
  if(i==2) {
    plot(ny,heataftot[,basemort,"abs","cfact.best","est"],main="",cex.main=1,
         type="n",bty="l",xlab=xlab,
       ylab=ylab,
       ylim=c(0,max(heataftotnoadaptdif[,basemort,"abs",,"ci.u"],na.rm=T)))
    box("plot",lty=1,"black")
    }
  
  # LOOP OVER ALTERNATIVE COUNTERFACTUAL TEMPERATURE VERSIONS AND SAVE LINEAR TREND ESTIMATES
  for (c in seq(cfactci)){
    
    fitlm <- lm(heataftot[,basemort,"abs",1+c,"est"] ~ ny)
    
    # SAVE LINEAR SLOPE
    attrmorttrend["Germany",i,"abs",1+c] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    attrmorttrendci["Germany",i,"abs",1+c,] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    attrmorttrendp["Germany",i,"abs",1+c] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit)) 
  }
  
  arrows(ny, heataftot[,basemort,"abs","cfact.best","ci.l"], ny, 
         heataftot[,basemort,"abs","cfact.best","ci.u"], col=collight[i],
         code=3, angle=90, length=0.02, lwd=0.7)
  points(ny,heataftot[,basemort,"abs","cfact.best","est"],col=col[i],pch=pch[i])

  # INCLUDE LINEAR REGRESSION LINE
  # TRY FITTING NON-LINEAR REGRESSION LINE 
  # NB: LINEAR BEST FIT BASED ON AIC AND BIC
  fitlm <- lm(heataftot[,basemort,"abs","cfact.best","est"] ~ ny)

  # SAVE LINEAR SLOPE
  attrmorttrend["Germany",i,"abs","cfact.best"] <- fitlm$coefficients[2]
  # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
  attrmorttrendci["Germany",i,"abs","cfact.best",] <- confint(fitlm,'ny',level=0.95)
  # AND CORRESPONDING P-VALUE
  attrmorttrendp["Germany",i,"abs","cfact.best"] <- summary(fitlm)$coefficients[2,4] 
  
  # COMPUTE ERROR BANDS (2*SE)
  pred = predict(fitlm, newdata = list(ny), se = TRUE)
  se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                               "lower" = fit-2*se.fit))
  # PLOT TREND LINE
  lines(ny,pred$fit,col=col[i],xpd=F,lty=1)
  
  # PLOT UNCERTAINTY BAND
  polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=collight[i],border=NA)
  
  # ADD LEGEND
  if (i==2){legend("top",c("with LE improvements","w/o LE improvements"),
                   lty=1,col=col, pch=pch, bty="n",cex=0.8)  }
}

# PLOT PANEL LABEL
text(ny[1],max(heataftotnoadaptdif[,basemort,"abs",,"ci.u"],na.rm=T),labels="a")

####################################
# RELATIVE DIFFERENCES FACT - CFACT
####################################

ylab <- expression('P'['CC']~'(%)')
ylim <- c(min(heataftot[,basemort,"rel","cfact.best","ci.l"],na.rm=T),max(heataftot[,basemort,"rel","cfact.best","est"]+20,na.rm=T))

# LOOP OVER SETUP WITH/WITHOUT LE IMPROVEMENTS
for (i in (c(2,1))){
  
  if (i==1){heataftot <- heataftotvardif} 
  else {heataftot <- heataftotnoadaptdif}
  
  if (i==2) {
    plot(ny,heataftot[,basemort,"rel","cfact.best","est"],main="",cex.main=1,
         type="n",bty="l",xlab=xlab,ylab=ylab,
       ylim=ylim)
    box("plot",lty=1,"black")
    }
  
  # LOOP OVER ALTERNATIVE COUNTERFACTUAL TEMPERATURE VERSIONS AND SAVE LINEAR TREND ESTIMATES
  for (c in seq(cfactci)){
    
    fitlm <- lm(heataftot[,basemort,"rel",1+c,"est"] ~ ny)
    
    # SAVE LINEAR SLOPE
    attrmorttrend["Germany",i,"rel",1+c] <- fitlm$coefficients[2]
    # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
    attrmorttrendci["Germany",i,"rel",1+c,] <- confint(fitlm,'ny',level=0.95)
    # AND CORRESPONDING P-VALUE
    attrmorttrendp["Germany",i,"rel",1+c] <- summary(fitlm)$coefficients[2,4] 
    
    # COMPUTE ERROR BANDS (2*SE)
    pred = predict(fitlm, newdata = list(ny), se = TRUE)
    se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                                 "lower" = fit-2*se.fit)) 
    
  }
  
  arrows(ny, heataftot[,basemort,"rel","cfact.best","ci.l"], ny, heataftot[,basemort,"rel","cfact.best","ci.u"], col=collight[i],
         code=3, angle=90, length=0.02, lwd=0.7)
  points(ny,heataftot[,basemort,"rel","cfact.best","est"],col=col[i],pch=pch[i])
  
  # INCLUDE LINEAR REGRESSION LINE
  # TRY FITTING NON-LINEAR REGRESSION LINE 
  # NB: LINEAR BEST FIT BASED ON AIC AND BIC
  fitlm <- lm(heataftot[,basemort,"rel","cfact.best","est"] ~ ny)

  # SAVE LINEAR SLOPE
  attrmorttrend["Germany",i,"rel","cfact.best"] <- fitlm$coefficients[2]
  # 95% CONFIDENCE INTERVALS OF SLOPE ESTIMATE
  attrmorttrendci["Germany",i,"rel","cfact.best",] <- confint(fitlm,'ny',level=0.95)
  # AND CORRESPONDING P-VALUE
  attrmorttrendp["Germany",i,"rel","cfact.best"] <- summary(fitlm)$coefficients[2,4]
  
  # COMPUTE ERROR BANDS (2*SE)
  pred = predict(fitlm, newdata = list(ny), se = TRUE)
  se_bands <- with(pred, cbind("upper" = fit+2*se.fit, 
                               "lower" = fit-2*se.fit))
  # PLOT TREND LINE
  lines(ny,pred$fit,col=col[i],xpd=F,lty=1)
  
  # PLOT UNCERTAINTY BAND
  polygon(c(ny,rev(ny)),c(se_bands[,"lower"],rev(se_bands[,"upper"])),col=collight[i],border=NA)

}

# PLOT PANEL LABEL
text(ny[1],ylim[2],labels="b")

dev.off()

#########################################
# PLOT LINEAR SLOPES BY CITY AND BOXPLOT
##########################################

pch <- c(1,5)

# TEXT SIZE FOR ANNOTATION OF SIGNFICANCE TESTS
cex <- c(1.2,0.8)
offset <- c(0.015,0.1)

ylab <- vector("list",2)
ylab[[1]] <-  expression(paste(Delta,'AF'['CC']~'(%/year)'))
ylab[[2]] <-  expression(paste(Delta,'P'['CC']~'(%/year)'))

plotname <- paste0("FigS7_",fext,".pdf")
pdf(plotname,width=2.5,height=4)
layout(matrix(1:2,ncol=1,byrow = T))
par(mar=c(2, 4.1, 1, 1.5), mgp=c(2.5, 1, 0),cex=0.75)

# LOOP OVER ABSOLUTE/RELATIVE ATTRIBUTABLE MEASURES
for (j in 1:2){

  # COMPUTE WILCOXON-RANK-TEST ON CITY-SPECIFIC TRENDS
  # TO DETERMINE WHETHER ACCOUNTING FOR OR NOT ACCOUNTING FOR LIFE EXPECTANCY GAINS MADE A DIFFERENCE
  wrt <- wilcox.test(as.vector(attrmorttrend[1:ncities,"adapt",j,]),as.vector(attrmorttrend[1:ncities,"noadapt",j,]),alternative="two.sided")
  labels <- wrtlabel(wrt$p.value)
  print(wrt$p.value)  
  
  if(j==1){ylim<-range(attrmorttrend[,,j,],na.rm=T)*1.3}else{ylim<-c(0,1.2)}
  plot(rep(1,ncities),attrmorttrend[1:ncities,1,j,"cfact.best"],type="n",xaxt="n",
       ylab=ylab[[j]],xlab="",xlim=c(0,3),ylim=ylim,cex.axis=1)
  axis(side=1,at=c(1,2),labels=c("with LE\nimprovements","w/o LE\nimprovements"),cex.axis=0.6)
  
  # PLOT BOXPLOTS FOR WITH/WITHOUT LE IMPROVEMENTS
  for (i in 1:2) {
    # BOXPLOT ACROSS ALL COUNTERFACTUAL TEMPERATURES
    boxplot(as.vector(attrmorttrend[1:ncities,i,j,]),add=T,at=i,border="black",col=collight[i],xaxt="n",yaxt="n",outpch=18)
    # POINTS FOR BEST ESTIMATES
    points(rep(i,ncities),attrmorttrend[1:ncities,i,j,"cfact.best"],col=col[i],pch=pch[i])
    points(i,attrmorttrend[ncities+1,i,j,"cfact.best"],pch=4,col="black",cex=1.7)
  }
  abline(h=0,lty=2,col="grey")
  
  # ADD SIGNIFICANCE LEVEL OF DIFFERENCES 
  yline <- ylim[2]-offset[j]
  lines(c(1,2),c(yline,yline),lty=1,col="black")
  lines(c(1,1),c(yline,yline-offset[j]/5),lty=1,col="black")
  lines(c(2,2),c(yline,yline-offset[j]/5),lty=1,col="black")
  text(1.5,yline+offset[j]/2,labels=labels,cex=cex[j])
  
  # PLOT PANEL LABEL
  text(0,ylim[2],labels=panellabels[j])
}

dev.off()

####
