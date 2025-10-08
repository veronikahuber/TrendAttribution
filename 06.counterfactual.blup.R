#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

###################################################################################
# 1) DERIVE FACTUAL AND COUNTERFACTUAL BLUPs 
# FROM LONGITUDINAL META-REGRESSION

# 2) PLOT FACTUAL AND COUNTERFACTUAL TEMPERATURE-MORTALITY ASSOCIATIONS BY CITY
# SUPPLEMENTARY FIGURES 2a,b AND 4
###################################################################################

######################################################
# FACTUAL AND COUNTERFACTUAL BLUPs
######################################################

# DERIVE FACTUAL BLUPs
blupf <- blup.mixmeta(modelfinal,vcov=T)
names(blupf) <- paste(cityinfo$city,cityinfo$period)

# GET BLUP RESIDUALS
blupres <- blup(modelfinal,vcov=T,type="residual")

#############################################################################
# PREDICT COUNTERFACTUAL FIXED EFFECT COEFFICIENTS FROM META-REGRESSION MODEL

# KEEP LIFE EXPECTANCY AT VALUE OF FIRST SUBPERIOD
cflife <- c(sapply(cities,function(city) rep(cityinfo[cityinfo$city==city & cityinfo$period==period[1],c("life.expectancy")],length(yearlist))))

# ASSEMBLE COUNTERFACTUAL META-PREDICTORS
datapredcf <- data.frame(year=cityinfo$year,
                       average.annual.temp=cityinfo$average.annual.temp,
                       life.expectancy=cflife,
                       sum.heat.days=cityinfo$sum.heat.days,
                       average.population.age=cityinfo$average.population.age)

predcf <- predict(modelfinal, datapredcf, vcov=T)

# SUM COUNTERFACTUAL FIXED EFFECT PREDICTIONS AND BLUP RESIDUALS

blupcf<-vector("list",length(predcf))

for (i in 1:length(predcf)) {
  blupcoef <- predcf[[i]]$fit+blupres[[i]]$blup
  blupvcov <- predcf[[i]]$vcov+blupres[[i]]$vcov
  
  blupcf[[i]] <- list(blupcoef,blupvcov)
  names(blupcf[[i]]) <- paste(list("blup","vcov"))
}

names(blupcf) <- paste(cityinfo$city,cityinfo$period)

# #####################################################################
# # CHECK FACTUAL BLUPs == FIXED EFFECT PREDICTION + BLUP RESIDUALS
# #####################################################################
# 
# # ASSEMBLE FACTUAL META-PREDICTORS
# datapredf <- data.frame(year=cityinfo$year,
#                        life.expectancy=cityinfo$life.expectancy,
#                        sum.heat.days=cityinfo$sum.heat.days)
# 
# predf <- predict(modelfinal, datapredf, vcov=T)
# 
# blupftest<-vector("list",length(predf))
# 
# for (i in 1:length(predf)) {
#   blupcoef <- predf[[i]]$fit+blupres[[i]]$blup-blupf[[i]]$blup
#   blupvcov <- predf[[i]]$vcov+blupres[[i]]$vcov-blupf[[i]]$vcov
#   
#   blupftest[[i]] <- list(blupcoef,blupvcov)
#   names(blupftest[[i]]) <- paste(list("blup","vcov"))
# }
# 
# names(blupftest) <- paste(cityinfo$city,cityinfo$period)

##################################
# OBTAIN MMTP AND MMT ESTIMATES
##################################

# GENERATE THE MATRICES FOR STORING THE RESULTS
mmtcityf <- mmtpcityf  <- mmtcitycf <- mmtpcitycf  <- array(NA,dim=c(length(dlist),length(yearlist)),
                                   dimnames=list(cities,unique(cityinfo$period)))

# DEFINE PERCENTILE RANGE FOR SEARCHING OF MMT
mmtprange <- 25:99

# DEFINE MINIMUM MORTALITY VALUES
# FOR FACTUAL AND COUNTERFACTUAL BLUPs

# LOOP OVER CITIES
for(i in seq(length(dlist))) {
  
  # PRINT CITY
  cat(cities[i],"")
  
  # EXTRACT DATA
  data <- dlist[[i]] 
  
  # SUBSET FOR SUMMER-ONLY
  datasum <- subset(data, month %in% sm)
  
  # SET FIXED KNOTS AND BOUND
  setknots <- knots[i,]

  for (p in seq(length(yearlist))){
    ysub <- yearlist[[p]]
    tmean <- subset(datasum,year %in% ysub)$tmeanf
    predvar <- quantile(tmean,mmtprange/100,na.rm=T)

    # REDEFINE THE FUNCTION 
    argvar <- list(x=predvar,fun=varfun,
                   knots=setknots)
    if(!is.null(vardegree)) argvar$degree <- vardegree
    bvar <- do.call(onebasis,argvar)
    
    # DEFAULT: USE FACTUAL BLUPs 
    mmtpcityf[i,p] <- (mmtprange)[which.min((bvar%*%blupf[[paste(cities[i],period[p])]]$blup))]
    mmtcityf[i,p] <- quantile(tmean,mmtpcityf[i,p]/100,na.rm=T)
    
    # DEFAULT: USE COUNTERFACTUAL BLUPs 
    mmtpcitycf[i,p] <- (mmtprange)[which.min((bvar%*%blupcf[[paste(cities[i],period[p])]]$blup))]
    mmtcitycf[i,p] <- quantile(tmean,mmtpcitycf[i,p]/100,na.rm=T)
    
  }
}

########################################
# PLOT FACTUAL AND COUNTERFACTUAL MMTs
########################################

filename <- paste0("MMT_",fext,".pdf")
pdf(filename,width=8,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)


# PLOT
for(i in seq(length(dlist))) {
  ylab <- "MMT (°C)"
  plot(yearpred, seq(yearpred), type="n", ylim=range(c(mmtcityf,mmtcitycf))*c(0.93,1.07), 
       ylab=ylab, xlab="", xaxt="n", bty="l", las=1, mgp=c(2.5,1,0), cex.axis=0.8,
       main=cities[i])
  
  # PLOTXAXIS LABELS
  axis(1,at=yearpred,labels=NA)
  
  text(yearpred, par("usr")[3]-0.1, 
       srt = 35, adj = 1, xpd = TRUE,
       labels = period, cex = 0.8)
  
  points(yearpred, mmtcityf[i,],  type="o", col="black", pch=19)
  points(yearpred, mmtcitycf[i,],  type="o", col="red", pch=19)
}

dev.off()

########################################################
# OVERALL CUMULATIVE ERF PLOTS BY PERIOD FOR EACH CITY
# AND SAVE RR AT 99TH PERCENTILE OF TEMPERATURE
#######################################################

# DEFINE TEMPERATURE PERCENTILES FOR COLD AND HEAT
perc <- 99

# PREDEFINE ARRAY TO STORE RR BY CITY AND PERIOD
RRf <- RRcf <- array(NA,dim=c(length(cities),length(yearlist),3),
            dimnames=list(cities,unique(cityinfo$period),c("est","ci.l","ci.u")))

# PREDEFINE VECTOR TO STORE
t99fixcity <- array(NA,dim=length(cities),dimnames=list(cities))

# PLOT PARAMETERS
xlab <- expression(paste("Temperature (",degree,"C)"))

# PLOT FOR FACTUAL AND COUNTERFACTUAL BLUPs

for (c in 1:2){
  
  if (c==1) {
    filename <- paste0("FigS2a_",fext,".pdf")
    blup <- blupf
    mmtcity <- mmtcityf
    ylim <- c(0.9,4)
  }else{
    filename <- paste0("FigS2b_",fext,".pdf")
    blup <- blupcf
    mmtcity <- mmtcitycf
    ylim <- c(0.9,4)
  }

  pdf(filename,width=7,height=9)
  layout(matrix(1:16,ncol=4,byrow=T))
  par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)
  
  # PLOT LEGEND
  plot(0,0,type="n",main="",yaxt="n",xaxt="n",ann=F,bty="n")
  legend("top",unique(cityinfo$period),lty=1,col=cl,cex=1,bty="n")
  
  for(i in seq(length(dlist))) {

    # EXTRACT DATA
    data <- dlist[[i]] 
    
    # SUBSET FOR SUMMER-ONLY
    datasum <- subset(data, month %in% sm)
    
    # SET FIXED KNOTS
    setknots <- knots[i,]

    # 99TH PERCENTILE FROM COOLEST PERIOD
    t99fixcity[i] <- quantile(subset(datasum,year %in% yearlist[[2]])$tmeanf,0.99,na.rm=T)
    
    for (p in seq(length(yearlist))){
      ysub <- yearlist[[p]]
      tmean <- subset(datasum,year %in% ysub)$tmeanf
      cen <- mmtcity[i,p]
  
      argvar <- list(x=tmean,fun=varfun,knots=setknots,Bound=range(tmean,na.rm=T))
      if(!is.null(vardegree)) argvar$degree <- vardegree
      bvar <- do.call(onebasis,argvar)
      
      # PREDICT FUNCTION FOR PLOTTING
      cp <- crosspred(bvar,coef=blup[[paste(cities[i],period[p])]]$blup,vcov=blup[[paste(cities[i],period[p])]]$vcov,
                      model.link="log",by=0.1,cen=cen)
      # DERIVE RRs
      pointcp <- crosspred(bvar,coef=blup[[paste(cities[i],period[p])]]$blup,vcov=blup[[paste(cities[i],period[p])]]$vcov,
                           model.link="log",at=t99fixcity[i],cen=cen)
      
      if (c==1){
        RRf[i,p,"est"] <- pointcp$allRRfit
        RRf[i,p,"ci.l"] <- pointcp$allRRlow
        RRf[i,p,"ci.u"] <- pointcp$allRRhigh
      } else{
        RRcf[i,p,"est"] <- pointcp$allRRfit
        RRcf[i,p,"ci.l"] <- pointcp$allRRlow
        RRcf[i,p,"ci.u"] <- pointcp$allRRhigh
      }
  
      if (p==1){
  
        plot(cp, ylim=ylim,xlim=c(min(datasum$tmeanf,na.rm=T)-1,max(datasum$tmeanf,na.rm=T)+1),xlab=xlab, ylab="RR",
             lab=c(6,5,7), lwd=1,main=cities[i],type="n",ci="n")
        
      }
      
      lines(cp, lwd=1.5, col=cl[p], ci="area", ci.arg=list(col=cllight[p]))
  
      abline(v=mmtcity[i,p], lty=3, lwd=0.9, col=cl[p])
      abline(v=t99fixcity[i], lty=2, col=grey(0.8))
    } 
  }
  
  dev.off()
}

############################
# PLOT RR
############################

ylab <-expression(RR~at~99^th~percentile)

pdf(paste0("FigS4_",fext,".pdf"),width=9,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)

# PLOT

for(i in seq(length(dlist))) {

  plot(yearpred, seq(yearpred), type="n", ylim=range(c(RRf[i,,],RRcf[i,,]))*c(0.93,1.07), 
       ylab=ylab, xlab="", xaxt="n", bty="l", las=1, mgp=c(2.5,1,0), cex.axis=0.8,
       main=cities[i])
  
  # PLOTXAXIS LABELS
  axis(1,at=yearpred,labels=NA)
  
  text(yearpred, par("usr")[3]-0.1, 
       srt = 35, adj = 1, xpd = TRUE,
       labels = period, cex = 0.8)
  
  arrows(yearpred, RRf[i,,2], yearpred, RRf[i,,3], col="grey",
         code=3, angle=90, length=0.05, lwd=2)
  points(yearpred, RRf[i,,1],  type="o", col="black", pch=19)
  
  arrows(yearpred, RRcf[i,,2], yearpred, RRcf[i,,3], col="pink",
         code=3, angle=90, length=0.05, lwd=2)
  points(yearpred, RRcf[i,,1],  type="o", col="red", pch=19)
  
  if (i==1) {legend("topright",
         c("with LE improvements","w/o LE improvements"),
         col=c("black","red"),pch=19, bty="n",cex=0.9)}
}

dev.off()

###
