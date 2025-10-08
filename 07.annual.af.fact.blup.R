#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

#################################################################################################
# COMPUTE ANNUAL HEAT-RELATED MORTALITY BASED ON FACTUAL AND COUNTERFACTUAL TEMPERATURE DATA

# USE FACTUAL BLUPs FROM LONGITUDINAL META-REGRESSION-MODEL ('WITH LE IMPROVEMENTS')
################################################################################################

#####################################################
# COMPUTE ANNUAL HEAT-ATTRIBUTABLE NUMBERS/FRACTIONS
# FOR FACTUAL AND COUNTERFACTUAL DATA

# USING OBSERVED OR AVERAGED DAILY MORTALITY DATA
#####################################################

###################################
# CREATE OBJECTS TO STORE RESULTS
###################################

# YEAR VECTOR
ny <- range(yearlist)[1]:range(yearlist)[2]

# NUMBERS OF DAYS PER WARM SEASON
ndwarm <- length(subset(dlist[[1]], month %in% sm & year==2020)$date)

# COUNTERFACTUAL TEMPERATURE NAMES
tmeancfvar <- paste0("tmeancf",1:9)

# HEAT ATTRIBUTABLE NUMBERS AND FRACTIONS BY CITY AND YEAR
heatanvar <- heatafvar <- array(NA,dim=c(length(cities),length(ny),2,10,3),
                          dimnames=list(cities,ny,c("obs","sim"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# CLIMATE CHANGE ATTRIBUTABLE NUMBERS AND FRACTIONS BY CITY AND YEAR
heatanvardif <- heatafvardif <- array(NA,dim=c(length(cities),length(ny),2,2,9,3),
                                dimnames=list(cities,ny,c("obs","sim"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# CLIMATE CHANGE ATTRIBUTABLE NUMBERS AND FRACTIONS BY CITY OVER TOTAL STUDY PERIOD (1993-2022)
meanheatanvardif <- meanheatafvardif <- array(NA,dim=c(length(cities),2,2,9,3),
                                      dimnames=list(cities,c("obs","sim"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# HEAT ATTRIBUTABLE NUMBERS AND FRACTIONS FOR ALL CITIES BY YEAR
heatantotvar <- heataftotvar <- array(NA,dim=c(length(ny),2,10,3),
                          dimnames=list(ny,c("obs","sim"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# CLIMATE CHANGE ATTRIBUTABLE NUMBERS AND FRACTIONS FOR ALL CITIES BY YEAR
heatantotvardif <- heataftotvardif <- array(NA,dim=c(length(ny),2,2,9,3),
                                      dimnames=list(ny,c("obs","sim"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# CLIMATE CHANGE ATTRIBUTABLE NUMBERS AND FRACTIONS FOR ALL CITIES OVER TOTAL STUDY PERIOD (1993-2022)
meanheatantotvardif <- meanheataftotvardif <- array(NA,dim=c(2,2,9,3),
                                            dimnames=list(c("obs","sim"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# ANNUAL DEATH COUNTS
deathtot <- array(NA,dim=c(length(cities),length(ny)),dimnames=list(cities,ny))
deathtotavg <- array(NA,dim=length(cities),dimnames=list(cities))

# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE MONTE CARLO SIMULATIONS OF ATTRIBUTABLE DEATHS
arraysimvar <- array(NA,dim=c(length(cities),length(ny),2,10,nsim),
                  dimnames=list(cities,ny,c("obs","sim"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),paste0("sim",seq(nsim))))

arraysimvardif <- array(NA,dim=c(length(cities),length(ny),2,2,9,nsim),
                     dimnames=list(cities,ny,c("obs","sim"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),paste0("sim",seq(nsim))))

#######################################################
# COMPUTE ATTRIBUTABLE NUMBERS
#######################################################

# LOOP OVER CITIES
for (i in seq(length(dlist))){
  
  city <- cities[i]
  cat("\n\n",city,"\n")
  
  # EXTRACT CITY-SPECIFIC DATA AND RESTRICT TO STUDY-PERIOD
  datafull <- subset(dlist[[i]],year %in% ny)

  # SUBSET TO SUMMER 
  data <- subset(datafull, month %in% sm)
  
  # DERIVE AVERAGED MORTALITY SERIES
  deathdoy <- tapply(data$death,rep(seq(ndwarm),length(ny)),mean,na.rm=T)
  while(any(isna <- is.na(deathdoy)))
    deathdoy[isna] <- rowMeans(Lag(deathdoy,c(-1,1)),na.rm=T)[isna]
  
  # TOTAL AVERAGE MORTALITY PER SUMMER
  deathtotavg[i] <- sum(deathdoy)
  
  # TOTAL OBSERVED MORTALITY PER SUMMER
  deathtot[i,] <- tapply(data$death,data$year,sum,na.rm=T)

  # COMPUTE FORWARD MOVING AVERAGE OF OBSERVED DEATHS ACROSS LAG PERIOD 
  avgdeath <- rowMeans(as.matrix(tsModel:::Lag(datafull$death,-seq(0,lag))))
  
  # AND EXTRACT SUMMER DATA
  avgdeathsum <- avgdeath[datafull$month %in% sm]
  
  # LOOP OVER PERIODS
  for (p in seq(length(yearlist))){
    
    # EXTRACT YEARS FOR SUBPERIOD AND SUBSET DATA
    ysub <- yearlist[[p]]
    nyper <- length(ysub)
    datasub <- subset(data,year %in% ysub)

    # DEFINE ARGVAR, COEF-VCOV, AND CENVEC
    # BASED ON FACTUAL TEMPERATURE
    # WITH FIXED INTERNAL KNOTS AND SUB-PERIOD SPECIFIC BOUNDARY KNOTS
    
    bound <- range(datasub$tmeanf,na.rm=T)
    argvar <- list(fun=varfun,knots=knots[i,],Bound=bound)
    if(!is.null(vardegree)) argvar$degree <- vardegree
    
    coef <- blupf[[paste(cities[i],period[p])]]$blup
    vcov <- blupf[[paste(cities[i],period[p])]]$vcov

    cen <- mmtcityf[i,p]
    cenvec <- do.call(onebasis,c(list(x=cen),argvar))

    #################################
    # USE FACTUAL TEMPERATURE DATA
    #################################
    
    # DERIVE TEMPERATURE BASES AND HEAT INDICATORS 
    bvar <- do.call(onebasis,c(list(x=datasub$tmeanf),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=F)
    ind <- datasub$tmeanf>cen
    
    # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
    # AND SUM HEAT ATTRIBUTABLE DEATH PER YEAR
    
    # A) USE AVERAGE DAILY MORTALITY - "MODEL"
    deathproj <- rep(deathdoy,nyper)
    an <- (1-exp(-bvarcen%*%coef))*deathproj
    # SET NEGATIVE ESTIMATES TO 0
    an[an < 0] <- 0
    # SUM PER YEAR 
    heatanvar[i,as.character(ysub),"sim","fact","est"] <- tapply(an[ind],datasub$year[ind],sum,na.rm=T)
    
    # B) USE OBSERVED DAILY MORTALITY - "OBSERVATION"
    # GET FORWARD AVERAGED DEATH SUB-PERIOD
    avgdeathsub <- avgdeathsum[data$year %in% ysub]
    # (1) AND FACTUAL TEMPERATURE DATA
    anobs <- (1-exp(-bvarcen%*%coef))*avgdeathsub
    # SET NEGATIVE ESTIMATES TO 0
    anobs[anobs < 0] <- 0
    # SUM PER YEAR 
    heatanvar[i,as.character(ysub),"obs","fact","est"] <- tapply(anobs[ind],datasub$year[ind],sum,na.rm=T)
    
    # COMPUTE EMPIRICAL CONFIDENCE INTERVALS
    
    set.seed(13041975)
    coefsim <- mvrnorm(nsim,coef,vcov)
    
    # LOOP ACROSS MONTE CARLO ITERATIONS
    for(s in seq(nsim)) {
      
      # FACTUAL TEMPERATURE DATA
      an <- (1-exp(-bvarcen%*%coefsim[s,]))*deathproj
      an[an < 0] <- 0
      arraysimvar[i,as.character(ysub),"sim","fact",s] <- tapply(an[ind],datasub$year[ind],sum,na.rm=T)
      
      anobs <- (1-exp(-bvarcen%*%coefsim[s,]))*avgdeathsub
      anobs[anobs < 0] <- 0
      arraysimvar[i,as.character(ysub),"obs","fact",s] <- tapply(anobs[ind],datasub$year[ind],sum,na.rm=T)

    }
    
    #########################################
    # USE COUNTERFACTUAL TEMPERATURE DATA
    #########################################
    
    # LOOP OVER COUNTERFACTUAL TEMPERATURE VERSIONS
    
    for (c in seq(tmeancfvar)){
      
      # DERIVE TEMPERATURE BASES AND HEAT INDICATORS
      tmeancf <- datasub[,tmeancfvar[c]]
      bvarcf <- do.call(onebasis,c(list(x=tmeancf),argvar))
      bvarcencf <- scale(bvarcf,center=cenvec,scale=F)
      indcf <- tmeancf>cen
      
      # A) USE AVERAGE DAILY MORTALITY - "MODEL" 
      ancf <- (1-exp(-bvarcencf%*%coef))*deathproj
      # SET NEGATIVE ESTIMATES TO 0
      ancf[ancf < 0] <- 0
      # SUM PER YEAR 
      heatanvar[i,as.character(ysub),"sim",c+1,"est"] <- tapply(ancf[indcf],datasub$year[indcf],sum,na.rm=T)
      
      # B) USE OBSERVED DAILY MORTALITY - "OBSERVATION"
      anobscf <- (1-exp(-bvarcencf%*%coef))*avgdeathsub
      # SET NEGATIVE ESTIMATES TO 0
      anobscf[anobscf < 0] <- 0
      # SUM PER YEAR 
      heatanvar[i,as.character(ysub),"obs",c+1,"est"] <- tapply(anobscf[indcf],datasub$year[indcf],sum,na.rm=T)
    
      # COMPUTE EMPIRICAL CONFIDENCE INTERVALS
      # LOOP ACROSS ITERATIONS
      for(s in seq(nsim)) {
        
        ancf <- (1-exp(-bvarcencf%*%coefsim[s,]))*deathproj
        ancf[ancf < 0] <- 0
        arraysimvar[i,as.character(ysub),"sim",c+1,s] <- tapply(ancf[indcf],datasub$year[indcf],sum,na.rm=T)
        
        anobscf <- (1-exp(-bvarcencf%*%coefsim[s,]))*avgdeathsub
        anobscf[anobscf < 0] <- 0
        arraysimvar[i,as.character(ysub),"obs",c+1,s] <- tapply(anobscf[indcf],datasub$year[indcf],sum,na.rm=T)
        
      }
    }
  }
 
  # SUM MOVING AVERAGE PER YEAR
  totavgdeath <- tapply(avgdeathsum,data$year,sum,na.rm=T)
  
  # COMPUTE YEARLY RESCALING FACTOR TO ALIGN FORWARD AVERAGED DEATH WITH TOTAL ANNUAL DEATHS
  rf <- deathtot[i,]/totavgdeath

  # RESCALE
  heatanvar[i,,"obs","fact","est"] <- heatanvar[i,,"obs","fact","est"]*rf
  arraysimvar[i,,"obs","fact",] <- arraysimvar[i,,"obs","fact",]*matrix(rep(rf,each=nsim),ncol=nsim,byrow=T)
  
  # LOOP OVER COUNTERFACTUAL TEMPERATURE VERSIONS
  for (c in seq(tmeancfvar)){
    heatanvar[i,,"obs",c+1,"est"] <- heatanvar[i,,"obs",c+1,"est"]*rf
    arraysimvar[i,,"obs",c+1,] <- arraysimvar[i,,"obs",c+1,]*matrix(rep(rf,each=nsim),ncol=nsim,byrow=T)
  }
}

#########################################################################
# COMPUTE DIFFERENCES FACTUAL VERSUS COUNTERFACTUAL ESTIMATES
########################################################################

# LOOP OVER COUNTERFACTUAL TEMPERATURE VERSIONS
for (c in seq(tmeancfvar)){
  # ABSOLUTE DIFFERENCES
  heatanvardif[,,,"abs",c,"est"] <- heatanvar[,,,"fact","est"] - heatanvar[,,,c+1,"est"]
  arraysimvardif[,,,"abs",c,] <- arraysimvar[,,,"fact",] - arraysimvar[,,,c+1,]
  
  # RELATIVE DIFFERENCES
  heatanvardif[,,,"rel",c,"est"] <- (heatanvar[,,,"fact","est"] - heatanvar[,,,c+1,"est"])/heatanvar[,,,"fact","est"]*100
  arraysimvardif[,,,"rel",c,] <- (arraysimvar[,,,"fact",] - arraysimvar[,,,c+1,])/arraysimvar[,,,"fact",]*100
}

###########################################################################
# COMPUTE 95% EMPIRICAL CONFIDENCE INTERVALS AND ATTRIBUTABLE FRACTIONS
###########################################################################

#################################
# ATTRIBUTABLE NUMBERS BY CITY 
#################################

# BY YEAR
heatanvar[,,,,"ci.l"] <- apply(arraysimvar,1:4,quantile,0.025,na.rm=T)
heatanvar[,,,,"ci.u"] <- apply(arraysimvar,1:4,quantile,0.975,na.rm=T)

heatanvardif[,,,,,"ci.l"] <- apply(arraysimvardif,1:5,quantile,0.025,na.rm=T)
heatanvardif[,,,,,"ci.u"] <- apply(arraysimvardif,1:5,quantile,0.975,na.rm=T)

# SUM/AVERAGE FOR FULL STUDY PERIOD 
meanheatanvardif[,,"abs",,"est"] <- apply(heatanvardif[,,,"abs",,"est"],c(1,3,4),sum,na.rm=T)
meanheatanvardif[,,"abs",,"ci.l"] <- apply(apply(arraysimvardif[,,,"abs",,],c(1,3:5),sum,na.rm=T),1:3,quantile,0.025,na.rm=T)
meanheatanvardif[,,"abs",,"ci.u"] <- apply(apply(arraysimvardif[,,,"abs",,],c(1,3:5),sum,na.rm=T),1:3,quantile,0.975,na.rm=T)

# LOOP OVER COUNTERFACTUAL TEMPERATURE VERSIONS
for (c in seq(tmeancfvar)){
  meanheatanvardif[,,"rel",c,"est"] <- (meanheatanvardif[,,"abs",c,"est"]/apply(heatanvar[,,,"fact","est"],c(1,3),sum,na.rm=T))*100
  meanheatanvardif[,,"rel",c,"ci.l"] <- (apply(apply(arraysimvardif[,,,"abs",c,],c(1,3,4),sum,na.rm=T)/apply(arraysimvar[,,,"fact",],c(1,3,4),sum,na.rm=T),1:2,quantile,0.025,na.rm=T))*100
  meanheatanvardif[,,"rel",c,"ci.u"] <- (apply(apply(arraysimvardif[,,,"abs",c,],c(1,3,4),sum,na.rm=T)/apply(arraysimvar[,,,"fact",],c(1,3,4),sum,na.rm=T),1:2,quantile,0.975,na.rm=T))*100
}

#####################################################
# DIVIDE ATTRIBUTABLE NUMBERS BY TOTAL ANNUAL DEATHS
#####################################################

# CREATE MATRIX OF TOTAL MORTALITY PER SUMMER FROM AVERAGED DEATHS
deathtotproj <- t(sapply(deathtotavg,function(x) rep(x,length(ny))))
colnames(deathtotproj) <- ny

for (j in seq(3)){
  heatafvar[,,"obs","fact",j] <- heatanvar[,,"obs","fact",j]/deathtot*100
  heatafvar[,,"sim","fact",j] <- heatanvar[,,"sim","fact",j]/deathtotproj*100
  
  for (c in seq(tmeancfvar)){
    heatafvar[,,"obs",c+1,j] <- heatanvar[,,"obs",c+1,j]/deathtot*100
    heatafvar[,,"sim",c+1,j] <- heatanvar[,,"sim",c+1,j]/deathtotproj*100
    
    heatafvardif[,,"obs","abs",c,j] <- heatanvardif[,,"obs","abs",c,j]/deathtot*100
    heatafvardif[,,"sim","abs",c,j] <- heatanvardif[,,"sim","abs",c,j]/deathtotproj*100
  }
}

meanheatafvardif[,"obs","abs",,] <- meanheatanvardif[,"obs","abs",,]/rowSums(deathtot)*100
meanheatafvardif[,"sim","abs",,] <- meanheatanvardif[,"sim","abs",,]/rowSums(deathtotproj)*100

# RELATIVE DIFFERENCES IN ATTRIBUTABLE FRACTION 
# EQUAL TO RELATIVE DIFFERENCES IN ATTRIBUTABLE NUMBERS
heatafvardif[,,,"rel",,] <- heatanvardif[,,,"rel",,]
meanheatafvardif[,,"rel",,] <- meanheatanvardif[,,"rel",,]

#######################
# SUMMED ACROSS CITIES
#######################

##########
# BY YEAR
##########

# SUM ABSOLUTE ESTIMATES
heatantotvar[,,,"est"] <- apply(heatanvar[,,,,"est"],2:4,sum,na.rm=F)
heatantotvar[,,,"ci.l"] <- apply(apply(arraysimvar,2:5,sum,na.rm=F),1:3,quantile,0.025,na.rm=T)
heatantotvar[,,,"ci.u"] <- apply(apply(arraysimvar,2:5,sum,na.rm=F),1:3,quantile,0.975,na.rm=T)

# SUM ABSOLUTE DIFFERENCES
heatantotvardif[,,"abs",,"est"] <- apply(heatanvardif[,,,"abs",,"est"],2:4,sum,na.rm=F)
heatantotvardif[,,"abs",,"ci.l"] <- apply(apply(arraysimvardif[,,,"abs",,],2:5,sum,na.rm=F),1:3,quantile,0.025,na.rm=T)
heatantotvardif[,,"abs",,"ci.u"] <- apply(apply(arraysimvardif[,,,"abs",,],2:5,sum,na.rm=F),1:3,quantile,0.975,na.rm=T)

# COMPUTE RELATIVE DIFFERENCES
for (c in seq(tmeancfvar)){
  heatantotvardif[,,"rel",c,"est"] <- (heatantotvar[,,"fact","est"] - heatantotvar[,,c+1,"est"])/heatantotvar[,,"fact","est"]*100
  heatantotvardif[,,"rel",c,"ci.l"] <- apply(apply(arraysimvardif[,,,"abs",c,],2:4,sum,na.rm=F)/apply(arraysimvar[,,,"fact",],2:4,sum,na.rm=F),1:2,quantile,0.025,na.rm=T)*100
  heatantotvardif[,,"rel",c,"ci.u"] <- apply(apply(arraysimvardif[,,,"abs",c,],2:4,sum,na.rm=F)/apply(arraysimvar[,,,"fact",],2:4,sum,na.rm=F),1:2,quantile,0.975,na.rm=T)*100
}

#########################
# FOR FULL STUDY PERIOD
#########################

# SUM ABSOLUTE DIFFERENCES
meanheatantotvardif[,"abs",,"est"] <- apply(heatantotvardif[,,"abs",,"est"],2:3,sum,na.rm=T)
meanheatantotvardif[,"abs",,"ci.l"] <- apply(apply(arraysimvardif[,,,"abs",,],3:5,sum,na.rm=T),1:2,quantile,0.025,na.rm=T)
meanheatantotvardif[,"abs",,"ci.u"] <- apply(apply(arraysimvardif[,,,"abs",,],3:5,sum,na.rm=T),1:2,quantile,0.975,na.rm=T)

# COMPUTE RELATIVE DIFFERENCES
meanheatantotvardif[,"rel",,"est"] <- meanheatantotvardif[,"abs",,"est"]/apply(heatantotvar[,,"fact","est"],2,sum,na.rm=T)*100
for (c in seq(tmeancfvar)){
  meanheatantotvardif[,"rel",c,"ci.l"] <- apply(apply(arraysimvardif[,,,"abs",c,],3:4,sum,na.rm=T)/apply(arraysimvar[,,,"fact",],c(3,4),sum,na.rm=T),1,quantile,0.025,na.rm=T)*100
  meanheatantotvardif[,"rel",c,"ci.u"] <- apply(apply(arraysimvardif[,,,"abs",c,],3:4,sum,na.rm=T)/apply(arraysimvar[,,,"fact",],c(3,4),sum,na.rm=T),1,quantile,0.975,na.rm=T)*100
}

#######################################################
# DIVIDE ATTRIBUTABLE NUMBERS BY TOTAL ANNUAL DEATHS
#######################################################

# BASED ON OBSERVED MORTALITY
totdeathtot <- colSums(deathtot)
heataftotvar[,"obs",,] <- heatantotvar[,"obs",,]/totdeathtot*100
heataftotvardif[,"obs","abs",,] <- heatantotvardif[,"obs","abs",,]/totdeathtot*100 
meanheataftotvardif["obs","abs",,] <- meanheatantotvardif["obs","abs",,]/sum(totdeathtot)*100
# BASED ON AVERAGED MORTALITY
totdeathtotproj <- colSums(deathtotproj)
heataftotvar[,"sim",,] <- heatantotvar[,"sim",,]/totdeathtotproj*100
heataftotvardif[,"sim","abs",,] <- heatantotvardif[,"sim","abs",,]/totdeathtotproj*100
meanheataftotvardif["sim","abs",,] <- meanheatantotvardif["sim","abs",,]/sum(totdeathtotproj)*100

# RELATIVE DIFFERENCES IN ATTRIBUTABLE FRACTION 
# EQUAL TO RELATIVE DIFFERENCES IN ATTRIBUTABLE NUMBERS
heataftotvardif[,,"rel",,] <- heatantotvardif[,,"rel",,]
meanheataftotvardif[,"rel",,] <- meanheatantotvardif[,"rel",,]
####
