#####################################################################################
# R code for the analysis in 

# Huber, V., Breitner-Busch, S., Feldbusch, H. et al. Improvements in life expectancy 
# mask rising trends in heat-related excess mortality attributable to climate change. 
# Nat Commun 16, 11632 (2025). https://doi.org/10.1038/s41467-025-66681-0
######################################################################################

##############################################################################################
# COMPUTE ANNUAL HEAT-RELATED MORTALITY BASED ON FACTUAL AND COUNTERFACTUAL TEMPERATURE DATA

# USE COUNTERFACTUAL BLUPs FROM LONGITUDINAL META-REGRESSION-MODEL ('W/O  LE IMPROVEMENTS')
##############################################################################################

#####################################################
# COMPUTE ANNUAL HEAT-ATTRIBUTABLE NUMBERS/FRACTIONS
# FOR FACTUAL AND COUNTERFACTUAL DATA

# USING OBSERVED OR AVERAGED DAILY MORTALITY DATA
#####################################################

###################################
# CREATE OBJECTS TO STORE RESULTS
###################################

# NUMBERS OF DAYS PER WARM SEASON
ndwarm <- length(subset(dlist[[1]], month %in% sm & year==2020)$date)

# HEAT ATTRIBUTABLE NUMBERS AND FRACTIONS BY CITY
heatannoadapt <- heatafnoadapt <- array(NA,dim=c(length(cities),length(ny),2,10,3),
                          dimnames=list(cities,ny,c("obs","sim"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# CLIMATE CHANGE ATTRIBUTABLE NUMBERS AND FRACTIONS BY CITY
heatannoadaptdif <- heatafnoadaptdif <- array(NA,dim=c(length(cities),length(ny),2,2,9,3),
                                dimnames=list(cities,ny,c("obs","sim"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# HEAT ATTRIBUTABLE NUMBERS AND FRACTIONS FOR ALL CITIES
heatantotnoadapt <- heataftotnoadapt <- array(NA,dim=c(length(ny),2,10,3),
                          dimnames=list(ny,c("obs","sim"),c("fact","cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

# CLIMATE CHANGE ATTRIBUTABLE NUMBERS AND FRACTIONS FOR ALL CITIES
heatantotnoadaptdif <- heataftotnoadaptdif <- array(NA,dim=c(length(ny),2,2,9,3),
                                      dimnames=list(ny,c("obs","sim"),c("abs","rel"),c("cfact.best",paste0("cfact.ci.",1:8)),c("est","ci.l","ci.u")))

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

    # COUNTERFACTUAL BLUPs
    coef <- blupcf[[paste(cities[i],period[p])]]$blup
    vcov <- blupcf[[paste(cities[i],period[p])]]$vcov
    
    cen <- mmtcitycf[i,p]
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
    an[an < 0] <- 0
    heatannoadapt[i,as.character(ysub),"sim","fact","est"] <- tapply(an[ind],datasub$year[ind],sum,na.rm=T)
    
    # B) USE OBSERVED DAILY MORTALITY - "OBSERVATION"
    # GET FORWARD AVERAGED DEATH SUB-PERIOD
    avgdeathsub <- avgdeathsum[data$year %in% ysub]
    # (1) AND FACTUAL TEMPERATURE DATA
    anobs <- (1-exp(-bvarcen%*%coef))*avgdeathsub
    anobs[anobs < 0] <- 0
    heatannoadapt[i,as.character(ysub),"obs","fact","est"] <- tapply(anobs[ind],datasub$year[ind],sum,na.rm=T)
    
    set.seed(13041975)
    coefsim <- mvrnorm(nsim,coef,vcov)
    
    # - LOOP ACROSS ITERATIONS
    for(s in seq(nsim)) {
      
      # USE AVERAGE DAILY MORTALITY - "MODEL" 
      an <- (1-exp(-bvarcen%*%coefsim[s,]))*deathproj
      an[an < 0] <- 0
      arraysimvar[i,as.character(ysub),"sim","fact",s] <- tapply(an[ind],datasub$year[ind],sum,na.rm=T)

      # USE OBSERVED DAILY MORTALITY - "OBSERVATION"
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

      # USE AVERAGE DAILY MORTALITY - "MODEL" 
      ancf <- (1-exp(-bvarcencf%*%coef))*deathproj
      ancf[ancf < 0] <- 0
      heatannoadapt[i,as.character(ysub),"sim",c+1,"est"] <- tapply(ancf[indcf],datasub$year[indcf],sum,na.rm=T)
      
      # USE OBSERVED DAILY MORTALITY - "OBSERVATION"
      anobscf <- (1-exp(-bvarcencf%*%coef))*avgdeathsub
      anobscf[anobscf < 0] <- 0
      heatannoadapt[i,as.character(ysub),"obs",c+1,"est"] <- tapply(anobscf[indcf],datasub$year[indcf],sum,na.rm=T)
      
      # LOOP ACROSS MONTE CARLO ITERATIONS 
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

  # RESCALE
  heatannoadapt[i,,"obs","fact","est"] <- heatannoadapt[i,,"obs","fact","est"]*rf
  arraysimvar[i,,"obs","fact",] <- arraysimvar[i,,"obs","fact",]*matrix(rep(rf,each=nsim),ncol=nsim,byrow=T)
  
  # LOOP OVER COUNTERFACTUAL TEMPERATURE VERSIONS
  for (c in seq(tmeancfvar)){
    heatannoadapt[i,,"obs",c+1,"est"] <- heatannoadapt[i,,"obs",c+1,"est"]*rf
    arraysimvar[i,,"obs",c+1,] <- arraysimvar[i,,"obs",c+1,]*matrix(rep(rf,each=nsim),ncol=nsim,byrow=T)
  }
}

#########################################################################
# COMPUTE DIFFERENCES FACTUAL VERSUS COUNTERFACTUAL ESTIMATES

# LOOP OVER COUNTERFACTUAL TEMPERATURE VERSIONS
for (c in seq(tmeancfvar)){
  
  # ABSOLUTE DIFFERENCES
  heatannoadaptdif[,,,"abs",c,"est"] <- heatannoadapt[,,,"fact","est"] - heatannoadapt[,,,c+1,"est"]
  arraysimvardif[,,,"abs",c,] <- arraysimvar[,,,"fact",] - arraysimvar[,,,c+1,]
  
  # RELATIVE DIFFERENCES
  heatannoadaptdif[,,,"rel",c,"est"] <- (heatannoadapt[,,,"fact","est"] - heatannoadapt[,,,c+1,"est"])/heatannoadapt[,,,"fact","est"]*100
  arraysimvardif[,,,"rel",c,] <- (arraysimvar[,,,"fact",] - arraysimvar[,,,c+1,])/arraysimvar[,,,"fact",]*100

}

###########################################################################
# COMPUTE 95% EMPIRICAL CONFIDENCE INTERVALS AND ATTRIBUTABLE FRACTIONS
###########################################################################

###########
# BY CITY
###########

heatannoadapt[,,,,"ci.l"] <- apply(arraysimvar,1:4,quantile,0.025,na.rm=T)
heatannoadapt[,,,,"ci.u"] <- apply(arraysimvar,1:4,quantile,0.975,na.rm=T)

heatannoadaptdif[,,,,,"ci.l"] <- apply(arraysimvardif,1:5,quantile,0.025,na.rm=T)
heatannoadaptdif[,,,,,"ci.u"] <- apply(arraysimvardif,1:5,quantile,0.975,na.rm=T)

# DIVIDE ATTRIBUTABLE NUMBERS BY TOTAL ANNUAL DEATHS

# CREATE MATRIX OF TOTAL MORTALITY PER SUMMER FROM AVERAGED DEATHS
deathtotproj <- t(sapply(deathtotavg,function(x) rep(x,length(ny))))
colnames(deathtotproj) <- ny

for (j in seq(3)){
  heatafnoadapt[,,"obs","fact",j] <- heatannoadapt[,,"obs","fact",j]/deathtot*100
  heatafnoadapt[,,"sim","fact",j] <- heatannoadapt[,,"sim","fact",j]/deathtotproj*100
  
  for (c in seq(tmeancfvar)){
    heatafnoadapt[,,"obs",c+1,j] <- heatannoadapt[,,"obs",c+1,j]/deathtot*100
    heatafnoadapt[,,"sim",c+1,j] <- heatannoadapt[,,"sim",c+1,j]/deathtotproj*100
    
    heatafnoadaptdif[,,"obs","abs",c,j] <- heatannoadaptdif[,,"obs","abs",c,j]/deathtot*100
    heatafnoadaptdif[,,"sim","abs",c,j] <- heatannoadaptdif[,,"sim","abs",c,j]/deathtotproj*100
  }
}

# RELATIVE DIFFERENCES IN ATTRIBUTABLE FRACTION 
# EQUAL TO RELATIVE DIFFERENCES IN ATTRIBUTABLE NUMBERS
heatafnoadaptdif[,,,"rel",,] <- heatannoadaptdif[,,,"rel",,]

#######################
# SUMMED ACROSS CITIES
#######################

# SUM ABSOLUTE ESTIMATES
heatantotnoadapt[,,,"est"] <- apply(heatannoadapt[,,,,"est"],2:4,sum,na.rm=F)
heatantotnoadapt[,,,"ci.l"] <- apply(apply(arraysimvar,2:5,sum,na.rm=F),1:3,quantile,0.025,na.rm=T)
heatantotnoadapt[,,,"ci.u"] <- apply(apply(arraysimvar,2:5,sum,na.rm=F),1:3,quantile,0.975,na.rm=T)

# SUM ABSOLUTE DIFFERENCES
heatantotnoadaptdif[,,"abs",,"est"] <- apply(heatannoadaptdif[,,,"abs",,"est"],2:4,sum,na.rm=F)
heatantotnoadaptdif[,,"abs",,"ci.l"] <- apply(apply(arraysimvardif[,,,"abs",,],2:5,sum,na.rm=F),1:3,quantile,0.025,na.rm=T)
heatantotnoadaptdif[,,"abs",,"ci.u"] <- apply(apply(arraysimvardif[,,,"abs",,],2:5,sum,na.rm=F),1:3,quantile,0.975,na.rm=T)

# COMPUTE RELATIVE DIFFERENCES
for (c in seq(tmeancfvar)){
  heatantotnoadaptdif[,,"rel",c,"est"] <- (heatantotnoadapt[,,"fact","est"] - heatantotnoadapt[,,c+1,"est"])/heatantotnoadapt[,,"fact","est"]*100
  heatantotnoadaptdif[,,"rel",c,"ci.l"] <- apply(apply(arraysimvardif[,,,"abs",c,],2:4,sum,na.rm=F)/apply(arraysimvar[,,,"fact",],2:4,sum,na.rm=F),1:2,quantile,0.025,na.rm=T)*100
  heatantotnoadaptdif[,,"rel",c,"ci.u"] <- apply(apply(arraysimvardif[,,,"abs",c,],2:4,sum,na.rm=F)/apply(arraysimvar[,,,"fact",],2:4,sum,na.rm=F),1:2,quantile,0.975,na.rm=T)*100
}

# DIVIDE ATTRIBUTABLE NUMBERS BY TOTAL ANNUAL DEATHS
# BASED ON OBSERVED MORTALITY
totdeathtot <- colSums(deathtot)
heataftotnoadapt[,"obs",,] <- heatantotnoadapt[,"obs",,]/totdeathtot*100
heataftotnoadaptdif[,"obs","abs",,] <- heatantotnoadaptdif[,"obs","abs",,]/totdeathtot*100 
# BASED ON AVERAGED MORTALITY
totdeathtotproj <- colSums(deathtotproj)
heataftotnoadapt[,"sim",,] <- heatantotnoadapt[,"sim",,]/totdeathtotproj*100
heataftotnoadaptdif[,"sim","abs",,] <- heatantotnoadaptdif[,"sim","abs",,]/totdeathtotproj*100

# RELATIVE DIFFERENCES IN ATTRIBUTABLE FRACTION 
# EQUAL TO RELATIVE DIFFERENCES IN ATTRIBUTABLE NUMBERS
heataftotnoadaptdif[,,"rel",,] <- heatantotnoadaptdif[,,"rel",,]

####

