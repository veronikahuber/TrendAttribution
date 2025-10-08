#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

###########################################################################
# COMPUTE FIRST STAGE COEFFICIENTS IN DIFFERENT TIME PERIODS
# AVERAGE META-PREDICTORS BY PERIOD FOR SECOND STAGE

# MODEL-SPECIFICATION AS IN VICEDO-CABRERA ET AL. 2021
###########################################################################

# LOAD THE PACKAGES
library(dlnm) ; library(mixmeta) ; library(splines) ; library(tsModel) ; library(MASS);
library(scales)

####################################################
# IMPORT ANNUAL META-PREDICTORS FROM INKAR DATABASE
####################################################

metadat <- read.csv("MetaPredictors/INKAR_extrapolated_1993_2022.csv")

# ADD DATA ON % UNEMPLOYED
metadat$perc.unemployed <- metadat$unemployed/metadat$population.corrected*100

#####################################
# RUN FIRST-STAGE MODELS BY PERIODS
#####################################

# DEFINE THE MODEL

# DEFINE WARM-SEASON MONTHS
sm <- 6:9

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"
vardegree <- NULL
varper <- c(50, 90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 10
lagnk <- 2

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 4

# DEGREE OF FREEDOM FOR TREND
dftrend <- 1

# DEFINE THE PERIODS AND
yearlist <- list(1993:1997,1998:2002,2003:2007,2008:2012,2013:2017,2018:2022)

# DEFINE DUMMY VARIABLE TO INDICATE PRE-/POST-IMPLEMENTATION OF HEAT EARLY WARNING SYSTEMS
hpp <- c(0,0,1,1,1,1)

# SENSITIVITY ANALYSIS: ALTERNATIVE PERIOD DEFINITIONS 
#yearlist <- list(1993:1998,1999:2004,2005:2010,2011:2016,2017:2022)
#hpp <- c(0,0,1,1,1)

# SENSITIVITY ANALYSIS: EXCLUDE THE YEAR 2003 
#yearlist <- list(1993:1997,1998:2002,2004:2007,2008:2012,2013:2017,2018:2022)

period <- sapply(yearlist, function(x) paste(range(x), collapse="-"))

# CREATE OBJECTS TO STORE THE RESULTS
tmeanperparlist <- tmeansumlist <- vector("list", length(cities))

# DEFINE PERCENTILES FOR AVERAGE DISTRIBUTION
per <- c(0:9/10, 1:99, 991:1000/10)/100

# CREATE OBJECT TO STORE KNOT LOCATIONS BY CITY
knots <- array(NA,dim=c(length(cities),length(varper)),
               dimnames=list(cities,paste0(varper,"%")))

# RUN THE LOOP
for(i in seq(length(cities))) {
  
  # PRINT CITY
  cat(cities[i],"")
  
  # GET DATA
  data <- dlist[[i]]
  
  # GET CITY-SPECIFIC META-PREDICTOR DATA
  metacity <- metadat[metadat$cityname==cities[i],
                      c("year","population.corrected","women.quota",
                        "residents.over.65","average.population.age",
                        "life.expectancy","rest.life.expectancy","foreign.residents",
                        "unemployed","household.income","population.density",
                        "perc.unemployed")]
  
  # SUBSET FOR SUMMER-ONLY
  datasum <- subset(data, month %in% sm)
  
  # DEFINE CROSSBASIS BASED ON FULL PERIOD TEMPERATURE DISTRIBUTION 
  knots[i,] <- quantile(datasum$tmeanf, varper/100, na.rm=T)
  
  # DEFINE CROSS-BASIS FOR TEMPERATURE
  if(!is.null(vardegree)) {
    argvar <- list(fun = varfun, degree=vardegree, knots = knots[i,])
  }else{
    argvar <- list(fun = varfun, knots = knots[i,])
  }
  arglag <- list(knots = logknots(lag, lagnk))
  cbtmean <- crossbasis(datasum$tmeanf, lag = lag, argvar = argvar, arglag = arglag,
                        group = datasum$year)
  
  # PERFORM MODEL BY PERIOD
  parlist <- lapply(yearlist, function(ysub) {
    model <- glm(death ~ cbtmean + dow + ns(yday, df = dfseas):factor(year) +
                   ns(date, df = round(length(unique(year)) / dftrend / 10)),
                 data=datasum, family=quasipoisson, subset=year %in% ysub)
    redpred <- crossreduce(cbtmean, model, cen=mean(datasum$tmeanf, na.rm=T))
    t(c(coef(redpred), vechMat(vcov(redpred))))
  })
  
  # DERIVE META-PREDICTORS OF TEMPERATURE
  avgtmeanper <- sapply(yearlist,function(ysub) mean(subset(data,year %in% ysub)$tmeanf,na.rm=T))
  rangetmeanper <- sapply(yearlist,function(ysub) diff(range(subset(data,year %in% ysub)$tmeanf,na.rm=T)))
  mstper <- sapply(yearlist,function(ysub) mean(subset(datasum,year %in% ysub)$tmeanf,na.rm=T))
  
  # COMPUTE NUMBER OF POTENTIAL HEAT ALERT DAYS PER PERIOD
  heatdayper <- sapply(yearlist,function(ysub) sum(subset(data,year %in% ysub)$eligible=="yes",na.rm=T))
  
  # AVERAGE SOCIO-ECONOMIC META-PREDICTORS BY PERIOD
  metacityper <- t(sapply(yearlist,function(ysub) colMeans(metacity[metacity$year %in% ysub,2:length(metacity)],na.rm=T)))

  # STORE PARAMETERS (COEF + VECTORIZED VCOV)
  par <- do.call(rbind, parlist)
  tmeanperpar <- data.frame(
    cities[i],
    period = period,
    year = sapply(yearlist, mean),
    avgtmeanper,
    rangetmeanper,
    mstper,
    heatdayper,
    metacityper,
    hpp,
    par,
    row.names=paste0(i, ".", seq(yearlist))
  )
  colnames(tmeanperpar)[1] <- "city"
  tmeanperparlist[[i]] <- tmeanperpar
  
  # CREATE COUNTRY-AVERAGE SUMMER TEMPERATURE DISTRIBUTION
  tmeansumlist[[i]] <- quantile(datasum$tmeanf, per, na.rm=T)
}

# RBIND COEF/VCOV TOGETHER IN DATAFRAMES
tmeanperpar <- do.call(rbind, tmeanperparlist)

# TEMPERATURE DISTRIBUTION (SUMMER ONLY)
avgtmeansum <- data.frame(perc=names(tmeansumlist[[1]]), 
                          tmean=apply(do.call(cbind, tmeansumlist), 1, mean))

# CREATE COUNTRY-AVERAGE TEMPERATURE TIME SERIES 
tmeancountry <- data.frame(date=data$date,year=data$year,month=data$month,tmeancountry=rowMeans(sapply(dlist,function(x) x$tmeanf),na.rm=T))

# COUNTRY-AVERAGE RANGE BY PERIOD
rangetmeancountryper <- sapply(yearlist,function(ysub) diff(range(subset(tmeancountry,year %in% ysub)$tmeancountry,na.rm=T)))

# COUNTRY-AVERAGE RANGE BY YEAR
rangetmeancountryyear <- sapply(1993:2022,function(ysub) diff(range(subset(tmeancountry,year %in% ysub)$tmeancountry,na.rm=T)))

# RESTRICT TO SUMMER MONTHS
tmeancountrysum <- subset(tmeancountry,month %in% sm)

####


