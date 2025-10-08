#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

################################################################################
# IMPORT AND PREPARE DATA FOR 15 MAJOR GERMAN CITIES
# 1993-2022
################################################################################

####################################################################################
# IMPORT MORTALITY, STATION CLIMATE DATA, AND ALERT DATA AND CONVERT TO LIST
####################################################################################

dat <- read.csv("MortalitySeries/MCC_Germany_1993_2022_onlyClimate.csv")

# CREATE CITY VECTOR
cities <- unique(dat$cityname)

# ADD DATE VARIABLES
dat$date <- as.Date(dat$date,"%m/%d/%Y")
dat$year <- as.integer(format(dat$date,"%Y"))
dat$month <- as.integer(format(dat$date,"%m"))
dat$day <- as.integer(format(dat$date,"%d"))
dat$yday <- as.POSIXlt(dat$date)$yday+1
dat$dow <- as.factor(weekdays(dat$date))

# MERGE WITH ALERT DATA
alertdat <- read.csv("MetaPredictors/GermanCitiesAlertDays1993To2022.csv")
dat$alert <- alertdat$hw.level
dat$eligible <- alertdat$hw.predict

# CONVERT TO LIST
dlist <- split(dat,as.factor(dat$cityname))
  
# REMOVE ORIGINAL INPUT DATA
rm(dat)
rm(alertdat)

####################################################
# ADD ERA5 DATA - FACTUAL AND COUNTERFACTUAL
####################################################

firstdate <-  "1993-01-01"
lastdate <- "2022-12-31"

# LOOP OVER CITIES
for (i in seq(dlist)) {
  
  # IMPORT ERA5 TEMPERATURE DATA
  tdat <- read.csv(paste0("TemperatureSeries/Detrended/TmeanERA5land",cities[i],"Counterfactual.csv"))
  tdat$date <- as.Date(tdat$date)
  
  # SELECT TIME PERIOD WITH MORTALITY COVERAGE
  tdat <- tdat[tdat$date >= firstdate & tdat$date <= lastdate,]
  
  # ADD ERA5-TEMPERATURE DATA TO DLIST

  dlist[[i]] <- cbind(dlist[[i]],tdat[,c("tmeanf",paste0("tmeancf",1:9))])

}

########################################
# PLOT WARM-SEASON MORTALITY BY YEAR
########################################

sm <- 6:9

pdf("WarmSeasonMortalityCity.pdf",width=8,height=9)
layout(matrix(1:15,ncol=3,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(3,1,0),las=1)

# LOOP OVER CITIES
for (i in seq(dlist)) {
  
  data <- dlist[[i]]
  subdata <- subset(data,data$month %in% sm) 
  
  mort <- tapply(subdata$death,subdata$year,mean,na.rm=T)
  
  barplot(mort,ylab="Warm season death count",xlab="",main=cities[i],ylim=c(0,max(mort)+2))
  
}

dev.off()

#############
