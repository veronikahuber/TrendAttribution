#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

##################################################################################
# COMPUTE COUNTERFACTUAL DAILY MEAN TEMPERATURES 
# BASED ON DETRENDING AGAINST 5-YEAR MOVING AVERAGE OF GLOBAL MEAN SURFACE TEMPERATURES

# DATA USED IN SENSITIVITY ANALYSIS
###################################################################################

####################################
# 1: IMPORT AND SMOOTH GMST DATA
####################################

# IMPORT MONTHLY SERIES OF HadCRUT5 GLOBAL MEAN SURFACE TEMPERATURE (GMST) (ANOMALY AGAINST 1961-1990)

gmst <- read.csv("TemperatureSeries/HadCRUT5/HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv")
gmst$Year <- as.integer(substr(gmst$Time,1,4))
gmst$Month <- as.integer(substr(gmst$Time,6,7))

# RE-SCALE AS DIFFERENCE TO PRE-INDUSTRIAL 1850-1900
preind <- colMeans(subset(gmst[,2:4],gmst$Year %in% 1850:1900))
gmst[,2:4] <- gmst[,2:4] - preind

# COMPUTE YEARLY AVERAGES
gmsty <- sapply(gmst[,2:4],function(x) tapply(x,gmst$Year,mean,na.rm=T))

gmsty <- data.frame(gmsty)
gmsty$Year <- unique(gmst$Year)
nby <- length(gmsty$Year)

# SMOOTH WITH 5-YEAR RUNNING MEAN
lag <- 4
gmsty4s <- matrix(NA,nrow=nby,ncol=3)

for (i in 1:3){
  
  # FORWARD MOVING AVERAGE
  avgt <- rowMeans(as.matrix(tsModel:::Lag(gmsty[,i],-seq(0,lag))))
  # RE-CENTRE WINDOWS
  gmsty4s[3:nby,i] <-avgt[1:(nby-2)]
  
  # FILL BEGINNING AND END-POINTS
  gmsty4s[1,i] <- mean(gmsty[1:(lag-1),i])
  gmsty4s[2,i] <- mean(gmsty[1:lag,i])
  gmsty4s[nby-1,i] <- mean(gmsty[(nby-lag+1):nby,i])
  gmsty4s[nby,i] <- mean(gmsty[(nby-lag+2):nby,i])
}

# FOR CHECKING: IMPORT DATA SMOOTHED WITH SINGULAR SPECTRUM ANALYSIS
dat <- read.csv("TemperatureSeries/HadCRUT5/HadCRUT.5.0.2.0.analysis.summary_series.global.ssa_smoothed.csv")
dat$Year <- as.integer(substr(dat$Time,1,4))

# RE-SCALE AS DIFFERENCE TO PRE-INDUSTRIAL 1850-1900 AND COMPUTE YEARLY AVERAGES
gmstssa <- dat[,5] - mean(subset(dat[,5],dat$Year %in% 1850:1900))
gmstssay <- tapply(gmstssa,dat$Year,mean,na.rm=T)

# PLOT GMST (RAW, 5-YEAR RUNNING MEAN, SSA SMOOTHED)

pdf("FigS10_gmsty.pdf",width=8,height=7)

plot(1850:2024,gmsty[gmsty$Year<=2024,1],xlab="Year",
     ylab="GMST anomaly (°C)", type="n",ylim=c(min(gmsty[,2]),max(gmsty[,3])))
polygon(c(1850:2024,rev(1850:2024)),
        c(gmsty[gmsty$Year<=2024,2],rev(gmsty[gmsty$Year<=2024,3])),
          col="lightblue",border=NA)
lines(1850:2024,gmsty[gmsty$Year<=2024,1],col="blue")
lines(1850:2024,gmsty4s[gmsty$Year<=2024,1],col="red",lty=1)
#lines(unique(dat$Year),gmstssay,col="black",lty=2)
legend("topleft",c("HadCRUT5","5-year running mean"),col=c("blue","red"),lty=1,bty="n")

dev.off()

#######################################################################
# 2: REGRESS LOCAL WARM-SEASON TEMPERATURE AGAINST SMOOTHED GMST AND
# 3: COMPUTE COUNTERFACTUAL TEMPERATURES SERIES
#######################################################################

cities <- c("Berlin",
               "Bremen",
               "Cologne",
               "Dortmund",
               "Dresden",
               "Duisburg",
               "Dusseldorf",
               "Essen",
               "Frankfurt",
               "Hamburg",
               "Hannover",
               "Leipzig",
               "Munich",
               "Nuremberg",
               "Stuttgart")

# DEFINE WARM-SEASON MONTHS
sm <- 6:9

# DEFINE STUDY PERIOD
stp <- 1950:2022

# SUBSET SMOOTHED GMST DATA TO STUDY PERIOD
gmstsub <- gmsty4s[gmsty$Year %in% stp,]

#############################
# PREPARE DAILY GMST VECTORS
#############################

# A: VERY SIMPLE APPROACH: USING 5-YEAR RUNNING MEAN YEARLY GMST ANOMALY

# GET DAILY DATE VECTORS CORRESPONDING TO STUDY PERIOD
dl <- as.Date(as.Date("1950-01-01"):as.Date("2022-12-31"))
date <- data.frame(matrix(NA,nrow=length(dl),ncol=3))
names(date) <- c("date","year","month")
date$date <- dl
date$year <- as.integer(format(date$date,"%Y"))
date$month <- as.integer(format(date$date,"%m"))

# PRE-DEFINE ARRAY
gmstdaily <- data.frame(matrix(NA,nrow=length(date$date),ncol=3))
names(gmstdaily) <- c("gmst.est","gmst.ci.l","gmst.ci.u")

# LOOP OVER YEARS
for (i in seq(stp)){
  ind <- (date$year==stp[i] & date$month %in% sm)
  
  # LOOP OVER GMST VERSION
  for (g in 1:3){
    gmstdaily[ind,g] <- rep(gmstsub[i,g],sum(ind==T))
  }
}

#############################################
# REGRESSION AND COUNTERFACTUAL COMPUTATION
#############################################

# PRE-DEFINE ARRAY TO STORE LINEAR REGRESSION COEFFICIENTS WITH 95% CI
beta <- array(NA,dim=c(length(cities),3,3),
              dimnames=list(cities,
                            c("beta.est","beta.ci.l","beta.ci.u"),
                            c("gmst.est","gmst.ci.l","gmst.ci.u")))

cf <- array(NA,dim=c(length(date$date),3,3),
              dimnames=list(date$date,
                            c("beta.est","beta.ci.l","beta.ci.u"),
                            c("gmst.est","gmst.ci.l","gmst.ci.u")))

# DEFINE COLORS 
collight <- c("lightblue","grey","pink")
col <- c("blue","black","red")

# INITIATE PLOT
pdf("FigS11_gmsty.pdf",width=9,height=7)
layout(matrix(1:15,ncol=5,byrow=T))
par(mar=c(4,3.8,3,1),mgp=c(2.5,1,0),las=1)

# LOOP OVER CITIES

for (i in seq(cities)){
  
  # IMPORT ERA5-LAND DATA
  path <- paste0("TemperatureSeries/ERA5_1950_2024/TmeanERA5land",cities[i],".csv")
  tdat <- read.csv(path)
  tdat$date <- as.Date(tdat$date, format="%Y-%m-%d")
  tdat$year <- as.integer(format(tdat$date,"%Y"))
  tdat$month <- as.integer(format(tdat$date,"%m"))
  
  # RESTRICT STUDY-PERIOD TO 1950-2022
  tdat <- subset(tdat,tdat$year %in% stp)

  # COMPUTE MONTHLY AVERAGE BY YEAR
  tdatm <- subset(tdat,tdat$month %in% sm)
  tdatm <- tapply(tdatm$tmean,tdatm$year,mean,na.rm=T)
  
  # LOOP OVER GMST VERSIONS
  # TO PLOT DATA POINTS
  for (g in 3:1){

    gmstdat <- gmstsub[,g]
    
    if (g==3){
      plot(gmstdat,tdatm,xlab="GMST anomaly (°C)", 
           ylab= "Warm-season temperature (°C)",
           main=cities[i],type="n",
           xlim=c(min(gmstsub[,2]),max(gmstsub[,3])),
           ylim=c(14,21))
    
    }
    
    points(gmstdat,tdatm,pch=19,col=collight[g])
  }
  
  # LOOP OVER GMST VERSIONS
  # TO PLOT LINEAR REGRESSION LINES
  
  for (g in 1:3){
    
    gmstdat <- gmstsub[,g]
    
    # ADD LINEAR REGRESSION LINE
    fitlm <- lm(tdatm ~ gmstdat)
    lines(gmstdat,fitlm$fitted.values,xpd=F,lty=1,col=col[g])
    
    # SAVE REGRESSION COEFFICIENTS WITH 95% CI
    beta[i,1,g] <- fitlm$coefficients[2]
    beta[i,2:3,g] <- confint(fitlm,'gmstdat',level=0.95)
    
    # COMPUTE COUNTERFACTUAL DAILY TEMPERATURES DURING WARM-SEASON MONTH
    
    # LOOP OVER REGRESSION COEFFICIENTS
    for (b in 1:3) {
      cf[,b,g] <- tdat$tmean - (beta[i,b,g]*gmstdaily[,g])
    }
    
  }

  if (i==1) {legend("topleft",c("Best estimate","2.5% CI","97.5% CI"),pch=19,col=collight,cex=0.8,bty="n")}
  
  # RE-SHAPE ARRAY TO DATAFRAME AND SAVE FACTUAL AND COUNTERFACTUAL TEMPERATURE DATA
  
  cfdf <- as.data.frame(cf)
  out <- cbind(tdat,cfdf)
  names(out) <- c("date","city","tmeanf","year","month",paste0("tmeancf",1:9))
  
  write.csv(out[,c("date","city","tmeanf",paste0("tmeancf",1:9))],
             paste0("TemperatureSeries/Detrended/TmeanERA5land",cities[i],"CounterfactualSensitivity.csv"),
            row.names=F)
  
}

dev.off()

#############



 