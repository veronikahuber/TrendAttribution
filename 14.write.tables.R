#################################################################################
# R code for the analysis in 

# Huber, V., et al. Life expectancy gains mask rising trends in heat-related
# excess mortality attributable to climate change. 
# Nature Communications, [to be added] (2025), doi: [to be added]

#################################################################################

##########################################
# WRITE OUTPUT AS TABLES
##########################################

################
# PREPARE DATA
################

#########################################################
# WARM-SEASON TEMPERATURE ATTRIBUTABLE TO CLIMATE CHANGE
# PERIOD AVERAGE, MIN, MAX
#########################################################

# CBIND CITY-ESTIMATES 
mstattrcity <- cbind(mstwarm[,"cfact.best"],apply(mstwarm[,2:9],1,min,na.rm=T),apply(mstwarm[,2:9],1,max,na.rm=T))

# CBIND COUNTRY-ESTIMATES
mstattrcountry <- c(mstwarmcountry["cfact.best"],min(mstwarmcountry[2:9]),max(mstwarmcountry[2:9]))

# RBIND AND ROUND TO 2 DIGITS
mstattr <- round(rbind(mstattrcity,mstattrcountry),digits=2)

# PREPARE FOR SAVING
mstattrout <- paste0(mstattr[,1]," (",
              mstattr[,2],", ",
              mstattr[,3],")")

#########################################################
# HEAT EXCESS MORTALITY ATTRIBUTABLE TO CLIMATE CHANGE
# PERIOD AVERAGE, MIN, MAX
#########################################################

################################
# ATTRIBUTABLE NUMBERS
################################

# BY CITY

# AVERAGE PER YEAR
ancityavg <- paste0(round(meanheatanvardif[,basemort,"abs","cfact.best","est"]/length(ny),digits=0)," (",
                    round(meanheatanvardif[,basemort,"abs","cfact.best","ci.l"]/length(ny),digits=0),", ",
                    round(meanheatanvardif[,basemort,"abs","cfact.best","ci.u"]/length(ny),digits=0),")")

# SUM OVER STUDY PERIOD
ancitysum <- paste0(round(meanheatanvardif[,basemort,"abs","cfact.best","est"],digits=0)," (",
                    round(meanheatanvardif[,basemort,"abs","cfact.best","ci.l"],digits=0),", ",
                    round(meanheatanvardif[,basemort,"abs","cfact.best","ci.u"],digits=0),")")

# MINIMUM ACROSS YEARS
ancitymin <-  matrix(data=NA,nrow=length(cities),ncol=2)
for (i in seq(length(cities))){
  ind <- heatanvardif[i,,basemort,"abs","cfact.best","est"] == min(heatanvardif[i,,basemort,"abs","cfact.best","est"],na.rm=T)
  # SET NA IN INDICATOR TO FALSE IF IT OCCURS
  ind[is.na(ind)] <- FALSE
  ancitymin[i,1] <- paste0(round(heatanvardif[i,ind,basemort,"abs","cfact.best","est"],digits=0)," (",
                      round(heatanvardif[i,ind,basemort,"abs","cfact.best","ci.l"],digits=0),", ",
                      round(heatanvardif[i,ind,basemort,"abs","cfact.best","ci.u"],digits=0),")")
  ancitymin[i,2] <- ny[ind]
}

# MAXIMUM ACROSS YEARS
ancitymax <- matrix(data=NA,nrow=length(cities),ncol=2)
for (i in seq(length(cities))){
  ind <- heatanvardif[i,,basemort,"abs","cfact.best","est"] == max(heatanvardif[i,,basemort,"abs","cfact.best","est"],na.rm=T)
  # SET NA IN INDICATOR TO FALSE IF IT OCCURS
  ind[is.na(ind)] <- FALSE
  ancitymax[i,1] <- paste0(round(heatanvardif[i,ind,basemort,"abs","cfact.best","est"],digits=0)," (",
                         round(heatanvardif[i,ind,basemort,"abs","cfact.best","ci.l"],digits=0),", ",
                         round(heatanvardif[i,ind,basemort,"abs","cfact.best","ci.u"],digits=0),")")
  ancitymax[i,2] <- ny[ind]
}

# ATTRIBUTABLE NUMBERS ALL CITIES TOGETHER
ancountryavg <- paste0(round(meanheatantotvardif[basemort,"abs","cfact.best","est"]/length(ny),digits=0)," (",
                       round(meanheatantotvardif[basemort,"abs","cfact.best","ci.l"]/length(ny),digits=0),", ",
                       round(meanheatantotvardif[basemort,"abs","cfact.best","ci.u"]/length(ny),digits=0),")")
ancountrysum <- paste0(round(meanheatantotvardif[basemort,"abs","cfact.best","est"],digits=0)," (",
                       round(meanheatantotvardif[basemort,"abs","cfact.best","ci.l"],digits=0),", ",
                       round(meanheatantotvardif[basemort,"abs","cfact.best","ci.u"],digits=0),")")
ind <- heatantotvardif[,basemort,"abs","cfact.best","est"] == min(heatantotvardif[,basemort,"abs","cfact.best","est"],na.rm=T)
 
ancountrymin <- c(paste0(round(heatantotvardif[ind,basemort,"abs","cfact.best","est"],digits=0)," (",
                       round(heatantotvardif[ind,basemort,"abs","cfact.best","ci.l"],digits=0),", ",
                       round(heatantotvardif[ind,basemort,"abs","cfact.best","ci.u"],digits=0),")"),
                  ny[ind])

ind <- heatantotvardif[,basemort,"abs","cfact.best","est"] == max(heatantotvardif[,basemort,"abs","cfact.best","est"],na.rm=T)
ancountrymax <- c(paste0(round(heatantotvardif[ind,basemort,"abs","cfact.best","est"],digits=0)," (",
                       round(heatantotvardif[ind,basemort,"abs","cfact.best","ci.l"],digits=0),", ",
                       round(heatantotvardif[ind,basemort,"abs","cfact.best","ci.u"],digits=0),")"),
                  ny[ind])

##########################
# ATTRIBUTABLE FRACTIONS 
##########################
# DEFINE DIGITS
digits <- 2
# BY CITY
# AVERAGE PER YEAR
afcityavg <- paste0(round(meanheatafvardif[,basemort,"abs","cfact.best","est"],digits=digits)," (",
                    round(meanheatafvardif[,basemort,"abs","cfact.best","ci.l"],digits=digits),", ",
                    round(meanheatafvardif[,basemort,"abs","cfact.best","ci.u"],digits=digits),")")
# MINIMUM ACROSS YEARS
afcitymin <- matrix(data=NA,nrow=length(cities),ncol=2)
for (i in seq(length(cities))){
  ind <- heatafvardif[i,,basemort,"abs","cfact.best","est"] == min(heatafvardif[i,,basemort,"abs","cfact.best","est"],na.rm=T)
  # SET NA IN INDICATOR TO FALSE IF IT OCCURS
  ind[is.na(ind)] <- FALSE
  afcitymin[i,1] <- paste0(round(heatafvardif[i,ind,basemort,"abs","cfact.best","est"],digits=digits)," (",
                         round(heatafvardif[i,ind,basemort,"abs","cfact.best","ci.l"],digits=digits),", ",
                         round(heatafvardif[i,ind,basemort,"abs","cfact.best","ci.u"],digits=digits),")")
  afcitymin[i,2]  <- ny[ind]
}
# MAXIMUM ACROSS YEARS
afcitymax <-  matrix(data=NA,nrow=length(cities),ncol=2)
for (i in seq(length(cities))){
  ind <- heatafvardif[i,,basemort,"abs","cfact.best","est"] == max(heatafvardif[i,,basemort,"abs","cfact.best","est"],na.rm=T)
  # SET NA IN INDICATOR TO FALSE IF IT OCCURS
  ind[is.na(ind)] <- FALSE
  afcitymax[i,1] <- paste0(round(heatafvardif[i,ind,basemort,"abs","cfact.best","est"],digits=digits)," (",
                         round(heatafvardif[i,ind,basemort,"abs","cfact.best","ci.l"],digits=digits),", ",
                         round(heatafvardif[i,ind,basemort,"abs","cfact.best","ci.u"],digits=digits),")")
  afcitymax[i,2] <- ny[ind]
}


# ATTRIBUTABLE FRACTIONS ALL CITIES TOGETHER
afcountryavg <- paste0(round(meanheataftotvardif[basemort,"abs","cfact.best","est"],digits=digits)," (",
                       round(meanheataftotvardif[basemort,"abs","cfact.best","ci.l"],digits=digits),", ",
                       round(meanheataftotvardif[basemort,"abs","cfact.best","ci.u"],digits=digits),")")
ind <- heataftotvardif[,basemort,"abs","cfact.best","est"] == min(heataftotvardif[,basemort,"abs","cfact.best","est"],na.rm=T)
# SET NA IN INDICATOR TO FALSE IF IT OCCURS
ind[is.na(ind)] <- FALSE
afcountrymin <- c(paste0(round(heataftotvardif[ind,basemort,"abs","cfact.best","est"],digits=digits)," (",
                       round(heataftotvardif[ind,basemort,"abs","cfact.best","ci.l"],digits=digits),", ",
                       round(heataftotvardif[ind,basemort,"abs","cfact.best","ci.u"],digits=digits),")"),
                  ny[ind])

ind <- heataftotvardif[,basemort,"abs","cfact.best","est"] == max(heataftotvardif[,basemort,"abs","cfact.best","est"],na.rm=T)

afcountrymax <- c(paste0(round(heataftotvardif[ind,basemort,"abs","cfact.best","est"],digits=digits)," (",
                       round(heataftotvardif[ind,basemort,"abs","cfact.best","ci.l"],digits=digits),", ",
                       round(heataftotvardif[ind,basemort,"abs","cfact.best","ci.u"],digits=digits),")"),
                  ny[ind])

##########################
# ATTRIBUTABLE PROPORTIONS 
##########################
digits <- 1
# BY CITY
# AVERAGE PER YEAR
pcityavg <- paste0(round(meanheatafvardif[,basemort,"rel","cfact.best","est"],digits=digits)," (",
                    round(meanheatafvardif[,basemort,"rel","cfact.best","ci.l"],digits=digits),", ",
                    round(meanheatafvardif[,basemort,"rel","cfact.best","ci.u"],digits=digits),")")
# MINIMUM ACROSS YEARS
pcitymin <- matrix(data=NA,nrow=length(cities),ncol=2)
for (i in seq(length(cities))){
  ind <- heatafvardif[i,,basemort,"rel","cfact.best","est"] == min(heatafvardif[i,,basemort,"rel","cfact.best","est"],na.rm=T)
  # SET NA TO FALSE
  ind[is.na(ind)] <- FALSE
  pcitymin[i,1] <- paste0(round(heatafvardif[i,ind,basemort,"rel","cfact.best","est"],digits=digits)," (",
                         round(heatafvardif[i,ind,basemort,"rel","cfact.best","ci.l"],digits=digits),", ",
                         round(heatafvardif[i,ind,basemort,"rel","cfact.best","ci.u"],digits=digits),")")
  pcitymin[i,1] <- ny[ind]
}

# MAXIMUM ACROSS YEARS
pcitymax <-  matrix(data=NA,nrow=length(cities),ncol=2)
for (i in seq(length(cities))){
  ind <- heatafvardif[i,,basemort,"rel","cfact.best","est"] == max(heatafvardif[i,,basemort,"rel","cfact.best","est"],na.rm=T)
  # SET NA TO FALSE
  ind[is.na(ind)] <- FALSE
  pcitymax[i,1] <- paste0(round(heatafvardif[i,ind,basemort,"rel","cfact.best","est"],digits=digits)," (",
                         round(heatafvardif[i,ind,basemort,"rel","cfact.best","ci.l"],digits=digits),", ",
                         round(heatafvardif[i,ind,basemort,"rel","cfact.best","ci.u"],digits=digits),")")
  pcitymax[i,2] <- ny[ind]
}


# ATTRIBUTABLE PROPORTIONS ALL CITIES TOGETHER
pcountryavg <- paste0(round(meanheataftotvardif[basemort,"rel","cfact.best","est"],digits=digits)," (",
                       round(meanheataftotvardif[basemort,"rel","cfact.best","ci.l"],digits=digits),", ",
                       round(meanheataftotvardif[basemort,"rel","cfact.best","ci.u"],digits=digits),")")
ind <- heataftotvardif[,basemort,"rel","cfact.best","est"] == min(heataftotvardif[,basemort,"rel","cfact.best","est"],na.rm=T)
pcountrymin <- c(paste0(round(heataftotvardif[ind,basemort,"rel","cfact.best","est"],digits=digits)," (",
                       round(heataftotvardif[ind,basemort,"rel","cfact.best","ci.l"],digits=digits),", ",
                       round(heataftotvardif[ind,basemort,"rel","cfact.best","ci.u"],digits=digits),")"),
                 ny[ind])

ind <- heataftotvardif[,basemort,"rel","cfact.best","est"] == max(heataftotvardif[,basemort,"rel","cfact.best","est"],na.rm=T)
pcountrymax <- c(paste0(round(heataftotvardif[ind,basemort,"rel","cfact.best","est"],digits=digits)," (",
                       round(heataftotvardif[ind,basemort,"rel","cfact.best","ci.l"],digits=digits),", ",
                       round(heataftotvardif[ind,basemort,"rel","cfact.best","ci.u"],digits=digits),")"),
                 ny[ind])

###################
# TABLE 1
###################

out <- cbind(mstattrout,
             c(ancitysum,ancountrysum),
             c(ancityavg,ancountryavg),
             c(afcityavg,afcountryavg),
             c(pcityavg,pcountryavg))         

write.table(out, paste0("Tables/Table1_",fext,".csv"),sep=";",
            col.names = c("Temp","ANsum","AN","AF","P"),
            row.names = c(cities,"All cities"))


#######################
# TABLE S6
#######################

out <- cbind(rbind(ancitymax,ancountrymax),
             rbind(afcitymax,afcountrymax),
             rbind(pcitymax,pcountrymax))         

write.table(out, paste0("Tables/TableS6_",fext,".csv"),sep=";",
            col.names = c("AN","year","AF","year","P","year"),
            row.names = c(cities,"All cities"))

############################################
# 2) TABLE S1
#############################################
#############################################

# DEATH SUM
deathsum <- sapply(dlist, function(x) sum(x$death[x$month %in% sm],na.rm=T))
deathsumcountry <- sum(deathsum)

# AVERAGE WARM-SEASON TEMPERATURE 
avgtemp <- round(rowMeans(mst[,,"fact"],na.rm=T),digits=1)
avgtempcountry <-round(mean(mstfcountry,na.rm=T),digits=1)

# AVERAGE POPULATION
avgpop <- round(sapply(cities,function(x) mean(metadat$population.corrected[metadat$cityname==x],na.rm=T)),digits=-3)
avgpopcountry <- sum(avgpop)

# AVERAGE SHARE IN GERMAN POPULATION
popGermany <- read.table("MetaPredictors/Population_Germany_1993_2020.csv",
                       sep=";",encoding="UTF-8",skip=6)

popshare <- round(avgpop/mean(popGermany$V2)*100,digits=1)
popsharecountry <- round(avgpopcountry/mean(popGermany$V2)*100,digits=1)

table <- cbind(c(deathsum,deathsumcountry),
                 c(avgtemp,avgtempcountry),
                 c(avgpop,avgpopcountry),
                 c(popshare,popsharecountry))
rownames(table)<- c(cities,"All cities")
colnames(table) <- c("death","temp","pop","sharepop")
write.table(table, paste0("Tables/TableS1_",fext,".csv"),sep=";")

############################################
############################################
# 2) TABLE 2
#############################################
#############################################

plevel <- function(p){
  levels <- c(0.001,0.01,0.05, 0.1)
  if (p < levels[1]){
    return(paste0("p<",levels[1]))
  }else if (p < levels[2]){
    return(paste0("p<",levels[2]))
  }else if (p < levels[3]){
    return(paste0("p<",levels[3]))
  }else if (p < levels[4]){
    return(paste0("p<",levels[4]))
  }else{
    return("p>0.1")
  }
}

slopetable <- matrix(NA,nrow=6,ncol=4)

# DEFINE NUMBER OF DIGITS
digits=2

#######################
# TEMPERATURE SLOPES
#######################

# FACTUAL TEMPERATURES
slopetable[1,1] <- paste0(round(temptrend["Germany","fact"]*10,digits=digits)," (",
                          round(temptrendci["Germany","fact","ci.l"]*10,digits=digits),", ",
                          round(temptrendci["Germany","fact","ci.u"]*10,digits=digits),")") 

slopetable[1,2] <- plevel(temptrendp["Germany","fact"])

# COUNTERFACTUAL TEMPERATURES
slopetable[2,1] <- paste0(round(temptrend["Germany","cfact.best"]*10,digits=digits)," (",
                               round(temptrendci["Germany","cfact.best","ci.l"]*10,digits=digits),", ",
                               round(temptrendci["Germany","cfact.best","ci.u"]*10,digits=digits),")") 

slopetable[2,2] <- plevel(temptrendp["Germany","cfact.best"])

########################################################
# HEAT-RELATED MORTALITY - WITH/WITHOUT CLIMATE CHANGE
#######################################################

for (i in c(1,2)){
  slopetable[2+i,1] <- paste0(round(morttrend["Germany",i,"adapt"]*10,digits=digits)," (",
                            round(morttrendci["Germany",i,"adapt","ci.l"]*10,digits=digits),", ",
                            round(morttrendci["Germany",i,"adapt","ci.u"]*10,digits=digits),")") 
  slopetable[2+i,2] <- plevel(morttrendp["Germany",i,"adapt"])
  slopetable[2+i,3] <- paste0(round(morttrend["Germany",i,"noadapt"]*10,digits=digits)," (",
                            round(morttrendci["Germany",i,"noadapt","ci.l"]*10,digits=digits),", ",
                            round(morttrendci["Germany",i,"noadapt","ci.u"]*10,digits=digits),")") 
  slopetable[2+i,4] <- plevel(morttrendp["Germany",i,"noadapt"])
}

########################################################
# ATTRIBUTABLE MORTALITY - ABSOLUTE AND RELATIVE
#######################################################

for (i in c(1,2)){
  slopetable[4+i,1] <- paste0(round(attrmorttrend["Germany","adapt",i,"cfact.best"]*10,digits=digits)," (",
                              round(attrmorttrendci["Germany","adapt",i,"cfact.best","ci.l"]*10,digits=digits),", ",
                              round(attrmorttrendci["Germany","adapt",i,"cfact.best","ci.u"]*10,digits=digits),")") 
  slopetable[4+i,2] <- plevel(attrmorttrendp["Germany","adapt",i,"cfact.best"])
  slopetable[4+i,3] <- paste0(round(attrmorttrend["Germany","noadapt",i,"cfact.best"]*10,digits=digits)," (",
                              round(attrmorttrendci["Germany","noadapt",i,"cfact.best","ci.l"]*10,digits=digits),", ",
                              round(attrmorttrendci["Germany","noadapt",i,"cfact.best","ci.u"]*10,digits=digits),")") 
  slopetable[4+i,4] <- plevel(attrmorttrendp["Germany","noadapt",i,"cfact.best"])
}

write.table(slopetable,paste0("Tables/Table2_",fext,".csv"),sep=";")

##################################################
# SUPPLMENTARY TABLE 4
# SLOPE ESTIMATES BASED ON ATTRIBUTABLE NUMBERS
##################################################

slopetable <- matrix(NA,nrow=4,ncol=4)

# DEFINE NUMBER OF DIGITS
digits=0

########################################################
# HEAT-RELATED DEATHS - WITH/WITHOUT CLIMATE CHANGE
#######################################################

for (i in c(1,2)){
  slopetable[i,1] <- paste0(round(antrend["Germany",i,"adapt"],digits=digits)," (",
                              round(antrendci["Germany",i,"adapt","ci.l"],digits=digits),", ",
                              round(antrendci["Germany",i,"adapt","ci.u"],digits=digits),")") 
  slopetable[i,2] <- plevel(antrendp["Germany",i,"adapt"])
  slopetable[i,3] <- paste0(round(antrend["Germany",i,"noadapt"],digits=digits)," (",
                              round(antrendci["Germany",i,"noadapt","ci.l"],digits=digits),", ",
                              round(antrendci["Germany",i,"noadapt","ci.u"],digits=digits),")") 
  slopetable[i,4] <- plevel(antrendp["Germany",i,"noadapt"])
}

########################################################
# ATTRIBUTABLE DEATHS - ABSOLUTE AND RELATIVE
#######################################################

for (i in c(1,2)){
  
  if (i==2) {digits <- 2}
  
  slopetable[2+i,1] <- paste0(round(attrantrend["Germany","adapt",i,"cfact.best"],digits=digits)," (",
                              round(attrantrendci["Germany","adapt",i,"cfact.best","ci.l"],digits=digits),", ",
                              round(attrantrendci["Germany","adapt",i,"cfact.best","ci.u"],digits=digits),")") 
  slopetable[2+i,2] <- plevel(attrantrendp["Germany","adapt",i,"cfact.best"])
  slopetable[2+i,3] <- paste0(round(attrantrend["Germany","noadapt",i,"cfact.best"],digits=digits)," (",
                              round(attrantrendci["Germany","noadapt",i,"cfact.best","ci.l"],digits=digits),", ",
                              round(attrantrendci["Germany","noadapt",i,"cfact.best","ci.u"],digits=digits),")") 
  slopetable[2+i,4] <- plevel(attrantrendp["Germany","noadapt",i,"cfact.best"])
}

write.table(slopetable,paste0("Tables/TableS4_",fext,".csv"),sep=";")

###########################################
# REFORMAT FOR SUPPLEMENTARY TABLE S5
###########################################

slopetable <- matrix(NA,nrow=1,ncol=8)

# DEFINE NUMBER OF DIGITS
digits=1

slopetable[1,1] <- paste0(round(morttrend["Germany",1,"adapt"]*10,digits=digits)," (",
                          round(morttrendci["Germany",1,"adapt","ci.l"]*10,digits=digits),", ",
                          round(morttrendci["Germany",1,"adapt","ci.u"]*10,digits=digits),") ", 
                          plevel(morttrendp["Germany",1,"adapt"])) 

slopetable[1,2] <- paste0(round(morttrend["Germany",2,"adapt"]*10,digits=digits)," (",
                          round(morttrendci["Germany",2,"adapt","ci.l"]*10,digits=digits),", ",
                          round(morttrendci["Germany",2,"adapt","ci.u"]*10,digits=digits),") ", 
                          plevel(morttrendp["Germany",2,"adapt"]))  

slopetable[1,3] <- paste0(round(attrmorttrend["Germany","adapt",1,"cfact.best"]*10,digits=digits)," (",
                          round(attrmorttrendci["Germany","adapt",1,"cfact.best","ci.l"]*10,digits=digits),", ",
                          round(attrmorttrendci["Germany","adapt",1,"cfact.best","ci.u"]*10,digits=digits),") ",
                          plevel(attrmorttrendp["Germany","adapt",1,"cfact.best"]))

slopetable[1,4] <- paste0(round(attrmorttrend["Germany","adapt",2,"cfact.best"]*10,digits=digits)," (",
                          round(attrmorttrendci["Germany","adapt",2,"cfact.best","ci.l"]*10,digits=digits),", ",
                          round(attrmorttrendci["Germany","adapt",2,"cfact.best","ci.u"]*10,digits=digits),") ",
                          plevel(attrmorttrendp["Germany","adapt",2,"cfact.best"]))

slopetable[1,5] <- paste0(round(morttrend["Germany",1,"noadapt"]*10,digits=digits)," (",
                          round(morttrendci["Germany",1,"noadapt","ci.l"]*10,digits=digits),", ",
                          round(morttrendci["Germany",1,"noadapt","ci.u"]*10,digits=digits),") ", 
                          plevel(morttrendp["Germany",1,"noadapt"])) 

slopetable[1,6] <- paste0(round(morttrend["Germany",2,"noadapt"]*10,digits=digits)," (",
                          round(morttrendci["Germany",2,"noadapt","ci.l"]*10,digits=digits),", ",
                          round(morttrendci["Germany",2,"noadapt","ci.u"]*10,digits=digits),") ", 
                          plevel(morttrendp["Germany",2,"noadapt"]))

slopetable[1,7] <- paste0(round(attrmorttrend["Germany","noadapt",1,"cfact.best"]*10,digits=digits)," (",
                          round(attrmorttrendci["Germany","noadapt",1,"cfact.best","ci.l"]*10,digits=digits),", ",
                          round(attrmorttrendci["Germany","noadapt",1,"cfact.best","ci.u"]*10,digits=digits),") ",
                          plevel(attrmorttrendp["Germany","noadapt",1,"cfact.best"]))

slopetable[1,8] <- paste0(round(attrmorttrend["Germany","noadapt",2,"cfact.best"]*10,digits=digits)," (",
                          round(attrmorttrendci["Germany","noadapt",2,"cfact.best","ci.l"]*10,digits=digits),", ",
                          round(attrmorttrendci["Germany","noadapt",2,"cfact.best","ci.u"]*10,digits=digits),") ",
                          plevel(attrmorttrendp["Germany","noadapt",2,"cfact.best"]))

write.table(slopetable,paste0("Tables/TableS5_",fext,".csv"),sep=";")

###########