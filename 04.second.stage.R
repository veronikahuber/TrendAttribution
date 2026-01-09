#####################################################################################
# R code for the analysis in 

# Huber, V., Breitner-Busch, S., Feldbusch, H. et al. Improvements in life expectancy 
# mask rising trends in heat-related excess mortality attributable to climate change. 
# Nat Commun 16, 11632 (2025). https://doi.org/10.1038/s41467-025-66681-0
######################################################################################

############################################
# SECOND-STAGE MIXED-EFFECT META-REGRESSION
############################################

library(corrplot)
library(lmtest)

# ESTRACT COEF/VCOV FROM FIRST-STAGE MODELS
coef <- as.matrix(tmeanperpar[,grep("b", names(tmeanperpar))])
vcov <- as.matrix(tmeanperpar[,grep("V", names(tmeanperpar))])

# CITY-SPECIFIC META-DATA
cityinfo <- tmeanperpar[,1:19,]

##########################################
# DEFINE VERSION FOR SENSITIVITY ANALYSES
##########################################

# DEFINE VERSION NAME TAG FOR PLOTS
vers <- c("default","gmsty","sim","5subperiods","without2003")

# CHOOSE VERSION
fext <- vers[1]

#############################################
# PLOT META-PREDICTORS CORRELATION MATRIX
#############################################

# RE-NAME CLIMATE VARIABLES
names(cityinfo)[4:7] <- c("average.annual.temp","range.temp","average.summer.temp","sum.heat.days")

# CHANGE POPULATION SIZE TO PER MIO.
cityinfo$population.corrected <- cityinfo$population.corrected/1000000

# PLOT MATRIX CORRELATION FOR META-PREDICTORS
M <- cor(cityinfo[,4:19])

testRes = cor.mtest(cityinfo[,4:19], conf.level = 0.95)

pdf(paste0("CorrelationMatrix1_",fext,".pdf"))
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.7, order = 'AOE', diag=FALSE)
dev.off()

# TAKE OUT HIGHLY-CORRELATED META-PREDICTORS
cityinfo = cityinfo[,!(names(cityinfo) %in% c("unemployed",
                                              #"average.population.age",
                                              "rest.life.expectancy",
                                              "hpp"))]

# PLOT SECOND MATRIX CORRELATION FOR SELECTED META-PREDICTORS
M <- cor(cityinfo[,4:16])

testRes = cor.mtest(cityinfo[,4:16], conf.level = 0.95)

pdf(paste0("CorrelationMatrix2_",fext,".pdf"))
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8,diag=FALSE)
dev.off()

######################################
# ADD ADDITIONAL META-PREDICTORS
######################################

# ADD INFO ON FEDERAL STATES
fstate <- c("Berlin","Bremen","NRW","NRW","Sachsen","NRW",
                     "NRW","NRW","Hessen","Hamburg","Niedersachsen",
                     "Sachsen","Bayern","Bayern","BW")
cityinfo$fstate <- rep(fstate,each=length(yearlist))

# ADD INFO ON EAST/WEST
eastwest <- c("E","W","W","W","E","W","W","W","W","W","W","E","W","W","W")
cityinfo$eastwest <- rep(eastwest,each=length(yearlist))
  
####################
# FIND BEST MODEL
####################

# RANDOM-INTERCEPT ONLY MODEL
model0 <- mixmeta(coef~1, vcov, data=cityinfo, method="ml",
                  random=~1|city, bscov="diag")

# ADD NESTED RANDOM EFFECT FOR EAST/WEST
model01 <- mixmeta(coef~1, vcov, data=cityinfo, method="ml",
                   random=~1|eastwest/city, bscov="diag")

# ADD NESTED RANDOM EFFECT FOR FEDERAL STATE
model02 <- mixmeta(coef~1, vcov, data=cityinfo, method="ml",
                   random=~1|fstate/city, bscov="diag")

# BASE MODEL WITH NATURAL SPLINE OF YEAR AS META-PREDICTOR, RANDOM INTERCEPT CITY
model03 <- mixmeta(coef~1+ns(year,knots=median(year)), vcov, data=cityinfo, method="ml",
                       random=~1|city, bscov="diag")
drop1(model03, test="Chisq")

# MODEL WITH RANDOM SLOPES FOR YEAR TERMS
model04 <- mixmeta(coef~1+ns(year,knots=median(year)), vcov, data=cityinfo, method="ml",
                        random=~1+ns(year,knots=median(year))|city, bscov="diag")
drop1(model04, test="Chisq")

# # # SINGLE PREDICTOR MODELS
model1 <- update(model0, .~+ns(year,knots=median(year))+average.annual.temp)
drop1(model1, test="Chisq")
model2 <- update(model0, .~+ns(year,knots=median(year))+range.temp)
drop1(model2, test="Chisq")
model3 <- update(model0, .~+ns(year,knots=median(year))+average.summer.temp)
drop1(model3, test="Chisq")
model4 <- update(model0, .~+ns(year,knots=median(year))+sum.heat.days)
drop1(model4, test="Chisq")
model5 <- update(model0, .~+ns(year,knots=median(year))+population.corrected)
drop1(model5, test="Chisq")
model6 <- update(model0, .~+ns(year,knots=median(year))+women.quota)
drop1(model6, test="Chisq")
model7 <- update(model0, .~+ns(year,knots=median(year))+residents.over.65)
drop1(model7, test="Chisq")
model8 <- update(model0, .~+ns(year,knots=median(year))+life.expectancy)
drop1(model8, test="Chisq")
model9 <- update(model0, .~+ns(year,knots=median(year))+foreign.residents)
drop1(model9, test="Chisq")
model10 <- update(model0, .~+ns(year,knots=median(year))+household.income)
drop1(model10, test="Chisq")
model11 <- update(model0, .~+ns(year,knots=median(year))+population.density)
drop1(model11, test="Chisq")
model12 <- update(model0, .~+ns(year,knots=median(year))+perc.unemployed)
drop1(model12, test="Chisq")
model13 <- update(model0, .~+ns(year,knots=median(year))+average.population.age)
drop1(model13, test="Chisq")

# # # COMPARE MODELS
AIC(model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13)
BIC(model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13)

# FULL MODEL (NOTE: POPULATION DATA TAKEN OUT TO AVOID CONVERGENCE ISSUES)
fullmodel <- update(model0, .~ns(year,knots=median(year))+average.annual.temp+range.temp+average.summer.temp+sum.heat.days+
                       women.quota+residents.over.65+life.expectancy+foreign.residents+household.income+
                        population.density+perc.unemployed+average.population.age)

# TEST (BACKWARD) (NOTE: NO VARIABLE IS IDENTIFIED AS SIGNIFICANT IN THE LRT)
drop1(fullmodel, test="Chisq")

# MODEL SELECTION (STEP FORWARD) (NOTE: POPULATION DATA TAKEN OUT TO AVOID CONVERGENCE ISSUES)
step(model0, .~ns(year,knots=median(year))+average.annual.temp+range.temp+average.summer.temp+sum.heat.days+
                  women.quota+residents.over.65+life.expectancy+foreign.residents+household.income+
                  population.density+perc.unemployed+average.population.age)

# FINAL MODEL (STEP FORWARD)
modelforward <- update(model0, .~ns(year,knots=median(year))+average.annual.temp+
                         sum.heat.days+life.expectancy+average.population.age)
drop1(modelforward, test="Chisq")

modelfinal <- modelforward

######

