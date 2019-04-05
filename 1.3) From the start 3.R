#------------------------------------------------------------------------------
# Running the stochastic scenarios for the North Sea model
#------------------------------------------------------------------------------

# Preliminary stuff

rm(list=ls())

##number of iterations
its<-100

setwd("~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea")

##the various rates of change:
rts<-as.list(seq(10, 50, 10))
##save them out so we can pull them back in in the next script:
save(rts, file="EWS output/rts.RData")

##the number of from the start runs of the model we impliment
from.the.start<-3
set.seed=from.the.start

source("R/paramNorthSeaModel.r")
source("R/calibration_funcs.r")
source("R/fmsyFuncs.r")
#source('~/Desktop/Work/Cod stock EWS/Julia code/3. Multispecies size spectrum model/multispecies_sizebasedmodels/NorthSea/R/fmsyFuncs.r', chdir = TRUE)

source("../R/SizeBasedModel.r")
source("../R/SelectivityFuncs.r")
source("../R/plots.r")
source("../R/Indicators.r")
source("../R/summaryFuncs.r")
source("../R/plots.r")

##new size returning function
source('~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/extract mean weight by spp function.R', chdir = TRUE)

#------------------------------------------------------------------------------

stochScen<-function(modelnum="11",
                    iterset="D",
                    fixedFeedingLevel = NA,
                    fixedPredationMortality = NA,
                    fishingOff=F, 
                    hold.year, 
                    rate.of.change, 
                    iterations, 
                    fmsy.to.zero, 
                    endYear,
                    reimpliment.fishing, 
                    from.the.start){

#----------------------------------------------------------
# Set up the stochastic projections
# Project 200 years at 1967 F level
# From 1967 to 2010 at 'true' Fs
# From 2011 to 2015 linear change in F to Fmsy
# 2016 to 2050 constant F according to scenario
##load to run in function:
  # rate.of.change=30
  # hold.year=2050
  # modelnum="11"
  # iterset="E"
  # fishingOff=F
  # fixedFeedingLevel = NA
  # fixedPredationMortality = NA
  # fmsy.to.zero=T
  # endYear<-2300
  # reimpliment.fishing=2200
  #iters<-iter<-iterations<-1

#load the models
load(file=paste("output/baseModel_",modelnum,".Rdata",sep=""))
load(file=paste("output/modelF0_",modelnum,".Rdata",sep=""))

# Load the real Fs and rescale them
load("input/Fmat.RData")
 
baseF <- diag(baseModel$param$Q)

Fmat<-Fmat[names(baseF),]

relativeF <- sweep(Fmat,1,baseF,"/")

# getSSFmsy(modelscen="baseModel",modelnum="32")

# Load Fmsy data - make sure you get right one 

#  hold at status quo, or multiplier of status quo, one species at a time
if (length(grep("D",iterset))==1) {
  fmsy<- relativeF[,"2010"]
  rate.of.change<-0
}

# all Fs to ICES Fmsy
if (length(grep("E",iterset))==1) {
fmsy<-(c(0.2,0.2,0.2,
0.25,
0.2,
0.2,
0.22,
0.2,
0.25,
0.3,
0.19,
0.3))/baseF
 }

# Use Fmult - relative

startYear <- 1767
endYear <- endYear
nYears <- endYear-startYear+1

# effort is a 2D array: time steps x species
# set up effort array, fill it in, copy it across
effort <- array(0,dim = c(nrow(baseModel$param$species), nYears), dimnames = list(species = baseModel$param$species$species, year = startYear:endYear))
# put in real data
effort[,as.character(1967:2010)] <- relativeF[,as.character(1967:2010)]
# set the historical data - recycling should be OK
effort[,as.character(startYear:1966)] <- relativeF[,as.character(1967)]
# effort[,as.character(startYear:1966)] <- 1
# effort[,as.character(startYear:1966)] <- relativeF[,as.character(1967)]
# set historical to almost unfished
# effort[,as.character(startYear:1966)] <- 0.05
# From 2015 we are at FMSY - again recyclying should be OK
effort[,as.character(2010:endYear)] <- relativeF[,"2010"]

##set to FMSY to zero?
if(fmsy.to.zero==T & rate.of.change!=0){
  effort[,as.character((hold.year+rate.of.change):endYear)]<-0
}
if(fmsy.to.zero==F & rate.of.change!=0){
  effort[,as.character((hold.year+rate.of.change):endYear)]<-fmsy
}

# Finally, set  rate of change between fishing efforts, based on rate.of.change which is the time it takes to do so
for (i in 1:nrow(baseModel$param$species)){
    effort[i,as.character(hold.year:(hold.year+rate.of.change))] <- seq(from = effort[i,as.character(hold.year)], to = effort[i,as.character(hold.year+rate.of.change)], length = rate.of.change+1)
}

##################################################################################################################
##impliment an increase in fishing pressure again after the decline? 
##as per julia's suggestion
##this will be the year that fishing starts again, and it will increase to fmsy
#reimpliment.fishing<-2200
if(is.numeric(reimpliment.fishing) & reimpliment.fishing<endYear){
  effort[,as.character((reimpliment.fishing+rate.of.change):endYear)]<-fmsy
  
  for (i in 1:nrow(baseModel$param$species)){
    effort[i,as.character(reimpliment.fishing:(reimpliment.fishing+rate.of.change))] <- seq(from = effort[i,as.character(reimpliment.fishing)], to = effort[i,as.character(reimpliment.fishing+rate.of.change)], length = rate.of.change+1)
  }
}else{cat("no reimplimentation of fishing, check inputs\n")}
##################################################################################################################
##plot out the effort data to check:
dd<-melt(effort)
#ggplot(subset(dd, year>1990), aes(x=year, y=value, col=species))+geom_line()

# Set up parameters based on base model but with correct number of time steps
projParam <- baseModel$param


projParam$tmax <- nYears
projModel <- Setup(projParam,ContinueCalculation=TRUE,initialcommunity=baseModel)
if(length(grep("unfished",iterset))==1) projModel <- Setup(projParam,ContinueCalculation=TRUE,initialcommunity=modelF0)

# projModel <- Setup(projParam) - Not starting from the calibrated model has a big effect!!
# Set the new effort (need to transpose)

projModel$effort <- t(effort)

# Option to run without fishing to get stochastic model baseline
 if (fishingOff==T) projModel$effort[,] <- 0


# Make some noise!
# sd_srr <- 0.7 # Lognormally distributed noise applied to RDD - sd is 0.7 (see Julia Blanchard email 28/07/12) 
# can we try species-specific values? 

sd_srr <- c(0.9,0.9,0.65,0.8,0.48,0.35,0.73,1.0,0.48,1.0,0.57,0.44) 
sd_rrPP <- 0 # applied to nPP (all sizes of background spectrum affected equally - multiplicative noise)
backgroundAutocorrelationFactor <- 0.5 # doesn't matter as background noise turned off

# Get the unfished summary - used for relative biomass
summaryF0 <- summarySpecies(modelF0, meanTimeSteps = 10,minw=0)

# base biomass
biomassF0 <- summaryF0$biomass
SSBF0 <- summaryF0$SSB


# Run a deterministic test projection - has effort set up worked?
test <- Project(projModel,fixedFeedingLevel=fixedFeedingLevel,fixedPredationMortality=fixedPredationMortality)

# plotResults(test)
 plotBioTime(test)
 save(test,file=paste("EWS output/deter_model",modelnum,iterset,rate.of.change,".RData",sep=""))
 tempSp <- summarySpeciesTime(test, baseBiomass = biomassF0, baseSSB = SSBF0, minw=10)
 save( tempSp,file=paste("EWS output/deter_popI_model",modelnum,iterset,rate.of.change,".RData",sep=""))
 # tempCom1 <- summaryCommunityTime(test, species = c(3,5:11), minw = 10, maxw = 10^5)
 tempCom1 <- summaryCommunityTime(test, species = c(3,5:12), minw = 10, maxw = 10^3)
 tempCom2 <- summaryCommunityTime(test, species = c(3,5:12), minw = 10, maxw = 10^4)
 tempCom3 <- summaryCommunityTime(test, species = c(3,5:12), minw = 10, maxw = 10^5)
 
 save(tempCom3,file=paste("EWS output/deter_comI3_model",modelnum,iterset,rate.of.change,".RData",sep=""))
 

# Important!
  iters <- iterations
# Empty output objects that will get filled up
opSp <- NULL
opCom1 <- NULL
opCom2 <- NULL
opCom3 <- NULL
opCom4 <- NULL


meant<-mint<-maxt<-array(NA, dim=c(12,length(projModel$w),iters))
meantG<-mintG<-maxtG<-array(NA, dim=c(12,length(projModel$w),iters))

# Loop over the iterations and store the projected models
for (iter in 1:iters){
    cat("\nRunning iter: ", iter, "of rate:",rate.of.change, "\n", "% of rate done:", (iter/iterations)*100)
    # Project the model with noise
    tempModel <- Project(projModel, sd_srr=sd_srr, sd_rrPP=sd_rrPP, backgroundAutocorrelationFactor=backgroundAutocorrelationFactor,  fixedFeedingLevel=fixedFeedingLevel,fixedPredationMortality=fixedPredationMortality)
    # plotResults(tempModel)
    # bt <- plotBioTime(tempModel)

    #{x<-model<-tempModel;baseBiomass = biomassF0; baseSSB = SSBF0; minw=10; maxw = max(x$w); minl=NULL; maxl=NULL; search.effort=0.01}
    # Get summary species results
    tempSp <- summarySpeciesTime_SIZE(tempModel, baseBiomass = biomassF0, baseSSB = SSBF0, minw=10)
    tempSp <- cbind(iter = iter, tempSp)
    opSp <- rbind(opSp,tempSp)

    # get summary community results 

    tempCom1 <- summaryCommunityTime(tempModel, species = c(3,5:12), minw = 10, maxw = 10^3)
    tempCom1 <- cbind(iter = iter, tempCom1)
    opCom1 <- rbind(opCom1,tempCom1)
    
    tempCom2 <- summaryCommunityTime(tempModel, species = c(3,5:12), minw = 10, maxw = 10^4)
    tempCom2 <- cbind(iter = iter, tempCom2)
    opCom2 <- rbind(opCom2,tempCom2)
    
    tempCom3 <- summaryCommunityTime(tempModel, species = c(3,5:12), minw = 10, maxw = 10^5)
    tempCom3 <- cbind(iter = iter, tempCom3)
    opCom3 <- rbind(opCom3,tempCom3)
    
    # get species size spectra and size-specific growth rates 
    
    meant[,,iter]<-apply(test$N[(1985-startYear):(1995-startYear),,],c(2,3),mean)
    maxt[,,iter]<-apply(test$N[(1985-startYear):(1995-startYear),,],c(2,3),max)
    mint[,,iter]<-apply(test$N[(1985-startYear):(1995-startYear),,],c(2,3),min)

    meantG[,,iter]<-apply(test$gg[(1985-startYear):(1995-startYear),,],c(2,3),mean)
    maxtG[,,iter]<-apply(test$gg[(1985-startYear):(1995-startYear),,],c(2,3),max)
    mintG[,,iter]<-apply(test$gg[(1985-startYear):(1995-startYear),,],c(2,3),min)

    
}

opCom1$iter<-opCom1$iter+200
opCom2$iter<-opCom2$iter+200
opCom3$iter<-opCom2$iter+200
opSp$iter<-opSp$iter+200
#------------------------------------------------------------------------------
save(opCom1,file=paste("EWS output/opCom1_",modelnum,iterset,rate.of.change,"fts=",from.the.start,".Rdata",sep=""))
save(opCom2,file=paste("EWS output/opCom2_",modelnum,iterset,rate.of.change,"fts=",from.the.start,".Rdata",sep=""))
save(opCom3,file=paste("EWS output/opCom3_",modelnum,iterset,rate.of.change,"fts=",from.the.start,".Rdata",sep=""))

##add in the rate of change
opSp$rate<-rate.of.change
opSp$hold.year<-hold.year
opSp$reimpliment.fishing<-reimpliment.fishing
opSp$from.the.start=from.the.start
save(opSp,file=paste("EWS output/opSp_",modelnum,iterset,rate.of.change,"fts=",from.the.start,".Rdata",sep=""))

stochGw<-list(meanG=meantG[,,iter],maxG=maxtG[,,iter],minG=mintG[,,iter])
stochNw<-list(meanN=meant[,,iter],maxN=maxt[,,iter],mint=mint[,,iter])

save(stochGw,file=paste("EWS output/stochBaseGG",modelnum,iterset,rate.of.change,"fts=",from.the.start,".RData",sep=""))
save(stochNw,file=paste("EWS output/stochBaseN",modelnum,iterset,rate.of.change,"fts=",from.the.start,".RData",sep=""))


# end stochScen function

 }

###################################################################################################################################################
# Model 11. Multispecies (Interactions) Model

# # model projection A -  hold at Fs from 2010
#  stochScen(modelnum="11",iterset="D",fishingOff=F,fixedFeedingLevel = NA,fixedPredationMortality = NA, hold.year=2040, rate.of.change=10, iterations=its)
#  
# model projection B -  ICES Fmsys 
library(parallel)

#mclapply(rts, function(rate, its){
#  stochScen(modelnum="11",iterset="E",fishingOff=F,fixedFeedingLevel = NA,fixedPredationMortality = NA, hold.year=2050, rate.of.change=rate, iterations=its)
#}, its=its,   mc.cores = 7)

lapply(rts, function(rate, its){
  stochScen(modelnum="11",
            iterset="E",
            fishingOff=F,
            fixedFeedingLevel = NA,
            fixedPredationMortality = NA, 
            hold.year=2040, 
            rate.of.change=rate, 
            iterations=its, 
            fmsy.to.zero=T, 
            endYear=2200,
            reimpliment.fishing=NA,
            from.the.start = from.the.start)
}, its=its)

stochScen(modelnum="11",
          iterset="D",
          fishingOff=F,
          fixedFeedingLevel = NA,
          fixedPredationMortality = NA, 
          hold.year=2040, 
          rate.of.change=10, 
          iterations=its, 
          fmsy.to.zero=T,
          endYear=2200,
          reimpliment.fishing=NA,
          from.the.start = from.the.start)