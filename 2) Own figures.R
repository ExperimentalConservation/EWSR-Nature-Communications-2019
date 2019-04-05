## own code for extracting data and making figures from julias outputs
rm(list=ls())
library(data.table)
library(deSolve)
library(viridis)
library(plyr)
library(ggplot2)
library(reshape2)
setwd("~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea")
source("../R/summaryFuncs.R")	
source("../R/Indicators.R")
source("R/CalcsPlotsPaper2.R")
source("../R/SizeBasedModel.r")
source("R/paramNorthSeaModel.r")
source("../R/SelectivityFuncs.r")
source("../R/plots.r")

##load in the output data:
##these specify which model
run="11"
iterset="E"

##laod the rates
load("EWS output/rts.RData");

##get the output files given the different rates of change, load them, and bind them into a data.table
dd.1<-rbindlist(lapply(rts, function(x){
  fts1<-get(load(paste("EWS output/opSp_",run,iterset,x,"fts=1",".RData",sep="")))
  fts2<-get(load(paste("EWS output/opSp_",run,iterset,x,"fts=2",".RData",sep="")))
  fts3<-get(load(paste("EWS output/opSp_",run,iterset,x,"fts=3",".RData",sep="")))
  return(rbind(fts1, fts2, fts3))
}))

##read in the data with no decrease in fishing pressures:
dd<-rbind(dd.1, 
          get(load('~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/opSp_11D0fts=1.Rdata')),
          get(load('~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/opSp_11D0fts=2.Rdata')),
          get(load('~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/opSp_11D0fts=3.Rdata')))

dd<-dd[dd$time>200,]

##add in years
yrange=1967:2200
dd$year<-yrange

########################################################################################################################################################
###For abundance (biomass)
  ##make similar figures to julia's, with multiple rates on a plot
  quantBIO<-melt(tapply(dd$biomass,list(dd$year,dd$species, dd$rate),quantile,0.5))#;quatSSB$year<-rownames(paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(quantBIO)<-c("year", "species","rate","biomass")
  
  l.quantBIO<-melt(tapply(dd$biomass,list(dd$year,dd$species, dd$rate),quantile,0.05))#[paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(l.quantBIO)<-c("year", "species","rate","lower")
  u.quantBIO<-melt(tapply(dd$biomass,list(dd$year,dd$species, dd$rate),quantile,0.95))#[paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(u.quantBIO)<-c("year", "species","rate","upper")
  
  quantBIO$lower<-l.quantBIO$lower
  quantBIO$upper<-u.quantBIO$upper
  
  #data for polygons
  pols<-data.frame(year=c(quantBIO$year,rev(quantBIO$year)), species=c(as.character(quantBIO$species), rev(as.character(quantBIO$species))), rate=c(quantBIO$rate, rev(quantBIO$rate)),y=c(quantBIO$lower, rev(quantBIO$upper)))
  
  
  quantBIO$rate<-round(1/quantBIO$rate, 3);quantBIO$rate[which(quantBIO$rate=="Inf")]<-0
  quantBIO$rate<-paste("Rate of change =",quantBIO$rate)
  
  ##ggplot of all spp with polygons
  pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/all spp model outcomes.pdf", width=12, height=7)
    p1<-ggplot(quantBIO, aes(x=year, y=biomass))+facet_wrap(~species, scales="free_y")+theme_bw()+geom_vline(xintercept=unique(dd$hold.year), col="black")
    #p1<-p1+geom_polygon(data=pols, aes(x=year, y=y, fill=factor(rate)), alpha=0.05)
    p1<-p1+geom_line(aes(col=factor(rate)), size=1)+scale_colour_manual(values = (viridis(7, option = "B")[1:6]), name="Rate of decline in fishing pressure")+ylab("Mean biomass")
    p1+xlab("Year")
  dev.off()
########################################################################################################################################################
## for body size
  ##make similar figures to julia's, with multiple rates on a plot
  dd$log.mean.size<-(dd$mean.size)
  quantSIZE<-melt(tapply(dd$log.mean.size,list(dd$year,dd$species, dd$rate),quantile,0.5))#;quatSSB$year<-rownames(paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(quantSIZE)<-c("year", "species","rate","mean.size")
  
  l.quantSIZE<-melt(tapply(dd$log.mean.size,list(dd$year,dd$species, dd$rate),quantile,0.05))#[paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(l.quantSIZE)<-c("year", "species","rate","lower")
  u.quantSIZE<-melt(tapply(dd$log.mean.size,list(dd$year,dd$species, dd$rate),quantile,0.95))#[paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(u.quantSIZE)<-c("year", "species","rate","upper")
  
  quantSIZE$lower<-l.quantSIZE$lower
  quantSIZE$upper<-u.quantSIZE$upper
  
  #data for polygons
  pols.size<-data.frame(year=c(quantSIZE$year,rev(quantSIZE$year)), species=c(as.character(quantSIZE$species), rev(as.character(quantSIZE$species))), rate=c(quantSIZE$rate, rev(quantSIZE$rate)),y=c(quantSIZE$lower, rev(quantSIZE$upper)))
  
  quantSIZE$rate<-round(1/quantSIZE$rate, 3);quantSIZE$rate[which(quantSIZE$rate=="Inf")]<-0
  quantSIZE$rate<-paste("Rate of change =",quantSIZE$rate)
  
  
  ##ggplot of all spp with polygons
  pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/all spp model outcomes - SIZE.pdf", width=12, height=7)
    p2<-ggplot(quantSIZE, aes(x=year, y=(mean.size)))+facet_wrap(~species, scales="free_y")+theme_bw()
    #p2<-p2+geom_polygon(data=pols.size, aes(x=year, y=y, fill=factor(rate)), alpha=0.2)
    p2<-p2+geom_line(aes(col=factor(rate)), size=1)+geom_vline(xintercept=unique(dd$hold.year), col="black", lty="solid")+scale_colour_manual(values = (viridis(7, option = "B")[1:6]), name="Rate of decline in fishing pressure")
    p2+ylab("Mean body size (g)")+xlab("Year")
  dev.off()
########################################################################################################################################################  
  ## for SD body size
  ##make similar figures to julia's, with multiple rates on a plot
  dd$log.sd.size<-(dd$sd.size)
  quantSIZEsd<-melt(tapply(dd$sd.size,list(dd$year,dd$species, dd$rate),quantile,0.5))#;quatSSB$year<-rownames(paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(quantSIZEsd)<-c("year", "species","rate","sd.size")
  
  l.quantSIZEsd<-melt(tapply(dd$sd.size,list(dd$year,dd$species, dd$rate),quantile,0.05))#[paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(l.quantSIZEsd)<-c("year", "species","rate","lower")
  u.quantSIZEsd<-melt(tapply(dd$sd.size,list(dd$year,dd$species, dd$rate),quantile,0.95))#[paste(c(yrange)),order(c(11,5,8,10,4,3,9,12,2,7,1,6))]
  names(u.quantSIZEsd)<-c("year", "species","rate","upper")
  
  quantSIZEsd$lower<-l.quantSIZEsd$lower
  quantSIZEsd$upper<-u.quantSIZEsd$upper
  
  #data for polygons
  pols.size.sd<-data.frame(year=c(quantSIZEsd$year,rev(quantSIZEsd$year)), species=c(as.character(quantSIZEsd$species), rev(as.character(quantSIZEsd$species))), rate=c(quantSIZEsd$rate, rev(quantSIZEsd$rate)),y=c(quantSIZEsd$lower, rev(quantSIZEsd$upper)))
  
  quantSIZEsd$rate<-round(1/quantSIZEsd$rate, 3);quantSIZEsd$rate[which(quantSIZEsd$rate=="Inf")]<-0
  quantSIZEsd$rate<-paste("Rate of change =",quantSIZEsd$rate)
  
  ##ggplot of all spp with polygons
  pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/all spp model outcomes - SIZE SD.pdf", width=12, height=7)
  p2<-ggplot(quantSIZEsd, aes(x=year, y=(sd.size)))+facet_wrap(~species, scales="free_y")+theme_bw()
  #p2<-p2+geom_polygon(data=pols.size, aes(x=year, y=y, fill=factor(rate)), alpha=0.2)
  p2<-p2+geom_line(aes(col=factor(rate)), size=1)+geom_vline(xintercept=unique(dd$hold.year), col="black", lty="solid")+scale_colour_manual(values = (viridis(7, option = "B")[1:6]), name="Rate of decline in fishing pressure")
  p2+ylab("SD body size")+xlab("Year")
  dev.off()
  ########################################################################################################################################################
save(dd, file="~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/dd.RData")
