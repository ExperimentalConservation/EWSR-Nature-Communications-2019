##############################################################################################################################
## ews analysis
rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape2)
library(rlist)
library(gtools)
library(gridExtra)
library(RColorBrewer)
library(ggthemes)
library(scales)
library(viridis)
##############################################################################################################################
##load in the data from the different sims, this is biomass etc data
load(file="~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/dd.RData")

##source in the gam smoothing file
source('~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Gam smoothing.R', chdir = TRUE)

##and the weighted EWS file from bioaXiv, we will later set the weightings to 1 , but the code is generalised:
source('~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Weighted composite ews.R', chdir = TRUE)

##function for rolling mean
rolling_mean <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- mean(x[1:i], na.rm=T);
  }    
  return(result);
}
##function for rolling sd
rolling_sd <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- sd(x[1:i], na.rm=T);
  }    
  return(result);
}
##############################################################################################################################
##look at it
head(dd)

##lets just take cod first of all
cod<-subset(dd, species=="Cod")
head(cod)
ggplot(cod, aes(x=year, y=biomass, group=iter, col=rate))+geom_line(alpha=0.2)+facet_wrap(~rate)+theme_bw()+geom_vline(xintercept=2010, col="maroon")+geom_vline(xintercept=2040, col="royalblue4")

##we need to calculate when the species have "recovered", how do we set this?
##by the mean pop size? with gam fitted to it?
quantNUM<-melt(tapply(dd$biomass,list(dd$year,dd$species, dd$rate),quantile,0.5))
names(quantNUM)<-c("year", "species","rate","biomass")
quantNUM<-subset(quantNUM, species=="Cod")

quantSIZE<-melt(tapply(dd$log.mean.size,list(dd$year,dd$species, dd$rate),quantile,0.5))
names(quantSIZE)<-c("year", "species","rate","mean.size")
quantSIZE<-subset(quantSIZE, species=="Cod")
quantNUM$size<-quantSIZE$mean.size

quantSIZEsd<-melt(tapply(dd$log.sd.size,list(dd$year,dd$species, dd$rate),quantile,0.5))
names(quantSIZEsd)<-c("year", "species","rate","sd.size")
quantSIZEsd<-subset(quantSIZEsd, species=="Cod")
quantNUM$sd.size<-quantSIZEsd$sd.size

quantNUM$rate<-round(1/quantNUM$rate, 3);quantNUM$rate[which(quantNUM$rate=="Inf")]<-0

p1<-ggplot(quantNUM, aes(x=year, y=biomass))+theme_bw()+geom_vline(xintercept=dd$hold.year[1], col="black", lty="dashed")+theme(legend.position="top")
p1<-p1+geom_line(aes(col=factor(rate)), size=1, alpha=1)+scale_fill_brewer(palette="Spectral", name = "Rate of change")
p1<-p1+scale_colour_manual(values = (viridis(7, option = "B")[1:6]), name="Rate of decline in fishing pressure")
p1<-p1+ylab("Mean biomass (tonnes)")+xlab("Year")+xlim(1967, 2200)+annotate("text", label="a.", x=1967, y=max(quantNUM$biomass)*0.95)
cols = rev(viridis(7, option = "B")[1:6])

p1.1<-ggplot(quantNUM, aes(x=year, y=size/1000))+theme_bw()+geom_vline(xintercept=dd$hold.year[1], col="black", lty="dashed")+theme(legend.position="top")
p1.1<-p1.1+geom_line(aes(col=factor(rate)), size=1, alpha=1)
p1.1<-p1.1+ylab("Mean body mass (kg)")+xlab("Year")+xlim(1967, 2200)+annotate("text", label="b.", x=1967, y=max(quantNUM$size/1000)*0.95)
p1.1<-p1.1+scale_colour_manual(values = (viridis(7, option = "B")[1:6]), name="Rate of decline in fishing pressure")

p1.3<-ggplot(quantNUM, aes(x=year, y=sd.size/1000))+theme_bw()+geom_vline(xintercept=dd$hold.year[1], col="black", lty="dashed")+theme(legend.position="top")
p1.3<-p1.3+geom_line(aes(col=factor(rate)), size=1, alpha=1)+annotate("text", label="c.", x=1967, y=max(quantNUM$sd.size/1000)*0.95)
p1.3<-p1.3+ylab("Mean stamdard deviation of body mass (kg)")+xlab("Year")+xlim(1967, 2200)#+geom_smooth(aes(group=rate, col=factor(rate)))
p1.3<-p1.3+scale_colour_manual(values = (viridis(7, option = "B")[1:6]), name="Rate of decline in fishing pressure")
  #scale_color_brewer(palette="Spectral", name = "Time frame over which fishing pressure declines")

#####lets just do cod:
dd1<-subset(dd, species=="Cod" & year>=2010)

##across the various rates of change, split and calcualte the inflection point for the mean population.
## we are going to do this by looking at the movement away from its historic baseline (between 2010 and 2040 where the fishing effort is fixed)
##

dd1.inc.rec<-mclapply(split(dd1, list(dd1$rate, dd1$iter)), function(x){
                #x<-split(dd1, list(dd1$rate, dd1$iter))[[50]]
                if(!is.numeric(x$reimpliment.fishing[1])){
                  breaks <- x$year[which(x$year >= x$hold.year[1] & x$year <= x$hold.year[1]+x$rate[1]*2)]
                  mse <- numeric(length(breaks))
                  for(i in 1:length(breaks)){
                    piecewise1 <- lm(x$biomass ~ (x$year < breaks[i]) + (x$year>=breaks[i]))
                    mse[i] <- summary(piecewise1)[6]
                  }
                  mse <- as.numeric(mse)
                  x$recovery.point<-breaks[which(mse==min(mse))]
                }else{
                  breaks.1 <- x$year[which(x$year >= x$hold.year[1] & x$year <= x$hold.year[1]+x$rate[1]*2)]
                  breaks.2 <- x$year[which(x$year >= round(x$reimpliment.fishing[1]*0.99) & x$year <= x$reimpliment.fishing[1]+x$rate[1]*2)]
                  all.breaks<-expand.grid(breaks.1, breaks.2);names(all.breaks)<-c("Break1", "Break2")
                  mse <- numeric(length(all.breaks$Break1))
                  for(i in 1:length(all.breaks$Break1)){
                    piecewise1 <- lm(x$biomass ~ (x$year < all.breaks$Break1[i]) + (x$year>=all.breaks$Break1[i] & x$year<all.breaks$Break2[i]))
                    mse[i] <- summary(piecewise1)[6]
                  }
                  mse <- as.numeric(mse)
                  both.breaks<-all.breaks[which(mse==min(mse)),]
                  x$recovery.point<-both.breaks$Break1
                  x$collapse.point<-both.breaks$Break2
                }
                #cat(breaks[which(mse==min(mse))], "\n")
                return(x)
              }, mc.cores=detectCores()-1)

##violin plot of recovery points:
plot.violin<-unique(rbindlist(dd1.inc.rec)[,c("iter", "rate", "recovery.point")])
plot.violin<-subset(plot.violin, rate!=0)
#cols = rev(colorRampPalette(brewer.pal(11,"Spectral"))(6)[2:6])
cols = rev(viridis(7, option = "B")[2:6])

#quantNUM$rate<-round(1/quantNUM$rate, 3);quantNUM$rate[which(quantNUM$rate=="Inf")]<-0
p1.2<-ggplot(plot.violin, aes(x=factor(round(1/rate, 3)), y=recovery.point, fill=factor(rate), group=factor(rate)), drop=T)
  p1.2<-p1.2+geom_dotplot(aes(col=factor(rate)),fill="white", binaxis='y', stackdir='center', binwidth = 0.1)+theme_bw()#geom_violin(alpha=0.7)+theme_bw()
p1.2<-p1.2+theme(legend.position="top")+ylab("Time of recovery")+xlab("Rate of decline in fishing pressure")
p1.2<-p1.2+scale_colour_manual(values=cols, name = "Rate of decline in fishing pressure")+annotate("text", label="d.", x=0.8, y=2093)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1.1)

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/cod dynamics with recovery times.pdf", width=13, height=4)
p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p1.1 + theme(legend.position="none"),
                               p1.3 + theme(legend.position="none"),
                               p1.2 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(6, 1))
dev.off()

##plot an example of broken stick:
x<-subset(dd1, rate==50 & iter==210)
breaks <- x$year[which(x$year >= x$hold.year[1] & x$year <= x$hold.year[1]+x$rate[1]*2)]
mse <- numeric(length(breaks))
for(i in 1:length(breaks)){
  piecewise1 <- lm(x$biomass ~ (x$year < breaks[i]) + (x$year>=breaks[i]))
  mse[i] <- summary(piecewise1)[6]
}

i.num<-which(unlist(mse)==min(unlist(mse)))
piecewise1 <- lm(x$biomass ~ (x$year < breaks[i.num]) + (x$year>=breaks[i.num]))
cfs <- as.data.frame(coef(summary(piecewise1)))$Estimate
#pols<-
  pp1<-ggplot(x, aes(x=year, y=biomass))+geom_line(col="NA")+ylab("Biomass")+xlab("Year")+theme_bw()
  pp1<-pp1+annotate("rect",xmin=2010, xmax=2040, ymin=0, ymax=2500000, alpha=0.2, fill="darkorchid1") 
  pp1<-pp1+annotate("rect", xmin=2040, xmax=x$hold.year[1]+i.num, ymin=0, ymax=2500000, alpha=0.3, fill="gold")
  pp1<-pp1+geom_segment(x=min(x$year), y=sum(cfs), xend=x$hold.year[1]+i.num, yend=sum(cfs))
  pp1<-pp1+geom_segment(x=x$hold.year[1]+i.num, y=cfs[1], xend=max(x$year), yend=cfs[1])
  pp1<-pp1+geom_vline(xintercept=x$hold.year[1]+i.num, lty="solid")+geom_vline(xintercept=x$hold.year[1], lty="dashed")
  pp1<-pp1+geom_vline(xintercept=x$hold.year[1]+x$rate[1], lty="solid", col="maroon")
  pp1<-pp1+geom_line()+xlim(2010, 2200)+scale_x_continuous(breaks=c(2010, 2040, 2100, 2200))
  
pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/broken stick.pdf", width=7, height=3)
  pp1
dev.off()

fin.res<-rbindlist(mclapply(dd1.inc.rec, function(y){
  
  ##make a dataframe of the data to calcualte ewss on. 
  ##sd.szie is fixed just so I dont have to change code in the function later
  ews.dat<-data.frame(year=y$year, biomass=y$biomass, mean.size=y$mean.size, sd.size=y$sd.size)
  #ews.dat<-subset(ews.dat, year>=2030)
  ##just look at the data from 2040 onwards, 
  #ews.dat<-subset(ews.dat, year>=y$hold.year[1]-10)
  
  ##list of ews to test independantly 
  metrics<-c("cv", "ar1", "mean.size", "sd.size")
  to.test.l<-list(NULL)
  for(jj in 1:length(metrics)){
  	#jj=1
	to.test.l[[jj]]<-split(combinations(n = length(metrics), r = jj, v = metrics, repeats.allowed = FALSE), seq(nrow(combinations(n = length(metrics), r = jj, v = metrics, repeats.allowed = FALSE))))
  }
  to.test<-unlist(to.test.l, recursive=FALSE)

  ##object to store results
  res<-NULL
  for(i in 1:length(to.test)){
    #i=10
    ##set the weighing to 1 - just required to make the code run, doesnt do anythnig
    W<-data.frame(inds=sort(unlist(to.test[i])), "wei"=1)
    ##run the EWS from clements & ozgul nat comes and save out:
    res[[i]]<-W_composite_ews(dat=ews.dat, indicators=sort(unlist(to.test[i])), weights=W, threshold=2, plotIt=F)
  }
  
  bind.res<-rbindlist(res)
  bind.res$species<-y$species[1]
  bind.res$iter<-y$iter[1]
  bind.res$rate<-y$rate[1]
  bind.res$hold.year<-y$hold.year[1]
  bind.res$recovery.point<-y$recovery.point[1]
  bind.res$biomass<-y$biomass[match(bind.res$time, y$year)]
  return(bind.res)
}, mc.cores=detectCores()-1))

##calculate the strength of ewss
fin.res$str<-(fin.res$metric.score-fin.res$rolling.mean)/fin.res$rolling.sd

##save out the ews output:
save(fin.res,file="~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/ews.res.RData")
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##then run again with varying time series length:
## ews analysis
rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape2)
library(rlist)
library(gtools)
library(gridExtra)
library(RColorBrewer)
##############################################################################################################################
##load in the data from the different sims, this is biomass etc data
load(file="~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/dd.RData")

##source in the gam smoothing file
source('~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Gam smoothing.R', chdir = TRUE)

##and the weighted EWS file from bioaXiv, we will later set the weightings to 1 , but the code is generalised:
source('~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Weighted composite ews.R', chdir = TRUE)

##function for rolling mean
rolling_mean <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- mean(x[1:i], na.rm=T);
  }    
  return(result);
}
##function for rolling sd
rolling_sd <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- sd(x[1:i], na.rm=T);
  }    
  return(result);
}
##############################################################################################################################
#####lets just do cod:
dd1<-subset(dd, species=="Cod" & year>=2010)

##across the various rates of change, split and calcualte the inflection point for the mean population.
## we are going to do this by looking at the movement away from its historic baseline (between 2010 and 2040 where the fishing effort is fixed)
##

dd1.inc.rec<-mclapply(split(dd1, list(dd1$rate, dd1$iter)), function(x){
  #x<-split(dd1, list(dd1$rate, dd1$iter))[[50]]
  if(!is.numeric(x$reimpliment.fishing[1])){
    breaks <- x$year[which(x$year >= x$hold.year[1] & x$year <= x$hold.year[1]+x$rate[1]*2)]
    mse <- numeric(length(breaks))
    for(i in 1:length(breaks)){
      piecewise1 <- lm(x$biomass ~ (x$year < breaks[i]) + (x$year>=breaks[i]))
      mse[i] <- summary(piecewise1)[6]
    }
    mse <- as.numeric(mse)
    x$recovery.point<-breaks[which(mse==min(mse))]
  }else{
    breaks.1 <- x$year[which(x$year >= x$hold.year[1] & x$year <= x$hold.year[1]+x$rate[1]*2)]
    breaks.2 <- x$year[which(x$year >= round(x$reimpliment.fishing[1]*0.99) & x$year <= x$reimpliment.fishing[1]+x$rate[1]*2)]
    all.breaks<-expand.grid(breaks.1, breaks.2);names(all.breaks)<-c("Break1", "Break2")
    mse <- numeric(length(all.breaks$Break1))
    for(i in 1:length(all.breaks$Break1)){
      piecewise1 <- lm(x$biomass ~ (x$year < all.breaks$Break1[i]) + (x$year>=all.breaks$Break1[i] & x$year<all.breaks$Break2[i]))
      mse[i] <- summary(piecewise1)[6]
    }
    mse <- as.numeric(mse)
    both.breaks<-all.breaks[which(mse==min(mse)),]
    x$recovery.point<-both.breaks$Break1
    x$collapse.point<-both.breaks$Break2
  }
  #cat(breaks[which(mse==min(mse))], "\n")
  return(x)
}, mc.cores=detectCores()-1)

fin.res.data.length<-rbindlist(mclapply(dd1.inc.rec, function(y){
  
  #y<-dd1.inc.rec[[18]]
  #plot(y$year,y$biomass, type="l")
  #abline(v=y$recovery.point)
  
  ##make a dataframe of the data to calcualte ewss on. 
  ##sd.szie is fixed just so I dont have to change code in the function later
  ews.dat<-data.frame(year=y$year, biomass=y$biomass, mean.size=y$mean.size, sd.size=y$sd.size)
  
  ##just look at the data from 2040 onwards, 
  #ews.dat<-subset(ews.dat, year>=y$hold.year[1]-10)
  
  ##list of ews to test independantly 
  metrics<-c("cv", "ar1", "mean.size", "sd.size")
  to.test.l<-list(NULL)
  for(jj in 1:length(metrics)){
    #jj=1
    to.test.l[[jj]]<-split(combinations(n = length(metrics), r = jj, v = metrics, repeats.allowed = FALSE), seq(nrow(combinations(n = length(metrics), r = jj, v = metrics, repeats.allowed = FALSE))))
  }
  to.test<-unlist(to.test.l, recursive=FALSE)
  
  ##object to store results
  res<-NULL
  ##counter
  counter<-0
  ##loop to increase length of data
  for(l in (2010:2039)){
    #l=2011
    ##subset the data to vary the length to analyse
    ews.dat.sub<-subset(ews.dat, year>=l)
    ##loop through all the metrics
    for(i in 1:length(to.test)){
      counter<-counter+1
      #i=3
      ##set the weighing to 1 - just required to make the code run, doesnt do anythnig
      W<-data.frame(inds=sort(unlist(to.test[i])), "wei"=1)
      ##run the EWS from clements & ozgul nat comes and save out:
      ews.sub.res<-W_composite_ews(dat=ews.dat.sub, indicators=sort(unlist(to.test[i])), weights=W, threshold=2.5, plotIt=F)
      ews.sub.res$data.length<-2040-ews.dat.sub$year[1]
      res[[counter]]<-ews.sub.res
    }
  }
  bind.res<-rbindlist(res)
  bind.res$species<-y$species[1]
  bind.res$iter<-y$iter[1]
  bind.res$rate<-y$rate[1]
  bind.res$hold.year<-y$hold.year[1]
  bind.res$recovery.point<-y$recovery.point[1]
  bind.res$biomass<-y$biomass[match(bind.res$time, y$year)]
  return(bind.res)
}, mc.cores=detectCores()-1))

##calculate the strength of ewss
fin.res.data.length$str<-(fin.res.data.length$metric.score-fin.res.data.length$rolling.mean)/fin.res.data.length$rolling.sd
fin.res.data.length<-as.data.table(fin.res.data.length)
##save out the ews output:
save(fin.res.data.length,file="~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/fin.res.data.length.RData")
