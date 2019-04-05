##############################################################################################################################
rm(list=ls())
library(grid)
library(ggplot2)
library(data.table)
library(reshape2)
library(rlist)
library(parallel)
library(gtools)
library(gridExtra)
library(viridis)
library(RColorBrewer)
##############################################################################################################################
##load in the ews output:
load("~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/ews.res.RData")
##############################################################################################################################
##strength of signal analysis?
head(fin.res, 3)
levels(fin.res$metric.code)<-c("AR(1)", "CV", "Mean size", "SD size", "AR(1) + CV", "AR(1) + Mean size", "AR(1) + SD size", "CV + mean size", "CV + SD size", "Mean size + SD size", "AR(1) + CV + mean size", "AR(1) + CV + SD size", "AR(1) + mean size + SD size", "CV + mean size + SD size", "AR(1) + CV + mean size + SD size")

##lets plot the intervals as with the biomass.
##first calcualte the quantiles
quantMetric<-melt(tapply(fin.res$metric.score,list(fin.res$time,fin.res$rate, fin.res$metric.code),quantile,0.5))
  names(quantMetric)<-c("year","rate","code","strength")
l.quantMetric<-melt(tapply(fin.res$metric.score,list(fin.res$time,fin.res$rate, fin.res$metric.code),quantile,0.05))
  names(l.quantMetric)<-c("year","rate","code","strength")
u.quantMetric<-melt(tapply(fin.res$metric.score,list(fin.res$time,fin.res$rate, fin.res$metric.code),quantile,0.95))
  names(u.quantMetric)<-c("year","rate","code","strength")
quantMetric$lower<-l.quantMetric$strength
quantMetric$upper<-u.quantMetric$strength

#data for polygons
pols.metric<-data.frame(year=c(quantMetric$year,rev(quantMetric$year)), rate=c(quantMetric$rate, rev(quantMetric$rate)), code=c(as.character(quantMetric$code), rev(as.character(quantMetric$code))), y=c(quantMetric$lower, rev(quantMetric$upper)))

##add in inc.size vector for plotting
quantMetric<-rbindlist(lapply(split(quantMetric, list(quantMetric$code)), function(HH){
  #HH<-split(quantMetric, list(quantMetric$code))[[10]]
  if(length(grep("size",HH$code[1], value = TRUE))>0){
    HH$inc.size<-T
  }else{HH$inc.size<-F}
  return(HH)
}))


quantMetric$fishing.end<-fin.res$hold.year[1]+quantMetric$rate
quantMetric$fishing.end[which(quantMetric$rate=="Rate of change = 0")]<-NA

quantMetric$rate<-round(1/quantMetric$rate, 3);quantMetric$rate[which(quantMetric$rate=="Inf")]<-0
quantMetric$rate<-paste("Rate of decline in fishing pressure =",quantMetric$rate)

##change so that there are 3 levels of line type (size only, abundance only, size and abundance) for plots
quantMetric$line.type<-ifelse(quantMetric$inc.size==FALSE, 0, 1)
quantMetric$size.only<-ifelse(match(quantMetric$code, c("Mean size", "SD size", "Mean size + SD size"))>0, 1, 0)
quantMetric$size.only[which(is.na(quantMetric$size.only))]<-0
quantMetric$line.type.plot<-quantMetric$line.type+quantMetric$size.only
line.types<-data.frame(code=c(1:3), line=c("dotdash", "dotted", "solid"))
quantMetric$line.type.plot<-as.factor(quantMetric$line.type.plot)
levels(quantMetric$line.type.plot)<-c("Abundance only", "Both", "Size only")
facet_labels<-data.frame(rate=unique(quantMetric$rate),
                         label=c("a.", "f.", "e.", "d.", "c.", "b."))
quantMetric$fishing.end[which(quantMetric$rate=="Rate of decline in fishing pressure = 0")]<-NA
black.vertical.lines<-data.frame(rate=unique(quantMetric$rate), 
                                 int=c(NA, rep(2040, length(unique(quantMetric$rate))-1)))

##ggplot of metric strength  with polygons
p2<-ggplot(quantMetric, aes(x=year, y=strength))+theme_bw()+xlim(min(quantMetric$year), 2100)
p2<-p2+geom_vline(data=black.vertical.lines, aes(xintercept=int), col="black", size=1)
p2<-p2+geom_line(aes(col=(code), lty=factor(line.type.plot)), size=1)+geom_hline(yintercept=2, size=1, lty="dashed")+theme(strip.background = element_rect(fill=c("mintcream")))
p2<-p2+facet_wrap(~rate, ncol=3)+xlab("Year")+ylab("Mean strength of signal")+theme(legend.position = "right")+geom_vline(aes(xintercept=fishing.end), col="maroon", size=1)
p2<-p2+scale_linetype_manual(values = c("solid", "dotted", "dashed"))+theme(legend.text = element_text(size=1))+theme(legend.key.size =  unit(0.2, "in"), legend.key.width = unit(2, "line"))
#p2<-p2+geom_text(data=facet_labels, aes(label=label), x = 2015, y = 9.4)
pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/strength of all ews over time.pdf", width=13, height=6)
  p2+labs(color='Metric', lty="Includes size?")+theme(legend.text=element_text(size=9))+theme(legend.direction = "vertical", legend.key=element_rect(size=0.01))
dev.off()

##re add info for analysis
quantMetric$fishing.end[which(quantMetric$rate=="Rate of decline in fishing pressure = 0")]<-2040

 ##for a select few EWS
sel.codes<-c("CV", "mean size", "AR(1)", "SD size")
sel.res<-quantMetric[which(quantMetric$code%in%sel.codes), ]

##ok so we can see that there is a big rise in the indicators before the dashed lines (when the fishing pressures hit their lowest point)
##however, what we want to know is if there is a big increase before the "recovery point" as determined by the broken stick regression

##pull out the control sims
null.mod<-as.data.table(subset(fin.res, rate==0))

##the other data that isnt from the null model (i.e. the recovery data):
non.null<-as.data.table(subset(fin.res, rate!=0))
non.null<-fin.res

##############################################################################################################################
##so lets look at the strength of those signals over time in the recovery vs over exploited:

##return the row with the highest strength for each metric in the period 2050-recovery
ews.strengths<-rbindlist(mclapply(split(fin.res, list(fin.res$rate, fin.res$iter, fin.res$metric.code)), function(k){
  #k<-split(fin.res, list(fin.res$rate, fin.res$iter, fin.res$metric.code))[[50]]
  if(k$rate[1]!=0){
    k<-subset(k, time>=k$hold.year[1] & time<recovery.point)
  }else{
    k<-subset(k, time>=k$hold.year[1] & time<2090)
  }
  return(k[which(k$str==max(k$str)),])
  #return(k)
}, mc.cores=detectCores()-1))

##so look at the mean strength of the signals across these:
mean.strengths<-tapply(ews.strengths$str, list(ews.strengths$rate, ews.strengths$metric.code), mean)

##calcuate standard error for them too to add on error bars:
se.strengths<-melt(tapply(ews.strengths$str, list(ews.strengths$rate, ews.strengths$metric.code), function(x) sd(x)/sqrt(length(x))))
mean.strengths<-melt(mean.strengths);names(mean.strengths)<-c("rate", "code", "mean.strength")
mean.strengths$se<-se.strengths$value

##add in a colour vector:
mean.strengths$col<-1
mean.strengths$col[which(mean.strengths$rate==0)]<-0

##plot all mean strengths
mean.strength.plots<-mean.strengths
mean.strength.plots$rate<-round(1/mean.strength.plots$rate, 3);mean.strength.plots$rate[which(mean.strength.plots$rate=="Inf")]<-0

pp2.all<-ggplot(mean.strength.plots, aes(x=factor(rate), y=mean.strength))+geom_bar(stat="identity", aes(fill=factor(col)), col="darkgrey")+facet_wrap(~code, ncol=3)+theme_bw()+geom_hline(yintercept=2)
pp2.all<-pp2.all+geom_errorbar(aes(ymin=mean.strength-se, ymax=mean.strength+se), width=.3)+theme(strip.background = element_rect(fill=c("mintcream")))
pp2.all<-pp2.all+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))+theme(legend.position="none")
pp2.all<-pp2.all+ylab("Mean maximum strength of signal")+xlab("Rate of decline in fishing pressure")#+ggtitle("Mean strength of signal over pre-recovery period")

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/mean max strength ALL metrics.pdf", width=9, height=12)
  pp2.all
dev.off()

##select out the best, max, and worst (strongest and weakest) signals:
ave.mean.strengths<-melt(tapply(subset(mean.strengths, rate!=0)$mean.strength, list(subset(mean.strengths, rate!=0)$code), mean))
mean.strengths$ave<-ave.mean.strengths$value[match(mean.strengths$code, ave.mean.strengths$indices)]
best<-subset(mean.strengths, code==mean.strengths$code[which(mean.strengths$ave==max(mean.strengths$ave))])
max<-subset(mean.strengths, code==mean.strengths$code[which(mean.strengths$mean.strength==max(mean.strengths$mean.strength))])

#mean.strengths$diff<-
no.change<-subset(mean.strengths, rate==0)
change<-melt(tapply(subset(mean.strengths, rate!=0)$mean.strength, subset(mean.strengths, rate!=0)$code, mean))
change$no.change<-no.change[order(no.change$code),]$mean.strength
change$diff<-change$value-change$no.change
biggest.diff.strength<-subset(change, diff==max(change$diff))$Var1

sel.mean.strengths<-subset(mean.strengths, code==biggest.diff.strength)
sel.mean.strengths<-unique(sel.mean.strengths)
sel.mean.strengths$rate<-round(1/sel.mean.strengths$rate, 3);sel.mean.strengths$rate[which(sel.mean.strengths$rate=="Inf")]<-0
##plot it
pp2<-ggplot(sel.mean.strengths, aes(x=factor(rate), y=mean.strength))+geom_bar(stat="identity", aes(fill=factor(col)), col="darkgrey")+facet_wrap(~code, ncol=3)+theme_bw()+geom_hline(yintercept=2)
pp2<-pp2+geom_errorbar(aes(ymin=mean.strength-se, ymax=mean.strength+se), width=.3)+theme(strip.background = element_rect(fill=c("mintcream")))
pp2<-pp2+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))+theme(legend.position="none")+annotate("text", x=0.8, y=2.7, label="a.")
pp2<-pp2+ylab("Mean maximum strength of signal")+xlab("Rate of decline in fishing pressure")#+ggtitle("Mean strength of signal over pre-recovery period")
pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/mean max strength.pdf", width=3*length(unique(sel.mean.strengths$code)), height=3)
  pp2
dev.off()

##############################################################################################################################
##so this is for a single EWS based on a threshold
thres.sigs<-2

##and the 95 quantile of the rate=0 strength of signals as a second threshold
rate.zero.95q<-quantile(subset(fin.res, rate==0)$metric.score, 0.95)

##add in whether the signal has passed a 2 sigma threshold:
null.mod$threshold.crossed.2sigma<-NA
null.mod$threshold.crossed.2sigma[which(null.mod$metric.score>(null.mod$rolling.mean+(thres.sigs*null.mod$rolling.sd)))]<-1
##do the same for the non-null
non.null$threshold.crossed.2sigma<-NA
non.null$threshold.crossed.2sigma[which(non.null$metric.score>(non.null$rolling.mean+(thres.sigs*non.null$rolling.sd)))]<-1

##add in whether the signal has passed a 95th quantile of zero rate of change threshold:
null.mod$threshold.crossed.95sigma<-NA
null.mod$threshold.crossed.95sigma[which(null.mod$metric.score>(null.mod$rolling.mean+(rate.zero.95q*null.mod$rolling.sd)))]<-1
##do the same for the non-null
non.null$threshold.crossed.95sigma<-NA
non.null$threshold.crossed.95sigma[which(non.null$metric.score>(non.null$rolling.mean+(rate.zero.95q*non.null$rolling.sd)))]<-1

##plot out an example for broken stick

##collapse the data to single earliest warning signal which is after the start of the release of fishing pressures, 
##and before the point of recovery 
ews.sigs<-rbindlist(mclapply(split(non.null, list(non.null$rate, non.null$iter, non.null$metric.code)), function(k){
  ##k<-split(non.null, list(non.null$rate, non.null$iter, non.null$metric.code))[[2000]]
 #k<-subset(non.null, metric.code=="AR1 + SD size" & iter==103 & rate==50)
  if(k$rate[1]!=0){
    k<-subset(k, time>=k$hold.year[1] & time<recovery.point)
  }else{
    k<-subset(k, time>=k$hold.year[1] & time<2090)
  }
  sub.k<-subset(k, threshold.crossed.2sigma==1)
  
  if(nrow(sub.k)==0){
    sub.k<-k[1,]
    sub.k$ews<-NA
  }else{
    sub.k<-sub.k[1,]
    sub.k$ews<-sub.k$time
  }
  return(sub.k)
}, mc.cores=detectCores()-1))

##calculate the proportion that show ewss at 2sigma prior to the recovery point
plot.props<-melt(tapply(ews.sigs$ews, list(ews.sigs$rate, ews.sigs$metric.code), function(x)length(na.omit(x))/length(x)))
names(plot.props)<-c("rate", "code", "proportion")
plot.props$col<-1
plot.props$col[which(plot.props$rate==0)]<-0

##plot all metrics
plot.props.all<-plot.props
plot.props.all$rate<-round(1/plot.props.all$rate, 3);plot.props.all$rate[which(plot.props.all$rate=="Inf")]<-0

pp1.all<-ggplot(plot.props.all, aes(x=factor(rate), y=proportion))+geom_bar(stat="identity", aes(fill=factor(col)), col="black")+facet_wrap(~code, ncol=3)
pp1.all<-pp1.all+theme_bw()#+ggtitle(bquote(list(atop("Proportion of time series showing at least one EWSR at a 2" ~sigma ~ "threshold"))))
pp1.all<-pp1.all+xlab("Rate of decline in fishing pressure")+ylab("Proportion showing EWSR")+theme(strip.background = element_rect(fill=c("mintcream")))
#pp1<-pp1+geom_hline(aes(yintercept=min(proportion)))

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/proportion showing signals ALL METRICS.pdf", width=9, height=12)
  pp1.all+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))+theme(legend.position="none")
dev.off()

##so the maximum proportion is 0.9, pretty good!
max(tapply(ews.sigs$ews, list(ews.sigs$rate, ews.sigs$metric.code), function(x)length(na.omit(x))/length(x)))

##find the 5 metrics with the highest mean proportion of signals:
ews.sigs.no.zero<-subset(ews.sigs, rate!=0)
props<-melt(tapply(ews.sigs.no.zero$ews, list(ews.sigs.no.zero$metric.code), function(x)length(na.omit(x))/length(x)))
names(props)<-c("Var1", "value")

##add the mean values (excluding 0 rate of change) to the full data set
plot.props$mean.prop<-props$value[match(plot.props$code, props$Var1)]
plot.props<-rbindlist(lapply(split(plot.props,list(plot.props$code)), function(h){
  #h<-split(plot.props,list(plot.props$code))[[1]]
  h$diff<-h$mean.prop[1]-h$proportion[which(h$rate==0)]
  return(h)
}))

##  
best.props<-props$Var1[which(props$value==max(props$value))]
biggest.diff<-unique(plot.props$code[which(plot.props$diff==max(plot.props$diff))])
#top.5<-c(as.character(biggest.diff), as.character(best.props))
top.5<-biggest.diff

plot.props$rate<-round(1/plot.props$rate, 3);plot.props$rate[which(plot.props$rate=="Inf")]<-0

##plot out the proportion of populations exhibiting EWSs
pp1<-ggplot(plot.props[which(plot.props$code%in%top.5),], aes(x=factor(rate), y=proportion))+geom_bar(stat="identity", aes(fill=factor(col)), col="black")+facet_wrap(~code, ncol=length(top.5))
pp1<-pp1+theme_bw()#+ggtitle(bquote(list(atop("Proportion of time series showing at least one EWSR at a 2" ~sigma ~ "threshold"))))
pp1<-pp1+xlab("Rate of decline in fishing pressure")+ylab("Proportion showing EWSR")+theme(strip.background = element_rect(fill=c("mintcream")))+annotate("text", x=0.8, y=0.95, label="b.")
#pp1<-pp1+geom_hline(aes(yintercept=min(proportion)))
pp1<-pp1+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))+theme(legend.position="none")+ylim(0,1)
pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/proportion showing signals.pdf", width=4*length(unique(plot.props[which(plot.props$code%in%top.5),]$code)), height=3)
  pp1
dev.off()

##distribution of timings of EWSs prior to the recovery of the system
thres.sigs.dist<-2
ews.sigs$years.prior.recovery<-ews.sigs$recovery.point-ews.sigs$ews
plot.sigs<-subset(ews.sigs, rate!=0 & metric.code==biggest.diff)
plot.sigs$rate<-round(1/plot.sigs$rate, 3);plot.sigs$rate[which(plot.sigs$rate=="Inf")]<-0

cols = rev(colorRampPalette(brewer.pal(11,"Spectral"))(6)[2:6])
vio1<-ggplot(plot.sigs, aes(x=factor(rate), y=years.prior.recovery, fill=factor(rate), group=rate), drop=T)
vio1<-vio1+geom_dotplot(aes(col=factor(rate)),fill="white", binaxis='y', stackdir='center', binwidth = 0.4)
vio1<-vio1+ylab("Years prior to recovery")+xlab("Rate of decline in fishing pressure")
vio1<-vio1+scale_fill_manual(values=cols, name = "Rate of decline in fishing pressure")
vio1<-vio1+theme_bw()+theme(legend.position="none")+annotate("text", x=0.8, y=44, label="c.")
vio1<-vio1+theme(strip.background = element_rect(fill=c("mintcream")))+facet_wrap(~metric.code)+scale_colour_manual(values = (viridis(7, option = "B")[2:7]), name="Rate of decline in fishing pressure")
#pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/distribution of signals prior to recovery.pdf", width=4, height=3)
  vio1+scale_colour_manual(values = (viridis(7, option = "B")[2:7]), name="Rate of decline in fishing pressure")
#dev.off()

#######Figure 3 for the MS
pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/Fig. 3.pdf", width=12, height=3)
  grid.arrange(pp2, pp1, vio1, ncol=3)
dev.off()
##############################################################################################################################
##consecutive EWSs analysis?

##set a sigma
threshold<-2

##calculate the longest consecutive run of signals
mean.ews_t<-rbindlist(mclapply(split(fin.res, list(fin.res$rate, fin.res$iter, fin.res$metric.code)), function(k, threshold){
        #k<-split(fin.res, list(fin.res$rate, fin.res$iter, fin.res$metric.code))[[2000]]
        #k<-subset(fin.res, metric.code=="AR1 + SD size" & iter==103 & rate==50)
        if(k$rate[1]!=0){
          k<-subset(k, time>=k$hold.year[1] & time<recovery.point)
        }else{
          k<-subset(k, time>=k$hold.year[1] & time<2090)
        }
        k$threshold.crossed<-0
        k$threshold.crossed[which(k$str>2)]<-1
        ##give each run of sigals (ews or not) a unique code to split them by
        k$seq.code<-0
        ##a counter
        counter<-0
        ##assuming there are more than 1 rows in the data:
        if(length(k$time)>1){
        for(i in 2:length(k$threshold.crossed)){
          #i=2
          if(k$threshold.crossed[i]==k$threshold.crossed[i-1]){
            k$seq.code[i]<-counter
          }else{
            counter<-counter+1
            k$seq.code[i]<-counter
          }
        }}else{
          k$seq.code<-counter
        }
        
        ##add in a vector for the timing of signals:
        k$sig.time<-0
        ##calculate longet run of positive signals:
        sig.k.l<-rbindlist(lapply(split(k, k$seq.code), function(o){
          #o<-split(k, k$seq.code)[[2]]
          if(o$threshold.crossed[1]==1){
            o$max.sig.length<-length(o$seq.code)
            if(length(o$seq.code)>=2)
              o$sig.time[2]<-1
            }else{
              o$max.sig.length<-0
            }
          return(o)
        }))
        ##add in maximum number of consecutive signals
        sig.k.l$max.sig.length<-max(sig.k.l$max.sig.length)
        
        ##calcualte the eariest >2 consecutive signals
        earliest.ews<-min(subset(sig.k.l, sig.time==1)$time)
        sig.k.l$earliest.sig<-earliest.ews
        
        
        
  return(sig.k.l[1,])
}, threshold=threshold, mc.cores=detectCores()-2))

##calcualte the mean number of consecutive signals
con.sigs<-melt(tapply(mean.ews_t$max.sig.length, list(mean.ews_t$metric.code, mean.ews_t$rate), mean));names(con.sigs)<-c("metric", "rate", "mean.con.sigs")
head(con.sigs)

##and the standard error of this
con.sigs.se<-melt(tapply(mean.ews_t$max.sig.length, list(mean.ews_t$metric.code, mean.ews_t$rate), function(x) sd(x)/sqrt(length(x))));names(con.sigs.se)<-c("metric", "rate", "se.con.sigs")
con.sigs$se<-con.sigs.se$se.con.sigs

##plot all of them
plot.con<-con.sigs
plot.con$rate<-round(1/plot.con$rate, 3);plot.con$rate[which(plot.con$rate=="Inf")]<-0
plot.con$col<-1
plot.con$col[which(plot.con$rate==0)]<-0

p5<-ggplot(plot.con, aes(x=factor(rate), y=mean.con.sigs))+geom_bar(stat="identity", aes(fill=factor(col)), col="darkgrey")
p5<-p5+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))
p5<-p5+theme_bw()+theme(legend.position = "none")
p5<-p5+geom_errorbar(aes(ymin=mean.con.sigs-se, ymax=mean.con.sigs+se), width=.3)
p5<-p5+xlab("Rate of decline in fishing pressure")+ylab("Mean number of consecutive EWSR")#+ggtitle(bquote(list("Mean number of consecutive EWSR\nat a 2 sigma  threshold")))
p5<-p5+facet_wrap(~metric, ncol=3)+theme(strip.background = element_rect(fill=c("mintcream")))

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/mean consecutive signals ALL METRICS.pdf", width=9, height=12)
  p5+theme(plot.margin=unit(c(1,0,0.3,0.3),"cm"))
dev.off()

##calcualte the mean per metric
head(con.sigs)
con.sigs.non.zero<-subset(con.sigs, rate!=0)
mean.across.rates<-melt(tapply(con.sigs.non.zero$mean.con.sigs, list(con.sigs.non.zero$metric), mean))
best.con.met<-mean.across.rates[which(mean.across.rates[,2]==max(mean.across.rates[,2])),1]

##mean across recovery treatments 
mean.line.plot<-mean.across.rates[which(mean.across.rates[,2]==max(mean.across.rates[,2])),2]

##subset out the best metric (highest mean across the recovery)
plot.best.con<-subset(con.sigs, metric==best.con.met)
plot.best.con$col<-1
plot.best.con$col[which(plot.best.con$rate==0)]<-0

plot.best.con$rate<-round(1/plot.best.con$rate, 3);plot.best.con$rate[which(plot.best.con$rate=="Inf")]<-0

p5<-ggplot(plot.best.con, aes(x=factor(rate), y=mean.con.sigs))+geom_bar(stat="identity", aes(fill=factor(col)), col="darkgrey")
p5<-p5+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))
p5<-p5+theme_bw()+theme(legend.position = "none")
p5<-p5+geom_errorbar(aes(ymin=mean.con.sigs-se, ymax=mean.con.sigs+se), width=.3)
#p5<-p5+geom_point(aes(col=factor(col)))
p5<-p5+xlab("Rate of decline in fishing pressure")+ylab("Mean number of consecutive EWSR")#+ggtitle(bquote(list("Mean number of consecutive EWSR\nat a 2 sigma  threshold")))
p5<-p5+facet_wrap(~metric)+theme(strip.background = element_rect(fill=c("mintcream")))
p5<-p5+geom_hline(yintercept = mean.line.plot, lty="dashed")+annotate("text", x=0.8, y=8.5, label="a.")

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/mean consecutive signals.pdf", width=4, height=3)
  p5
dev.off()

###then calculate the proportion showing EWSR with 2 consecutive signals:
best.con.sigs<-subset(mean.ews_t, metric.code==best.con.met)
best.con.sigs$threshold.con<-ifelse(best.con.sigs$max.sig.length<2, 0, 1)

plot.prop.con.sigs<-rbindlist(lapply(split(best.con.sigs, best.con.sigs$rate), function(jk){
  #jk<-split(no.col.con, best.con.sigs$rate)[[1]]
  return(data.frame(rate=jk$rate[1], proportion=sum(jk$threshold.con)/length(jk$threshold.con)))
}))

plot.prop.con.sigs$rate<-round(1/plot.prop.con.sigs$rate, 3);plot.prop.con.sigs$rate[which(plot.prop.con.sigs$rate=="Inf")]<-0
plot.prop.con.sigs$col<-1
plot.prop.con.sigs$col[which(plot.prop.con.sigs$rate==0)]<-0
plot.prop.con.sigs$metric<-best.con.sigs$metric.code[1]

##plot out the proportion of populations exhibiting EWSs
p5.prop<-ggplot(plot.prop.con.sigs, aes(x=factor(rate), y=proportion))+facet_wrap(~metric)+geom_bar(stat="identity", aes(fill=factor(col)), col="black")
#p5.prop<-p5.prop+geom_point(aes(col=factor(col)))
p5.prop<-p5.prop+theme_bw()#+ggtitle(bquote(list(atop("Proportion of time series showing at least one EWSR at a 2" ~sigma ~ "threshold"))))
p5.prop<-p5.prop+xlab("Rate of decline in fishing pressure")+ylab("Proportion of time series showing\nEWSR in at least 2 consecutive years")+theme(strip.background = element_rect(fill=c("mintcream")))
p5.prop<-p5.prop+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))+theme(legend.position="none")+ylim(0,1)+annotate("text", x=0.8, y=0.95, label="b.")

##############################################################################################################################

##calculate the proportion that show ewss at 2sigma prior to the recovery point
mean.ews_t$earliest.sig[which(mean.ews_t$earliest.sig=="Inf")]<-NA
plot.props.con<-melt(tapply(mean.ews_t$earliest.sig, list(mean.ews_t$rate, mean.ews_t$metric.code), function(x)length(na.omit(x))/length(x)))
names(plot.props.con)<-c("rate", "code", "proportion")
plot.props.con$col<-1
plot.props.con$col[which(plot.props$rate==0)]<-0

##so the maximum proportion is 0.86
max(tapply(mean.ews_t$earliest.sig, list(mean.ews_t$rate, mean.ews_t$metric.code), function(x)length(na.omit(x))/length(x)))

##find the 5 metrics with the highest mean proportion of signals:
ews.sigs.no.zero<-subset(mean.ews_t, rate!=0)
props<-melt(tapply(ews.sigs.no.zero$earliest.sig, list(ews.sigs.no.zero$metric.code), function(x)length(na.omit(x))/length(x)))
names(props)<-c("Var1", "value")

##add the mean values (excluding 0 rate of change) to the full data set
plot.props.con$mean.prop<-props$value[match(plot.props.con$code, props$Var1)]
#names(plot.props.con)<-c("rate", "code", "proportion")
plot.props<-rbindlist(lapply(split(plot.props.con,list(plot.props.con$code)), function(h){
  #h<-split(plot.props.con,list(plot.props.con$code))[[1]]
  h$diff<-h$mean.prop[1]-h$proportion[which(h$rate==0)]
  return(h)
}))

##  
best.props<-props$Var1[which(props$value==max(props$value))]
biggest.diff<-unique(plot.props$code[which(plot.props$diff==max(plot.props$diff))])
#top.5<-c(as.character(biggest.diff), as.character(best.props))
#top.5<-c(as.character(biggest.diff))
top.5<-best.con.met

plot.props$rate<-round(1/plot.props$rate, 3);plot.props$rate[which(plot.props$rate=="Inf")]<-0
plot.props$col<-1
plot.props$col[which(plot.props$rate==0)]<-0
##plot out the proportion of populations exhibiting EWSs
pp1<-ggplot(plot.props[which(plot.props$code%in%top.5),], aes(x=factor(rate), y=proportion))+geom_bar(stat="identity", aes(fill=factor(col)), col="black")+facet_wrap(~code, ncol=length(top.5))
pp1<-pp1+theme_bw()#+ggtitle(bquote(list(atop("Proportion of time series showing at least one EWSR at a 2" ~sigma ~ "threshold"))))
pp1<-pp1+xlab("Rate of decline in fishing pressure")+ylab("Proportion showing EWSR")+theme(strip.background = element_rect(fill=c("mintcream")))
#pp1<-pp1+geom_hline(aes(yintercept=min(proportion)))

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/proportion showing consecutive signals.pdf", width=4*length(unique(plot.props[which(plot.props$code%in%top.5),]$code)), height=3)
  pp1+scale_fill_manual(values=alpha(c("darkorange", "deepskyblue"), 1))+theme(legend.position="none")+ylim(0, 1)
dev.off()

#####do follow up analyses with the 2 consecutive signals:

##violin plot of consecutive sigs ews 
vio.con<-subset(mean.ews_t, metric.code==plot.best.con$metric[1] & rate!=0 & earliest.sig!="Inf")
vio.con$years.prior.recovery<-vio.con$recovery.point-vio.con$earliest.sig
vio.con$rate<-round(1/vio.con$rate, 3);vio.con$rate[which(vio.con$rate=="Inf")]<-0

cols = rev(colorRampPalette(brewer.pal(11,"Spectral"))(6)[2:6])
p5.vio<-ggplot(vio.con, aes(x=factor(rate), y=years.prior.recovery))
p5.vio<-p5.vio+geom_dotplot(aes(col=factor(rate)),fill="white", binaxis='y', stackdir='center', binwidth = 0.4)
p5.vio<-p5.vio+scale_fill_manual(values=rev(cols), name = "Rate of decline in fishing pressure")
p5.vio<-p5.vio+theme_bw()+theme(legend.position = "none")+annotate("text", x=0.75, y=43, label="c.")+scale_colour_manual(values = (viridis(7, option = "B")[2:7]), name="Rate of decline in fishing pressure")
p5.vio<-p5.vio+ylab("Years prior to recovery")+xlab("Rate of decline in fishing pressure")+facet_wrap(~metric.code)+theme(strip.background = element_rect(fill=c("mintcream")))

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/Fig. 4.pdf", width=12, height=3)
  grid.arrange(p5, p5.prop, p5.vio, ncol=3)
dev.off()

################################################################################################################
################################################################################################################
###ROC
################################################################################################################

##thresholds for ROC curves
thresholds<-c(seq(0.01, 6, 0.2), 2)

##the data to look across
fin.res
##split it by metric
auc.split<-split(fin.res, list(fin.res$metric.code))

auc.curves<-rbindlist(mclapply(auc.split, function(m){
  #m<-auc.split[[1]]
  ##subset the data by collapse and non collapse
  coll<-subset(m, rate!=0)
  non.coll<-subset(m, rate==0)

  ##subset to only look at data between the ceccasation of fishing and the recovery 
  coll<-subset(coll, time>2040 & time < recovery.point)
  non.coll<-subset(non.coll, time>2040 & time < 2090)
  
  results<-NULL
  for(v in 1:length(thresholds)){
    #v=10
    ##split the recovery data by iteration and rate
    split.col<-split(coll, list(coll$iter, coll$rate))
    ##calculate whether each time series which recovers shows an EWS at a given threshold
    res.coll<-unlist(lapply(split.col, function(z){
      #z<-split.col[[1]]
      if(length(which(z$str>thresholds[v]))>0){return(1)}else{return(0)}
    }))
    ##cauluate the proportion that do:
    TP.rat<-sum(res.coll)/length(res.coll)

    ##split the non-recovery data by iteration and rate
    split.non.col<-split(non.coll, list(non.coll$iter, non.coll$rate))
    ##calculate whether each time series WHICH DOESNT RECOVER shows an EWS at a given threshold
    res.non.coll<-unlist(lapply(split.non.col, function(zz){
      #zz<-split.non.col[[1]]
      if(length(which(zz$str>thresholds[v]))>0){return(0)}else{return(1)}
    }))
    ##cauluate the proportion that do:
    FP.rat<-1-(sum(res.non.coll)/length(res.non.coll))
    #############################################################################################
    ###now lets do it for 2 consecutive signals too:
    #############################################################################################
    ##calculate whether each time series which recovers shows an EWS at a given threshold
    res.coll.con<-unlist(lapply(split.col, function(z){
      #z<-split.col[[1]]
      if(min(diff(which(z$str>thresholds[v])))==1){return(1)}else{return(0)}
    }))
    ##cauluate the rroportion that do:
    TP.rat.con<-sum(res.coll.con)/length(res.coll.con)

    ##calculate whether each time series WHICH DOESNT RECOVER shows an EWS at a given threshold
    res.non.coll.con<-unlist(lapply(split.non.col, function(zz){
      #zz<-split.non.col[[1]]
      if(min(diff(which(zz$str>thresholds[v])))==1){return(0)}else{return(1)}
    }))
    ##cauluate the proportion that do:
    FP.rat.con<-1-(sum(res.non.coll.con)/length(res.non.coll.con))

    results[[v]]<-data.frame(metric=coll$metric.code[1],
               threshold=thresholds[v],
               TP.rat=TP.rat,
               FP.rat=FP.rat,
               con.TP.rat=TP.rat.con,
               con.FP.rat=FP.rat.con)
  }

  fin.results<-rbindlist(results)

  return(fin.results)
  }, mc.cores=7))

##add in inc.size vector for plotting
auc.curves<-rbindlist(lapply(split(auc.curves, list(auc.curves$metric)), function(HH){
  #HH<-split(quantMetric, list(quantMetric$code))[[10]]
  if(length(grep("size",HH$metric[1], value = TRUE))>0){
    HH$inc.size<-T
  }else{HH$inc.size<-F}
  return(HH)
}))

auc.curves

auc1<-ggplot(auc.curves, aes(x=FP.rat, y=TP.rat, col=metric))+geom_line(aes(linetype=inc.size))+geom_abline(slope=1, intercept=0)+theme_bw()+ggtitle("a. Single signal")+ylab("TPR")+xlab("FPR")+labs(col="Metric", linetype="Includes size")
auc2<-ggplot(auc.curves, aes(x=con.FP.rat, y=con.TP.rat, col=metric))+geom_line(aes(linetype=inc.size))+geom_abline(slope=1, intercept=0)+theme_bw()+ggtitle("b. Two consecutive signals")+ylab("TPR")+xlab("FPR")

sig2<-subset(auc.curves, threshold ==2)

auc1<-auc1+geom_point(data=sig2, aes(x=FP.rat, y=TP.rat, col=metric))
auc2<-auc2+geom_point(data=sig2, aes(x=con.FP.rat, y=con.TP.rat, col=metric))

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)
  }

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/AUC.pdf", width=11, height=6)
   grid_arrange_shared_legend(auc1, auc2, ncol=2, nrow = 1)
dev.off()


