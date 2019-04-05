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
library(RColorBrewer)
##############################################################################################################################
##load in the ews output:
load("~/Desktop/Work/Recovery EWS/R code/sizeBasedStuff/NorthSea/EWS output/fin.res.data.length.RData")
fin.res.data.length<-as.data.table(fin.res.data.length)
##############################################################################################################################
##strength of signal analysis?
fin.res.data.length
levels(fin.res.data.length$metric.code)<-c("AR(1)", "CV", "Mean size", "SD size", "AR(1) + CV", "AR(1) + Mean size", "AR(1) + SD size", "CV + mean size", "CV + SD size", "Mean size + SD size", "AR(1) + CV + mean size", "AR(1) + CV + SD size", "AR(1) + mean size + SD size", "CV + mean size + SD size", "AR(1) + CV + mean size + SD size")

##############################################################################################################################
##only look at the following metric:
#fin.res.data.length<-subset(fin.res.data.length, metric.code=="AR1 + SD size")

##make a key vector:
fin.res.data.length[,key_var:=paste(metric.code, iter, rate, data.length)]

##vector for splitting the data.table up in mclapply
s <- split(seq(nrow(fin.res.data.length)), fin.res.data.length$key_var)

##
result <- rbindlist(lapply(s, function(.indx){ 
  ##.indx<-s[[1]]
  subdat <- fin.res.data.length[.indx, ] 
  thresh<-2
  ##for single signals
  if(subdat$rate[1]!=0){
    subdat<-subset(subdat, time>=subdat$hold.year[1] & time<recovery.point)
    TP<-if(length(which(subdat$str>thresh))>0){1}else{0}
  }else{
    subdat<-subset(subdat, time>=subdat$hold.year[1] & time<2090)
    FP<-if(length(which(subdat$str>thresh))>0){1}else{0}
  }

  ##for consecutive signals
  if(subdat$rate[1]!=0){
    subdat<-subset(subdat, time>=subdat$hold.year[1] & time<recovery.point)
    if(length(which(subdat$str>thresh))>0){
      TP.con<-if(min(diff(which(subdat$str>thresh)))==1){1}else{0}
    }else{TP.con<-0}
  }else{
    subdat<-subset(subdat, time>=subdat$hold.year[1] & time<2090)
    if(length(which(subdat$str>thresh))>0){
      FP.con<-if(min(diff(which(subdat$str>thresh)))==1){1}else{0}
    }else{FP.con<-0}
  }
  ##return results:  
  return(data.frame(metric.code=subdat$metric.code[1],
             rate=subdat$rate[1],
             data.length=subdat$data.length[1],
             iter=subdat$iter[1],
             #total.data.length=subdat$data.length[1]+(min(subdat$time[which(subdat$str>thresh)])-subdat$hold.year[1]),
             TP=if(subdat$rate[1]!=0){TP}else{NA},
             FP=if(subdat$rate[1]==0){FP}else{NA},
             TP.con=if(subdat$rate[1]!=0){TP.con}else{NA},
             FP.con=if(subdat$rate[1]==0){FP.con}else{NA}))
}))

split.res<-split(result, list(result$data.length, result$metric.code))

dd<-rbindlist(lapply(split.res, function(m){
  #m<-split.res[[1]]
  non.rec<-subset(m, rate==0)
  rec<-subset(m, rate!=0)
  FP.ratio<-sum(non.rec$FP)/length(non.rec$FP)
  FP.con.ratio<-sum(non.rec$FP.con)/length(non.rec$FP.con)
  TP.ratio<-sum(rec$TP)/length(rec$TP)
  TP.con.ratio<-sum(rec$TP.con)/length(rec$TP.con)
  
  return(data.frame(metric.code=m$metric.code[1],
             data.length=m$data.length[1],
             FP.ratio=FP.ratio,
             FP.con.ratio=FP.con.ratio,
             TP.ratio=TP.ratio,
             TP.con.ratio=TP.con.ratio))
}))
ddd<-melt(dd, id=c(1,2))

levels(ddd$variable)<-c("False positive ratio\n(single signal)", 
                        "False positive ratio\n(consecutive signals)", 
                        "True positive ratio\n(single signal)",
                        "True positive ratio\n(consecutive signals)")

ddd$col<-"Single signal"
ddd$col[which(ddd$variable=="False positive ratio\n(consecutive signals)")]<-"Consecutive signals"
ddd$col[which(ddd$variable=="True positive ratio\n(consecutive signals)")]<-"Consecutive signals"

ddd$sig<-"True positive"
ddd$sig[which(ddd$variable=="False positive ratio\n(single signal)")]<-"False positive"
ddd$sig[which(ddd$variable=="False positive ratio\n(consecutive signals)")]<-"False positive"

single.metric<-subset(ddd, metric.code=="AR(1) + SD size")
pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/data length.pdf", width=7, height=4)
  p1<-ggplot(single.metric, aes(x=data.length, y=value, group=interaction(sig, col), linetype=sig, col=factor(col)))+geom_line()+theme_bw()
  p1<-p1+ylab("Proportion")+xlab("Training data available prior to 2040 (years)")+facet_wrap(~metric.code)
  p1<-p1+scale_colour_manual(values=alpha(c("darkorange", "deepskyblue"), 1))
  p1<-p1+theme(strip.background = element_rect(fill=c("mintcream")))+labs(color='Signal type', linetype="")
  p1
dev.off()

##then for all metrics:

pdf("~/Desktop/Work/Recovery EWS/R code/Various fishing pressures/Plots/data length ALL.pdf", width=9, height=12)
p1<-ggplot(ddd, aes(x=data.length, y=value, group=interaction(sig, col), linetype=sig, col=factor(col)))+geom_line()+theme_bw()
p1<-p1+ylab("Proportion")+xlab("Training data available prior to 2040 (years)")+facet_wrap(~metric.code)
p1<-p1+scale_colour_manual(values=alpha(c("darkorange", "deepskyblue"), 1))
p1<-p1+theme(strip.background = element_rect(fill=c("mintcream")))+labs(color='Signal type', linetype="")
p1+facet_wrap(~metric.code, ncol=3)
dev.off()