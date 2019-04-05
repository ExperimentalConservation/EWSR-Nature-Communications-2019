###########################################################################################
library(mgcv)
eps <- 1e-7 ## finite difference interval

##functino to calculate the slopes, derivatives and ci's of the data, based on the Gam method from Burthe et al.:

is.there = function(x=0, L.CI, U.CI){

		pos<- ifelse(x<U.CI, 1, -1) 

		negs<-ifelse(x>L.CI, -1, 1)

		return(pos+negs)}

													##round(length(timeseries)/4)	

gam_smoothing<-function(years, timeseries,knots){

		if(length(which(timeseries<=0))==0){

		gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian(link="log"))}else{

		gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian)}

		time.series.fit<-predict(gam1, newdata=data.frame(years=years), type="response")

		

		X0<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")

		X1<-predict(gam1, newdata=data.frame(years=years+eps), type= "lpmatrix")

		Xi<-(X1-X0)/eps

		df <- Xi%*%coef(gam1)              ## ith smooth derivative 

		df.sd <- rowSums(Xi%*%gam1$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5

		#plot(years,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)))##plot 'em

		#lines(years,df+2*df.sd,lty=2);lines(years,df-2*df.sd,lty=2)

		splines<-data.frame(years=years,deriv=df,U.CI=df+2*df.sd, L.CI=df-2*df.sd)	

		splines$sign<-is.there(0, splines$L.CI, splines$U.CI)/2

		splines$fit<-time.series.fit

		return(splines)}

###########################################################################################
