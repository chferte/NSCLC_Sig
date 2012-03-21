## Charles Ferté
## Sage Bionetworks
## February 14th 2012

# load the libraries
library(survcomp)
library(xtable)

# load Charles data
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/")


## evaluate the performance of my survival predictor 
x1  <-"my predictor object for dataset1 (cox model)"
x2 <- "my predictor object for dataset2 (cox model)" # normally lenght(datset 1)==lenght(dataset2) for comparing the index

y1 <-  "a vector of event times (i.e: the time of the survival we want to predict) in dataset1"
z1 <- "a vector of event occurence indicators (i.e: the event status of the survival we want to predict) in dataset1"

y2 <-  "a vector of event times (i.e: the time of the survival we want to predict) in dataset2"
z2 <- "a vector of event occurence indicators (i.e: the event status of the survival we want to predict) in dataset2"


##  evaluating the performance of the predcitor using concordance index
## to dysplay the standard error and the lower/upper bound of the confidence intervall (here, alpha =.05)

library(survcomp)
cindex1 <- concordance.index(x=x1, surv.time=y1,  surv.event=z1, method="noether", na.rm=TRUE, alpha= .05)
cindex2 <- concordance.index(x=x2, surv.time=y2,  surv.event=z2, method="noether", na.rm=TRUE,alpha= .05)


cindex1$c.index # returns the cindex1
cindex1$se #standard error of the estimate
cindex1$lower # lower CI bound of cindex1
cindex1$upper # upper CI bound of cindex1

cindex2$c.index # returns the cindex2
cindex2$se #standard error of the estimate
cindex2$lower # lower CI bound of cindex2
cindex2$upper # upper CI bound of cindex2

## to compare two concordance index (Student t test for dependant samples)
cindex.comp(cindex1, cindex2)

## to display the different concordance index in a forest plot setting
forestplot(labeltext=rbind(c("predictor 1","   "),c("predictor 2","   ")),mean=c(cindex1$c.index,cindex2$c.index),lower=c(cindex1$lower,cindex2$lower,.4),upper=c(cindex1$upper,cindex2$upper),is.summary=FALSE,xlab="Concordance Index",zero=0.5,graphwidth = unit(2,"inches"),xlog=F,col=meta.colors(box="royalblue", line="darkblue", zero="darkred"),boxsize=.5,x.ticks=seq(0,1,.1))

## to compute and dysplay the Time dependent AUCs (AUC(t))
library(risksetROC)
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/Perf_Eval_SurvMod/")
pdf("AUC_ROC.pdf")
AUC_COX_1 <- risksetAUC(Stime=y1,status=z1, marker=x1, method="Cox", tmax=200, col="green", lwd=1,lty=1,main="Time dependent AUC for the predictors x1 and x2 in dataset 1")                       
AUC_COX_2 <- risksetAUC(Stime=y2,status=z2, marker=x2, method="Cox", span=.4, tmax=200, plot=FALSE)
lines(AUC_COX_2$utimes, AUC_COX_2$AUC,col="blue", lty=1, lwd=1)
dev.off()





