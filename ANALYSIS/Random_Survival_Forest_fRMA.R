### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

### script for running modelling prediction

setwd("~/FELLOW/cferte/MATRIX_RESP_OBJECTS/")
load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")
load("y_TS.Rdata")
load("y_VS.Rdata")
load("y_VS2.Rdata")
load("y_OS_TS.Rdata")
load("y_OS_VS.Rdata")
load("y_OS_VS2.Rdata")





############################################################################################################################
## ### univariate anbalysis
############################################################################################################################
library(randomSurvivalForest)
library(survival)

mySurv <-Surv(y_OS_TS[,1],y_OS_TS[,2])

myCoxFun <- function(x){
  summary(coxph(mySurv ~ as.numeric(x)))$logtest["pvalue"]
}


myResults <- apply(MATRIX_TS, 1, myCoxFun)
hist(myResults)
summary(myResults)
summary(p.adjust(myResults,method="BH"))
blah <- myResults<.1
blah1 <- p.adjust(myResults,method="BH")<.05
MATRIX_TSb <- t(MATRIX_TS[blah,])
MATRIX_VSb <- t(MATRIX_VS[blah,])

all(rownames(y_OS_TS)==rownames(MATRIX_TSb))
all(rownames(y_OS_VS)==rownames(MATRIX_VSb))

MATFINAL <- as.data.frame(cbind(MATRIX_TSb,y_OS_TS))
MATVAL <- as.data.frame(cbind(MATRIX_VSb,y_OS_VS))

colnames(MATFINAL)[!(colnames(MATFINAL) %in% c("time", "status"))] <- paste("v", colnames(MATFINAL)[!(colnames(MATFINAL) %in% c("time", "status"))], sep="")
colnames(MATVAL)[!(colnames(MATVAL) %in% c("time", "status"))] <- paste("v", colnames(MATVAL)[!(colnames(MATVAL) %in% c("time", "status"))], sep="")



############################################################################################################################
### run random forest survival
############################################################################################################################

setwd("/home/cferte/FELLOW/cferte/RSF/")
# fitBase <- rsf(Survrsf(time, status) ~ ., data = MATFINAL)
#save(fitBase,file="fitBase_RSF.Rdata")
load("fitBase_RSF.Rdata")
names(fitBase)
summary(fitBase$err.rate)
pred_VS <-predict(fitBase,MATVAL)
names(pred_VS)
summary(pred_VS$err.rate)
save(pred_VS,file="pred_VS_RSF.Rdata")

