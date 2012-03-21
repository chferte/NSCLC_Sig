## charles ferté


### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

############################################################################################################################
### script for running modelling prediction
############################################################################################################################

PATH <- "/Volumes"

setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/",sep=""))

load("TS_CLIN.Rdata")
load("VS_CLIN.Rdata")
load("VS2_CLIN.Rdata")
load("VS3_CLIN.Rdata")

TOTAL <- rbind(DirClinF,GirClinF)
TOTAL <- rbind(TOTAL,HouClinF)
TOTAL <- as.data.frame(rbind(TOTAL,ZhuClinF))

TOTAL$OS <- as.numeric(TOTAL$MONTHS_TO_LAST_CONTACT_OR_DEATH)
TOTAL$OSC <- as.numeric(TOTAL$VITAL_STATUS)
MATRIX <- t(TOTAL)
MATRIX <- MATRIX[-c(1,11,12,13,14,20),]

############################################################################################################################
###### cox ph  analysis
############################################################################################################################
library(randomSurvivalForest)
library(survival)

mySurv <-Surv(TOTAL$OS,TOTAL$OSC)

myCoxFun <- function(x){
  summary(coxph(mySurv ~ as.numeric(x)))$logtest["pvalue"]
}

colnames(TOTAL)
plot(survfit(Surv(TOTAL$OS,TOTAL$OSC) ~ TOTAL$P_Stage), col=1:length(TOTAL$P_Stage),main="OS according to pStage in 464 NSCLC", xlab="months")

summary(coxph(Surv(TOTAL$OS, TOTAL$OSC) ~ TOTAL$P_Stage+TOTAL$SMOKING+TOTAL$Histology+TOTAL$H_Grade+TOTAL$RACE+TOTAL$Age+TOTAL$GENDER+TOTAL$SITE ,method="breslow",data=TOTAL))

myResults <- apply(MATRIX, 1, myCoxFun)
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
