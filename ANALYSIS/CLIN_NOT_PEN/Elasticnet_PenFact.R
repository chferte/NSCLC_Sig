##### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### March, 2 2012

### script for describing the NSCLC datasets
options(stringsAsFactors=FALSE)

#load the libraries
library(metaGEO)
library(snm)
library(lattice)
library(affy)
library(survival)
library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)
library(survival)
library(risksetROC)
library(caret)
library(survcomp)

## run Elasticnet to predict survival
## NOT penalizing the clinical covariates

## load the DirClinF & ZhuClinF clinicalcovariates file
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/")
load("VS_CLIN.Rdata")
load("TS_CLIN.Rdata")

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/Zhu_CEL/")
load("metaGEO_VS_rmSITE.Rdata")
load("metaGEO_VS.Rdata")
MATRIX_VS <- dat
MATRIX_VS1 <- new.dat

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/Dir_CEL")
load("metaGEO_TS_rmSITE.Rdata")
load("metaGEO_TS.Rdata")
MATRIX_TS <- dat
MATRIX_TS1 <- new.dat

rm(dat,new.dat)

## integrating the clinical covariates without NAs
## transform the clinical covariates into matrix
tmp <- DirClinF[,c("Histology","P_Stage","GENDER", "Age")]
table(is.na(DirClinF$Histology))
table(is.na(DirClinF$P_Stage))
table(is.na(DirClinF$GENDER))
table(is.na(DirClinF$Age))

tmp1 <- ZhuClinF[,c("Histology","P_Stage","GENDER", "Age")]
table(is.na(ZhuClinF$Histology))
table(is.na(ZhuClinF$P_Stage))
table(is.na(ZhuClinF$GENDER))
table(is.na(ZhuClinF$Age))

tmp$Histology <- as.numeric(as.factor(tmp$Histology))
tmp$P_Stage <- as.numeric(as.factor(tmp$P_Stage))
tmp$GENDER <- as.numeric(as.factor(tmp$GENDER))

tmp1$Histology <- as.numeric(as.factor(tmp1$Histology))
tmp1$P_Stage <- as.numeric(as.factor(tmp1$P_Stage))
tmp1$GENDER <- as.numeric(as.factor(tmp1$GENDER))

tmp <- t(tmp)
tmp1 <- t(tmp1)

MATRIX_TS <- rbind(MATRIX_TS,tmp)
MATRIX_VS <- rbind(MATRIX_VS,tmp1)
MATRIX_TS1 <- rbind(MATRIX_TS1,tmp)
MATRIX_VS1 <- rbind(MATRIX_VS1,tmp1)

##############################################################################################
#####      using sample to to cross validation on dat
##############################################################################################

C <- DirClinF
C1 <- DirClinF

## set the variables
x <- t(MATRIX_TS)
x1 <- t(MATRIX_VS)
y <-  Surv(C$MONTHS_TO_LAST_CONTACT_OR_DEATH,C$VITAL_STATUS)
y1 <-  Surv(C1$MONTHS_TO_LAST_CONTACT_OR_DEATH,C1$VITAL_STATUS)

## build the Grid alpha/lambda
alphas <- exp(-1*seq(0,10,1))
lambdas <- exp(seq(-4,3,1))
GRID <- expand.grid(.family="cox",.alpha=alphas,.lambda=lambdas)

CI_TS <- c()
CI_VS <- c()
for(i in 1:dim(GRID)[1]){
  fit  <- try(glmnet(x,y,family="cox",alpha=GRID$.alpha[i],lambda=GRID$.lambda[i],penalty.factor=c(rep(1,12161),rep(0,4))))
  
  if( class(fit) == "try-error" ){
    CI_TS <- c(CI_TS, NA)
    CI_VS <- c(CI_VS, NA)
    
  } else{
    y_E_TS <- predict(fit,x,type="link")
    CI_TS <- c(CI_TS,survConcordance(Surv(y[,1],y[,2])~y_E_TS)$concordance)
    
    y_E_VS <- predict(fit,x1,type="link")
    CI_VS <- c(CI_VS,survConcordance(Surv(y1[,1],y1[,2])~y_E_VS)$concordance)
  }
}

##

CI_TOTAL <- cbind(GRID$.alpha,GRID$.lambda,CI_TS,CI_VS)
rownames(CI_TOTAL)<- c(1:88)
colnames(CI_TOTAL)[1]<-"alpha"
colnames(CI_TOTAL)[2]<-"lambda"
CI_TOTAL_dat<-as.data.frame(CI_TOTAL)

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/ANALYSIS/CLIN_NOT_PEN/")
save(CI_TOTAL_dat,file="CI_dat.Rdata")


##############################################################################################
#####      using sample to to cross validation on new.dat
##############################################################################################


C <- DirClinF
C1 <- DirClinF

## set the variables
x <- t(MATRIX_TS1)
x1 <- t(MATRIX_VS1)
y <-  Surv(C$MONTHS_TO_LAST_CONTACT_OR_DEATH,C$VITAL_STATUS)
y1 <-  Surv(C1$MONTHS_TO_LAST_CONTACT_OR_DEATH,C1$VITAL_STATUS)

## build the Grid alpha/lambda
alphas <- exp(-1*seq(0,10,1))
lambdas <- exp(seq(-4,3,1))
GRID <- expand.grid(.family="cox",.alpha=alphas,.lambda=lambdas)

CI_TS <- c()
CI_VS <- c()
for(i in 1:dim(GRID)[1]){
  fit  <- try(glmnet(x,y,family="cox",alpha=GRID$.alpha[i],lambda=GRID$.lambda[i],penalty.factor=c(rep(1,12161),rep(0,4))))
  
  if( class(fit) == "try-error" ){
    CI_TS <- c(CI_TS, NA)
    CI_VS <- c(CI_VS, NA)
    
  } else{
    y_E_TS <- predict(fit,x,type="link")
    CI_TS <- c(CI_TS,survConcordance(Surv(y[,1],y[,2])~y_E_TS)$concordance)
    
    y_E_VS <- predict(fit,x1,type="link")
    CI_VS <- c(CI_VS,survConcordance(Surv(y1[,1],y1[,2])~y_E_VS)$concordance)
  }
}

##

CI_TOTAL <- cbind(GRID$.alpha,GRID$.lambda,CI_TS,CI_VS)
rownames(CI_TOTAL)<- c(1:88)
colnames(CI_TOTAL)[1]<-"alpha"
colnames(CI_TOTAL)[2]<-"lambda"
CI_TOTAL_new_dat<-as.data.frame(CI_TOTAL)

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/ANALYSIS/CLIN_NOT_PEN/")
save(CI_TOTAL_new_dat,file="CI_new_dat.Rdata")


#################################################################################################################################################################################
############################################################ Visualizing the results of the concordance index   ###########################################################
#################################################################################################################################################################################

#################################################################################################################################################################################
############################################################ Visualizing the results of the concordance index   ###########################################################
#################################################################################################################################################################################
SNM <- CI_TOTAL_new_dat
SNM <- SNM[!is.na(SNM$CI_TS),]
SNM <- SNM[SNM$CI_TS>.51,]
SNM <- SNM[order(SNM$CI_TS),]

US <- CI_TOTAL_dat
US <- US[!is.na(US$CI_TS),]
US <- US[US$CI_TS>.51,]
US <- US[order(US$CI_TS),]

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/ANALYSIS/CLIN_NOT_PEN/")

library(lattice)
pdf("Performance_Supervised.pdf")
xyplot(SNM$CI_TS +SNM$CI_VS ~ c(1:33),type="l",auto.key=list(columns=2,space="bottom",points=FALSE,lines=TRUE),ylim=c(.3,1), 
       main="performance of the prediction using
       SUPERVISED normalization method (SNM)",xlab="combinations of (alpha,lambda)",ylab="concordance index")
dev.off()

pdf(("Performance_Unsupervised.pdf"))
xyplot(US$CI_TS +US$CI_VS ~ c(1:31),type="l",auto.key=list(columns=2),ylim=c(.3,1), 
       main="performance of the prediction using
       UNSUPERVISED normalization method (US)",xlab="combinations of (alpha,lambda)",ylab="concordance index")
dev.off()

pdf("Performance_VS.pdf")
xyplot(SNM$CI_VS +US$CI_VS ~ c(1:33),type="l",auto.key=list(columns=2),ylim=c(.3,1), 
       main="performance of the prediction of the VALIDATION SET
       comparing SNM vs US methods ",xlab="combinations of (alpha,lambda)",ylab="concordance index")
dev.off()

pdf("Performance_TS.pdf")
xyplot(SNM$CI_TS +US$CI_TS ~ c(1:33),type="l",auto.key=list(columns=2),ylim=c(.3,1), 
       main="performance of the prediction of the TRAINING SET
       comparing SNM vs US methods",xlab="combinations of (alpha,lambda)",ylab="concordance index")
dev.off()

