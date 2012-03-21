##### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### March, 2 2012

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

## load the DirClinF & ZhuClinF clinicalcovariates file
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/")
load("VS_CLIN.Rdata")
load("TS_CLIN.Rdata")

## load the dataset normalized by snm where site is not removed 
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/Zhu_CEL/")
fits <- runWorkflow(".")
dat <- exprs(fits$hgu133a)
save(dat, file="metaGEO_VS.Rdata")

## load the dataset normalized and where SITE is removed:
bio.var <- model.matrix(~ ZhuClinF$Histology)
adj.var <- model.matrix(~ ZhuClinF$SCANBATCH)
snm.fit <- snm(dat, 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)

new.dat <- snm.fit$norm.dat


# make the colnames of newdat=colnames(dat)
colnames(new.dat) <- colnames(dat)
save(new.dat,file="metaGEO_VS_rmSITE.Rdata")

## check the dat and new.dat colnames  match with DirClinF rownames
identical(colnames(dat),rownames(ZhuClinF))

## run Elasticnet to predict survival
##load the data

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


##############################################################################################
#####      using sample to to cross validation on new.dat
##############################################################################################

C <- DirClinF
C1 <- ZhuClinF

## set the variables
x <- t(MATRIX_TS1)
x1 <- t(MATRIX_VS1)
y <-  Surv(C[,12],C[,11])
y1 <-  Surv(C1[,12],C1[,11])

## build the Grid alpha/lambda
alphas <- exp(-1*seq(0,10,1))
lambdas <- exp(seq(-4,3,1))
GRID <- expand.grid(.family="cox",.alpha=alphas,.lambda=lambdas)

CI_TS <- c()
CI_VS <- c()
for(i in 1:dim(GRID)[1]){
  fit  <- try(glmnet(x,y,family="cox",alpha=GRID$.alpha[i],lambda=GRID$.lambda[i]))
  
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

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/")
save(CI_TOTAL_new_dat,file="CI_new_dat.Rdata")


##############################################################################################
#####      using sample to to cross validation on dat
##############################################################################################

x <- t(MATRIX_TS)
x1 <- t(MATRIX_VS)
C <- DirClinF[tmp1,]
C1 <- DirClinF[tmp2,]

## set the variables

y <-  Surv(C[,12],C[,11])
y1 <-  Surv(C1[,12],C1[,11])

## build the Grid alpha/lambda
alphas <- exp(-1*seq(0,10,1))
lambdas <- exp(seq(-4,3,1))
GRID <- expand.grid(.family="cox",.alpha=alphas,.lambda=lambdas)

CI_TS <- c()
CI_VS <- c()
for(i in 1:dim(GRID)[1]){
  fit  <- try(glmnet(x,y,family="cox",alpha=GRID$.alpha[i],lambda=GRID$.lambda[i]))
  
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

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/")
save(CI_TOTAL_dat,file="CI_dat.Rdata")

#################################################################################################################################################################################
############################################################ Visualizing the results of the concordance index   ###########################################################
#################################################################################################################################################################################

## load the dat and new dat objects


setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/")
load("CI_dat.Rdata")
load("CI_new_dat.Rdata")

identical(CI_TOTAL_dat,CI_TOTAL_new_dat)


## load the DirClinF clinicalcovariates file
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/")
load("VS_CLIN.Rdata")
load("TS_CLIN.Rdata")



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

setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/SNM/")

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

