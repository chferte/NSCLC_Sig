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


## load the dataset normalized by snm where site is not removed 
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/Hou_CEL/")
fits <- runWorkflow(".")
dat <- exprs(fits$hgu133plus2)
save(dat, file="metaGEO_VS2.Rdata")

## load the DirClinF clinicalcovariates file
load("VS2_CLIN.Rdata")

## load the dataset normalized and where SITE is removed:
DirClinF$H_Grade[is.na(DirClinF$H_Grade)] <- 2
ids <- which(!is.na(DirClinF$H_Grade))
bio.var <- model.matrix(~ DirClinF$H_Grade)
adj.var <- model.matrix(~ DirClinF$SITE)
snm.fit <- snm(dat, 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)

new.dat <- snm.fit$norm.dat

# make the colnames of newdat=colnames(dat)
colnames(new.dat) <- colnames(dat)

## check the dat and new.dat colnames  match with DirClinF rownames
identical(colnames(dat),rownames(DirClinF))

## save the products into dat and new.dat
save(dat, file="metaGEO_TS.Rdata")
save(new.dat,file="metaGEO_TS_rmSITE.Rdata")




## run Elasticnet to predict survival





##############################################################################################
#####      using sample to to cross validation on new.dat
##############################################################################################
tmp <- c(1:299)
tmp1 <- sample(tmp,199)
tmp2 <- setdiff(tmp,tmp1)
x <- new.dat[,tmp1]
x1 <- new.dat[,tmp2]
C <- DirClinF[tmp1,]
C1 <- DirClinF[tmp2,]

## set the variables
x <- t(x)
x1 <- t(x1)
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


save(CI_TOTAL_dat,file="CI_new_dat.Rdata")


##############################################################################################
#####      using sample to to cross validation on dat
##############################################################################################

x <- dat[,tmp1]
x1 <- dat[,tmp2]
C <- DirClinF[tmp1,]
C1 <- DirClinF[tmp2,]

## set the variables
x <- t(x)
x1 <- t(x1)
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


save(CI_TOTAL_dat,file="CI_dat.Rdata")

#################################################################################################################################################################################
############################################################ Visualizing the results of the concordance index   ###########################################################
#################################################################################################################################################################################

## load the dat and new dat objects
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/Dir_CEL/")
load("metaGEO_TS.Rdata")
load("metaGEO_TS_rmSITE.Rdata")
## load the DirClinF clinicalcovariates file
load("TS_CLIN.Rdata")

# split the TS and the VS in 2 sets 
tmp <- c(1:299)
tmp1 <- sample(tmp,199)
tmp2 <- setdiff(tmp,tmp1)
x <- dat[,tmp1]
x1 <- dat[,tmp2]
C <- DirClinF[tmp1,]
C1 <- DirClinF[tmp2,]

## set the variables
x <- t(x)
x1 <- t(x1)
y <-  Surv(C[,12],C[,11])
y1 <-  Surv(C1[,12],C1[,11])

## run glmnet cv.glmnet
cv.fit <- cv.glmnet(x,y, family = "cox", maxit = 1000)















#################################################################################################################################################################################
############################################################ Visualizing the results of the concordance index   ###########################################################
#################################################################################################################################################################################
TND <- CI_TOTAL_new.data
TND <- TND[!is.na(TND$CI_TS),]
TND <- TND[order(TND$CI_TS),]

TD <- CI_TOTAL_dat
TD <- TD[!is.na(TD$CI_TS),]
TD <- TD[order(TD$CI_TS),]

library(lattice)
xyplot(TND$CI_TS +TND$CI_VS ~ c(1:53),type="l",auto.key=list(columns=2),ylim=c(.3,1), 
       main="snm normalized data removing the SITE",xlab="combinations of (alpha,lambda)",ylab="concordance index")


xyplot(TD$CI_TS +TD$CI_VS ~ c(1:53),type="l",auto.key=list(columns=2),ylim=c(.3,1), 
       main="snm normalized data NOT removing the SITE",xlab="combinations of (alpha,lambda)",ylab="concordance index")


xyplot(TND$CI_VS +TD$CI_VS ~ c(1:53),type="l",auto.key=list(columns=2),ylim=c(.3,1), 
       main="snm normalized data removing/NOT removing the SITE (in the validation set only) ",xlab="combinations of (alpha,lambda)",ylab="concordance index")


xyplot(TND$CI_TS +TD$CI_TS ~ c(1:53),type="l",auto.key=list(columns=2),ylim=c(.3,1), 
       main="snm normalized data removing/NOT removing the SITE (in the training set only) ",xlab="combinations of (alpha,lambda)",ylab="concordance index")
