## charles Ferté
## Sage Bionetworks
## March 6th, 2012

library(survival)
options(stringsAsFactors=FALSE)

# load the DirClin
load("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/TS_CLIN.Rdata")
dat <- DirClinF
dats <- dat[,c("MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
dat <- dat[,c("P_Stage","Age","GENDER","RACE","SMOKING","SITE","H_Grade")]
dat$P_Stage <- as.factor(dat$P_Stage)
dat$GENDER <- as.factor(dat$GENDER)
dat$RACE <- as.factor(dat$RACE)
dat$SMOKING <- as.factor(dat$SMOKING)
dat$SITE <- as.factor(dat$SITE)
dat$H_Grade <- as.factor(dat$H_Grade)


## make a training set dat1 and a validation set dat2
dat$tmp <- 1
tmp <- c(1:299)
set.seed(1225)
CI_TS <- c()
CI_VS <- c()

f_K_fold <- function(Nobs,K=10){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

for(j in 1:1){

  DAT <- f_K_fold(299)
  for(s in c(1:10)){
  dat1 <- dat[DAT[[s]]$train,]
dat2 <- dat[DAT[[s]]$test,]
dats1 <- dats[DAT[[s]]$train,]
dats2 <- dats[DAT[[s]]$test,]


GRID <- expand.grid(a1=c(1,0),a2=c(1,0),a3=c(1,0),a4=c(1,0),a5=c(1,0),a6=c(1,0),a7=c(1,0))
GRID <- GRID[-128,]
GRID$tmp <- 1

a <- c()
b <- c()
TS <- c()
VS <- c()
for(i in 1:127){
  adat1 <- data.frame(dat1[,which(GRID[i,]==1)])
  for(k in 1: dim(adat1)[2]){
  adat1 <- adat1[!is.na(adat1[,k]),]
  }
  adat2 <- data.frame(dat2[,which(GRID[i,]==1)])
  for(k in 1: dim(adat2)[2]){
    adat2 <- adat2[!is.na(adat2[,k]),]
  }
 adat1$tmp <- NULL
adat2$tmp <- NULL
  
  blah1 <- dats1[rownames(adat1),]
  blah1 <- Surv(blah1$MONTHS_TO_LAST_CONTACT_OR_DEATH,blah1$VITAL_STATUS)
  
  blah2 <- dats2[rownames(adat2),]
  blah2 <- Surv(blah2$MONTHS_TO_LAST_CONTACT_OR_DEATH,blah2$VITAL_STATUS)
  
  aFit <- coxph(blah1~.,data=adat1)
  a <- survConcordance(blah1~predict(aFit),adat1)

  #coef <- matrix(aFit$coefficients)
  #mm <- model.matrix(~ . - 1, data=adat2)
  #predNew <- as.numeric(mm %*% coef)
  #predict(aFit, newdata=mm)
  b <- survConcordance(blah2~predict(aFit,newdata=adat2,reference="sample"),adat2)
  
  
  TS <- c(TS,a$concordance[names(a$concordance)=="concordant"])
VS <- c(VS,b$concordance[names(b$concordance)=="concordant"])
  }
CI_TS <- cbind(CI_TS,TS) 
CI_VS <- cbind(CI_VS,VS)

}
  }

rownames(CI_VS) <- c(1:127)
CI_VS <- t(CI_VS)

medianVS <- apply(CI_VS,2,median)

CI_VS <- CI_VS[,order(medianVS)]

boxplot.matrix(CI_VS,use.cols=TRUE,col="aquamarine4",main="Performance of OS prediction in the Dir dataset 
(clinical covariates only)", ylab="concordance index",xlab="127 different combinations of clinical covariates (no molecular features)
cox model, cross-validation (nfolds=10)
    NB:  the dashed blue line represent today's gold standard= P Stage only",xaxt='n',ylim=c(.4,1))

## looking for the BEST model
GRID$tmp <- NULL
colnames(GRID) <- c("P_Stage","Age","GENDER","RACE","SMOKING","SITE","H_Grade")
COL <- colnames(CI_VS)[127]
colnames(GRID)[which(GRID[COL,]==1)]
BEST <- CI_VS[,127]

# looking for the CI of the gold standard= using p_Stage only as predictor
p <- which(GRID$P_Stage==1 & GRID$Age==0 & GRID$GENDER==0 & GRID$RACE==0 & GRID$SMOKING==0 & GRID$SITE==0 & GRID$H_Grade)
medianVS[names(medianVS)==p]
P_STAGE_ONLY <-CI_VS[,p] 
abline(v=p,col="royalblue",lwd=1.5,lty=5)

## compare the best model found with the gold standard= P_STAGE_ONLY
TMP <- cbind(P_STAGE_ONLY,BEST)
boxplot(TMP,xlab= paste("NB: the BEST model takes in account:  ",paste(colnames(GRID)[which(GRID[COL,]==1)],collapse=", "),sep=""), ylab="concordance index",
        col=c("royalblue","aquamarine4"),xaxt='n')
TEST <- t.test(BEST,P_STAGE_ONLY)
title("Comparison of model performance: 
  BEST model vs. P Stage only")

print("variables included in the model for maximizing the concordance index:")
colnames(GRID)[which(GRID[COL,]==1)]
summary(CI_VS[,127])

## the model that includes only p stage and smoking is close to the best one in terms of performance 
GRID[colnames(CI_VS)[123],]
summary(CI_VS[,123])
