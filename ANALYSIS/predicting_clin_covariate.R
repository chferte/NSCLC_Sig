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
for(j in 1:100){
tmp1 <- sample(299,60)
dat1 <- dat[tmp1,]
dat2 <- dat[setdiff(tmp,tmp1),]
dats1 <- dats[tmp1,]
dats2 <- dats[setdiff(tmp,tmp1),]


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

rownames(CI_VS) <- c(1:127)
CI_VS <- t(CI_VS)
meanVS <- apply(CI_VS,2,median)

CI_VS <- CI_VS[,order(meanVS)]
boxplot.matrix(CI_VS,use.cols=TRUE,col="aquamarine4",main="survival prediction using clinical covariates only (coxph model)
               cross-validation (nfolds=5, iterations=100)", ylab="concordance index")





TOTAL <- as.data.frame(cbind(c(1:j),CI_TS,CI_VS))
t <- dim(TOTAL)[1]
library(lattice)
TOTAL <- TOTAL[order(TOTAL$CI_VS),]
xyplot(TOTAL$CI_TS + TOTAL$CI_VS~c(1:t),auto.key=list(columns=2),
       xlab="combinations of clinical covariates",main="predicting survival in Dir using clin covariates only (cross validation)",
       ylab="concordance index", type="l")



print("max CI_VS")
print(max(CI_VS))


print("variables included in the model for maximizing CI_TS:")
colnames(dat1[,which(GRID[which(CI_TS==max(CI_TS)),]==1)])

print("variables included in the model for maximizing CI_VS:")
colnames(dat2[,which(GRID[which(CI_VS==max(CI_VS)),]==1)])


plot(CI_TS,ylim=c(0,1), main="Cooncordance index in the Training set")
abline(h=.5,col="red",lty=2)

summary(coxph(blah2~.,data=dat2[,v]))$concordance

#########

mySurv <- function(x) {summary(coxph(Surv(dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH,dat1$VITAL_STATUS)~x))$logtest[3]}
mySurv(dat1$P_Stage)
mySurv(dat1$Age)
mySurv(dat1$GENDER)
mySurv(dat1$H_Grade)
mySurv(dat1$RACE)
mySurv(dat1$SMOKING)
mySurv(dat1$SITE)




