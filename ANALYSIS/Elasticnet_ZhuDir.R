### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

### script for running modelling prediction

#load the packages
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

  
# point the directory (choose method among = RMA, GCRMA, MAS5, dCHIP, metaGEO, fRMA or barcode)
method= "dCHIP"
PATH <- "/home/cferte/FELLOW/cferte/NSCLC_MA/ALLDATA_NORMALIZED/"
setwd(PATH)
load("AllData_MAS5.Rdata")

## load the matrix and response files
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/MAS5/")


load("y_TS.Rdata")
load("y_VS.Rdata")
load("y_OS_TS.Rdata")
load("y_OS_VS.Rdata")

T_CLIN <- rbind(DirClinF,ZhuClinF)
T_CLIN <- rbind(T_CLIN,HouClinF)
T_CLIN <- rbind(T_CLIN,GirClinF)
T_CLIN <- as.data.frame(T_CLIN)
T_CLIN$Study <- substr(T_CLIN$LABORATORY_BATCH,1,3)


AllData_MAS5 <- exprs(AllData_MAS5)
colnames(AllData_MAS5) <- sub(" ","_",colnames(AllData_MAS5))
tmp1 <- intersect(colnames(AllData_MAS5),names(y_VS))
MATRIX_VS <- AllData_MAS5[,tmp1]
tmp2 <- intersect(colnames(AllData_MAS5),names(y_TS))
MATRIX_TS <- AllData_MAS5[,tmp2]


##############################################################################################
### rescale the VS according to the TS and call the new p . n matrix YSCALED
##############################################################################################
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
	sd.y <- apply(Y, 1, sd)
	Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x
	Y.adj
}

X <- MATRIX_TS     # this my p x n training set
mean_x  <- apply(X,1,mean)
sd_x <- apply(X,1,sd)

MATRIX_TS_S <- MATRIX_TS
MATRIX_VS_S <- normalize_to_X(mean_x,sd_x,MATRIX_VS)


############################################################################################################################
############################################################################################################################
######### start ElasticNet --   COX proportional model
############################################################################################################################
############################################################################################################################

x <- t(MATRIX_TS)
y <-  Surv(y_OS_TS[,1],y_OS_TS[,2])


CI_TS <- c()
CI_VS <- c()
CI_VS2 <- c()
CI_TS_S <- c()
CI_VS_S <- c()
CI_VS2_S <- c()




alphas <- exp(-1*seq(0,10,1))
lambdas <- exp(seq(-4,3,1))


GRID <- expand.grid(.family="cox",.alpha=alphas,.lambda=lambdas)



for(i in 1:dim(GRID)[1]){
  fit  <- try(glmnet(x,y,family="cox",alpha=GRID$.alpha[i],lambda=GRID$.lambda[i]))

  if( class(fit) == "try-error" ){
    CI_TS <- c(CI_TS, NA)
    CI_VS <- c(CI_VS, NA)
    CI_VS2 <- c(CI_VS2, NA)
    CI_TS_S <- c(CI_TS_S, NA)
    CI_VS_S <- c(CI_VS_S, NA)
    CI_VS2_S <- c(CI_VS2_S, NA)
  } else{
y_E_TS <- predict(fit,x,type="link")
CI_TS <- rbind(CI_TS,concordance.index(y_E_TS,y_OS_TS[,1],y_OS_TS[,2],na.rm=T,method="noether"))
  
y_E_VS <- predict(fit,t(MATRIX_VS),type="link")
CI_VS <- rbind(CI_VS,concordance.index(y_E_VS,y_OS_VS[,1],y_OS_VS[,2],na.rm=T,method="noether"))


y_E_TS_S <- predict(fit,x,type="link")
CI_TS_S <- rbind(CI_TS_S,concordance.index(y_E_TS_S,y_OS_TS[,1],y_OS_TS[,2],na.rm=T,method="noether"))

y_E_VS_S <- predict(fit,t(MATRIX_VS_S),type="link")
CI_VS_S <- rbind(CI_VS_S,concordance.index(y_E_VS_S,y_OS_VS[,1],y_OS_VS[,2],na.rm=T,method="noether"))


    }
}

summary(as.numeric(CI_TS[names(CI_TS)=="c.index"][is.na(CI_TS[names(CI_TS)=="c.index"])==F]))
summary(as.numeric(CI_VS[names(CI_VS)=="c.index"][is.na(CI_VS[names(CI_VS)=="c.index"])==F]))

summary(as.numeric(CI_TS_S[names(CI_TS_S)=="c.index"][is.na(CI_TS_S[names(CI_TS_S)=="c.index"])==F]))
summary(as.numeric(CI_VS_S[names(CI_VS_S)=="c.index"][is.na(CI_VS_S[names(CI_VS_S)=="c.index"])==F]))

GRID[which(CI_TS==max(as.numeric(CI_TS[names(CI_TS)=="c.index"][is.na(CI_TS[names(CI_TS)=="c.index"])==F]))),]
GRID[which(CI_VS==max(as.numeric(CI_VS[names(CI_VS)=="c.index"][is.na(CI_VS[names(CI_VS)=="c.index"])==F]))),]


CI_TOTAL <- cbind(GRID$.alpha,GRID$.lambda,CI_TS,CI_VS)
rownames(CI_TOTAL)<- c(1:88)
colnames(CI_TOTAL)[1]<-"alpha"
colnames(CI_TOTAL)[2]<-"lambda"
CI_TOTAL<-as.data.frame(CI_TOTAL)
CI_TOTAL$GAL <- ifelse(is.na(CI_TOTAL$CI_TS),"blue","red")
CI_TOTAL$method <- method

setwd("~/FELLOW/cferte/NSCLC_MA/ANALYSIS/results_elasticnet_grid/")

png(paste("GRID_",method,".png",sep=""))
plot(log(as.numeric(CI_TOTAL$alpha)),log(as.numeric(CI_TOTAL$lambda)), col=CI_TOTAL$GAL, main=paste("alpha lambda GRID for",method),xlab="log(alpha)",ylab="log(lambda)",pch=20)
dev.off()

tmp <- paste("CI_",method,sep="")
assign(tmp,CI_TOTAL)
save(list=paste("CI_",method,sep=""), file=paste("CI_",method,".Rdata",sep=""))
     
############################################################################################################################




############################################################################################################################
## start Elastic Net  ---  BINOMIAL model (3 YEARS OS)
############################################################################################################################
############################################################################################################################
setwd("/home/cferte/FELLOW/cferte/RIDGE/")

library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)



#CROSS VALIDATION
ROC_Total <- NA
CV <- NA 
LAMBDA_1se <- NA
#png("ROC_binomial_CV_DIR.png")
alphas <- c(.01,.1,.2,.3,.4,.6,.7,.8,.85,.9,.95,.99)
 par(mfrow=c(4,3))
 for(alpha in alphas){
  CV <- cv.glmnet(t(MATRIX_TS),as.numeric(y_TS),maxit=5000,family="binomial",alpha=alpha)
  #plot(CV)    
 #}
#dev.off()
  LAMBDA_1se <- c(LAMBDA_1se,CV$lambda.1se)
  fit <- glmnet(t(MATRIX_TS),as.numeric(y_TS),alpha=alpha, lambda=CV$lambda.1se,family="binomial")
y_E <- predict(fit,t(MATRIX_TS),type="response")
 prediction.obj <- prediction(y_E,y_TS) 
  ROC <- performance(prediction.obj,"tpr","fpr")
plot(ROC)
  ROC <- performance(prediction.obj,"auc")
ROC_Total <- c(ROC_Total, ROC@y.values[[1]])
  }
#dev.off()
print(alphas)
print(LAMBDA_1se)
print(ROC_Total)

## VALIDATION into Zhu
ROC_Total <- NA
CV <- NA 
LAMBDA_1se <- NA
png("ROC_Binomial_Validation_Zhu.png")
alphas <- c(.01,.1,.2,.3,.4,.6,.7,.8,.85,.9,.95,.99)
 par(mfrow=c(4,3))
 for(alpha in alphas){
  CV <- cv.glmnet(t(MATRIX_TS),as.numeric(y_TS),nfolds=5,family="binomial",alpha=alpha)
   LAMBDA_1se <- c(LAMBDA_1se,CV$lambda.1se)    
  fit <- glmnet(t(MATRIX_TS),as.numeric(y_TS),alpha=alpha, lambda=CV$lambda.1se,family="binomial")
y_E <- predict(fit,t(MATRIX_VS),type="response")
 prediction.obj <- prediction(y_E,y_VS) 
  ROC <- performance(prediction.obj,"tpr","fpr")
plot(ROC)
  ROC <- performance(prediction.obj,"auc")
ROC_Total <- c(ROC_Total, ROC@y.values[[1]])
  }
dev.off()
print(alphas)
print(LAMBDA_1se)
print(ROC_Total)
#

## VALIDATION into Hou
ROC_Total <- NA
CV <- NA 
LAMBDA_1se <- NA
png("ROC_Binomial_Validation_Hou.png")
alphas <- c(.01,.1,.2,.3,.4,.6,.7,.8,.85,.9,.95,.99)
 par(mfrow=c(4,3))
 for(alpha in alphas){
  CV <- cv.glmnet(t(MATRIX_TS),as.numeric(y_TS),nfolds=5,family="binomial",alpha=alpha)
   LAMBDA_1se <- c(LAMBDA_1se,CV$lambda.1se)    
  fit <- glmnet(t(MATRIX_TS),as.numeric(y_TS),alpha=alpha, lambda=CV$lambda.1se,family="binomial")
y_E <- predict(fit,t(MATRIX_VS2),type="response")
 prediction.obj <- prediction(y_E,y_VS2) 
  ROC <- performance(prediction.obj,"tpr","fpr")
plot(ROC)
  ROC <- performance(prediction.obj,"auc")
ROC_Total <- c(ROC_Total, ROC@y.values[[1]])
  }
dev.off()
print(alphas)
print(LAMBDA_1se)
print(ROC_Total)
#


paste ("Working with the following alpha values:",alphas, sep= " ")
paste("AUCs of the prediction applied in the validation set: ", ROC_Total, sep="")


setwd("/home/cferte/FELLOW/cferte/LASSO/")
      
############################################################################################################################
## start Elastic Net  ---  BINOMIAL model (3 YEARS OS) on a 50 times loop
############################################################################################################################
############################################################################################################################
setwd("/home/cferte/FELLOW/cferte/RIDGE/")

library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)


#CROSS VALIDATION
ROC_Total <- NA
CV <- NA 
LAMBDA_1se <- NA
y_E <- NA
N <- NA
#png("ROC_binomial_CV_DIR.png")
alphas <- c(.01)
FIT <- NA
for(i in c(1:50)) {
  tmp <- c(1:307)
tmp1 <- sample(tmp,150)
setE1 <- DirExp1[,tmp1]
setC1 <- DirClin1[,tmp1]
SET1 <- rbind(setE1,setC1)
y <- SET1[22278,]
SET1 <- SET1[-22280,]
  SET1 <- SET1[-22279,]
  SET1 <- SET1[-22278,]
  x <- t(SET1)
  CV <- cv.glmnet(x,as.numeric(y),maxit=2000,family="binomial",alpha=.01)
  LAMBDA_1se <- c(LAMBDA_1se,CV$lambda.1se)
  fit <- glmnet(x,as.numeric(y),alpha=.01, lambda=CV$lambda.1se,family="binomial")
FIT <- c(FIT,fit$beta)
  y_E <- c(y_E,predict(fit,x,type="response"))
 N <- c(N,i) }
save(FIT,file="FIT_50_Boostrap.Rdata")
    save(y_E,file="yE_binomial_50_bootstrap.Rdata")
 prediction.obj <- prediction(y_E,y) 
  ROC <- performance(prediction.obj,"tpr","fpr")
plot(ROC)
  ROC <- performance(prediction.obj,"auc")
ROC_Total <- c(ROC_Total, ROC@y.values[[1]])
  }
#dev.off()
print(alphas)
print(LAMBDA_1se)
print(ROC_Total)

