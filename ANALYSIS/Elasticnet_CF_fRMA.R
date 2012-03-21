### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

### script for running modelling prediction

#load affy
library(affy)
library(survival)

# set Dir CEL files directory
setwd("/home/cferte/FELLOW/cferte/FRMA_BARCODE/")

#read the files
load("Dir_frma_29dec11.Robject")
load("Zhu_frma_29dec11.Robject")
load("Hou_frma_30dec11.Robject")
load("Gir_frma_30dec11.Robject")

DirExp <- exprs(Dirfrma)
ZhuExp <- exprs(Zhufrma)
HouExp <- exprs(Houfrma)
GirExp <- exprs(Girfrma)

# set Dir CEL files directory
setwd("/home/cferte/FELLOW/cferte/CLIN_DATA_FILES/")

#read the files
DirClin <- read.csv("DirectorData.csv")
ZhuClin <- read.csv("ZhuData.csv")
GirClin <- read.csv("GirardData.csv")
HouClin <- read.csv("HouData.csv")

rownames(DirClin) <- DirClin$CEL_ID
rownames(ZhuClin) <- ZhuClin$CEL_ID
rownames(GirClin) <- GirClin$CEL_ID
rownames(HouClin) <- HouClin$CEL_ID

## change ALive -> 0 & dead -> 1 in DirClin ZhuClin
DirClin$VITAL_STATUS<- ifelse(DirClin$VITAL_STATUS=="Dead",1,0)
ZhuClin$VITAL_STATUS<- ifelse(ZhuClin$VITAL_STATUS=="Dead",1,0)
HouClin$VITAL_STATUS<- ifelse(HouClin$VITAL_STATUS=="Dead",1,0)
GirClin$VITAL_STATUS<- ifelse(GirClin$VITAL_STATUS=="Dead",1,0)

# transform into numeric:
HouClin$MONTHS_TO_LAST_CONTACT_OR_DEATH <- as.numeric(HouClin$MONTHS_TO_LAST_CONTACT_OR_DEATH)


# create a new variable named Three years OS (probability of overall survival)
DirClin$THREE_YEAR_OS <- ifelse(DirClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & DirClin$VITAL_STATUS==0,1,0)
ZhuClin$THREE_YEAR_OS <- ifelse(ZhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & ZhuClin$VITAL_STATUS==0,1,0)
GirClin$THREE_YEAR_OS <- ifelse(GirClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & GirClin$VITAL_STATUS==0,1,0)
HouClin$THREE_YEAR_OS <- ifelse(HouClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & HouClin$VITAL_STATUS==0,1,0)


## in the Dir dataset: 
## rename the moffit CEL files that are differently named in the Exp file and in the Clin files
## (change all the " " and the "." into "_")

## do it in the Dir Clin file
rownames(DirClin) <- sub(" ","_",rownames(DirClin))
rownames(DirClin) <- paste(rownames(DirClin),".CEL",sep="")

DirClin$CEL_ID <- rownames(DirClin)

## do it in the Dir Exp file
colnames(DirExp) <- sub("Moff.","Moff_",colnames(DirExp))

## do it in the clin Zhu file
rownames(ZhuClin) <-  paste(rownames(ZhuClin),".CEL",sep="")



# select in DirClin and ZhuClin only the patients without any adjuvant treatment:
DirClin1  <- DirClin[DirClin$totalAdjuvant==0,]
ZhuClin1  <- ZhuClin[ZhuClin$totalAdjuvant==0,]
GirClin1  <- GirClin[GirClin$totalAdjuvant==0,]
HouClin1  <- HouClin[HouClin$totalAdjuvant==0,]

#create a matrix with only the CEL_ID and the 3YOS
DirClin1 <- t(DirClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")])
ZhuClin1 <- t(ZhuClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")])
GirClin1 <- t(GirClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")])
HouClin1 <- t(HouClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")])



####################################################################
### match the columns of Exp files with the rows of the Clin1 files
####################################################################

#### match the DIR file ####

# two colnames of DirClin1 were inappropriately written and do not match with DirExp ones:
colnames(DirClin1)[126] <- colnames(DirExp)[200]
colnames(DirClin1)[181]  <- colnames(DirExp)[257]

# create DirExp1 where are only the samples that match with the samles in DirClin1
tmp1 <- intersect(colnames(DirClin1),colnames(DirExp))
DirExp1 <- DirExp[,tmp1]
DirClin1 <- DirClin1[,tmp1]
rm(tmp1)

#check if the names are the same between DirClin1 and DirExp1
table(colnames(DirClin1)==colnames(DirExp1))


#### match the ZHU file ####

# create ZhuExp1 where are only the samples that match with the samles in ZhuClin1
tmp1 <- intersect(colnames(ZhuClin1),colnames(ZhuExp))
ZhuExp1 <- ZhuExp[,tmp1]
ZhuClin2 <- ZhuClin1[,tmp1]
rm(tmp1)
#check if the names are the same between ZhuClin1 and ZhuExp1
table(colnames(ZhuClin2)==colnames(ZhuExp1))


#### match the HOU file ####

## check the Hou file: if names are the same between Exp and Clin file
colnames(HouClin1) <- paste(colnames(HouClin1),".CEL.gz",sep="")
tmp1 <- match(colnames(HouClin1),colnames(HouExp))
HouExp1 <- HouExp[,tmp1]
rm(tmp1)
# check if the names are the same between ZhuClin1 and ZhuExp1
table(colnames(HouClin1)==colnames(HouExp1))

#### match the GIR file ####

## check the Gir file: if names are the same between Exp and Clin file
colnames(GirClin1) <- paste(colnames(GirClin1),".CEL.gz",sep="")
tmp1 <- match(colnames(GirClin1),colnames(GirExp))
GirExp1 <- GirExp[,tmp1]
rm(tmp1)
# check if the names are the same between ZhuClin1 and ZhuExp1
table(colnames(GirClin1)==colnames(GirExp1))

##############################################################################################
## keep the Gir & Hou Probes that are common with the Dir & Zhu Probes
##############################################################################################
tmp1 <- rownames(GirExp1)[match(rownames(DirExp1),rownames(GirExp1))]
tmp1 <- tmp1[!is.na(tmp1)]
length(tmp1)
head(tmp1)
GirExp1 <- GirExp1[tmp1,]
DirExp1 <- DirExp1[tmp1,]
ZhuExp1 <- ZhuExp1[tmp1,]
HouExp1 <- HouExp1[tmp1,]
rm(tmp1)
# check if the names are the same between DirExp1 and GirExp1
table(rownames(DirExp1)==rownames(GirExp1))
table(rownames(ZhuExp1)==rownames(GirExp1))
table(rownames(HouExp1)==rownames(GirExp1))

##############################################################################################
#####      using sample to to cross validation
##############################################################################################
tmp <- c(1:307)
tmp1 <- sample(tmp,154)
tmp2 <- setdiff(tmp,tmp1)
setE1 <- DirExp1[,tmp1]
setE2 <- DirExp1[,tmp2]
setC1 <- DirClin1[,tmp1]
setC2 <- DirClin1[,tmp2]
dim(setE1) 
dim(setC1)
dim(setE2)
dim(setC2)
table(colnames(setC1)==colnames(setC2))
table(colnames(setE1)==colnames(setC1))
table(colnames(setE2)==colnames(setC2))

SET1 <- rbind(setE1,setC1)
SET2 <- rbind(setE2,setC2)
rm(tmp,tmp1,tmp2,setE1,setE2,setC1,setC2)

y_SET1 <- SET1[22278,]
y_OS_SET1 <- t(SET1[c(22279,22280),])
colnames(y_OS_SET1) <-c("time","status")
SET1  <- SET1[-22280,]
SET1  <- SET1[-22279,]
SET1  <- SET1[-22278,]
dim(SET1)
length(y_SET1)
dim(y_OS_SET1)

y_SET2 <- SET2[22278,]
y_OS_SET2 <- t(SET2[c(22279,22280),])
colnames(y_OS_SET2) <-c("time","status")
SET2  <- SET2[-22280,]
SET2  <- SET2[-22279,]
SET2  <- SET2[-22278,]
dim(SET2)
length(y_SET2)
dim(y_OS_SET2)


##############################################################################################
#####      Combining the data Clin and Exp in each Dataset
##############################################################################################

#combine the DirExp1 and the DirClin1
MATRIX_TS <- rbind(DirExp1,DirClin1)
y_TS <- MATRIX_TS[22278,]
y_OS_TS <- t(MATRIX_TS[c(22279,22280),])
colnames(y_OS_TS) <-c("time","status")
MATRIX_TS  <- MATRIX_TS[-22280,]
MATRIX_TS  <- MATRIX_TS[-22279,]
MATRIX_TS  <- MATRIX_TS[-22278,]
dim(MATRIX_TS)
length(y_TS)
dim(y_OS_TS)

##combine the ZhuExp1 and the ZhuClin1
MATRIX_VS <- rbind(ZhuExp1,ZhuClin1)
y_VS <- MATRIX_VS[22278,]
y_OS_VS <- t(MATRIX_VS[c(22279,22280),])
colnames(y_OS_VS) <-c("time","status")
MATRIX_VS  <- MATRIX_VS[-22280,]
MATRIX_VS  <- MATRIX_VS[-22279,]
MATRIX_VS  <- MATRIX_VS[-22280,]
dim(MATRIX_VS)
length(y_VS)
dim(y_OS_VS)

##combine the HouExp1 and the HouClin1
MATRIX_VS2 <- rbind(HouExp1,HouClin1)
y_VS2 <- MATRIX_VS2[22278,]
y_OS_VS2 <- t(MATRIX_VS2[c(22279,22280),])
colnames(y_OS_VS2) <-c("time","status")
MATRIX_VS2  <- MATRIX_VS2[-22280,]
MATRIX_VS2  <- MATRIX_VS2[-22279,]
MATRIX_VS2  <- MATRIX_VS2[-22278,]
dim(MATRIX_VS2)
length(y_VS2)
dim(y_OS_VS2)

##combine the GirExp1 and the GirClin1
MATRIX_VS3 <- rbind(GirExp1,GirClin1)
y_VS3 <- MATRIX_VS3[22278,]
y_OS_VS3 <- t(MATRIX_VS3[c(22279,22280),])
colnames(y_OS_VS3) <-c("time","status")
MATRIX_VS3  <- MATRIX_VS3[-22280,]
MATRIX_VS3  <- MATRIX_VS3[-22279,]
MATRIX_VS3  <- MATRIX_VS3[-22278,]
dim(MATRIX_VS3)
length(y_VS3)
dim(y_OS_VS3)



##############################################################################################
### rescale the Zhu according to the Dir and call the new p . n matrix YSCALED
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
MATRIX_VS_S <- normalize_to_X(mean_x,sd_x,MATRIX_VS)
MATRIX_VS2_S <- normalize_to_X(mean_x,sd_x,MATRIX_VS2)
MATRIX_VS3_S <- normalize_to_X(mean_x,sd_x,MATRIX_VS3)






############################################################################################################################
## ### check survival curves in TS & VS
############################################################################################################################
Sur1 <- as.data.frame(y_OS_TS)
Sur2 <- as.data.frame(y_OS_VS)
Sur3 <- as.data.frame(y_OS_VS2)
Sur4 <- as.data.frame(y_OS_VS3)
Sur1$Dataset <- "Dir"
Sur2$Dataset <- "Zhu"
Sur3$Dataset <- "Hou"
Sur4$Dataset <- "Gir"
tot1 <- rbind(Sur1,Sur2)
tot1 <- rbind(tot1,Sur3)
tot1 <- rbind(tot1,Sur4)

png("Diff_Surv_datasets.png")
plot(survfit(Surv(tot1$time,tot1$status)~tot1$Dataset))
survfit(Surv(tot1$time,tot1$status)~tot1$Dataset)
survdiff(Surv(tot1$time, tot1$status)~factor(tot1$Dataset),rho=0)
dev.off()

rm(Sur1,Sur2,Sur3,Sur4)

############################################################################################################################
############################################################################################################################
######### start ElasticNet --   COX proportional model
############################################################################################################################
############################################################################################################################
library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)
library(survival)
library(risksetROC)

x <- t(MATRIX_TS)
y <-  Surv(y_OS_TS[,1], y_OS_TS[,2])

LAMBDA_1se <- c()
LAMBDA_min <- c()
#png("cox_cv_set1_frma.png")
#alpha <- .05
alphas <-c(.01,.1,.2,.3,.4,.6,.7,.8,.85,.9,.95,.99)
par(mfrow=c(4,3))

for(alpha in alphas)
  {
  CV <- cv.glmnet(x,y,maxit=1000,family="cox",alpha=alphas)
  LAMBDA_1se <- c(LAMBDA_1se,CV$lambda.1se)
  LAMBDA_min <- c(LAMBDA_min,CV$lambda.min)
}
  fit <- glmnet(x,y,alpha=alphas, lambda=CV$lambda.min,family="cox")  
  print(fit)

yE_TS <- predict (fit,t(MATRIX_TS),type="link")
yE_VS <- predict (fit,t(MATRIX_VS),type="link")
yE_VS2 <- predict (fit,t(MATRIX_VS2),type="link")
yE_VS3 <- predict (fit,t(MATRIX_VS3),type="link")

AUC_COX_TS <- risksetAUC(Stime=y_OS_TS[,1],status=y_OS_TS[,2], marker=yE_TS, method="Cox", tmax=200, col="red", lwd=3,lty=2,main="AUC over the time (TS VS VS2 VS3)")
AUC_COX_VS <- risksetAUC(Stime=y_OS_VS[,1],status=y_OS_VS[,2], marker=yE_VS, method="Cox", span=.4, tmax=200, plot=FALSE)
AUC_COX_VS2 <- risksetAUC(Stime=y_OS_VS2[,1],status=y_OS_VS2[,2], marker=yE_VS2, method="Cox", span=.4, tmax=200, plot=FALSE)
AUC_COX_VS3 <- risksetAUC(Stime=y_OS_VS3[,1],status=y_OS_VS3[,2], marker=yE_VS3, method="Cox", span=.4, tmax=200, plot=FALSE)

lines(AUC_COX_VS$utimes, AUC_COX_VS$AUC,col="green", lty=1, lwd=3)
lines(AUC_COX_VS2$utimes, AUC_COX_VS2$AUC,col="darkblue", lty=3,lwd=3)
lines(AUC_COX_resp_e$utimes, AUC_COX_link_e$AUC,col="green", lty=4,lwd=3)
}
print(alphas)
print(LAMBDA_min)


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

