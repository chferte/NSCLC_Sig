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
load("Dir_Barcode_29dec11.Robject")
load("Zhu_Barcode_29dec11.Robject")
load("Hou_Barcode_30dec11.Robject")
load("Gir_Barcode_30dec11.Robject")

DirExp <- Dir_Barcode
ZhuExp <- Zhu_Barcode
HouExp <- Hou_Barcode
GirExp <- Gir_Barcode

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
tmp1 <- sample(tmp,150)
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

y_SET1 <- SET1[22216,]
y_OS_SET1 <- t(SET1[c(22217,22218),])
colnames(y_OS_SET1) <-c("time","status")
SET1  <- SET1[-22216,]
SET1  <- SET1[-22217,]
SET1  <- SET1[-22218,]
dim(SET1)
length(y_SET1)
dim(y_OS_SET1)

y_SET2 <- SET2[22216,]
y_OS_SET2 <- t(SET2[c(22217,22218),])
colnames(y_OS_SET2) <-c("time","status")
SET2  <- SET2[-22216,]
SET2  <- SET2[-22217,]
SET2  <- SET2[-22218,]
dim(SET2)
length(y_SET2)
dim(y_OS_SET2)


##############################################################################################
#####      Combining the data Clin and Exp in each Dataset
##############################################################################################

#combine the DirExp1 and the DirClin1
MATRIX_TS <- rbind(DirExp1,DirClin1)
y_TS <- MATRIX_TS[22216,]
y_OS_TS <- t(MATRIX_TS[c(22217,22218),])
colnames(y_OS_TS) <-c("time","status")
MATRIX_TS  <- MATRIX_TS[-22218,]
MATRIX_TS  <- MATRIX_TS[-22217,]
MATRIX_TS  <- MATRIX_TS[-22216,]
dim(MATRIX_TS)
length(y_TS)
dim(y_OS_TS)

##combine the ZhuExp1 and the ZhuClin1
MATRIX_VS <- rbind(ZhuExp1,ZhuClin1)
y_VS <- MATRIX_VS[22216,]
y_OS_VS <- t(MATRIX_VS[c(22217,22218),])
colnames(y_OS_VS) <-c("time","status")
MATRIX_VS  <- MATRIX_VS[-22218,]
MATRIX_VS  <- MATRIX_VS[-22217,]
MATRIX_VS  <- MATRIX_VS[-22216,]
dim(MATRIX_VS)
length(y_VS)
dim(y_OS_VS)

##combine the HouExp1 and the HouClin1
MATRIX_VS2 <- rbind(HouExp1,HouClin1)
y_VS2 <- MATRIX_VS2[22216,]
y_OS_VS2 <- t(MATRIX_VS2[c(22217,22218),])
colnames(y_OS_VS2) <-c("time","status")
MATRIX_VS2  <- MATRIX_VS2[-22218,]
MATRIX_VS2  <- MATRIX_VS2[-22217,]
MATRIX_VS2  <- MATRIX_VS2[-22216,]
dim(MATRIX_VS2)
length(y_VS2)
dim(y_OS_VS2)

##combine the GirExp1 and the GirClin1
MATRIX_VS3 <- rbind(GirExp1,GirClin1)
y_VS3 <- MATRIX_VS3[22216,]
y_OS_VS3 <- t(MATRIX_VS3[c(22217,22218),])
colnames(y_OS_VS3) <-c("time","status")
MATRIX_VS3  <- MATRIX_VS3[-22218,]
MATRIX_VS3  <- MATRIX_VS3[-22217,]
MATRIX_VS3  <- MATRIX_VS3[-22216,]
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
rm(X)





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

#png("Diff_Surv_datasets.png")
plot(survfit(Surv(tot1$time,tot1$status)~tot1$Dataset))
survfit(Surv(tot1$time,tot1$status)~tot1$Dataset)
survdiff(Surv(tot1$time, tot1$status)~factor(tot1$Dataset),rho=0)
#dev.off()

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


#MATRIX_TS <- SET1
#y_OS_TS <- y_OS_SET1
#MATRIX_VS <- SET2
#y_OS_VS <- y_OS_SET2
LAMBDA_1se <- c()
LAMBDA_min <- c()
#png("cox_cv_set1_frma.png")
#alpha <- .05
alphas <- c(0,.2,.4,.6,.8,1)
#par(mfrow=c(2,3))
x <- t(MATRIX_TS)
y <-  Surv(y_OS_TS[,1], y_OS_TS[,2])
for(alpha in alphas)
  {
  CV <- cv.glmnet(x,y,maxit=5000,family="cox",alpha=alpha)
  LAMBDA_1se <- c(LAMBDA_1se,CV$lambda.1se)
  #print(LAMBDA_1se)
  LAMBDA_min <- c(LAMBDA_min,CV$lambda.min)
  #print(LAMBDA_min)
#plot(CV)
}
  
  #dev.off()

fit <- glmnet(x,y,alpha=alphas, lambda=CV$lambda.min,family="cox")    
  #Coeff <- coef(fit,s=CV$lambda.min)
  #Active_index <- which(Coeff != 0)
  #Active_coefficients <- Coeff[Active_index]
  #for (f in colnames(mTS)[Active_index]) { print(f); print(summary(coxph(Surv(y_OS_TS[,1], y_OS_TS[,2]) ~ mTS[,f]))$coefficients)                                       
   #                                        plot(survfit(y ~ 1))                                      
    #                                             }
  #}
#dev.off()
#print(CV)
print(fit)
#Z <- y_OS_VS
yE_TS <- predict (fit,t(MATRIX_TS),type="link")
yE_VS <- predict (fit,t(MATRIX_VS),type="link")
yE_VS2 <- predict (fit,t(MATRIX_VS2),type="link")
yE_VS3 <- predict (fit,t(MATRIX_VS3),type="link")

#png("AUC_time.png")
AUC_COX_TS <- risksetAUC(Stime=y_OS_TS[,1],status=y_OS_TS[,2], marker=yE_TS, method="Cox", tmax=200, col="red", lwd=3,lty=2,main="AUC over the time (TS VS VS2 VS3)")
AUC_COX_VS <- risksetAUC(Stime=y_OS_VS[,1],status=y_OS_VS[,2], marker=yE_VS, method="Cox", span=.4, tmax=200, plot=FALSE)
AUC_COX_VS2 <- risksetAUC(Stime=y_OS_VS2[,1],status=y_OS_VS2[,2], marker=yE_VS2, method="Cox", span=.4, tmax=200, plot=FALSE)
AUC_COX_VS3 <- risksetAUC(Stime=y_OS_VS3[,1],status=y_OS_VS3[,2], marker=yE_VS3, method="Cox", span=.4, tmax=200, plot=FALSE)

lines(AUC_COX_VS$utimes, AUC_COX_VS$AUC,col="green", lty=1, lwd=3)
lines(AUC_COX_VS2$utimes, AUC_COX_VS2$AUC,col="darkblue", lty=3,lwd=3)
lines(AUC_COX_resp_e$utimes, AUC_COX_link_e$AUC,col="green", lty=4,lwd=3)
#dev.off()

#yETS <- predict (fit,t(MATRIX_VS),type="response",s=CV$lambda.min)
#blah <- ifelse (yE<=median(yETS),"low-risk of Death","high-risk of death")
#Z <- as.data.frame(Z) 
#Z$blah<- blah
#survfit(Surv(Z[,1], Z[,2])~factor(Z[,3]))
#survdiff(Surv(Z[,1], Z[,2])~factor(Z[,3]),rho=0)
#png("KM_fit_set2.png")
#plot(survfit(Surv(Z[,1], Z[,2])~factor(Z[,3])))
#dev.off()


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
png("ROC_binomial_CV_DIR.png")
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
dev.off()
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
## start LASSO on MATRIX  BINOMIAL (3 YEARS OS)
############################################################################################################################
############################################################################################################################

library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)

png("binomial_cv_Dir_frma.png")
CV <- cv.glmnet(t(MATRIX_TS),as.numeric(y_TS),nfolds=5,family="binomial")
plot(CV)
CV
dev.off()

fit <- glmnet(t(MATRIX_TS),as.numeric(y_TS),alpha=1, family="binomial")
y_E <- t(MATRIX_VS) %*% fit$beta + matrix(rep(fit$a0,each=41),nrow=41)
y_E <- sign(y_E)
y_E <- as.matrix(y_E)
y_E <- ifelse(y_E==-1,0,1)






ACC_y <- c()
SAMPLE <- dim(y_E)[1]
for(i in 1:100) {
ACC_y <- c(ACC_y,(SAMPLE-sum(xor(y_E[,i],y_VS)))/SAMPLE)  
}
max(ACC_y)
png(file="Accuracy_Binomial_prediction_Zhu_frma.png")
plot(ACC_y,pch=20)
dev.off()
dim(MATRIX_TS)
length(y_TS)

                    ############################################################################################################################
############################################################################################################################
######### start Ridge on MATRIX  ---   COX proportional model
############################################################################################################################
############################################################################################################################
setwd("/home/cferte/FELLOW/cferte/RIDGE/")
library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)

## cross validation in the Training set (TS)
png("cox_cv_Dir_frma.png")
CV <- cv.glmnet(t(MATRIX_TS),alpha=0.1,y_OS_TS,nfolds=5,family="cox")
plot(CV)
CV
dev.off()

fit <- glmnet(t(MATRIX_TS),y_OS_TS,alpha=0.1, lambda=CV$lambda.1se,family="cox")
### y_E <- t(MATRIX_VS) %*% fit$beta
y_E <- predict(fit,t(MATRIX_TS),type="response")


tmp1 <- which(y_OS_VS[,2]==1)
Y_REF <- y_OS_VS[tmp1,]
Y_REF <- Y_REF[,-2]
RMSE_y <- c()
R2_y <- c()
for(i in 1:dim(y_E)[2]) { 
RMSE_y <- c(RMSE_y,  sum((Y_REF - y_E[tmp1,i])^2))
R2_y <- c(R2_y,cor(Y_REF,y_E[tmp1,i]))
}
png("Accuracy_COX_RMSE_into_Zhu_frma.png")
plot(RMSE_y)
dev.off()
png("Accuracy_COX_R2_into_Zhu_frma.png")
plot(R2_y)
dev.off()


min(RMSE_y)
which(RMSE_y==min(RMSE_y))
