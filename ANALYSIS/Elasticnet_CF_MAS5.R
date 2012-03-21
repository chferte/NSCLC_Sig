### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

A <- survfit(Surv(DirClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,DirClin$VITAL_STATUS)~DirClin$totalAdjuvant)
B <- survfit(Surv(ZhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH,ZhuClin$VITAL_STATUS)~ZhuClin$totalAdjuvant)
C <- survfit(Surv(DirClin1$MONTHS_TO_LAST_CONTACT_OR_DEATH,DirClin1$VITAL_STATUS)~1)
D <- survfit(Surv(ZhuClin1$MONTHS_TO_LAST_CONTACT_OR_DEATH,ZhuClin1$VITAL_STATUS)~1)

### script for running modelling prediction

#load affy
library(affy)
library(survival)

# set Dir CEL files directory
setwd("/Volumes/cferte/FELLOW/cferte/MAS5_NORM/")

#read the files
DirExp <- exprs(readExpressionSet("Dir_MAS5_28dec11.txt"))
ZhuExp <- exprs(readExpressionSet("Zhu_MAS5_28dec11.txt"))

# set Dir CEL files directory
setwd("/Volumes/cferte/FELLOW/cferte/CLIN_DATA_FILES/")

#read the files
DirClin <- read.csv("DirectorData.csv")
ZhuClin <- read.csv("ZhuData.csv")

rownames(DirClin) <- DirClin$CEL_ID
rownames(ZhuClin) <- ZhuClin$CEL_ID

## change ALive -> 0 & dead -> 1 in DirClin ZhuClin
DirClin$VITAL_STATUS<- ifelse(DirClin$VITAL_STATUS=="Dead",1,0)
ZhuClin$VITAL_STATUS<- ifelse(ZhuClin$VITAL_STATUS=="Dead",1,0)

# create a new variable named Three years OS (probability of overall survival)
DirClin$THREE_YEAR_OS <- ifelse(DirClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & DirClin$VITAL_STATUS==0,1,0)
ZhuClin$THREE_YEAR_OS <- ifelse(ZhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & ZhuClin$VITAL_STATUS==0,1,0)

## rename the moffit CEL files that are differently named in the Exp file and in the Clin files
##(change all the " " and the "." into "_")

## do it in the Dir Clin file
rownames(DirClin) <- sub(" ","_",rownames(DirClin))
rownames(DirClin) <- paste(rownames(DirClin),".CEL",sep="")

DirClin$CEL_ID <- rownames(DirClin)

## do it in the Dir Exp file
colnames(DirExp) <- sub("Moff.","Moff_",colnames(DirExp))

## do it in the clin Zhu file
rownames(ZhuClin) <-  paste(rownames(ZhuClin),".CEL.gz",sep="")

# select in DirClin and ZhuClin only the patients without any adjuvant treatment:
DirClin1  <- DirClin[DirClin$totalAdjuvant==0,]
ZhuClin1  <- ZhuClin[ZhuClin$totalAdjuvant==0,]




#create a mtrix with only the CEL_ID and the 3YOS
DirClin1 <- t(DirClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")])
ZhuClin1 <- t(ZhuClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")])


# two colnames of DirClin1 were inappropriately written and do not match with DirExp ones:
colnames(DirClin1)[126] <- colnames(DirExp)[200]
colnames(DirClin1)[181]  <- colnames(DirExp)[257]

# create DirExp1 where are only the samples that match with the samles in DirClin1
tmp1 <- match(colnames(DirClin1),colnames(DirExp))
DirExp1 <- DirExp[,tmp1]
rm(tmp1)
#check if the names are the same between DirClin1 and DirExp1
table(colnames(DirClin1)==colnames(DirExp1))

# create ZhuExp1 where are only the samples that match with the samles in ZhuClin1
tmp1 <- match(colnames(ZhuClin1),colnames(ZhuExp))
ZhuExp1 <- ZhuExp[,tmp1]
rm(tmp1)
#check if the names are the same between ZhuClin1 and ZhuExp1
table(colnames(ZhuClin1)==colnames(ZhuExp1))

#combine the DirExp1 and the DirClin1
MATRIX_TS <- rbind(DirExp1,DirClin1)
y_TS <- MATRIX_TS[22284,]
y_OS_TS <- t(MATRIX_TS[c(22285,22286),])
colnames(y_OS_TS) <-c("time","status")
MATRIX_TS  <- MATRIX_TS[-22286,]
MATRIX_TS  <- MATRIX_TS[-22285,]
MATRIX_TS  <- MATRIX_TS[-22284,]
dim(MATRIX_TS)
length(y_TS)
dim(y_OS_TS)

##combine the ZhuExp1 and the ZhuClin1
MATRIX_VS <- rbind(ZhuExp1,ZhuClin1)
y_VS <- MATRIX_VS[22284,]
y_OS_VS <- t(MATRIX_VS[c(22285,22286),])
colnames(y_OS_VS) <-c("time","status")
MATRIX_VS  <- MATRIX_VS[-22286,]
MATRIX_VS  <- MATRIX_VS[-22285,]
MATRIX_VS  <- MATRIX_VS[-22284,]
dim(MATRIX_VS)
length(y_VS)
dim(y_OS_VS)

### check survival curves in TS & VS
png("OS_TS.png")
OS1 <- y_OS_TS[,1]
OSS1 <- y_OS_TS[,2]
plot(survfit(Surv(OS1, OSS1)~1) ,main="Dir OS",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS1, OSS1)~1)
dev.off()

png("OS_VS.png")
OS2 <- y_OS_VS[,1]
OSS2 <- y_OS_VS[,2]
plot(survfit(Surv(OS2, OSS2)~1) ,main="Dir OS",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS2, OSS2)~1)
dev.off()



############################################################################################################################
############################################################################################################################
######### start Lasso on MATRIX  ---   COX proportional model
############################################################################################################################
############################################################################################################################

library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)

fit <- glmnet(t(MATRIX_TS),y_OS_TS,alpha=1, family="cox")
y_E <- t(MATRIX_VS) %*% fit$beta

tmp1 <- which(y_OS_VS[,2]==1)
Y_REF <- y_OS_VS[tmp1,]
Y_REF <- Y_REF[,-2]
RMSE_y <- c()
R2_y <- c()
for(i in 1:dim(y_E)[2]) { 
RMSE_y <- c(RMSE_y,  sum((Y_REF - y_E[tmp1,i])^2))
R2_y <- c(R2_y,cor(Y_REF,y_E[tmp1,i]))
}
pdf("Accuracy RMSE.pdf")
plot(RMSE_y)
dev.off()
pdf("Accuracy R2.pdf")
plot(R2_y)
dev.off()


min(RMSE_y)
which(RMSE_y==min(RMSE_y))

CV <- cv.glmnet(t(MATRIX_TS),y_OS_TS,nfolds=5,family="cox")
plot(cv_fit)
cv_fit

############################################################################################################################
## start Lasso on MATRIX  BINOMIAL (3 YEARS OS)
############################################################################################################################
############################################################################################################################

library(Biobase)
library(MASS)
library(glmnet)
library(corpcor)
library(ROCR)
library(synapseClient)

fit <- glmnet(t(MATRIX_TS),as.numeric(y_TS),alpha=1, family="binomial")
y_E <- t(MATRIX_VS) %*% fit$beta + matrix(rep(fit$a0,each=41),nrow=41)
y_E <- sign(y_E)
y_E <- as.matrix(y_E)
y_E <- ifelse(y_E==-1,0,1)

CV <- cv.glmnet(t(MATRIX_TS),as.numeric(y_TS),nfolds=5,family="binomial")
plot(CV)
CV




ACC_y <- c()
SAMPLE <- dim(y_E)[1]
for(i in 1:100) {
ACC_y <- c(ACC_y,(SAMPLE-sum(xor(y_E[,i],y_VS)))/SAMPLE)  
}
max(ACC_y)
png(file="Accuracy_prediction_Zhu.png")
plot(ACC_y,pch=20)
dev.off()
dim(MATRIX_TS)
length(y_TS)
