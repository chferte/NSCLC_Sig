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
method= "barcode"
PATH <- "/home/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/"
setwd(paste(PATH,method,sep=""))


## load the matrix and response files
load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")
load("y_TS.Rdata")
load("y_VS.Rdata")
load("y_VS2.Rdata")
load("y_OS_TS.Rdata")
load("y_OS_VS.Rdata")
load("y_OS_VS2.Rdata")


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
MATRIX_VS2_S <- normalize_to_X(mean_x,sd_x,MATRIX_VS2)

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
CI_TS <- c(CI_TS,concordance.index(y_E_TS,y_OS_TS[,1],y_OS_TS[,2],na.rm=T,method="noether")[1])
  
y_E_VS <- predict(fit,t(MATRIX_VS),type="link")
CI_VS <- c(CI_VS,concordance.index(y_E_VS,y_OS_VS[,1],y_OS_VS[,2],na.rm=T,method="noether")[1])

y_E_VS2 <- predict(fit,t(MATRIX_VS2),type="link")
CI_VS2 <- c(CI_VS2,concordance.index(y_E_VS2,y_OS_VS2[,1],y_OS_VS2[,2],na.rm=T,method="noether")[1])

y_E_TS_S <- predict(fit,x,type="link")
CI_TS_S <- c(CI_TS_S,concordance.index(y_E_TS_S,y_OS_TS[,1],y_OS_TS[,2],na.rm=T,method="noether")[1])

y_E_VS_S <- predict(fit,t(MATRIX_VS_S),type="link")
CI_VS_S <- c(CI_VS_S,concordance.index(y_E_VS_S,y_OS_VS[,1],y_OS_VS[,2],na.rm=T,method="noether")[1])

y_E_VS2_S <- predict(fit,t(MATRIX_VS2_S),type="link")
CI_VS2_S <- c(CI_VS2_S,concordance.index(y_E_VS2_S,y_OS_VS2[,1],y_OS_VS2[,2],na.rm=T,method="noether")[1])
    }
}

summary(as.numeric(CI_TS[names(CI_TS)=="c.index"][is.na(CI_TS[names(CI_TS)=="c.index"])==F]))
summary(as.numeric(CI_VS[names(CI_VS)=="c.index"][is.na(CI_VS[names(CI_VS)=="c.index"])==F]))
summary(as.numeric(CI_VS2[names(CI_VS2)=="c.index"][is.na(CI_VS2[names(CI_VS2)=="c.index"])==F]))
summary(as.numeric(CI_TS_S[names(CI_TS_S)=="c.index"][is.na(CI_TS_S[names(CI_TS_S)=="c.index"])==F]))
summary(as.numeric(CI_VS_S[names(CI_VS_S)=="c.index"][is.na(CI_VS_S[names(CI_VS_S)=="c.index"])==F]))
summary(as.numeric(CI_VS2_S[names(CI_VS2_S)=="c.index"][is.na(CI_VS2_S[names(CI_VS2_S)=="c.index"])==F]))


GRID[which(CI_TS==max(as.numeric(CI_TS[names(CI_TS)=="c.index"][is.na(CI_TS[names(CI_TS)=="c.index"])==F]))),]
GRID[which(CI_VS==max(as.numeric(CI_VS[names(CI_VS)=="c.index"][is.na(CI_VS[names(CI_VS)=="c.index"])==F]))),]
GRID[which(CI_VS2==max(as.numeric(CI_VS2[names(CI_VS2)=="c.index"][is.na(CI_VS2[names(CI_VS2)=="c.index"])==F]))),]

CI_TOTAL <- cbind(GRID$.alpha,GRID$.lambda,CI_TS,CI_VS,CI_VS2)
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
save(list=paste("CI_",method,sep=""), file=paste("CI_",method,".Rdata",sep="")
     
  