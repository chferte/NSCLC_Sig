### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### February, 19th 2012

### script for evaluating the correlation of HR between the different dataset

#
library(survival)

#point the data
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/RMA/")
load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")
load("MATRIX_VS4.Rdata")

load("y_TS.Rdata")
load("y_VS.Rdata")
load("y_VS2.Rdata")
load("y_VS4.Rdata")
load("y_OS_TS.Rdata")
load("y_OS_VS.Rdata")
load("y_OS_VS2.Rdata")
load("y_OS_VS4.Rdata")

############################################################################################################################
## ### univariate anbalysis
############################################################################################################################

library(survival)

y_OS_TS <- as.data.frame(y_OS_TS)
y_OS_TS[,1] <- as.numeric(y_OS_TS[,1])
y_OS_TS[,2] <- as.numeric(y_OS_TS[,2])

y_OS_VS <- as.data.frame(y_OS_VS)
y_OS_VS[,1] <- as.numeric(y_OS_VS[,1])
y_OS_VS[,2] <- as.numeric(y_OS_VS[,2])

y_OS_VS2 <- as.data.frame(y_OS_VS2)
y_OS_VS2[,1] <- as.numeric(y_OS_VS2[,1])
y_OS_VS2[,2] <- as.numeric(y_OS_VS2[,2])

y_OS_VS4 <- as.data.frame(y_OS_VS4)
y_OS_VS4[,1] <- as.numeric(y_OS_VS4[,1])
y_OS_VS4[,2] <- as.numeric(y_OS_VS4[,2])


MATRIX_TS <- apply(MATRIX_TS,2,as.numeric)
MATRIX_VS <- apply(MATRIX_VS,2,as.numeric)
MATRIX_VS2 <- apply(MATRIX_VS2,2,as.numeric)
MATRIX_VS4 <- apply(MATRIX_VS4,2,as.numeric)

myCoxHR <- function(x){
  summary(coxph(mySurv ~ as.numeric(x)))$coefficients[2]
}

myCoxP <- function(x){
  summary(coxph(mySurv ~ as.numeric(x)))$coefficients[5]
}

mySurv <-Surv(y_OS_TS[,1],y_OS_TS[,2])
HR_TS <- apply(MATRIX_TS[1:5000,], 1, myCoxHR)

mySurv <-Surv(y_OS_VS[,1],y_OS_VS[,2])
HR_VS <- apply(MATRIX_VS[1:5000,], 1, myCoxHR)

mySurv <-Surv(y_OS_VS2[,1],y_OS_VS2[,2])
HR_VS2 <- apply(MATRIX_VS2[1:5000,], 1, myCoxHR)

mySurv <-Surv(y_OS_VS4[,1],y_OS_VS4[,2])
HR_VS4 <- apply(MATRIX_VS4[1:5000,], 1, myCoxHR)

##############################################
mySurv <-Surv(y_OS_TS[,1],y_OS_TS[,2])
P_TS <- apply(MATRIX_TS, 1, myCoxP)
P_TS <- p.adjust(P_TS, method = "BH", n = length(P_TS))

mySurv <-Surv(y_OS_VS[,1],y_OS_VS[,2])
P_VS <- apply(MATRIX_VS, 1, myCoxP)
P_VS <- p.adjust(P_VS, method = "BH", n = length(P_VS))

mySurv <-Surv(y_OS_VS2[,1],y_OS_VS2[,2])
P_VS2 <- apply(MATRIX_VS2, 1, myCoxP)
P_VS2 <- p.adjust(P_VS2, method = "BH", n = length(P_VS2))

mySurv <-Surv(y_OS_VS4[,1],y_OS_VS4[,2])
P_VS4 <- apply(MATRIX_VS4, 1, myCoxP)
P_VS4 <- p.adjust(P_VS4, method = "BH", n = length(P_VS4))

png("density_probes_univariate_BH.png")
par(mfrow=c(2,2))
hist(P_TS,main="TS",xlim=c(0,1),breaks=100)
abline(v=.05,lty=2,col="red")
hist(P_VS,main="VS",xlim=c(0,1),breaks=100)
abline(v=.05,lty=2,col="red")
hist(P_VS2,main="VS2",xlim=c(0,1),breaks=100)
abline(v=.05,lty=2,col="red")
hist(P_VS4,main="VS4",xlim=c(0,1),breaks=100)
abline(v=.05,lty=2,col="red")
title(main= "density of the probe sets associated with survival across datasets", outer=T)
dev.off()

setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/SNM/HR_CORRELATIONS/")

png("HR_TS_VS.png")
par(mfrow=c(2,3))

plot(log(HR_TS),log(HR_VS),xlab="log(HR) TS",ylab="log(HR) VS")
abline(h=0,col="Blue",lty=2)
abline(v=0,col="Blue",lty=2)

plot(log(HR_TS),log(HR_VS2),xlab="log(HR) TS",ylab="log(HR) VS2")
abline(h=0,col="Blue",lty=2)
abline(v=0,col="Blue",lty=2)

plot(log(HR_TS),log(HR_VS4),xlab="log(HR) TS",ylab="log(HR) VS4")
abline(h=0,col="Blue",lty=2)
abline(v=0,col="Blue",lty=2)

plot(log(HR_VS),log(HR_VS2),xlab="log(HR) VS",ylab="log(HR) VS2")
abline(h=0,col="Blue",lty=2)
abline(v=0,col="Blue",lty=2)

plot(log(HR_VS),log(HR_VS4),xlab="log(HR) VS",ylab="log(HR) VS4")
abline(h=0,col="Blue",lty=2)
abline(v=0,col="Blue",lty=2)

plot(log(HR_VS2),log(HR_VS4),xlab="log(HR) VS2",ylab="log(HR) VS4")
abline(h=0,col="Blue",lty=2)
abline(v=0,col="Blue",lty=2)

title(main="distribution of the Hazard ratios across datasets", outer=T)
dev.off()

##########
# set the working directory for plotting the COR COR
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/SNM/CORCOR/")

## charles' method
#VS <- MATRIX_VS
#TS <- MATRIX_TS
#CORCOR <- c()
#for(k in 1:50){
# a <- 100*(k-1)+1
#  b <- k*100
#COR_TS <- c()
#COR_VS <- c()
#  for(i in a:b){
#COR_TS <- rbind(COR_TS,cor(TS[i,],t(TS)))
#COR_VS <- rbind(COR_VS,cor(VS[i,],t(VS)))

#}
#for(r in 1:100) {CORCOR <- c(CORCOR,cor(COR_VS[r,],COR_TS[r,]))}
#print(k)
#  }
  

## altrenative quicker and "all probe sets" method (Charles & In Sock)
TS <- MATRIX_TS
VS <- MATRIX_VS
COROICO <- c()
CORFUNTS <- function(x){cor(x,t(TS))}
CORFUNVS <- function(y){cor(y,t(VS))}
for(i in 1:22215){
  a <- as.vector(CORFUNTS(TS[i,]))
  b <- as.vector(CORFUNVS(VS[i,]))
  COROICO <- c(COROICO,cor(a,b))
  print(i)
}
png("CORCOR_22K_TS_VS.png")
plot(COROICO,xlab="probes 1:22215",ylab="correlation(COR_TS*COR_VS)",ylim=c(-1,1))
dev.off()
save(COROICO,file="CORCOR_22K_TS_VS.Rdata")
#

TS <- MATRIX_TS
VS2 <- MATRIX_VS2
COROICO <- c()
CORFUNTS <- function(x){cor(x,t(TS))}
CORFUNVS2 <- function(y){cor(y,t(VS2))}
for(i in 1:22215){
  a <- as.vector(CORFUNTS(TS[i,]))
  b <- as.vector(CORFUNVS2(VS2[i,]))
  COROICO <- c(COROICO,cor(a,b))
  print(i)
}
png("CORCOR_22K_TS_VS2.png")
plot(COROICO,xlab="probes 1:22215",ylab="correlation(COR_TS*COR_VS2)",ylim=c(-1,1))
dev.off()
save(COROICO,file="CORCOR_22K_TS_VS2.Rdata")
#


VS <- MATRIX_VS
VS2 <- MATRIX_VS2
COROICO <- c()
CORFUNVS <- function(x){cor(x,t(VS))}
CORFUNVS2 <- function(y){cor(y,t(VS2))}
for(i in 1:22215){
  a <- as.vector(CORFUNVS(VS[i,]))
  b <- as.vector(CORFUNVS2(VS2[i,]))
  COROICO <- c(COROICO,cor(a,b))
  print(i)
}
png("CORCOR_22K_VS_VS2.png")
plot(COROICO,xlab="probes 1:22215",ylab="correlation(COR_VS*COR_VS2)",ylim=c(-1,1))
dev.off()
#


VS <- MATRIX_VS
VS4 <- MATRIX_VS4
COROICO <- c()
CORFUNVS <- function(x){cor(x,t(VS))}
CORFUNVS4 <- function(y){cor(y,t(VS4))}
for(i in 1:22215){
  a <- as.vector(CORFUNVS(VS[i,]))
  b <- as.vector(CORFUNVS4(VS4[i,]))
  COROICO <- c(COROICO,cor(a,b))
  print(i)
}
png("CORCOR_22K_VS_VS4.png")
plot(COROICO,xlab="probes 1:22215",ylab="correlation(COR_VS*COR_VS4)",ylim=c(-1,1))
dev.off()
save(COROICO,file="CORCOR_22K_VS_VS4.Rdata")
#


VS2 <- MATRIX_VS2
VS4 <- MATRIX_VS4
COROICO <- c()
CORFUNVS2 <- function(x){cor(x,t(VS2))}
CORFUNVS4 <- function(y){cor(y,t(VS4))}
for(i in 1:22215){
  a <- as.vector(CORFUNVS2(VS2[i,]))
  b <- as.vector(CORFUNVS4(VS4[i,]))
  COROICO <- c(COROICO,cor(a,b))
  print(i)
}
png("CORCOR_22K_VS2_VS4.png")
plot(COROICO,xlab="probes 1:22215",ylab="correlation(COR_VS2*COR_VS4)",ylim=c(-1,1))
dev.off()
save(COROICO,file="CORCOR_22K_VS2_VS4.Rdata")
#


TS <- MATRIX_TS
VS4 <- MATRIX_VS4
COROICO <- c()
CORFUNTS <- function(x){cor(x,t(TS))}
CORFUNVS4 <- function(y){cor(y,t(VS4))}
for(i in 1:22215){
  a <- as.vector(CORFUNTS(TS[i,]))
  b <- as.vector(CORFUNVS4(VS4[i,]))
  COROICO <- c(COROICO,cor(a,b))
  print(i)
}
png("CORCOR_22K_TS_VS4.png")
plot(COROICO,xlab="probes 1:22215",ylab="correlation(COR_TS*COR_VS4)",ylim=c(-1,1))
dev.off()
save(COROICO,file="CORCOR_22K_TS_VS4.Rdata")
##
##
##




## plot the densityplots
library(lattice)
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/SNM/CORCOR/")

load("CORCOR_22K_TS_VS.Rdata")
TS_VS <- COROICO
load("CORCOR_22K_TS_VS2.Rdata")
TS_VS2 <- COROICO
load("CORCOR_22K_TS_VS4.Rdata")
TS_VS4 <- COROICO
load("CORCOR_22K_VS_VS2.Rdata")
VS_VS2 <- COROICO
load("CORCOR_22K_VS_VS4.Rdata")
VS_VS4 <- COROICO
load("CORCOR_22K_VS2_VS4.Rdata")
VS2_VS4 <- COROICO
rm(COROICO)

obj <- data.frame( cors = c(TS_VS,TS_VS2,TS_VS4,VS_VS2,VS_VS4,VS2_VS4),
  			 studies = rep(c("Dir_Zhu","Dir_Hou","Dir_TCGA","Zhu_Hou","Hou_TCGA","Hou_TCGA"),each=22215))

png("dens_plot_ALLCORS.png")
densityplot(obj$cors, groups=obj$studies, auto.key=list(space="right"), plot.points=FALSE, lwd=3, xlab="Correlations")
dev.off()

## plot the first gene
png("COR_TS_VS.png")
plot(COR_TS[1,],COR_VS[1,])
dev.off()

png("COR_VS_VS2.png")
plot(COR_VS[1,],COR_VS2[1,])
dev.off()

png("COR_TS_VS2.png")
plot(COR_TS[1,],COR_VS2[1,])
dev.off()

png("COR_VS_VS4.png")
plot(COR_VS[1,],COR_VS4[1,])
dev.off()

png("COR_TS_VS4.png")
plot(COR_TS[1,],COR_VS4[1,])
dev.off()

png("COR_VS2_VS4.png")
plot(COR_VS2[1,],COR_VS4[1,])
dev.off()
#
library(lattice)
densityplot(cor(COR_VS2[1,],COR_VS4[1,]))
