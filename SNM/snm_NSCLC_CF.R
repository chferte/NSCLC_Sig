### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### March, 2 2012

#load the libraries
library(metaGEO)
library(snm)
library(lattice)


## load the data:
method <- "RMA"
PATH <- "/Volumes"
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/",method,sep=""))

#read the files
load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")
load("MATRIX_VS3.Rdata")
load("MATRIX_VS4.Rdata")

# read the clin files
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/",sep=""))
load("TS_CLIN.Rdata")
load("VS_CLIN.Rdata")
load("VS2_CLIN.Rdata")
load("VS3_CLIN.Rdata")
load("VS4_CLIN.Rdata")

## transform the H_Grdae in VS4_CLIN
LuscClinF$P_Stage <- substr(LuscClinF$P_Stage,7,nchar(LuscClinF$P_Stage))
LuscClinF$P_Stage <- sub("III","3",LuscClinF$P_Stage)
LuscClinF$P_Stage <- sub("II","2",LuscClinF$P_Stage)
LuscClinF$P_Stage <- sub("I","1",LuscClinF$P_Stage)

## process SVDs to get the data
Y <- MATRIX_VS
CLIN <- ZhuClinF
s <- fs(Y)

par(mfrow=c(3,3))

# histology
CB <- CLIN$Histology
CBN <- "Histology"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))

# pathological stage
CB <- CLIN$P_Stage
CBN <- "p-Stage"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))

# Gender
CB <- CLIN$GENDER
CBN <- "Gender"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))

# Grade
CB <- CLIN$H_Grade
CBN <- "Grade"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))

# SMOKING
CB <- CLIN$SMOKING
CBN <- "SMOKING"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))

# SITE
CB <- CLIN$SITE
CBN <- "SITE"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))
xyplot(s$v[,2] ~ s$v[,1], groups=CLIN$SITE)


# LABORATORY BATCH
CB <- CLIN$LABORATORY_BATCH
CBN <- "Laboratory Batch (Date)"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))

# VITAL STATUS
CB <- CLIN$VITAL_STATUS
CBN <- "Events (Deaths)"
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN,sep=""))




## survival in alive patients
CB <- CLIN$MONTHS_TO_LAST_CONTACT_OR_DEATH[CLIN$VITAL_STATUS==0]
CBN <- "survival in alive patients"
plot(density(CB), xlab="survival in alivae patients (months)", main="density of the estimated follow-up of the patients", col= "black", lwd=2)


### implement brig's snm method to my data

bio.var <- DirClinF[,c("VITAL_STATUS","H_Grade")]
bio.var$H_Grade <- as.factor(bio.var$H_Grade)
adj.var <- DirClinF[,c("SITE","LABORATORY_BATCH")]
int.var <- data.frame(array=factor(1:dim(Y)[2]))
options("contrasts")
X <- model.matrix(~bio.var$H_Grade + bio.var$VITAL_STATUS,bio.var)
Z <- model.matrix(~adj.var$SITE,adj.var)
int.var <- int.var


u <- fs(Y)


blah <- snm(raw.dat=Y, bio.var=X, adj.var=Z, int.var=int.var)
ks.test(snm.fit$pval[true.nulls], "punif")
hist(snm.fit$pval[true.nulls])





SITE.effect <- snm:::sim.probe.specific(data, adj.var$batches, 0.3, list(func=rnorm,params=c(mean=0,sd=0.3)))
age.effect <- snm:::sim.probe.specific(data, adj.var$age, 0.2, list(func=rnorm, params=c(mean=1,sd=0.1)))
M <- data + group.effect + batches.effect + age.effect
array.effect <- snm:::sim.intensity.dep(M, int.var$array, 2, list(func=rnorm, params=c(mean=0,sd=2)))
E <- matrix(rnorm(length(data),0,0.25), nr=nrow(data), nc=ncol(data))
Y <- M + array.effect + E
true.nulls <- which(group.effect[,1] == group.effect[,26])
plot2 <- function(x,y, ...) {
  plot(Y[,c(x,y)], ...);
  abline(0,1,col="red",lwd=3,lty=2)
}
png(file="samples_33_vs_44.png")
plot2(44,33,xlab="Sample 44", ylab="Sample 33")
dev.off()
u <- fs(Y)
X <- model.matrix(~bio.var$groups)
Z <- model.matrix(~adj.var$age + adj.var$batches)
int.var <- data.frame(array=factor(1:ncol(Y)))
snm.fit <- snm(Y, bio.var=X, adj.var=Z, int.var=int.var, rm.adj=TRUE)
ks.test(snm.fit$pval[true.nulls], "punif")
hist(snm.fit$pval[true.nulls])