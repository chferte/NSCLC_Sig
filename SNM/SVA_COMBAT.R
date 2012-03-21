### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

### script for running combat in a file:
library(sva)
library(pamr)
library(limma)

## Obtain phenotypic data
PATH <- "/home"
method <-  "RMA"
folder1 <- "/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/"
setwd(paste(PATH,folder1,method,sep=""))

load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")
load("MATRIX_VS3.Rdata")

TOTAL <- cbind(MATRIX_TS,MATRIX_VS)
TOTAL <- cbind(TOTAL,MATRIX_VS2)
TOTAL <- cbind(TOTAL,MATRIX_VS3)
TOTAL <- apply(TOTAL,2,as.numeric)

folder2  <- "/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/"
setwd(paste(PATH,folder2,sep=""))
load("TS_CLIN.Rdata")
load("VS_CLIN.Rdata")
load("VS2_CLIN.Rdata")
load("VS3_CLIN.Rdata")

TCLIN <- rbind(DirClinF,ZhuClinF)
TCLIN <- rbind(TCLIN,GirClinF)
TCLIN <- rbind(TCLIN,HouClinF)

tmp1 <- match(colnames(TOTAL),rownames(TCLIN))
TCLIN <- TCLIN[tmp1,]

## check if colnames(TOTAL)=rownames(TCLIN)
all(rownames(TCLIN)==colnames(TOTAL))

TCLIN <- as.data.frame(TCLIN)
TCLIN$Study <- substr(TCLIN$LABORATORY_BATCH,1,3)

# estimate the number of latent factors
mod <- model.matrix(~as.factor(Study),data=TCLIN)
mod0 <- model.matrix(~1,data=TCLIN)

n.sv = num.sv(TOTAL,mod)
n.sv

# estimate the surrogate variables
svobj = sva(TOTAL,mod,mod0,n.sv=n.sv)


## adjust for surrogate variable using the f.pvalue function
# estimate the differntially expressed probes without adjusting to surrogate variables.
pValues = f.pvalue(TOTAL,mod,mod0)
qValues = p.adjust(pValues,method="BH")
# ...and adjusting to surrogate variables:
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(TOTAL,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

## Adjusting for surrogate variables using the limma package

## Applying the ComBat function to adjust for known batches= study ! variable of interest = 3 year OS
batch <- TCLIN$Study
mod <- model.matrix(~as.factor(THREE_YEAR_OS),data=TCLIN)
TOTAL_CB = ComBat(dat=TOTAL, batch=batch, mod=mod, numCovs=NULL, par.prior=FALSE,prior.plot=FALSE)
pValuesCombat = f.pvalue(TOTAL_CB,mod,mod0)
qValuesCombat = p.adjust(pValuesCombat,method="BH")
table(qValuesCombat<.05)
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/")
save(TOTAL_CB, file="TOTAL_CB.Rdata")
save(TOTAL, file="TOTAL.Rdata")
save(TCLIN, file="TCLIN.Rdata")


## separate the matrix according to the different studies again in order to perform then the analysis on it:

tmp1 <- match(rownames(DirClinF),colnames(TOTAL_CB))
MATRIX_TS_CB <- TOTAL_CB[,tmp1]
dim(MATRIX_TS_CB)

tmp1 <- match(rownames(ZhuClinF),colnames(TOTAL_CB))
MATRIX_VS_CB <- TOTAL_CB[,tmp1]
dim(MATRIX_VS_CB)

tmp1 <- match(rownames(HouClinF),colnames(TOTAL_CB))
MATRIX_VS2_CB <- TOTAL_CB[,tmp1]
dim(MATRIX_VS2_CB)

tmp1 <- match(rownames(GirClinF),colnames(TOTAL_CB))
MATRIX_VS3_CB <- TOTAL_CB[,tmp1]
dim(MATRIX_VS3_CB)


save(MATRIX_TS_CB,file="MATRIX_TS_CB.Rdata")
save(MATRIX_VS_CB,file="MATRIX_VS_CB.Rdata")
save(MATRIX_VS2_CB,file="MATRIX_VS2_CB.Rdata")
save(MATRIX_VS3_CB,file="MATRIX_VS3_CB.Rdata")