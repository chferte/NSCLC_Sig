### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### February, 15 2012

### script for plotting the eigengenes of the normalized expression according to each dataset

### script for describing the NSCLC datasets
options(stringsAsFactors=FALSE)

#load affy
library(affy)
library(corpcor)
library(lattice)

# point the directory &choose method among = RMA, GCRMA, MAS5, dCHIP, metaGEO, fRMA or barcode

# Set the normalization method ("RMA","GCRMA","fRMA", "MAS5", "barcode", "dCHIP", "metaGEO")
#for(method in c("RMA","GCRMA", "MAS5", "dCHIP", "metaGEO")) 
 #   {

  method <- "RMA"
  PATH <- "/Volumes"
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/",method,sep=""))
  
#read the files
load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")
load("MATRIX_VS3.Rdata")
load("MATRIX_VS4.Rdata")


## merge the H133A (Zhu and Dir) and call it DirZhu
tmp1 <- match(rownames(MATRIX_TS),rownames(MATRIX_VS))

TSVS <- cbind(MATRIX_TS, MATRIX_VS[tmp1,])
tmp1 <- match(rownames(TSVS),rownames(MATRIX_VS4))
TSVS <- cbind(TSVS, MATRIX_VS4[tmp1,])
  
## merge the H133plus2 (Hou and Gir) and call it HouGir
tmp2 <- match(rownames(MATRIX_VS2),rownames(MATRIX_VS3))
VS23 <- cbind(MATRIX_VS2,MATRIX_VS3[tmp2,])
  
## Since we will process the PCA on the commmon probes only (H133A), let's merge DirZhu and HouGir and LUSC and call it TOTAL
tmp3 <- rownames(VS23)[match(rownames(TSVS),rownames(VS23))]
tmp3 <- tmp3[!is.na(tmp3)]
TOTAL <- cbind(TSVS[tmp3,],VS23[tmp3,])
  
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

## merge the clin fil into TOT_CLIN
DirClinF$RN <- rownames(DirClinF)
ZhuClinF$RN <- rownames(ZhuClinF)
HouClinF$RN <- rownames(HouClinF)
LuscClinF$RN <- rownames(LuscClinF)
GirClinF$RN <- rownames(GirClinF)

TOTCLIN <- merge(DirClinF,LuscClinF,all=T)
TOTCLIN <- merge(TOTCLIN,HouClinF,all=T)
TOTCLIN <- merge(TOTCLIN,ZhuClinF,all=T)
TOTCLIN <- merge(TOTCLIN,GirClinF,all=T)
rownames(TOTCLIN) <-TOTCLIN$RN

## arrange TOTAL colnames and TOTCLIN rownames in the same order
tmp <- intersect(colnames(TOTAL),rownames(TOTCLIN))
TOTCLIN <- TOTCLIN[tmp,]
TOTAL <- TOTAL[,tmp]
identical(colnames(TOTAL),rownames(TOTCLIN))
  
### Start SVD 
data <- as.matrix(MATRIX_TS) # set the data we want to explain 
CB <- DirClinF$SMOKING  # set CB as the Clinico-Biological variable of interest
CBN <- "Smoking Status"          # name of the CB
print(method)
print(table(CB))
dim(data)


setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/PCA/VAR_ADJ",sep=""))

s <- fast.svd(data - rowMeans(data))
#png(paste(file="Var_Exp_",method,".png",sep=""))
plot(s$d^2/sum(s$d^2), main=paste("Percent Variability Explained ",method,sep=""), ylab="% Variability Explained", xlab="")
#dev.off()

#png(file=paste("first_SVD_ ",method,".png",sep=""))
boxplot(split(s$v[, 1], CB), main=paste("First SVD by ",CBN," (",method,")",sep=""))
#dev.off()

#png(file=paste("second_SVD_ ",method,".png",sep=""))
boxplot(split(s$v[, 2], CB), main=paste("Second SVD by ",CBN," (",method,")",sep=""))
#dev.off()

#png(file=paste("SVD1_SVD2_ ",method,".png",sep=""))
print(xyplot(s$v[,2] ~ s$v[,1], groups=CB, xlab="SVD1" , ylab="SVD2", pch= 20 ,cex=1, 
             main = paste("SVD1 vs SVD2 colored by :  ",CBN,sep=""), 
             auto.key=list(col=1:length(unique(CB)), points=FALSE,columns=length((unique(CB)))),col=1:length(unique(CB))))
#dev.off()

#
##  
}

#






