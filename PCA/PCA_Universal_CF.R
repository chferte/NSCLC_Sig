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
for(method in c("RMA","GCRMA", "MAS5", "dCHIP", "metaGEO")) 
    {

  PATH <- "/home"
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/",method,sep=""))
print(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/",method,sep=""))
  
#read the files
load("MATRIX_TS.Rdata")
load("MATRIX_VS.Rdata")
load("MATRIX_VS2.Rdata")
load("MATRIX_VS3.Rdata")
load("MATRIX_VS4.Rdata")

#rename the colnames of the datests including the name of the datset
colnames(MATRIX_TS)  <- paste("DIR_",colnames(MATRIX_TS),sep="")
colnames(MATRIX_VS2)  <- paste("HOU_",colnames(MATRIX_VS2),sep="")
colnames(MATRIX_VS)  <- paste("ZHU_",colnames(MATRIX_VS),sep="")
colnames(MATRIX_VS3)  <- paste("GIR_",colnames(MATRIX_VS3),sep="")
colnames(MATRIX_VS4)  <- paste("LUS_",colnames(MATRIX_VS4),sep="")

## Dir  Zhu and Lusc are H133A whereas Hou and Gir are H133plus2. 
## We will process the PCA on the commmon probes only (H133A)

## first merge the H133A then the H133plus2

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

### Start SVD in TOTAL
## SVD scripts adapted form Brig & Brian


blah <- as.matrix(TOTAL)
study <- substr(colnames(TOTAL),1,3)
print(method)
print(table(study))
dim(TOTAL)

setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/PCA/",method,"_PCA",sep=""))

s <- fast.svd(blah - rowMeans(blah))
png(paste(file="Var_Exp_",method,".png",sep=""))
plot(s$d^2/sum(s$d^2), main=paste("Percent Variability Explained ",method,sep=""), ylab="% Variability Explained", xlab="")
dev.off()

png(file=paste("first_SVD_ ",method,".png",sep=""))
boxplot(split(s$v[, 1], study), main=paste("First SVD by batch (study)",method,sep=""))
dev.off()

png(file=paste("second_SVD_ ",method,".png",sep=""))
boxplot(split(s$v[, 2], study), main=paste("Second SVD by batch (study) ",method,sep=""))
dev.off()

png(file=paste("SVD1_SVD2_ ",method,".png",sep=""))
print(xyplot(s$v[,2] ~ s$v[,1], groups=study, xlab="SVD1
                     Dir= Blue   Hou= Red   TCGA-LUSC= Gold   Zhu= Green   Gir= Black" , ylab="SVD2", pch= 20 ,cex=1, 
             main = paste("SVD1 vs SVD2 colored by datasets normalized separately (",method,")",sep=""), col=c("royalblue","black","brown1","gold","aquamarine4")))
dev.off()

#
##  
}

#






