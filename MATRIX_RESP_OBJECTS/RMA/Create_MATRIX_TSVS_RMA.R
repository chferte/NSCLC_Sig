### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

### script for creating objects derivating from the normalization files

#load affy & survival
library(affy)
library(survival)

# Set the normalization method ("RMA","GCRMA","fRMA", "MAS5", "barcode", "dCHIP", "metaGEO")
NORM <- "RMA"
PATH <- "/home"
Folder1  <- "/cferte/FELLOW/cferte/NSCLC_MA/"

# set the expression file directory
setwd(paste(PATH,Folder1,"DIR_NORMALIZED/",sep=""))
load(paste("Dir_",NORM,".Rdata",sep=""))
setwd(paste(PATH,Folder1,"ZHU_NORMALIZED/",sep=""))
load(paste("Zhu_",NORM,".Rdata",sep=""))
setwd(paste(PATH,Folder1,"GIR_NORMALIZED/",sep=""))
load(paste("Gir_",NORM,".Rdata",sep=""))
setwd(paste(PATH,Folder1,"HOU_NORMALIZED/",sep=""))
load(paste("Hou_",NORM,".Rdata",sep=""))


DirExp <- exprs(Dir_RMA)
ZhuExp <- exprs(Zhu_RMA)
HouExp <- exprs(Hou_RMA)
GirExp <- exprs(Gir_RMA)

# set Dir CEL files directory
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/",sep=""))

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




# transpose the clin matrix to make the sample names = colnames (ie==expression files)
DirClin1 <- t(DirClin1)
ZhuClin1 <- t(ZhuClin1)
GirClin1 <- t(GirClin1)
HouClin1 <- t(HouClin1)





####################################################################
### match the columns of Exp files with the rows of the Clin1 files
####################################################################

#### match the DIR file ####

# two colnames of DirClin1 were inappropriately written and do not match with DirExp ones:
colnames(DirClin1) <- sub("__","_",colnames(DirClin1)) 

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



##############################################################################################
######     deleting all the control probles (AFF...)
##############################################################################################
DirExp1 <- DirExp1[-grep("AFFX", rownames(DirExp1)),]
ZhuExp1 <- ZhuExp1[-grep("AFFX", rownames(ZhuExp1)),]
HouExp1 <- HouExp1[-grep("AFFX", rownames(HouExp1)),]
GirExp1 <- GirExp1[-grep("AFFX", rownames(GirExp1)),]

# check if the names are the same between DirExp1 and GirExp1
table(rownames(DirExp1)==rownames(GirExp1))
table(rownames(ZhuExp1)==rownames(GirExp1))
table(rownames(HouExp1)==rownames(GirExp1))

##############################################################################################
#####      Combining the data Clin and Exp in each Dataset
##############################################################################################

DirClinF <- t(DirClin1)
ZhuClinF <- t(ZhuClin1)
HouClinF <- t(HouClin1)
GirClinF <- t(GirClin1)


DirClinF$SITE.1 <- NULL
HouClinF$SITE <- NA
HouClinF$WARNING <- NA
tmp1 <- match(colnames(DirClinF),colnames(HouClinF))
HouClinF <- HouClinF[,tmp1]

#create a matrix with only the CEL_ID and the 3YOS
DirClin1 <- DirClin1[c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS"),]
ZhuClin1 <- ZhuClin1[c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS"),]
GirClin1 <- GirClin1[c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS"),]
HouClin1 <- HouClin1[c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS"),]


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

############################################################################################################################
## ### create objects and save them in the MATRIX_RES_OBJECTS directory:
############################################################################################################################
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/",NORM,sep=""))

save(MATRIX_TS,file="MATRIX_TS.Rdata")
save(MATRIX_VS,file="MATRIX_VS.Rdata")
save(MATRIX_VS2,file="MATRIX_VS2.Rdata")
save(MATRIX_VS3,file="MATRIX_VS3.Rdata")

save(y_TS,file="y_TS.Rdata")
save(y_VS,file="y_VS.Rdata")
save(y_VS2,file="y_VS2.Rdata")
save(y_VS3,file="y_VS3.Rdata")

save(y_OS_TS,file="y_OS_TS.Rdata")
save(y_OS_VS,file="y_OS_VS.Rdata")
save(y_OS_VS2,file="y_OS_VS2.Rdata")
save(y_OS_VS3,file="y_OS_VS3.Rdata")

folder2  <- "/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/"
setwd(paste(PATH,folder2,sep=""))

save(DirClinF,file="TS_CLIN.Rdata")
save(ZhuClinF,file="VS_CLIN.Rdata")
save(HouClinF,file="VS2_CLIN.Rdata")
save(GirClinF,file="VS3_CLIN.Rdata")
