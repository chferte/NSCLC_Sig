### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 6th 2012

### script for creating objects derivating from the normalization files

### script for describing the NSCLC datasets
options(stringsAsFactors=FALSE)

#load affy & survival
library(affy)
library(survival)
library(metaGEO)
library(corpcor)

# Set the normalization method ("RMA","GCRMA","fRMA", "MAS5", "barcode", "dCHIP", "metaGEO")
for(NORM in c("RMA","GCRMA", "MAS5", "dCHIP", "metaGEO")) 
    {
PATH <- "/Volumes"
Folder1  <- "/cferte/FELLOW/cferte/NSCLC_MA/"

# set the expression file directory
setwd(paste(PATH,Folder1,"Dir_NORMALIZED/",sep=""))
load(paste("Dir_",NORM,".Rdata",sep=""))
setwd(paste(PATH,Folder1,"Zhu_NORMALIZED/",sep=""))
load(paste("Zhu_",NORM,".Rdata",sep=""))
setwd(paste(PATH,Folder1,"Gir_NORMALIZED/",sep=""))
load(paste("Gir_",NORM,".Rdata",sep=""))
setwd(paste(PATH,Folder1,"Hou_NORMALIZED/",sep=""))
load(paste("Hou_",NORM,".Rdata",sep=""))
setwd(paste(PATH,Folder1,"LUSC_NORMALIZED/",sep=""))
load(paste("LUSC_",NORM,".Rdata",sep=""))

ifelse(NORM=="metaGEO",assign("DirExp",exprs(get(paste("Dir_",NORM,sep=""))$hgu133a)),assign("DirExp",exprs(get(paste("Dir_",NORM,sep="")))))
ifelse(NORM=="metaGEO",assign("ZhuExp",exprs(get(paste("Zhu_",NORM,sep=""))$hgu133a)),assign("ZhuExp",exprs(get(paste("Zhu_",NORM,sep="")))))
ifelse(NORM=="metaGEO",assign("HouExp",exprs(get(paste("Hou_",NORM,sep=""))$hgu133plus2)),assign("HouExp",exprs(get(paste("Hou_",NORM,sep="")))))
ifelse(NORM=="metaGEO",assign("GirExp",exprs(get(paste("Gir_",NORM,sep=""))$hgu133plus2)),assign("GirExp",exprs(get(paste("Gir_",NORM,sep="")))))
ifelse(NORM=="metaGEO",assign("LuscExp",exprs(get(paste("LUSC_",NORM,sep=""))$hthgu133a)),assign("LuscExp",exprs(get(paste("LUSC_",NORM,sep="")))))

dim(DirExp)
dim(ZhuExp)
dim(HouExp)
dim(GirExp)
dim(LuscExp)

# set the CLIN_DATA files directory
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/",sep=""))

#read the files
DirClin <- read.csv("DirectorData.csv")
ZhuClin <- read.csv("ZhuData.csv")
GirClin <- read.csv("GirardData.csv")
HouClin <- read.csv("HouData.csv")
LuscClin <- read.csv("LUSCData.csv")

rownames(DirClin) <- DirClin$CEL_ID
rownames(ZhuClin) <- ZhuClin$CEL_ID
rownames(GirClin) <- GirClin$CEL_ID
rownames(HouClin) <- HouClin$CEL_ID
rownames(LuscClin) <- LuscClin$CEL_ID

## change ALive -> 0 & Dead -> 1 in DirClin ZhuClin
DirClin$VITAL_STATUS<- ifelse(DirClin$VITAL_STATUS=="Dead",1,0)
ZhuClin$VITAL_STATUS<- ifelse(ZhuClin$VITAL_STATUS=="Dead",1,0)
HouClin$VITAL_STATUS<- ifelse(HouClin$VITAL_STATUS=="Dead",1,0)
GirClin$VITAL_STATUS<- ifelse(GirClin$VITAL_STATUS=="Dead",1,0)
LuscClin$VITAL_STATUS <- ifelse(LuscClin$VITAL_STATUS=="Dead",1,0)

# transform into numeric:
HouClin$MONTHS_TO_LAST_CONTACT_OR_DEATH <-as.numeric(sub(",",".",HouClin$MONTHS_TO_LAST_CONTACT_OR_DEATH))
HouClin$Age <-as.numeric(sub(",",".",HouClin$Age))

# create a new variable named Three years OS (probability of overall survival)
DirClin$THREE_YEAR_OS <- ifelse(DirClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & DirClin$VITAL_STATUS==0,1,0)
ZhuClin$THREE_YEAR_OS <- ifelse(ZhuClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & ZhuClin$VITAL_STATUS==0,1,0)
GirClin$THREE_YEAR_OS <- ifelse(GirClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & GirClin$VITAL_STATUS==0,1,0)
HouClin$THREE_YEAR_OS <- ifelse(HouClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & HouClin$VITAL_STATUS==0,1,0)
LuscClin$THREE_YEAR_OS <- ifelse(LuscClin$MONTHS_TO_LAST_CONTACT_OR_DEATH>36 & LuscClin$VITAL_STATUS==0,1,0)

## in the Dir dataset: 
## rename the moffit CEL files that are differently named in the Exp file and in the Clin files
## (change all the " " and the "." into "_")
rownames(DirClin) <- sub(" ","_",rownames(DirClin))
rownames(DirClin) <- paste(rownames(DirClin),".CEL",sep="")

# two colnames of DirClin1 were inappropriately written and do not match with DirExp ones:
rownames(DirClin) <- sub("__","_",rownames(DirClin))
DirClin$CEL_ID <- rownames(DirClin)

## do it in the Dir Exp file
colnames(DirExp) <- sub("Moff.","Moff_",colnames(DirExp))

## in the Zhu Clin dataset ...some CEL_ID were inappropriately written.
rownames(ZhuClin) <-  paste(rownames(ZhuClin),".CEL",sep="")

# select in DirClin and ZhuClin only the patients without any adjuvant treatment:
DirClin1  <- DirClin[DirClin$totalAdjuvant==0,]
ZhuClin1  <- ZhuClin[ZhuClin$totalAdjuvant==0,]
GirClin1  <- GirClin[GirClin$totalAdjuvant==0,]
HouClin1  <- HouClin[HouClin$totalAdjuvant==0,]
LuscClin1  <- LuscClin[LuscClin$totalAdjuvant==0,]

## get rid of all the patients with stage 3B
DirClin1  <- DirClin1[DirClin1$P_Stage!="3B",]
ZhuClin1  <- ZhuClin1[ZhuClin1$P_Stage!="3B",]
HouClin1  <- HouClin1[HouClin1$P_Stage!="3B",]
GirClin1  <- GirClin1[GirClin1$P_Stage!="3B",]
LuscClin1  <- LuscClin1[LuscClin1$P_Stage!="3B",]

## get rid of all the patients with stage 4
DirClin1  <- DirClin1[DirClin1$P_Stage!="4",]
ZhuClin1  <- ZhuClin1[ZhuClin1$P_Stage!="4",]
HouClin1  <- HouClin1[HouClin1$P_Stage!="4",]
GirClin1  <- GirClin1[GirClin1$P_Stage!="4",]
LuscClin1  <- LuscClin1[LuscClin1$P_Stage!="4",]

####################################################################
### match the columns of Exp files with the rows of the Clin1 files
####################################################################

#### match the DIR file ####

# create DirExp1 where are only the samples that match with the samles in DirClin1
tmp1 <- intersect(rownames(DirClin1),colnames(DirExp))
DirExp1 <- DirExp[,tmp1]
DirClin1 <- DirClin1[tmp1,]
rm(tmp1)
#check if the names are the same between DirClin1 and DirExp1
table(rownames(DirClin1)==colnames(DirExp1))


#### match the ZHU file ####

# create ZhuExp1 where are only the samples that match with the samles in ZhuClin1
tmp1 <- intersect(rownames(ZhuClin1),colnames(ZhuExp))
ZhuExp1 <- ZhuExp[,tmp1]
ZhuClin1 <- ZhuClin1[tmp1,]
rm(tmp1)
#check if the names are the same between ZhuClin1 and ZhuExp1
table(rownames(ZhuClin1)==colnames(ZhuExp1))


#### match the HOU file ####

## check the Hou file: if names are the same between Exp and Clin file
rownames(HouClin1) <- paste(rownames(HouClin1),".CEL.gz",sep="")
tmp1 <- intersect(rownames(HouClin1),colnames(HouExp))
HouExp1 <- HouExp[,tmp1]
HouClin1 <- HouClin1[tmp1,] 
rm(tmp1)
# check if the names are the same between ZhuClin1 and ZhuExp1
table(rownames(HouClin1)==colnames(HouExp1))

#### match the GIR file ####

## check the Gir file: if names are the same between Exp and Clin file
rownames(GirClin1) <- paste(rownames(GirClin1),".CEL.gz",sep="")
tmp1 <- intersect(rownames(GirClin1),colnames(GirExp))
GirExp1 <- GirExp[,tmp1]
GirClin1 <- GirClin1[tmp1,]
rm(tmp1)
# check if the names are the same between ZhuClin1 and ZhuExp1
table(rownames(GirClin1)==colnames(GirExp1))

#### match the LUSC file ####
tmp1 <- intersect(rownames(LuscClin1), colnames(LuscExp))
LuscExp1 <- LuscExp[,tmp1]
LuscClin1 <- LuscClin1[tmp1,]
rm(tmp1)
table(rownames(LuscClin1)==colnames(LuscExp1))


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

tmp1 <- match(rownames(DirExp1),rownames(LuscExp1))
LuscExp1 <- LuscExp1[tmp1,]



##############################################################################################
######     deleting all the control probles (AFF...)
##############################################################################################
ifelse(NORM=="metaGEO",DirExp1,assign("DirExp1",DirExp1[-grep("AFFX", rownames(DirExp1)),]))
ifelse(NORM=="metaGEO",ZhuExp1,assign("ZhuExp1",ZhuExp1[-grep("AFFX", rownames(ZhuExp1)),]))
ifelse(NORM=="metaGEO",HouExp1,assign("HouExp1",HouExp1[-grep("AFFX", rownames(HouExp1)),]))
ifelse(NORM=="metaGEO",GirExp1,assign("GirExp1",GirExp1[-grep("AFFX", rownames(GirExp1)),]))
ifelse(NORM=="metaGEO",LuscExp1,assign("LuscExp1",LuscExp1[-grep("AFFX", rownames(LuscExp1)),]))

# check if the rownames are the same between DirExp1 HouExp1, LuscExp1, ZhuExp1 and GirExp1
identical(rownames(DirExp1),rownames(LuscExp1))
identical(rownames(GirExp1),rownames(LuscExp1))
identical(rownames(ZhuExp1),rownames(LuscExp1))
identical(rownames(HouExp1),rownames(LuscExp1))

##############################################################################################
#####      Combining the data Clin and Exp in each Dataset
##############################################################################################

DirClinF <- DirClin1
ZhuClinF <- ZhuClin1
HouClinF <- HouClin1
GirClinF <- GirClin1
LuscClinF <- LuscClin1

DirClinF$SITE.1 <- NULL
HouClinF$SITE <- NA
HouClinF$WARNING <- NA
tmp1 <- match(colnames(DirClinF),colnames(HouClinF))
HouClinF <- HouClinF[,tmp1]
tmp1 <- match(colnames(DirClinF),colnames(LuscClinF))
LuscClinF <- LuscClinF[,tmp1]

# check if the colnames are the same between DirExp1 HouExp1, LuscExp1, ZhuExp1 and GirExp1
identical(colnames(DirClinF),colnames(GirClinF))
identical(colnames(GirClinF),colnames(ZhuClinF))
identical(colnames(ZhuClinF),colnames(HouClinF))
identical(colnames(HouClinF),colnames(LuscClinF))

#create a dataframe with only the CEL_ID, the OS, vital status and the 3YOS
DirClin1 <- DirClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
ZhuClin1 <- ZhuClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
GirClin1 <- GirClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
HouClin1 <- HouClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
LuscClin1 <- LuscClin1[,c("THREE_YEAR_OS","MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]

#combine the DirExp1 and the DirClin1
MATRIX_TS <- DirExp1
y_TS <- DirClin1[,c("THREE_YEAR_OS")]
y_OS_TS <- DirClin1[,c("MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
colnames(y_OS_TS) <-c("time","status")
dim(MATRIX_TS)
length(y_TS)
dim(y_OS_TS)

##combine the ZhuExp1 and the ZhuClin1
MATRIX_VS <- ZhuExp1
y_VS <- ZhuClin1[,c("THREE_YEAR_OS")]
y_OS_VS <- ZhuClin1[,c("MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
colnames(y_OS_VS) <-c("time","status")
dim(MATRIX_VS)
length(y_VS)
dim(y_OS_VS)

##combine the HouExp1 and the HouClin1
MATRIX_VS2 <- HouExp1
y_VS2 <- HouClin1[,c("THREE_YEAR_OS")]
y_OS_VS2 <- HouClin1[,c("MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
colnames(y_OS_VS2) <-c("time","status")
dim(MATRIX_VS2)
length(y_VS2)
dim(y_OS_VS2)

##combine the GirExp1 and the GirClin1
MATRIX_VS3 <- GirExp1
y_VS3 <- GirClin1[,c("THREE_YEAR_OS")]
y_OS_VS3 <- GirClin1[,c("MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
colnames(y_OS_VS3) <-c("time","status")
dim(MATRIX_VS3)
length(y_VS3)
dim(y_OS_VS3)

##combine the GirExp1 and the GirClin1
MATRIX_VS4 <- LuscExp1
y_VS4 <- LuscClin1[,c("THREE_YEAR_OS")]
y_OS_VS4 <- LuscClin1[,c("MONTHS_TO_LAST_CONTACT_OR_DEATH","VITAL_STATUS")]
colnames(y_OS_VS4) <-c("time","status")
dim(MATRIX_VS4)
length(y_VS4)
dim(y_OS_VS4)



############################################################################################################################
## ### create objects and save them in the MATRIX_RES_OBJECTS directory:
############################################################################################################################
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/",NORM,sep=""))

save(MATRIX_TS,file="MATRIX_TS.Rdata")
save(MATRIX_VS,file="MATRIX_VS.Rdata")
save(MATRIX_VS2,file="MATRIX_VS2.Rdata")
save(MATRIX_VS3,file="MATRIX_VS3.Rdata")
save(MATRIX_VS4,file="MATRIX_VS4.Rdata")

save(y_TS,file="y_TS.Rdata")
save(y_VS,file="y_VS.Rdata")
save(y_VS2,file="y_VS2.Rdata")
save(y_VS3,file="y_VS3.Rdata")
save(y_VS4,file="y_VS4.Rdata")

save(y_OS_TS,file="y_OS_TS.Rdata")
save(y_OS_VS,file="y_OS_VS.Rdata")
save(y_OS_VS2,file="y_OS_VS2.Rdata")
save(y_OS_VS3,file="y_OS_VS3.Rdata")
save(y_OS_VS4,file="y_OS_VS4.Rdata")

folder2  <- "/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/"
setwd(paste(PATH,folder2,sep=""))

save(DirClinF,file="TS_CLIN.Rdata")
save(ZhuClinF,file="VS_CLIN.Rdata")
save(HouClinF,file="VS2_CLIN.Rdata")
save(GirClinF,file="VS3_CLIN.Rdata")
save(LuscClinF,file="VS4_CLIN.Rdata")

print(NORM)
#
}

#####
