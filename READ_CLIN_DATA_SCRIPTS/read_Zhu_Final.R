## POINT TO DATA DIRECTORY
dataPath <- "/home/cferte/NSCLC_Prognostic_Signature/Zhu_GSE14814/series_matrix/"
dataPath2 <- "/home/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/"


## READ IN THE TAB DELIMITED FILES
dat1 <- read.delim2(file = paste(dataPath, "Zhu_2010_JCO_microarray_pts_clin_info-1.txt", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
rownames(dat1) <- dat1$arrayname
dat2 <- read.delim2(file = paste(dataPath, "GSE14814_series_matrix.txt", sep=""), na.strings=c("NA", "", " "), skip=30, header=F, nrows=2, dec=",", as.is=T)
dat2 <- dat2[, -1]
Zhu43_in_Dir <- read.csv(file = paste(dataPath2,"Zhu43_in_Dir.csv",sep=""), na.strings=c("NA", "", " "), header=T, dec=",", as.is=T)


## CREATE GSE VARIABLE
CEL_ID <- as.character((dat2[2,]))
names(CEL_ID) <- sub("patient_sample_", "", dat2[1,])

dat1$CEL_ID <- NA
dat1[names(CEL_ID), "CEL_ID"] <- CEL_ID
dat1$CEL_ID[is.na(dat1$CEL_ID)] <- dat1$arrayname[is.na(dat1$CEL_ID)]


## TRANSFORM OS.Time (years) into MONTHS_TO_LAST_CONTACT_OR_DEATH (months)
dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH <- dat1$OS.time * 12

## CREATE NEW VARIABLE NAMED VITAL_STATUS
dat1$VITAL_STATUS  <- dat1$OS.status

## create new columns "Histology"
dat1$Histology  <- dat1$Histology.type
dat1$Histology[dat1$Histology=="SQCC"] <- "SCC"
dat1$Histology[dat1$Histology=="LCUC"] <- "LCC"

## create a new variable H_Grade
dat1$H_Grade  <- NA
tmp1 <- match(Zhu43_in_Dir$CEL_ID,dat1$CEL_ID)
dat1$H_Grade[tmp1] <-Zhu43_in_Dir$H_Grade 

## create new columns "GENDER"
dat1$GENDER  <- dat1$Sex

## create new columns "Smoking"
dat1$SMOKING  <- NA
dat1$SMOKING[tmp1] <- Zhu43_in_Dir$SMOKING

## create new columns "Race"
dat1$RACE  <- NA
dat1$RACE[tmp1] <-Zhu43_in_Dir$RACE

## create new columns "TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT"
dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT  <- NA
dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT[tmp1]  <-Zhu43_in_Dir$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT

## CREATE A COLUMN ADJUVANT TREATMENT
dat1$totalAdjuvant[dat1$Post.Surgical.Treatment=="OBS"] <- 0
dat1$totalAdjuvant[dat1$Post.Surgical.Treatment=="ACT"] <- 1

## create new columns "REC_STATUS"
dat1$REC_STATUS  <- NA
dat1$REC_STATUS[tmp1] <- Zhu43_in_Dir$REC_STATUS

## create new columns "SITE"
dat1$SITE  <- NA
dat1$SITE[tmp1] <- Zhu43_in_Dir$SITE

## create new columns ""LABORATORY_BATCH"
dat1$LABORATORY_BATCH  <- NA
dat1$LABORATORY_BATCH  <- paste("Zhu",dat1$LABORATORY_BATCH,sep="")

##CREATE NEW COLUMNS: WARNING
dat1$WARNING  <- NA
dat1$WARNING[tmp1] <- Zhu43_in_Dir$WARNING

##CREATE NEW COLUMNS: Surg_Type
dat1$Surg_Type  <- NA

## CREATE NEW VARIABLE NAMED P_Stage
dat1$P_Stage <- dat1$Stage

## CREATE new table...ZhuData
ZhuData  <- dat1[,c("CEL_ID", "Histology", "H_Grade","P_Stage","GENDER", "Age", "SMOKING", "RACE", "Surg_Type","totalAdjuvant", "VITAL_STATUS", "MONTHS_TO_LAST_CONTACT_OR_DEATH", "REC_STATUS" ,"TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT","SITE","LABORATORY_BATCH", "WARNING")]

## define CEL path in order to scan the date of production of the cel files
celPath  <- "/home/cferte/NSCLC_Prognostic_Signature/Zhu_GSE14814/CEL_files/compressed_Zhu_CEL_files/"
library(affy)
fns  <- list.celfiles(full.names=T,path= celPath)
blah <- ReadAffy(filenames=fns)

## scan the date of the CEL files and attribute a number to each BatchScan
ScanDate  <- blah@protocolData@data
tmpNames <- strsplit(rownames(ScanDate), ".", fixed=T)
ScanDate$CEL_ID <- unlist(lapply(tmpNames, "[[", 1))

ScanDate$ScanDate <- substr(ScanDate$ScanDate,1,8)
ScanDate$ScanDate  <- as.Date(ScanDate$ScanDate, "%m/%d/%y")
ScanDate$SCANBATCH  <- paste("Zhu", as.numeric(factor(ScanDate$ScanDate)), sep="_")

##merge the tables
ZhuData <- merge(x=ZhuData,y=ScanDate,by="CEL_ID", all.x=T, all.y=F)

## name OS and PSTAGE:
OS <- ZhuData$MONTHS_TO_LAST_CONTACT_OR_DEATH
BATCH  <- substr(ZhuData$LABORATORY_BATCH,4,5)
PSTAGE  <- ZhuData$P_Stage
SCANBATCH  <- ZhuData$SCANBATCH

## plot OS according to pSTAGE
library(survival)
OSs <- ifelse(ZhuData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ PSTAGE), col=1:length(unique(PSTAGE)),main="OS according to pSTAGE",xlab="months",ylab="Overall Survival (%)")

## plot OS according to SCANBATCH
library(survival)
OSs <- ifelse(ZhuData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ SCANBATCH), col=1:length(unique(SCANBATCH)),main="OS according to SCANBATCH",xlab="months",ylab="Overall Survival (%)")

## create and paste the new table girard data in  the fellow/cferte folder
setwd(dataPath2)
write.table(ZhuData, file="ZhuData.csv", sep=",",row.names=FALSE,quote=FALSE)

             
             