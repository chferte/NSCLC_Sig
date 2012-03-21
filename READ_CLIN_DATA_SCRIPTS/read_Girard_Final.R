### Charles Ferté
### Sage Bionetworks,
### Seattle,WA

## POINT TO DATA DIRECTORY
dataPath <- "/home/cferte/NSCLC_Prognostic_Signature/Girard_GSE31908/series_matrix/"

## READ IN THE TAB DELIMITED FILES
dat1 <- read.delim2(file = paste(dataPath, "GSE31908-GPL570_series_matrix_CORRECTED.txt", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")

## these are the different columns i want to create
FCOL <- c("CEL_ID", "Histology", "H_Grade","P_Stage","GENDER", "Age", "SMOKING", "RACE", "Surg_Type","totalAdjuvant", "VITAL_STATUS", "MONTHS_TO_LAST_CONTACT_OR_DEATH", "REC_STATUS" ,"TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT","SITE","LABORATORY_BATCH")
length(FCOL)

## rename sample ID (dat3) into CEL_ID 
colnames(dat1)[2]  <- FCOL[1]
colnames(dat1)[8]  <- FCOL[2]
dat1$Histology[dat1$Histology=="Lung adenocarcinoma"] <- "ADC"

## arrange column GENDER
dat1$GENDER[dat1$GENDER=="gender: Female"] <- "Female"
dat1$GENDER[dat1$GENDER=="gender: Male"] <- "Male"

##arrange column Age
dat1$Age <-NA

##arrange column Surg_Type
dat1$Surg_Type <-NA

##arrange column SITE
dat1$SITE <-NA

##arrange column Surg_Type
dat1$WARNING <-NA

##arrange column RACE
dat1$RACE <-substr(dat1$RACE,7,(nchar(dat1$RACE)))
dat1$RACE[dat1$RACE %in% c("Not Reported","Unknown")] <- NA


## remove column 26 as it is NULL
dat1 <- dat1[,-26]
dat1 <- dat1[,-3:-7]


## arrange EGFR STATUS COLLUMN
dat1$EGFR_STATUS <-substr(dat1$EGFR_STATUS,13,(nchar(dat1$EGFR_STATUS)))

## arrange ADJUVANT CHEMO RX INTO totalAdjuvant
dat1$ADJUVANT_CHEMO <-substr(dat1$ADJUVANT_CHEMO,17,(nchar(dat1$ADJUVANT_CHEMO)))
dat1$ADJUVANT_RX <-substr(dat1$ADJUVANT_RX,14,(nchar(dat1$ADJUVANT_RX)))
dat1$totalAdjuvant <- ifelse(dat1$ADJUVANT_CHEMO=="No" & dat1$ADJUVANT_RX=="No","0","1")

## arrange H_Grade column
dat1$H_Grade  <- substr(dat1$H_Grade,19,nchar(dat1$H_Grade))
dat1$H_Grade[dat1$H_Grade=="WELL DIFFERENTIATED"]  <- "1"
dat1$H_Grade[dat1$H_Grade=="Moderate Differentiation"]  <- "2"


## arrange P_Stage column
dat1$pN <- substr(dat1$pN,21,22)
dat1$pT <- substr(dat1$pT,21,22)
dat1$P_Stage <- paste(dat1$pT,dat1$pN,sep="")
dat1$P_Stage[dat1$P_Stage=="T1N0"] <- "1A"
dat1$P_Stage[dat1$P_Stage=="T2N0"] <- "1B"
dat1$P_Stage[dat1$P_Stage=="T2N1"] <- "2B"

## create new columns SMOKING
dat1$SMOKING <- substr(dat1$SMOKING,18,nchar(dat1$SMOKING))
dat1$SMOKING[dat1$SMOKING=="Smoked in the past"] <- "PAST"
dat1$SMOKING[dat1$SMOKING=="Never smoked"] <- "NEVER"

## create new columns LABORATORY BATCH
dat1$LABORATORY_BATCH <- substr(dat1$LABORATORY_BATCH,19,nchar(dat1$LABORATORY_BATCH))
dat1$LABORATORY_BATCH <- paste("Gir",dat1$LABORATORY_BATCH,sep="")


## arrange  MONTHS_TO_LAST_CONTACT_OR_DEATH (months)
dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH <- substr(dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH,34,nchar(dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH))
dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH <- as.numeric(dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH)

## create new columns "TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT"
dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT  <- NA
dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE  <- as.character(dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE)
dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE <- substr(dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE,30,nchar(dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE))
dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT  <- as.character(dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT)
dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT <- substr(dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT,37,nchar(dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT))
dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE  <- as.numeric(dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE)
dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT  <- as.numeric(dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT)
dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- dat1$MONTHS_TO_LAST_CLINICAL_ASSESSMENT
dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT[is.na(dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT)] <- dat1$MONTHS_TO_FIRST_PROGRESSION_RELPASE[is.na(dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT)]

## create new columns "REC_STATUS"
dat1$REC_STATUS <- substr(dat1$REC_STATUS,30,nchar(dat1$REC_STATUS))

## create new columns "VITAL_STATUS"
dat1$VITAL_STATUS <- substr(dat1$VITAL_STATUS,14,nchar(dat1$VITAL_STATUS))

## CREATE new table...GirardData
GirardData  <- dat1[,c("CEL_ID", "Histology", "H_Grade","P_Stage","GENDER", "Age", "SMOKING", "RACE", "Surg_Type","totalAdjuvant", "VITAL_STATUS", "MONTHS_TO_LAST_CONTACT_OR_DEATH", "REC_STATUS" ,"TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT","SITE","LABORATORY_BATCH", "WARNING")]
rm(dat1)

## define CEL path in order to scan the date of production of the cel files
celPath  <- "/home/cferte/NSCLC_Prognostic_Signature/GEO_FTP/GSE31908/"
library(affy)
fns  <- list.celfiles(full.names=T,path= celPath)
blah <- ReadAffy(filenames=fns)

## scan the date of the CEL files and attribute a number to each BatchScan
ScanDate  <- blah@protocolData@data
ScanDate$CEL_ID <- rownames(ScanDate)
ScanDate$CEL_ID <- substr(ScanDate$CEL_ID,1,nchar(ScanDate$CEL_ID)-7)
ScanDate$ScanDate <- substr(ScanDate$ScanDate,1,8)
ScanDate$ScanDate  <- as.Date(ScanDate$ScanDate, "%m/%d/%y")
ScanDate$SCANBATCH  <- paste("Gir", as.numeric(factor(ScanDate$ScanDate)), sep="_")

##merge the tables
GirardData <- merge(x=GirardData,y=ScanDate,by="CEL_ID", all.x=T, all.y=F)

## name OS BATCH SCANBATCH and PSTAGE:
OS <- GirardData$MONTHS_TO_LAST_CONTACT_OR_DEATH
BATCH  <- substr(GirardData$LABORATORY_BATCH,4,5)
PSTAGE  <- GirardData$P_Stage
SCANBATCH  <- GirardData$SCANBATCH

## plot OS according to pSTAGE
library(survival)
OSs <- ifelse(GirardData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ PSTAGE), col=1:length(unique(PSTAGE)),main="OS according to pSTAGE",xlab="months",ylab="Overall Survival (%)")

## plot OS according to BATCH
library(survival)
OSs <- ifelse(GirardData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ BATCH), col=1:length(unique(BATCH)),main="OS according to BATCH",xlab="months",ylab="Overall Survival (%)")

## plot OS according to SCANBATCH
library(survival)
OSs <- ifelse(GirardData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ SCANBATCH), col=1:length(unique(SCANBATCH)),main="OS according to SCANBATCH",xlab="months",ylab="Overall Survival (%)")

## create and paste the new table girard data in  the fellow/cferte folder
setwd("/home/cferte/FELLOW/cferte/")
write.table(GirardData, file="GirardData.csv", sep=",",row.names=FALSE,quote=FALSE)
