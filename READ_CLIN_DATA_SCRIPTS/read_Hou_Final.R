## POINT TO DATA DIRECTORY
dataPath <- "/home/cferte/NSCLC_Prognostic_Signature/Hou_GSE19188/publication/"


## READ IN THE TAB DELIMITED FILES
dat1 <- read.delim2(file = paste(dataPath, "GSE19188_clin_info.txt", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")



## keep only the tumor samples
dat1  <-  dat1[dat1$TissueType=="tumor",]

## keep only the samples with OS available
dat1  <-  dat1[dat1$Overall.Survival..month. !="Not available",]

## create variable totalAdjuvant = 0 for all patients since none received adj chemotherapy/neoadjuvant chemotherapy
dat1$totalAdjuvant <- "0"

##RENAME COLUMN GSM NUMBER...INTO CEL_ID
dat1$CEL_ID  <- dat1$Array.GSM.IDs
### Charles Ferté
### Sage Bionetworks
### Seattle, WA


##CREATE NEW COLUMNS: Histology (here all adc)
dat1$Histology <- dat1$Cell.Type

##CREATE NEW COLUMNS: H_Grade
dat1$H_Grade <- NA

##CREATE NEW COLUMNS: P_Stage according to the values of the Stage Column
dat1$P_Stage  <- NULL
dat1$P_Stage [dat1$Stage=="IA"]  <- "1A"
dat1$P_Stage [dat1$Stage=="IB"]  <- "1B"
dat1$P_Stage [dat1$Stage=="IIA"]  <- "2A"
dat1$P_Stage [dat1$Stage=="IIB"]  <- "2B"
dat1$P_Stage [dat1$Stage=="IIIA"]  <- "3A"
dat1$P_Stage [dat1$Stage=="IV"]  <- "4"


##CREATE NEW COLUMNS:Age
dat1$Age  <-dat1$AgeatDignosis

##CREATE NEW COLUMNS: Surg_Type (Pneumonecctomy 1 lobectomy 0)
dat1$Surg_Type  <- NA

##CREATE NEW COLUMNS: SMOKING (CURRENT PAST NEVER NA)
dat1$SMOKING  <-  NA

##ARRANGE COLUMN RACE (1= caucasian /white, 2 unknown, 3 other)
dat1$RACE  <- NA
dat1$RACE[dat1$Ethinic == "1"] <- "White"
dat1$RACE[dat1$Ethinic == "2"] <- NA
dat1$RACE[dat1$Ethinic == "3"] <- "Other"

##CREATE NEW COLUMNS: REC_STATUS (event status death 1 alive 0) and OS (time to last clinical evaluation or death)
dat1$REC_STATUS  <- NA

##CREATE NEW COLUMNS: TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT
dat1$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- NA

##CREATE NEW COLUMNS: VITAL_STATUS (event status death 1 alive 0)
dat1$VITAL_STATUS  <- dat1$Status..0.alive.1.deceased.
dat1$VITAL_STATUS[dat1$VITAL_STATUS==1] <- "Dead"
dat1$VITAL_STATUS[dat1$VITAL_STATUS==0] <- "Alive"

##CREATE NEW COLUMNS: "MONTHS_TO_LAST_CONTACT_OR_DEATH"
dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH <- dat1$Overall.Survival..month.

## CREATE NEW COLUMN LABORATORY_BATCH
dat1$LABORATORY_BATCH  <- "HouNA"

## CREATE NEW COLUMN GENDER
dat1$GENDER <- dat1$Gender
dat1$GENDER[dat1$GENDER=="F"] <- "Female"
dat1$GENDER[dat1$GENDER=="M"] <- "Male"

## CREATE new table...DirectorData
HouData  <- dat1[,c("CEL_ID", "Histology", "H_Grade","P_Stage","GENDER", "Age", "SMOKING", "RACE", "Surg_Type","totalAdjuvant", "VITAL_STATUS", "MONTHS_TO_LAST_CONTACT_OR_DEATH", "REC_STATUS" ,"TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT","LABORATORY_BATCH")]
View(HouData)

## define CEL path in order to scan the date of production of the cel files
celPath  <- "/home/cferte/NSCLC_Prognostic_Signature/GEO_FTP/GSE19188/"

library(affy)
fns  <- list.celfiles(full.names=T,path= celPath)
blah <- ReadAffy(filenames=fns)

## scan the date of the CEL files and attribute a number to each BatchScan
ScanDate  <- blah@protocolData@data
ScanDate$CEL_ID <- rownames(ScanDate)
ScanDate$CEL_ID <- substr(ScanDate$CEL_ID,1,nchar(ScanDate$CEL_ID)-7)
ScanDate$ScanDate <- substr(ScanDate$ScanDate,1,8)
ScanDate$ScanDate  <- as.Date(ScanDate$ScanDate, "%m/%d/%y")
ScanDate$SCANBATCH  <- paste("Hou", as.numeric(factor(ScanDate$ScanDate)), sep="_")

##merge the tables
HouData <- merge(x=HouData,y=ScanDate,by="CEL_ID", all.x=T, all.y=F)

## name OS BATCH SCANBATCH and PSTAGE:
OS <- HouData$MONTHS_TO_LAST_CONTACT_OR_DEATH
BATCH  <- substr(HouData$LABORATORY_BATCH,4,5)
PSTAGE  <- HouData$P_Stage
SCANBATCH  <- HouData$SCANBATCH

## plot OS according to pSTAGE
library(survival)
OSs <- ifelse(HouData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ PSTAGE), col=1:length(unique(PSTAGE)),main="OS according to pSTAGE",xlab="months",ylab="Overall Survival (%)")

## plot OS according to BATCH
library(survival)
OSs <- ifelse(HouData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ BATCH), col=1:length(unique(BATCH)),main="OS according to BATCH",xlab="months",ylab="Overall Survival (%)")

## plot OS according to SCANBATCH
library(survival)
OSs <- ifelse(HouData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ SCANBATCH), col=1:length(unique(SCANBATCH)),main="OS according to SCANBATCH",xlab="months",ylab="Overall Survival (%)")

## create and paste the new table girard data in  the fellow/cferte folder
setwd("/home/cferte/FELLOW/cferte/")
write.table(HouData, file="HouData.csv", sep=",",row.names=FALSE,quote=TRUE)
