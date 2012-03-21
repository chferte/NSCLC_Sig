## POINT TO DATA DIRECTORY
dataPath <- "/Volumes/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Mitra_GSE9971/publication/"


## READ IN THE TAB DELIMITED FILES
dat1 <- read.delim2(file = paste(dataPath, "lung cancer raw data_Mitra.txt", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
dat2 <- read.delim2(file = paste(dataPath, "correspondance_GSM_SampleID_Mitra.txt", sep=""), na.strings=c("NA", "", " "), dec=",",header=F, as.is=T)

## arrange dat2
dat2 <- dat2[,-1]
dat2 <- t(dat2)
dat2 <- data.frame(dat2)
colnames(dat2)  <- c("sample_ID","CEL_ID")
dat2$GEO.R.NR <- ifelse(substr(dat2$sample_ID,1,3)=="Rec",paste("R",substr(dat2$sample_ID,17,19),sep=""),paste("NR",substr(dat2$sample_ID,21,23),sep=""))
dat1 <- dat1[-29,]
dat1 <- dat1[-28,]

##MERGE THE TWO TABLES dat1 dat2
datN <- merge(x=dat1, y=dat2, by="GEO.R.NR", all.x=T, all.y=T)

## remove dat1,dat2
rm(dat1)
rm(dat2)

## create variable totalAdjuvant = 0 for all patients since none received adj chemotherapy/neoadjuvant chemotherapy
datN$totalAdjuvant <- "0"

# removre the three pts with neoadjuvant chemotherapy:UH0101-38, UH0202-02, UH0107-09



## Arrange COLUMNS: Histology (here all adc)
datN$Histology[datN$Histology %in% c("Adeno","adeno")] <- "ADC"
datN$Histology[datN$Histology =="Squamous"] <- "SCC"
datN$Histology[datN$Histology =="other"] <- NA

##CREATE NEW COLUMNS: H_Grade
datN$H_Grade <- NA

##CREATE NEW COLUMNS: P_Stage according to the values of the Stage Column
datN$P_Stage  <- NULL
datN$P_Stage [datN$STAGE=="IA"]  <- "1A"
datN$P_Stage [datN$STAGE=="IB"]  <- "1B"
datN$P_Stage [datN$STAGE=="IIA"]  <- "2A"
datN$P_Stage [datN$STAGE=="IIB"]  <- "2B"
datN$P_Stage [datN$STAGE=="IIIA"]  <- "3A"
datN$P_Stage [datN$STAGE=="IIIB"]  <- "3B"

##CREATE NEW COLUMNS:Age
datN$Age  <-datN$Age.at.Op

##CREATE NEW COLUMNS: Surg_Type (Pneumonecctomy 1 lobectomy 0)
datN$Surg_Type  <- NA

##CREATE NEW COLUMNS: SMOKING (CURRENT PAST NEVER NA)
datN$SMOKING  <-  NA
datN$SMOKING [datN$SMOKING.HX=="non smoker"]  <- "NEVER"
datN$SMOKING [datN$SMOKING.HX %in% c("past history","Brief smoking history")]  <- "PAST"
datN$SMOKING [datN$SMOKING.HX=="Smoker quit"]  <- "CURRENT"

##ARRANGE COLUMN RACE (transform Unknown into NA)
datN$RACE[datN$RACE=="U"] <- NA

##CREATE NEW COLUMNS: REC_STATUS (event status death 1 alive 0) and OS (time to last clinical evaluation or death)
datN$REC_STATUS <- 0
datN$REC_STATUS [!is.na(datN$Date.of.Recurrence)] <- 1 

##CREATE NEW COLUMNS: VITAL_STATUS (event status death 1 alive 0)
datN$VITAL_STATUS  <- 0
datN$VITAL_STATUS [!is.na(datN$Date.of.death)] <- 1

## consider operation date as a date of forma standard in r
datN$Operation.date <- as.Date(datN$Operation.date,format='%d/%m/%y')

## consider last time visit as a date of forma standard in r
datN$Last.visit.time <- as.Date(datN$Last.visit.time,format='%d/%m/%y')

## transform error in the date of reccurence column (row 5)
datN$Date.of.Recurrence[5]  <- "15/06/04"
datN$Date.of.Recurrence <- as.Date(datN$Date.of.Recurrence,format='%d/%m/%y')

# transform error in the date of death (row 22)
datN$Date.of.death[22]  <- "15/11/02"
datN$Date.of.death <- as.Date(datN$Date.of.death,format='%d/%m/%y')

##CREATE NEW COLUMNS: TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT
datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- datN$Date.of.Recurrence
datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT[is.na(datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT)] <- datN$Last.visit.time[is.na(datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT)]
datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- difftime(datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT, datN$Operation.date, units="days")
datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT *12/365.25

##CREATE NEW COLUMNS: "MONTHS_TO_LAST_CONTACT_OR_DEATH"
datN$MONTHS_TO_LAST_CONTACT_OR_DEATH  <- datN$Date.of.death
datN$MONTHS_TO_LAST_CONTACT_OR_DEATH[is.na(datN$Date.of.death)] <- datN$Last.visit.time[is.na(datN$Date.of.death)]
datN$MONTHS_TO_LAST_CONTACT_OR_DEATH <- difftime(datN$MONTHS_TO_LAST_CONTACT_OR_DEATH, datN$Operation.date, units="days")
datN$MONTHS_TO_LAST_CONTACT_OR_DEATH <- datN$MONTHS_TO_LAST_CONTACT_OR_DEATH*12/365.25

## CREATE NEW COLUMN LABORATORY_BATCH
datN$LABORATORY_BATCH  <- "MitNA"

## CREATE NEW COLUMN GENDER
datN$GENDER <- datN$SEX

## CREATE new table...DirectorData
MitraData  <- datN[,c("CEL_ID", "Histology", "H_Grade","P_Stage","GENDER", "Age", "SMOKING", "RACE", "Surg_Type","totalAdjuvant", "VITAL_STATUS", "MONTHS_TO_LAST_CONTACT_OR_DEATH", "REC_STATUS" ,"TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT","LABORATORY_BATCH")]
View(MitraData)
setwd("/volumes/cferte/FELLOW/cferte/")
write.table(MitraData, file="MitraData.csv", sep=",",row.names=FALSE,quote=FALSE)

## name OS and PSTAGE:
OS <- MitraData$MONTHS_TO_LAST_CONTACT_OR_DEATH
BATCH  <- substr(MitraData$LABORATORY_BATCH,4,5)
PSTAGE  <- MitraData$P_Stage

## plot OS according to pSTAGE
library(survival)
OSs <- ifelse(MitraData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ PSTAGE), col=1:length(unique(PSTAGE)),main="OS according to pSTAGE",xlab="months",ylab="Overall Survival (%)")

## plot OS according to BATCH
library(survival)
OSs <- ifelse(MitraData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ BATCH), col=1:length(unique(BATCH)),main="OS according to BATCH",xlab="months",ylab="Overall Survival (%)")
