## POINT TO DATA DIRECTORY
dataPath <- "/home/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Director_challenge/series_matrix/directors_suppl_files/"

## READ IN THE TAB DELIMITED FILES
dat1 <- read.delim2(file = paste(dataPath, "DCLungStudy_Clinical_Covariates_with_Hgrade_Try.txt", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
dat2 <- read.delim2(file = paste(dataPath, "UpdatedTherapiesfor_DCLungStudy_Clinical_Covariates_with_Hgrade.txt", sep=""), na.strings=c("NA", "", " "), dec=",",header=T, as.is=T)

## WE KNOW THAT THERE ARE SOME ROWS COLUMNS THAT DO NOT BELONG IN DC STUDY
dat1 <- dat1[ dat1$IN_DC_STUDY==1, ]

## WE KNOW THAT THERE ARE DUPLICATES IN DAT2 - SO ONLY LOOK AT UNIQUE
dat2 <- unique(dat2)

##MERGE THE TWO TABLES dat1 dat2
datNew <- merge(x=dat1, y=dat2, by="DC_STUDY_ID", all.x=T, all.y=F)

## CORRECT CHEMO
datNew$correctADJUVANT_CHEMO <- datNew$revisedADJUVANT_CHEMO
datNew$correctADJUVANT_CHEMO[ is.na(datNew$correctADJUVANT_CHEMO) ] <- datNew$ADJUVANT_CHEMO[ is.na(datNew$correctADJUVANT_CHEMO) ]

## CORRECT RT
datNew$correctADJUVANT_RT <- datNew$revisedADJUVANT_RT
datNew$correctADJUVANT_RT[is.na(datNew$correctADJUVANT_RT)] <- datNew$ADJUVANT_RT[is.na(datNew$correctADJUVANT_RT)]

## DELETE patinets with "unknown" CHEMO treatment status
datNew  <- datNew[datNew$correctADJUVANT_CHEMO != "Unknown",]

## DELETE patinets with "unknown" RT treatment status
datNew  <- datNew[datNew$correctADJUVANT_RT != "Unknown",]

## AVOID BOTH ADJUVANT THERAPIES
datNew$totalAdjuvant <- ifelse( datNew$correctADJUVANT_CHEMO=="No" & datNew$correctADJUVANT_RT=="No", 0, 1)

## KEEP ONLY PATIENTS (ie: the ones with OS data)
datNew  <- datNew [datNew$MONTHS_TO_LAST_CONTACT_OR_DEATH != "na",]
datNew <- datNew[!is.na(datNew$MONTHS_TO_LAST_CONTACT_OR_DEATH), ]

##CREATE NEW COLUMN: CEL_ID
datNew$CEL_ID <- datNew$MICROARRAY

##CREATE NEW COLUMNS: Histology (here all adc)
datNew$Histology <- "ADC"

##CREATE NEW COLUMNS: H_Grade
datNew$Histologic.grade <- sub("POORLY DIFFERENTIATED","3",datNew$Histologic.grade)
datNew$Histologic.grade <- sub("WELL DIFFERENTIATED","1",datNew$Histologic.grade)
datNew$Histologic.grade <- sub("Moderate Differentiation","2",datNew$Histologic.grade)
datNew$H_Grade  <- datNew$Histologic.grade

##CREATE NEW COLUMNS: P_Stage and remove those with unknown p stage status
datNew$P_Stage  <- NULL
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N0. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T1. ACCORDING TO AJCC CRITERIA"]  <- "1A"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N0. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T2. ACCORDING TO AJCC CRITERIA"]  <- "1B"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N1. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T1. ACCORDING TO AJCC CRITERIA"]  <- "2A"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N1. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T2. ACCORDING TO AJCC CRITERIA"]  <- "2B"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N0. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T3. ACCORDING TO AJCC CRITERIA"]  <- "2B"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N1. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T3. ACCORDING TO AJCC CRITERIA"]  <- "3A"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N2. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T1. ACCORDING TO AJCC CRITERIA"]  <- "3A"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N2. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T2. ACCORDING TO AJCC CRITERIA"]  <- "3A"
datNew$P_Stage [datNew$PATHOLOGIC_N_STAGE=="N2. ACCORDING TO AJCC CRITERIA" & datNew$PATHOLOGIC_T_STAGE=="T3. ACCORDING TO AJCC CRITERIA"]  <- "3A"
datNew$P_Stage [datNew$PATHOLOGIC_T_STAGE=="T4. ACCORDING TO AJCC CRITERIA"]  <- "3B"
datNew <- datNew[!is.na(datNew$P_Stage), ]

##CREATE NEW COLUMNS:Age
datNew$Age  <-datNew$AGE_AT_DIAGNOSIS

##CREATE NEW COLUMNS: Surg_Type (Pneumonecctomy 1 lobectomy 0)
datNew$Surg_Type  <- NA
colnames(datNew)

##CREATE NEW COLUMNS: SMOKING (CURRENT PAST NEVER NA)
datNew$SMOKING  <-  NULL
datNew$SMOKING[datNew$SMOKING_HISTORY== "Smoked in the past"]  <- "PAST"
datNew$SMOKING[datNew$SMOKING_HISTORY== "Currently smoking"]  <- "CURRENT"
datNew$SMOKING[datNew$SMOKING_HISTORY== "Never smoked"]  <- "NEVER"
datNew$SMOKING[datNew$SMOKING_HISTORY== "Unknown"]  <- NA

##ARRANGE COLUMN RACE (grouping Not reported and Unknown and missing values as NA)
datNew$RACE[datNew$RACE %in% c("Unknown(99)" ,"Not Reported(98)")] <- NA
datNew$RACE[datNew$RACE=="Asian(05)"] <- "Asian"
datNew$RACE[datNew$RACE %in% c("Black or African American(03)")] <- "Black/Af_Am"
datNew$RACE[datNew$RACE=="White(01)"] <- "White"
datNew$RACE[datNew$RACE=="Native Hawaiian or Other Pacific Islander(04)"] <- "Pacific"

##CREATE NEW COLUMNS: REC_Status (event status death 1 alive 0) and OS (time to last clinical evaluation or death)
datNew$REC_STATUS  <- datNew$FIRST_PROGRESSION_OR_RELAPSE
datNew$REC_STATUS[datNew$REC_STATUS=="Unknown"]  <- NA
table(datNew$REC_STATUS)

##CREATE NEW COLUMNS: TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT
datNew$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- NA
datNew$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT[datNew$REC_STATUS=="Yes" & !is.na(datNew$REC_STATUS)]  <- datNew$MONTHS_TO_FIRST_PROGRESSION[datNew$REC_STATUS=="Yes" & !is.na(datNew$REC_STATUS)]
datNew$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT[datNew$REC_STATUS=="No" & !is.na(datNew$REC_STATUS)]  <-datNew$MTHS_TO_LAST_CLINICAL_ASSESSMENT[datNew$REC_STATUS=="No" & !is.na(datNew$REC_STATUS)]

## ARRANGE COLUMN LABORATORY_BATCH
datNew$LABORATORY_BATCH  <- paste("Dir",datNew$LABORATORY_BATCH, sep="")

## CREATE new table...DirectorData
DirectorData  <- datNew[,c("CEL_ID", "Histology", "H_Grade","P_Stage","GENDER", "Age", "SMOKING", "RACE", "Surg_Type","totalAdjuvant", "VITAL_STATUS", "MONTHS_TO_LAST_CONTACT_OR_DEATH", "REC_STATUS" ,"TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT","SITE","LABORATORY_BATCH", "SITE","WARNING")]

## define CEL path in order to scan the date of production of the cel files
celPath  <- "/home/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Director_challenge/CEL_Files/caArray_jacob-00182_files-AFFYMETRIX_CEL/"
library(affy)
fns  <- list.celfiles(full.names=T,path= celPath)
blah <- ReadAffy(filenames=fns)

## scan the date of the CEL files and attribute a number to each BatchScan
ScanDate  <- blah@protocolData@data
ScanDate$CEL_ID <- rownames(ScanDate)
ScanDate$CEL_ID  <- substr(ScanDate$CEL_ID,1,nchar(ScanDate$CEL_ID)-4)
ScanDate$ScanDate <- substr(ScanDate$ScanDate,1,8)
ScanDate$ScanDate  <- as.Date(ScanDate$ScanDate, "%m/%d/%y")
ScanDate$SCANBATCH  <- paste("Dir", as.numeric(factor(ScanDate$ScanDate)), sep="_")

##merge the tables
DirectorData <- merge(x=DirectorData,y=ScanDate,by="CEL_ID", all.x=T, all.y=F, as.is=T)

## Separate the 43 samples form DirectorData that are actually belon,ging to the JBR10 (ie: the Zhu dataset) 
dataPath2 <- "/home/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Zhu_GSE14814/series_matrix/"
datZhu <- read.delim2(file = paste(dataPath2, "Zhu_2010_JCO_microarray_pts_clin_info-1.txt", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
rownames(datZhu) <- datZhu$arrayname
tmp1 <- match(DirectorData$CEL_ID,datZhu$arrayname)
datZhu_in_Dir <- DirectorData[!is.na(tmp1),]
DirectorData <- DirectorData[is.na(tmp1),]
rm(tmp1,datZhu,dataPath2)

## create and paste the new table Directordata in  the fellow/cferte folder
setwd("/home/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/")
write.table(DirectorData, file="DirectorData.csv", sep=",",row.names=FALSE,quote=FALSE)
write.table(datZhu_in_Dir, file="Zhu43_in_Dir.csv", sep=",",row.names=FALSE,quote=FALSE)

## name OS BATCH SCANBATCH and PSTAGE:
OS <- DirectorData$MONTHS_TO_LAST_CONTACT_OR_DEATH
BATCH  <- substr(DirectorData$LABORATORY_BATCH,4,5)
PSTAGE  <- DirectorData$P_Stage
SCANBATCH  <- DirectorData$SCANBATCH

## plot OS according to pSTAGE
library(survival)
OSs <- ifelse(DirectorData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ PSTAGE), col=1:length(unique(PSTAGE)),main="OS according to pSTAGE",xlab="months",ylab="Overall Survival (%)")

## plot OS according to BATCH
library(survival)
OSs <- ifelse(DirectorData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ BATCH), col=1:length(unique(BATCH)),main="OS according to BATCH",xlab="months",ylab="Overall Survival (%)")

## plot OS according to SCANBATCH
library(survival)
OSs <- ifelse(DirectorData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ SCANBATCH), col=1:length(unique(SCANBATCH)),main="OS according to SCANBATCH",xlab="months",ylab="Overall Survival (%)")

## plot OS according to Adjuvant Chemotherapy
library(survival)
OSs <- ifelse(DirectorData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ DirectorData$totalAdjuvant), col=1:length(unique(DirectorData$totalAdjuvant)),main="OS according to Adjuvant therapy",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS, OSs) ~ DirectorData$totalAdjuvant)

table(DirectorData$P_Stage, DirectorData$totalAdjuvant)

## plot OS according to Adjuvant Chemotherapy in P_Stage >2A
library(survival)
OSs <- ifelse(DirectorData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS[DirectorData$P_Stage %in% c("3A","3B")], OSs[DirectorData$P_Stage %in% c("3A","3B")]) ~ DirectorData$totalAdjuvant[DirectorData$P_Stage %in% c("3A","3B")]), col=1:length(unique(DirectorData$totalAdjuvant)),main="OS according to Adjuvant therapy",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS[DirectorData$P_Stage %in% c("2A","2B","3A","3B")], OSs[DirectorData$P_Stage %in% c("2A","2B","3A","3B")]) ~ DirectorData$totalAdjuvant[DirectorData$P_Stage %in% c("2A","2B","3A","3B")])
