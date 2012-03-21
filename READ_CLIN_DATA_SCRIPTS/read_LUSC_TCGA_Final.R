### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### February, 19th 2012

## reading the clin and exp data from LUSC TCGA

### script for avoiding dataframe/matrix
options(stringsAsFactors=FALSE)

# select PATH
PATH <- "/home"

# load library
library(synapseClient)
library(affy)
library(rSCR)

## POINT TO DATA DIRECTORY

# synapse login
synapseLogin()

unlink("~/FELLOW/cferte/temp_synapse",recursive=TRUE)
setwd("~/FELLOW/cferte/")
dir.create("./temp_synapse")

for(i in myFiles) { 
  file.symlink(i,"./temp_synapse")
  
}

#download the phenotype data from the level 1 SCC NSCLC clinical data and stoire them into LUSC_C
layers <- synapseQuery(paste('select * from layer where layer.tcgaLevel == "Level_1" and layer.parentId == "',datasets[grep("Lung Squamous Cell Carcinoma TCGA" ,datasets$dataset.name), 'dataset.id'],'"',sep=""))
layers$layer.platform
clinicalLayers <- synapseQuery(paste('select * from layer where layer.type == "C" and layer.parentId == "',datasets[grep("Lung Squamous Cell Carcinoma TCGA",datasets$dataset.name), 'dataset.id'],'"',sep=""))
## find the command line !!!
entity_14482 <- loadEntity(14482)
entity_14482$files
LUSC_C <- read.delim(paste(entity_14482$cacheDir,"/","clinical_patient_public_lusc.txt",sep=""), na.strings=c("NA", "", " ","[Not Applicable]"), as.is=T)
rownames(LUSC_C) <- LUSC_C$bcr_patient_barcode

# retreive the correspondance between the names of the CEL files and the barcode number and store it into tmp3
tmp <- handleTcgaClinical(4602)
tmp3 <- tmp
tmp3 <- tmp3[,c("bcr_sample_barcode","broad.mit.edu_LUSC.HT_HG-U133A-Array.Data.File","vital_status","bcr_aliquot_barcode")]
tmp3$bcr_patient_barcode <- substr(tmp3$bcr_sample_barcode,1,12)
tmp3 <- tmp3[!is.na(tmp3[,2]),]

# check that all the CEL files in tmp3 are the one that are in CELNAMES
table(tmp3[,2] %in% CELNAMES)

# match the barcode in tmp3 and LUSC_C
LUSC_C <- LUSC_C[match(tmp3$bcr_patient_barcode,LUSC_C$bcr_patient_barcode),]

# merge tmp3 and LUSC_C 
LUSC_C <- unique(merge(tmp3,LUSC_C))

# delete the duplicated barcodes in LUSC_C
LUSC_C <- LUSC_C[-c(grep(LUSC_C$bcr_patient_barcode[which(duplicated(LUSC_C$bcr_patient_barcode)==T)],LUSC_C$bcr_patient_barcode)),]

# keep only the ones with stage IA to IIIA 
LUSC_C <- LUSC_C[LUSC_C$tumor_stage %in% c("Stage IA","Stage IB","Stage IIA","Stage IIB","Stage IIIA"),]

## get rid of the patients with chemotherapy
grep("drug_name",colnames(tmp))
grep("regimen_indication",colnames(tmp))
grep("bcr_patient_barcode",colnames(tmp))
grep("bcr_drug_uuid",colnames(tmp))
grep("therapy_type",colnames(tmp))
chemo <- tmp[,c(3,34,36,44,46,47)]
chemo <- unique(chemo)
chemo <- chemo[!is.na(chemo$bcr_drug_uuid),]
LUSC_C <- LUSC_C[LUSC_C$bcr_patient_barcode %in% c(setdiff(LUSC_C$bcr_patient_barcode,chemo$bcr_patient_barcode)),]

LUSC_C <- LUSC_C[LUSC_C$pretreatment_history=="NO",]


## get rid of the patients with radiotherapy
colnames(tmp)[grep("radiation",colnames(tmp))]
radio <- tmp[,c("bcr_patient_barcode","bcr_radiation_uuid")]
radio <- unique(radio)
radio <- radio[!is.na(radio$bcr_radiation_uuid),]
LUSC_C <- LUSC_C[LUSC_C$bcr_patient_barcode %in% c(setdiff(LUSC_C$bcr_patient_barcode,radio$bcr_patient_barcode)),]

rm(radio,chemo)

datNew <- LUSC_C

# get rid  of the patients with other than R0 resection
datNew <- datNew[datNew$residual_tumor=="R0",]

##CREATE NEW COLUMN: CEL_ID
datNew$CEL_ID <- datNew[,4]

##CREATE NEW COLUMNS: Histology (here all adc)
datNew$Histology <- "SCC"

##CREATE NEW COLUMNS: H_Grade (NA in TCGA...)
datNew$H_Grade  <- NA

##CREATE NEW COLUMNS: P_Stage and remove those with unknown p stage status
datNew$P_Stage  <- datNew$tumor_stage
datNew <- datNew[!is.na(datNew$P_Stage), ]
datNew$P_Stage <- substr(datNew$P_Stage,7,nchar(datNew$P_Stage))
datNew$P_Stage <- sub("III","3",datNew$P_Stage)
datNew$P_Stage <- sub("II","2",datNew$P_Stage)
datNew$P_Stage <- sub("I","1",datNew$P_Stage)

##CREATE NEW COLUMNS:Age
datNew$Age  <-datNew$age_at_initial_pathologic_diagnosis

##CREATE NEW COLUMNS: Surg_Type (Pneumonecctomy 1 lobectomy 0)
datNew$Surg_Type  <- NA
colnames(datNew)

##CREATE NEW COLUMNS: SMOKING (CURRENT PAST NEVER NA)
datNew$SMOKING  <-  datNew$tobacco_smoking_history_indicator
datNew$SMOKING[datNew$SMOKING %in% c("Current reformed smoker for < or = 15 years","Current reformed smoker for > 15 years")]  <- "PAST"
datNew$SMOKING[datNew$SMOKING == "Current smoker"]  <- "CURRENT"
datNew$SMOKING[datNew$SMOKING == "Lifelong Non-smoker"]  <- "NEVER"
datNew$SMOKING[datNew$SMOKING == "[Not Available]"]  <- NA

##ARRANGE COLUMN RACE (grouping Not reported and Unknown and missing values as NA)
datNew$RACE <- NA

## KEEP ONLY PATIENTS  with OS data)
datNew$days_to_last_known_alive <- as.numeric(datNew$days_to_last_known_alive)*12/365.25
datNew$MONTHS_TO_LAST_CONTACT_OR_DEATH <- datNew$days_to_last_known_alive
datNew <- datNew[!is.na(datNew$MONTHS_TO_LAST_CONTACT_OR_DEATH), ]
datNew$VITAL_STATUS <- datNew$vital_status
datNew$VITAL_STATUS<- ifelse(datNew$vital_status=="LIVING","Alive","Dead")
datNew <- datNew[!is.na(datNew$VITAL_STATUS), ]

## ARRANGE COLUMN LABORATORY_BATCH
datNew$LABORATORY_BATCH  <- substr(datNew$bcr_aliquot_barcode,22,25)

## ARRANGE COLUMN SITE
datNew$SITE <- substr(datNew$bcr_sample_barcode,6,7)
 

## GENDER
datNew$GENDER <- datNew$gender

datNew$REC_STATUS <- NA
datNew$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- NA
datNew$WARNING <- NA
datNew$Surg_Type <- NA
datNew$totalAdjuvant <- 0

## CREATE new table...LUSCData
LUSCData  <- datNew[,c("CEL_ID", "Histology", "H_Grade","P_Stage","GENDER", "Age", "SMOKING", "RACE", "Surg_Type","totalAdjuvant", "VITAL_STATUS", "MONTHS_TO_LAST_CONTACT_OR_DEATH", "REC_STATUS" ,"TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT","SITE","LABORATORY_BATCH","WARNING")]
colnames(LUSCData)


## scan the date of the CEL files and attribute a number to each BatchScan
ScanDate  <- LUSC_E@phenoData@data
ScanDate$CEL_ID <- rownames(ScanDate)
ScanDate$CEL_ID  <- substr(ScanDate$CEL_ID,1,nchar(ScanDate$CEL_ID)-4)
ScanDate$SCANBATCH  <- paste("LUSC", as.numeric(factor(ScanDate$ScanDate)), sep="_")
ScanDate$sample <- NULL
ScanDate$ScanDate <- NA

##merge the tables
LUSCData <- merge(x=LUSCData,y=ScanDate,by="CEL_ID", all.x=T, all.y=F, as.is=T)

## create and paste the new table Directordata in  the fellow/cferte folder
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/",sep=""))
write.table(LUSCData, file="LUSCData.csv", sep=",",row.names=FALSE,quote=FALSE)

## name OS BATCH SCANBATCH and PSTAGE:
OS <- LUSCData$MONTHS_TO_LAST_CONTACT_OR_DEATH
BATCH  <- LUSCData$LABORATORY_BATCH
PSTAGE  <- LUSCData$P_Stage
SCANBATCH  <- LUSCData$SCANBATCH
SITE <- LUSCData$SITE

## plot OS according to pSTAGE
library(survival)
OSs <- ifelse(LUSCData$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".", OS))
plot(survfit(Surv(OS, OSs) ~ PSTAGE), col=1:length(unique(PSTAGE)),main="OS according to pSTAGE",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS, OSs) ~ PSTAGE)

## plot OS according to BATCH
library(survival)
plot(survfit(Surv(OS, OSs) ~ BATCH), col=1:length(unique(BATCH)),main="OS according to BATCH",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS, OSs) ~ BATCH)

## plot OS according to SITE
library(survival)
plot(survfit(Surv(OS, OSs) ~ SITE), col=1:length(unique(SITE)),main="OS according to SITE",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS, OSs) ~ SITE)
