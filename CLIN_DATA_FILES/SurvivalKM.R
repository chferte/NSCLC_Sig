## POINT TO DATA DIRECTORY
dataPath <- "/volumes/cferte/NSCLC_Prognostic_Signature/Director_challenge/series_matrix/directors_suppl_files/"


## READ IN THE TAB DELIMITED FILES
dat1 <- read.delim2(file = paste(dataPath, "DCLungStudy_Clinical_Covariates_with_Hgrade.txt", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
dat2 <- read.delim2(file = paste(dataPath, "UpdatedTherapiesfor_DCLungStudy_Clinical_Covariates_with_Hgrade.txt", sep=""), na.strings=c("NA", "", " "), dec=",",header=T, as.is=T)

## WE KNOW THAT THERE ARE SOME COLUMNS THAT DO NOT BELONG IN DC STUDY
dat1 <- dat1[ dat1$IN_DC_STUDY==1, ]

## WE KNOW THAT THERE ARE DUPLICATES IN DAT2 - SO ONLY LOOK AT UNIQUE
dat2 <- unique(dat2)

datNew <- merge(x=dat1, y=dat2, by="DC_STUDY_ID", all.x=T, all.y=F)

## CORRECT CHEMO
datNew$correctADJUVANT_CHEMO <- datNew$revisedADJUVANT_CHEMO
datNew$correctADJUVANT_CHEMO[ is.na(datNew$correctADJUVANT_CHEMO) ] <- datNew$ADJUVANT_CHEMO[ is.na(datNew$correctADJUVANT_CHEMO) ]

## CORRECT RT
datNew$correctADJUVANT_RT <- datNew$revisedADJUVANT_RT
datNew$correctADJUVANT_RT[is.na(datNew$correctADJUVANT_RT)] <- datNew$ADJUVANT_RT[is.na(datNew$correctADJUVANT_RT)]


## AVOID BOTH ADJUVANT THERAPIES
datNew$totalAdjuvant <- ifelse( datNew$correctADJUVANT_CHEMO=="No" & datNew$correctADJUVANT_RT=="No", 0, 1)

## RELOCATE the COLUMN "MICROARRAY" to be the first Column of a new file called DirectorData
DirectorData <- datNew[, c("MICROARRAY", colnames(datNew)[colnames(datNew) != "MICROARRAY"])]

## name OS and BATCH:
OS <- DirectorData$MONTHS_TO_LAST_CONTACT_OR_DEATH
OS <- as.numeric(sub(",", ".", OS))
boxplot(BATCH,OS)


library(survival)
OSs <- ifelse(dat1$VITAL_STATUS=="Alive", 0,1)
plot(survfit(Surv(OS, OSs) ~ BATCH), col=1:length(unique(BATCH)))

table(OS)
OS
