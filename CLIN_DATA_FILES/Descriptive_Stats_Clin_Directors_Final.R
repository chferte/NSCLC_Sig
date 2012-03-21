### Charles Ferté
### Sage Bionetworks
### Seattle, WA

### script for describing the NSCLC datasets

## POINT TO DATA DIRECTORY
dataPath <- "/volumes/cferte/FELLOW/cferte/"

## READ IN THE TAB DELIMITED FILES
dat1 <- read.csv(file = paste(dataPath, "DirectorData.csv", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
dat2 <- read.csv(file = paste(dataPath, "ZhuData.csv", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
dat4 <- read.csv(file = paste(dataPath, "GirardData.csv", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")
dat3 <- read.csv(file = paste(dataPath, "HouData.csv", sep=""), header=T, na.strings=c("NA", "", " "), as.is=T,dec=",")

## lump together batches according to the closer dates in scanbatch in dat1
table(dat1$ScanDate)
dat1$COR_BATCH [dat1$ScanDate=="2005-06-03"]<- "1"
dat1$COR_BATCH [dat1$ScanDate %in% c("2004-11-19","2004-10-08","2004-11-09","2004-11-10","2004-11-16")]<- "2"
dat1$COR_BATCH [dat1$ScanDate %in% c("2004-07-14","2004-07-16","2004-07-01","2004-07-02","2004-07-09","2004-07-13")]<- "3"
dat1$COR_BATCH [dat1$ScanDate %in% c("2004-05-05","2004-05-07","2004-05-12","2004-05-14","2004-05-18","2004-05-20","2004-05-25","2004-06-08","2004-06-09","2004-06-16")]<- "4"
dat1$COR_BATCH [dat1$ScanDate %in% c("2004-03-31","2004-04-02","2004-04-06","2004-04-08","2004-04-13")]<- "5"
dat1$COR_BATCH [dat1$ScanDate %in% c("2004-01-28","2004-01-30","2004-02-04")]<- "6"
table(dat1$COR_BATCH)

## merge the tables and create datN
datN <- merge(dat1,dat2,all=T)
datN <- merge(datN,dat4,all=T)
datN <- merge(datN,dat3,all=T)
datN <- as.data.frame(datN)


## cut off the 3B and IV p Stage
datN <- datN[datN$P_Stage %in% c("1A","1B","2A","2B","3A"), ]

## create a new variable "DATASET"
datN$DATASET <- substr(datN$LABORATORY_BATCH,1,3)

datN$Age <- as.numeric(datN$Age)
datN$MONTHS_TO_LAST_CONTACT_OR_DEATH <- as.numeric(datN$MONTHS_TO_LAST_CONTACT_OR_DEATH)
datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- as.numeric(datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT)


### descriptive statistics

# item name ,item number, nvalid, mean, sd, 
# median, mad, min, max, skew, kurtosis, se


## HISTOLOGY
Total <- table(datN$Histology)
HistD <- table(datN$Histology,datN$DATASET)
HISTOLOGY  <- as.data.frame(cbind(HistD,Total))
chisq.test(HistD)
fisher.test(HistD,hybrid=T)

## Histological Grade
Total <- table(datN$H_Grade)
HGradeD <- table(datN$H_Grade, datN$DATASET)
H_GRADE <- as.data.frame(cbind(HGradeD,Total))
rm(Total,HGradeD)
chisq.test(H_GRADE)

## p Stage
Total <- table(datN$P_Stage)
PSTAGED <- table(datN$P_Stage,datN$DATASET)
PSTAGE <- as.data.frame(cbind(PSTAGED,Total))
rm(Total,PSTAGED)
chisq.test(PSTAGE)

## GENDER
Total <- table(datN$GENDER)
GENDERD <- table(datN$GENDER,datN$DATASET)
NGENDER <- as.data.frame(cbind(GENDERD,Total))
rm(Total,GENDERD)
chisq.test(NGENDER)

## Survival OS according to the 4 datasets
library(survival)
datN$VITAL_STATUS <- ifelse(datN$VITAL_STATUS %in% c("Alive"," Alive"), 0,1)
OSs <- ifelse(datN$VITAL_STATUS %in% c("Alive"," Alive"), 0,1)
OS <- as.numeric(sub(",", ".",datN$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ datN$DATASET), col=1:length(unique(datN$DATASET)),main="OS according to the 4 Datasets",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~datN$DATASET)
summary(survfit(Surv(OS,OSs)~DATASET), times=36)


## Survival OS according to the adjuvant status
library(survival)
OSs <- ifelse(datN$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".",datN$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ datN$totalAdjuvant), col=1:length(unique(datN$totalAdjuvant)),main="OS according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~datN$totalAdjuvant)
summary(survfit(Surv(OS,OSs)~totalAdjuvant), times=36)

## Survival OS according to the P_STAGE
library(survival)
OSs <- ifelse(datN$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".",datN$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ datN$P_Stage), col=1:length(unique(datN$P_Stage)),main="OS according to the Pathological Stage",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~datN$P_Stage)
summary(survfit(Surv(OS,OSs)~PSTAGE), times=36)

par(mfrow= c(1,2))
plot(survfit(Surv(OS[datN$totalAdjuvant==0], OSs[datN$totalAdjuvant==0]) ~ datN$P_Stage[datN$totalAdjuvant==0]), col=1:length(unique(datN$P_Stage)),main="OS ~ P-Stage (OBS) ",xlab="months",ylab="Overall Survival (%)")
plot(survfit(Surv(OS[datN$totalAdjuvant==1], OSs[datN$totalAdjuvant==1]) ~ datN$P_Stage[datN$totalAdjuvant==1]), col=1:length(unique(datN$P_Stage)),main="OS ~ P-Stage (ACT) ",xlab="months",ylab="Overall Survival (%)")

## compute and plot OS according to the adjuvant status and each dataset
par(mfrow=c(2,2))
plot(survfit(Surv(OS[datN$DATASET=="Dir"],OSs[datN$DATASET=="Dir"])~datN$totalAdjuvant[datN$DATASET=="Dir"]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Director's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
plot(survfit(Surv(OS[datN$DATASET=="Zhu"],OSs[datN$DATASET=="Zhu"])~datN$totalAdjuvant[datN$DATASET=="Zhu"]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Zhu's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
plot(survfit(Surv(OS[datN$DATASET=="Hou"],OSs[datN$DATASET=="Hou"])~datN$totalAdjuvant[datN$DATASET=="Hou"]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Hou's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
plot(survfit(Surv(OS[datN$DATASET=="Gir"],OSs[datN$DATASET=="Gir"])~datN$totalAdjuvant[datN$DATASET=="Gir"]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Girard's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")

survfit(Surv(OS[datN$DATASET=="Dir"],OSs[datN$DATASET=="Dir"])~datN$totalAdjuvant[datN$DATASET=="Dir"]) 
survfit(Surv(OS[datN$DATASET=="Zhu"],OSs[datN$DATASET=="Zhu"])~datN$totalAdjuvant[datN$DATASET=="Zhu"]) 
survfit(Surv(OS[datN$DATASET=="Hou"],OSs[datN$DATASET=="Hou"])~datN$totalAdjuvant[datN$DATASET=="Hou"])
survfit(Surv(OS[datN$DATASET=="Gir"],OSs[datN$DATASET=="Gir"])~datN$totalAdjuvant[datN$DATASET=="Gir"])

## Survival OS according to Age65 (>65 = 1 ; <65=0)
library(survival)
datN$Age65 <- ifelse(datN$Age<65,0,1)
OSs <- ifelse(datN$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".",datN$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ datN$Age65), col=1:length(unique(datN$Age65)),main="OS according to Age > 65 y",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~datN$Age65)
summary(survfit(Surv(OS,OSs)~datN$Age65), times=36)

## Survival OS according to Age65 (>65 = 1 ; <65=0) and to the ACT status in the director's challenge'
plot(survfit(Surv(OS[datN$DATASET=="Dir" & datN$Age65==1],OSs[datN$DATASET=="Dir" & datN$Age65==1])~datN$totalAdjuvant[datN$DATASET=="Dir" & datN$Age65==1]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Director's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
plot(survfit(Surv(OS[datN$DATASET=="Dir" & datN$Age65==0],OSs[datN$DATASET=="Dir" & datN$Age65==0])~datN$totalAdjuvant[datN$DATASET=="Dir" & datN$Age65==0]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Director's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS[datN$DATASET=="Dir"& datN$Age65==1],OSs[datN$DATASET=="Dir"& datN$Age65==1])~datN$totalAdjuvant[datN$DATASET=="Dir"& datN$Age65==1])
survfit(Surv(OS[datN$DATASET=="Dir"& datN$Age65==0],OSs[datN$DATASET=="Dir"& datN$Age65==0])~datN$totalAdjuvant[datN$DATASET=="Dir"& datN$Age65==0])


## compute and plot OS according to the laboratory batch in each dataset
plot(survfit(Surv(OS[datN$DATASET=="Dir"],OSs[datN$DATASET=="Dir"])~datN$LABORATORY_BATCH[datN$DATASET=="Dir"]), col=1:length(unique(datN$LABORATORY_BATCH)), main="OS in the Director's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS[datN$DATASET=="Dir"],OSs[datN$DATASET=="Dir"])~datN$LABORATORY_BATCH[datN$DATASET=="Dir"])

plot(survfit(Surv(OS[datN$DATASET=="Zhu"],OSs[datN$DATASET=="Zhu"])~datN$totalAdjuvant[datN$DATASET=="Zhu"]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Zhu's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
plot(survfit(Surv(OS[datN$DATASET=="Hou"],OSs[datN$DATASET=="Hou"])~datN$totalAdjuvant[datN$DATASET=="Hou"]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Hou's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
plot(survfit(Surv(OS[datN$DATASET=="Gir"],OSs[datN$DATASET=="Gir"])~datN$totalAdjuvant[datN$DATASET=="Gir"]), col=1:length(unique(datN$totalAdjuvant)), main="OS in the Girard's dataset according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")


survfit(Surv(OS[datN$DATASET=="Zhu"],OSs[datN$DATASET=="Zhu"])~datN$totalAdjuvant[datN$DATASET=="Zhu"]) 
survfit(Surv(OS[datN$DATASET=="Hou"],OSs[datN$DATASET=="Hou"])~datN$totalAdjuvant[datN$DATASET=="Hou"])
survfit(Surv(OS[datN$DATASET=="Gir"],OSs[datN$DATASET=="Gir"])~datN$totalAdjuvant[datN$DATASET=="Gir"])

## evaluate if ACT/OBS is correlated with any covariates
table(datN$Histology,datN$totalAdjuvant)
table(datN$H_Grade,datN$totalAdjuvant)
table(datN$P_Stage,datN$totalAdjuvant)
table(datN$GENDER,datN$totalAdjuvant)
table(datN$Age65,datN$totalAdjuvant)
table(datN$SMOKING,datN$totalAdjuvant)
table(datN$RACE,datN$totalAdjuvant)


chisq.test(datN$Histology,datN$totalAdjuvant)
chisq.test(datN$H_Grade,datN$totalAdjuvant)
chisq.test(datN$P_Stage,datN$totalAdjuvant)
chisq.test(datN$GENDER,datN$totalAdjuvant)
chisq.test(datN$Age65,datN$totalAdjuvant)
chisq.test(datN$SMOKING,datN$totalAdjuvant)
chisq.test(datN$RACE,datN$totalAdjuvant)

barplot(table(datN$Histology,datN$totalAdjuvant), main="Histology")
barplot(table(datN$H_Grade,datN$totalAdjuvant),main="Histological Grade")
barplot(table(datN$P_Stage,datN$totalAdjuvant),main="pathological stage")
barplot(table(datN$GENDER,datN$totalAdjuvant), main="Gender")
barplot(table(datN$Age65,datN$totalAdjuvant),main="Age >65)")
barplot(table(datN$SMOKING,datN$totalAdjuvant),main="smoking status")
barplot(table(datN$RACE,datN$totalAdjuvant),main="Race")

par(mfrow=c(1,2))
barplot(table(datN$Histology[datN$totalAdjuvant==0]),main="Histology (OBS)")
barplot(table(datN$Histology[datN$totalAdjuvant==1]),main="Histology (ACT)")

par(mfrow=c(1,2))
barplot(table(datN$H_Grade[datN$totalAdjuvant==0]),main="Histological Grade (OBS)")
barplot(table(datN$H_Grade[datN$totalAdjuvant==1]),main="Histological Grade (ACT)")

par(mfrow=c(1,2))
barplot(table(datN$P_Stage[datN$totalAdjuvant==0]),main="Pathological stage(OBS)")
barplot(table(datN$P_Stage[datN$totalAdjuvant==1]),main="Pathological stage(ACT)")

par(mfrow=c(1,2))
barplot(table(datN$GENDER[datN$totalAdjuvant==0]),main="Gender (OBS)")
barplot(table(datN$GENDER[datN$totalAdjuvant==1]),main="Gender (ACT)")

par(mfrow=c(1,2))
barplot(table(datN$Age65[datN$totalAdjuvant==0]),main="Age >65 (OBS)")
barplot(table(datN$Age65[datN$totalAdjuvant==1]),main="Age >65 (ACT)")

par(mfrow=c(1,2))
barplot(table(datN$SMOKING[datN$totalAdjuvant==0]),main="SMOKING (OBS)")
barplot(table(datN$SMOKING[datN$totalAdjuvant==1]),main="SMOKING (ACT)")

## Survival OS according to the CORBATCH in dat1
library(survival)
dat1$VITAL_STATUS <- ifelse(dat1$VITAL_STATUS %in% c("Alive"," Alive"), 0,1)
OSs <- ifelse(dat1$VITAL_STATUS %in% c("Alive"," Alive"), 0,1)
OS <- as.numeric(sub(",", ".",dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ dat1$COR_BATCH), col=1:length(unique(dat1$COR_BATCH)),main="OS according to the CORBATCH",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~dat1$COR_BATCH)
summary(survfit(Surv(OS,OSs)~dat1$COR_BATCH), times=36)


## Survival OS according to the SITES in dat1
library(survival)
dat1$VITAL_STATUS <- ifelse(dat1$VITAL_STATUS %in% c("Alive"," Alive"), 0,1)
OSs <- ifelse(dat1$VITAL_STATUS %in% c("Alive"," Alive"), 0,1)
OS <- as.numeric(sub(",", ".",dat1$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ dat1$SITE), col=1:length(unique(dat1$SITE)),main="OS according to the SITES",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~dat1$SITE)
summary(survfit(Surv(OS,OSs)~dat1$SITE), times=36)
