### Charles Ferté
### Sage Bionetworks
### Seattle, WA

### script for describing the NSCLC datasets
options(stringsAsFactors=FALSE)

## POINT TO DATA DIRECTORY
PATH <- "/Volumes"
folder1 <- "/cferte/FELLOW/cferte/NSCLC_MA/CLIN_DATA_FILES/"
setwd(paste(PATH,folder1,sep=""))
dir()

load("TS_CLIN.Rdata")
load("VS_CLIN.Rdata")
load("VS2_CLIN.Rdata")
load("VS3_CLIN.Rdata")
load("VS4_CLIN.Rdata")

LuscClinF$LABORATORY_BATCH <- paste("Lus",LuscClinF$LABORATORY_BATCH,sep="")

## transform the H_Grdae in VS4_CLIN
LuscClinF$P_Stage <- substr(LuscClinF$P_Stage,7,nchar(LuscClinF$P_Stage))
LuscClinF$P_Stage <- sub("III","3",LuscClinF$P_Stage)
LuscClinF$P_Stage <- sub("II","2",LuscClinF$P_Stage)
LuscClinF$P_Stage <- sub("I","1",LuscClinF$P_Stage)

## trabnsform the Gendre in Lus
LuscClinF$GENDER <- ifelse(LuscClinF$GENDER=="FEMALE","Female","Male")

## merge the tables and create datN
datN <- merge(DirClinF,ZhuClinF,all=T)
datN <- merge(datN,HouClinF,all=T)
#datN <- merge(datN,GirClinF,all=T)
datN <- merge(datN,LuscClinF,all=T)


datN <- apply(datN,2,sub,pattern=",",replacement=".",fixed=T)
datN <- as.data.frame(datN)


## convert in as.numeric all the relevant variables
datN$H_Grade <- as.numeric(datN$H_Grade)
datN$Age <- as.numeric(datN$Age)
datN$VITAL_STATUS <- as.numeric(datN$VITAL_STATUS)
datN$MONTHS_TO_LAST_CONTACT_OR_DEATH <- as.numeric(datN$MONTHS_TO_LAST_CONTACT_OR_DEATH)
datN$REC_STATUS <- as.numeric(datN$REC_STATUS)
datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT <- as.nuemric(datN$TIME_TO_PROGRESSION_OR_LAST_CLIN_ASSESSMENT)
datN$THREE_YEAR_OS <- datN$THREE_YEAR_OS

datN$DATASET <- substr(datN$LABORATORY_BATCH,1,3)

### descriptive statistics

# item name ,item number, nvalid, mean, sd, 
# median, mad, min, max, skew, kurtosis, se

### AGE
summary(datN$Age)
tapply(datN$Age,datN$DATASET,summary)
boxplot(datN$Age~datN$DATASET, ylab="age at diagnosis", xlab="datasets")
library(lattice)
densityplot(datN$Age,groups=datN$DATASET,auto.key==list(col=4),main="age distribution according to the different datasets")
## HISTOLOGY
TOTAL <- table(datN$Histology)
HistD <- table(datN$Histology,datN$DATASET)
print(HistD)
HISTOLOGY  <- as.data.frame(cbind(HistD,Total))
chisq.test(HistD)

## Histological Grade
Total <- table(datN$H_Grade)
HGradeD <- table(datN$H_Grade, datN$DATASET)
print(HGradeD)
H_GRADE <- as.data.frame(cbind(HGradeD,Total))
rm(Total,HGradeD)
chisq.test(H_GRADE)

## p Stage
Total <- table(datN$P_Stage)
PSTAGED <- table(datN$P_Stage,datN$DATASET)
print(PSTAGED)
PSTAGE <- as.data.frame(cbind(PSTAGED,Total))
rm(Total,PSTAGED)
chisq.test(PSTAGE)

## GENDER
Total <- table(datN$GENDER)
GENDERD <- table(datN$GENDER,datN$DATASET)
print(GENDERD)
NGENDER <- as.data.frame(cbind(GENDERD,Total))
rm(Total,GENDERD)
chisq.test(NGENDER)

## RACE
RACED <- table(datN$RACE,datN$DATASET)
print(RACED)

## SMOKING
SMOKED <-table(datN$SMOKING,datN$DATASET) 
print(SMOKED)

## Survival OS according to the 5 datasets
library(survival)


OSs <- datN$VITAL_STATUS
OS <- as.numeric(sub(",", ".",datN$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ datN$DATASET),col=c("royalblue","black","brown1","aquamarine4"),main="OS according to the 4 Datasets",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~datN$DATASET)
survdiff(Surv(OS,OSs)~datN$DATASET)
summary(survfit(Surv(OS,OSs)~datN$DATASET), times=36)
summary(survfit(Surv(OS,OSs)~datN$DATASET), times=60)

## Survival OS according to the adjuvant status
library(survival)
OSs <- ifelse(datN$VITAL_STATUS=="Alive", 0,1)
OS <- as.numeric(sub(",", ".",datN$MONTHS_TO_LAST_CONTACT_OR_DEATH))
plot(survfit(Surv(OS, OSs) ~ datN$DATASET), col=1:length(unique(datN$DATASET)),main="OS according to the adjuvant status",xlab="months",ylab="Overall Survival (%)")
survfit(Surv(OS,OSs)~datN$totalAdjuvant)
summary(survfit(Surv(OS,OSs)~datN$DATASET), times=36)
summary(survfit(Surv(OS,OSs)~datN$DATASET), times=60)

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
