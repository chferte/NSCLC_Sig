# Girard NSCLC Dataset Normalization rma

# Charles Ferté
# Sage Bionetworks
# Seattle, WA

# load affy
library(affy)

# this datset is separated into two different parts
# 1/ GPL96_HGU133A 17 samples
# 2/ GPL570_HGU133plus2 50 samples

#1/ Normaliszation of the first dataset (GPL96_HGU133A 17 samples)
# set Girard GPL96_U133A CEL files directory
setwd("/Volumes/cferte/NSCLC_Prognostic_Signature/Girard_GSE31908/CEL_Files/GPL96_U133A/")

#read CEL files as Affy objects
rawdata <- ReadAffy()

# Extract perfect match intensities and store it into a new "PM"" file 
PM <- probes(rawdata,which="pm")

# Retrieve the Gene ID from the first column of the PM file and save it into a "AffyInfo" new file
AffyInfo <- dimnames(PM)[[1]]

# looking for the number of digits following each probe name within AffyInfo data set and save it into "cutpos file". Gives a -1 if no digit  
cutpos <- regexpr("\\d+$",AffyInfo,perl=T)

# extract the digits following the probe names in AffyInfo and save these IDs as "AffyID"
AffyID <- substr(AffyInfo,1,cutpos-1)

# take numeric objects from the "AffyInfo" dataset, which are the probe IDs (most are 1-11) and save them as "probe" file
probe <- as.numeric(substr(AffyInfo,cutpos,nchar(AffyInfo)))

# Raw data background corrects probe intensity values using RMA method 
data.bcg <- bg.correct(rawdata,method="rma")

# quantile normalization of the perfect match probe intensities
data.bcg.q <- normalize.AffyBatch.quantiles(data.bcg,type="pmonly")

# extract the perfect match intensities of the data.bcg.q file into a new file called "pm.bcg.q"
pm.bcg.q <- probes (data.bcg.q,which="pm")

# combine "AffyID","Probe", and "pm.bcg.q" into a new file called "Girard_GPL96_Nrma", which contains the  normalized PM data
Girard_GPL96_Nrma  <- cbind (AffyID, probe,pm.bcg.q)

# to have the expression measure in an excel file readable format, save the "Girard_GPL96_Nrma" file as a .csv file
write.table(Girard_GPL96_Nrma, file="Girard_GPL96_Nrma.csv",sep=",",row.names=FALSE,quote=FALSE)

#2/ Normaliszation of the first dataset (GPL570_HGU133plus2 50 samples)
# set Girard GPL570_HGU133plus2 CEL files directory
setwd("/Volumes/cferte/NSCLC_Prognostic_Signature/Girard_GSE31908/CEL_Files/GPL570_HGU133plus2/")

#read CEL files as Affy objects
rawdata <- ReadAffy()

# Extract perfect match intensities and store it into a new "PM"" file 
PM <- probes(rawdata,which="pm")

# Retrieve the Gene ID from the first column of the PM file and save it into a "AffyInfo" new file
AffyInfo <- dimnames(PM)[[1]]

# looking for the number of digits following each probe name within AffyInfo data set and save it into "cutpos file". Gives a -1 if no digit  
cutpos <- regexpr("\\d+$",AffyInfo,perl=T)

# extract the digits following the probe names in AffyInfo and save these IDs as "AffyID"
AffyID <- substr(AffyInfo,1,cutpos-1)

# take numeric objects from the "AffyInfo" dataset, which are the probe IDs (most are 1-11) and save them as "probe" file
probe <- as.numeric(substr(AffyInfo,cutpos,nchar(AffyInfo)))

# Raw data background corrects probe intensity values using RMA method 
data.bcg <- bg.correct(rawdata,method="rma")

# quantile normalization of the perfect match probe intensities
data.bcg.q <- normalize.AffyBatch.quantiles(data.bcg,type="pmonly")

# extract the perfect match intensities of the data.bcg.q file into a new file called "pm.bcg.q"
pm.bcg.q <- probes (data.bcg.q,which="pm")

# combine "AffyID","Probe", and "pm.bcg.q" into a new file called "Girard_GPL570_Nrma", which contains the  normalized PM data
Girard_GPL570_Nrma  <- cbind (AffyID, probe,pm.bcg.q)

# to have the expression measure in an excel file readable format, save the "Girard_GPL570_Nrma" file as a .csv file
write.table(Girard_GPL570_Nrma, file="Girard_GPL570_Nrma.csv",sep=",",row.names=FALSE,quote=FALSE)
