# Zhu_Nrma

# Charles Ferté

# load affy
library(affy)

# set Zhu CEL files directory
setwd("/Volumes/cferte/NSCLC_Prognostic_Signature/Zhu_2/CEL_files/Zhu40_CELfiles/")

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

# combine "AffyID","Probe", and "pm.bcg.q" into a new file called "Zhu_Nrma", which contains the  normalized PM data
Zhu_Nrma  <- cbind (AffyID, probe,pm.bcg.q)

# to have athe expression measure in an excel file readable format, save the "Zhu_Nrma" file as a .csv file
write.table(Zhu_Nrma, file="Zhu.Nrma.csv",sep=",",row.names=FALSE,quote=FALSE)
