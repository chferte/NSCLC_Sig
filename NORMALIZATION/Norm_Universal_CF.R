### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### 2012-01-30

### script for normalizing a file using the following methods: RMA, GCRMA, MAS5, dCHIP, Barcode, fRMA and metaGEO

# load affy
library(affy)
library(gcrma)
library(frma)
library(hgu133plus2frmavecs)
library(hgu133afrmavecs)
library(metaGEO)
library(snm)

## select dataset among Dir,Zhu, Gir,Hou and LUSC
dataset <- "Gir"

## where are the CEL files

Gir_CEL <- "/home/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Girard_GSE31908/CEL_Files/GPL570_HGU133plus2/"
Dir_CEL <- "/home/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Director_challenge/CEL_Files/caArray_jacob-00182_files-AFFYMETRIX_CEL/"
Zhu_CEL <- "/home/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Zhu_GSE14814/CEL_files/compressed_Zhu_CEL_files/"
Hou_CEL <- "/home/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Hou_GSE19188/CEL_files/GSE19188/"
LUSC_CEL <- "/home/cferte/temp_synapse/"

datapath1  <-get(paste(dataset,"_CEL",sep=""))

## where should be stored the normalized files
datapath2  <- paste("/home/cferte/FELLOW/cferte/NSCLC_MA/",dataset,"_NORMALIZED/",sep="")

# set the CEL files directory
setwd(datapath1)

#read CEL files as Affy objects
rawdata <- ReadAffy()


## set the chosen dataset as directory
setwd(datapath2)

## boxplot them all !
png(paste("boxplot_",dataset,"raw_data_expr.png",sep=""))
boxplot(rawdata)
dev.off()

## COMPUTE RMA & make a file
tmp <- paste(dataset,"_RMA",sep="")
assign(tmp,rma(rawdata))
save(list=paste(dataset,"_RMA",sep=""), file=paste(dataset,"_RMA",".Rdata",sep=""))

## COMPUTE GCRMA & make a file
tmp <- paste(dataset,"_GCRMA",sep="")
assign(tmp,gcrma(rawdata))
save(list=paste(dataset,"_GCRMA",sep=""), file=paste(dataset,"_GCRMA",".Rdata",sep=""))

## COMPUTE dCHIP & make a file
tmp <- paste(dataset,"_dCHIP",sep="")
assign(tmp,expresso(rawdata, normalize.method = "invariantset", bg.correct = FALSE, pmcorrect.method = "pmonly",summary.method = "liwong"))
save(list=paste(dataset,"_dCHIP",sep=""), file=paste(dataset,"_dCHIP",".Rdata",sep=""))

## COMPUTE METAGEO & make a file
setwd(datapath1)
tmp <- paste(dataset,"_metaGEO",sep="")
assign(tmp,runWorkflow("."))
setwd(datapath2)
save(list=paste(dataset,"_metaGEO",sep=""), file=paste(dataset,"_metaGEO",".Rdata",sep=""))

## COMPUTE MAS5.0 & make a file
tmp <- paste(dataset,"_MAS5",sep="")
assign(tmp,mas5(rawdata))
save(list=paste(dataset,"_MAS5",sep=""), file=paste(dataset,"_MAS5",".Rdata",sep=""))

## COMPUTE fRMA & make a file
tmp <- paste(dataset,"_fRMA",sep="")
assign(tmp,frma(rawdata, summarize = "random_effect"))
save(list=paste(dataset,"_fRMA",sep=""), file=paste(dataset,"_fRMA",".Rdata",sep=""))

## COMPUTE BARCODE & make a file
tmp1 <- paste(dataset,"_barcode",sep="")
assign(tmp1,barcode(get(tmp)))
save(list=paste(dataset,"_barcode",sep=""), file=paste(dataset,"_barcode",".Rdata",sep=""))
rm(tmp1)


# END# 