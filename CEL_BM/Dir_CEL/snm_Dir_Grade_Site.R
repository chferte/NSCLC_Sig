### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### March, 2 2012

# code to normalize the Dir using snm using site as adj var and removing or not the grade
# then to compare the different the performance of the models

#load the libraries
library(metaGEO)
library(snm)
library(lattice)

## set the CEL files directory
PATH <- "/home"
setwd(paste(PATH,"/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/Dir_CEL",sep=""))
load("TS_CLIN.Rdata")

#Run this code
fits <- runWorkflow("/home/cferte/FELLOW/cferte/NSCLC_MA/CEL_BM/Dir_CEL/",workflow="snm")
dat <- exprs(fits$hgu133a)

ids <- which(!is.na(DirClinF$H_Grade))
bio.var <- model.matrix(~ DirClinF$H_Grade[ids])
adj.var <- model.matrix(~ DirClinF$SITE[ids])
snm.fit <- snm(dat[,ids], 
               bio.var=bio.var, 
               adj.var=adj.var, 
               rm.adj=TRUE)

new.dat <- snm.fit$norm.dat

# Note we lost a few samples as they did not have grade

# now do your cross validation on the new.dat and dat.