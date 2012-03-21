### Charles Ferté
### Sage Bionetworks
### Seattle, WA
### January, 9th 2012

### script for plotting the eigengenes of the normalized expression according to each dataset

#load affy
library(affy)

# set Dir CEL files directory
setwd("/home/cferte/NSCLC_Prognostic_Signature/Director_challenge/CEL_file/")
Dir <- ReadAffy()
Dir_pm1000 <- pm(Dir)[1:1000,]

# set Dir CEL files directory
setwd("/home/cferte/NSCLC_Prognostic_Signature/Zhu_GSE14814/CEL_files/compressed_Zhu_CEL_files/")
Zhu <- ReadAffy()
Zhu_pm1000 <- pm(Zhu)[1:1000,]

# cbin the two tables
DirZhu <- cbind( Dir_pm1000, Zhu_pm1000)
dim(DirZhu)
setwd("/home/cferte/FELLOW/cferte/")
A <- cor(DirZhu,DirZhu)
dim(A)