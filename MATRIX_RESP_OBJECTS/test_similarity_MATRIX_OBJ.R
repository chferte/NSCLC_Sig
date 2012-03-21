## charles ferté
## sage bionetworks
## 29 fev 2012



# tester si les Gir_NORM de GCRMA et de MAS5 sont les memes...
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/Gir_NORMALIZED/")
dir()
load("Gir_GCRMA.Rdata")
load("Gir_MAS5.Rdata")
load("Gir_metaGEO.Rdata")
load("Gir_RMA.Rdata")
load("Gir_dCHIP.Rdata")

a <- exprs(Gir_GCRMA)
b <- exprs(Gir_MAS5)
c <- exprs(Gir_RMA)
d <- exprs(Gir_metaGEO$hgu133plus2)
e <- exprs(Gir_dCHIP)

identical(a,b)
identical(a,c)
identical(a,d)
identical(a,e)

dim(a)
dim(b)
dim(c)
dim(d)
dim(e)


# tester si les MATRIX_VS3 de GCRMA et de MAS5 sont les memes...
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/GCRMA/")
load("MATRIX_VS4.Rdata")
aa <- MATRIX_VS4
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/MAS5/")
load("MATRIX_VS4.Rdata")
bb <- MATRIX_VS4
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/RMA/")
load("MATRIX_VS4.Rdata")
cc <- MATRIX_VS4
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/metaGEO/")
load("MATRIX_VS4.Rdata")
dd <- MATRIX_VS4
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/MATRIX_RESP_OBJECTS/dCHIP/")
load("MATRIX_VS4.Rdata")
ee <- MATRIX_VS4
identical(aa,bb)
identical(aa,cc)
identical(aa,dd)
identical(aa,ee)

dim(aa)
dim(bb)
dim(cc)
dim(dd)
dim(ee)

