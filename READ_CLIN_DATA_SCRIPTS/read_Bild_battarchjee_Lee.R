## read bild / battachardee


setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Lu_GSE6253/publication/")
dat3 <- read.delim2(file="Lu_Clin_Data.txt",header=T, na.strings=c("NA", "", " "), skip=2,
                    as.is=T,dec=",")
colnames(dat3)
table(dat3$Institute.Author)
Bild <- dat3[dat3$Institute.Author=="Bild et al.",]
Bhattacharje <- dat3[dat3$Institute.Author=="Bhattacharjee et al.",]

## browse by histology and event:
table(Bild$Cancer.Death,Bild$Pathologya)

## read Lee et al,
setwd("/Volumes/cferte/FELLOW/cferte/NSCLC_MA/NSCLC_Prognostic_Signature/Lee_GSE8894/publication/")
Lee <- read.delim2(file="s138_clinic_info_111115_toFerte_updated200704.txt",header=T, na.strings=c("NA", "", " "),
                    as.is=T,dec=",")
LeeObs <- Lee[Lee$Adjuvant_Cx=="No",]
table(LeeObs$Histology)
table(LeeObs$status.for.event..survival.,LeeObs$Histology)
