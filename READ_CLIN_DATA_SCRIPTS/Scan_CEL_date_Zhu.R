library(affy)

celPath <- "/Volumes/cferte/NSCLC_Prognostic_Signature/Zhu_14814/CEL_files/compressed_Zhu_CEL_files/"

setwd(celPath)

dir()

abatch <- ReadAffy()

fns <- list.celfiles(path=celPath, listGzipped=TRUE, full.names=T)

fns

blah <- ReadAffy(filenames=fns)
@blah@protocolData@data
#blah@protocolData$ScanDate