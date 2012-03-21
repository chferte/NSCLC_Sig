datacache <- new.env(hash=TRUE, parent=emptyenv())

mg.hgu133plus2 <- function() showQCData("mg.hgu133plus2", datacache)
mg.hgu133plus2_dbconn <- function() dbconn(datacache)
mg.hgu133plus2_dbfile <- function() dbfile(datacache)
mg.hgu133plus2_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
mg.hgu133plus2_dbInfo <- function() dbInfo(datacache)

mg.hgu133plus2ORGANISM <- "Homo sapiens"

.onLoad <- function(libname, pkgname)
{
    require("methods", quietly=TRUE)
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "mg.hgu133plus2.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    txdb <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"ChipDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, txdb, envir=ns)
    namespaceExport(ns, dbNewname)
        
    ## Create the AnnObj instances
    ann_objs <- createAnnObjs.SchemaChoice("HUMANCHIP_DB", "mg.hgu133plus2", "chip mg.hgu133plus2", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("mg.hgu133plus2.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(mg.hgu133plus2_dbconn())
}

