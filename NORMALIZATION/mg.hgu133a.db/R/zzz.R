datacache <- new.env(hash=TRUE, parent=emptyenv())

mg.hgu133a <- function() showQCData("mg.hgu133a", datacache)
mg.hgu133a_dbconn <- function() dbconn(datacache)
mg.hgu133a_dbfile <- function() dbfile(datacache)
mg.hgu133a_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
mg.hgu133a_dbInfo <- function() dbInfo(datacache)

mg.hgu133aORGANISM <- "Homo sapiens"

.onLoad <- function(libname, pkgname)
{
    require("methods", quietly=TRUE)
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "mg.hgu133a.sqlite", package=pkgname, lib.loc=libname)
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
    ann_objs <- createAnnObjs.SchemaChoice("HUMANCHIP_DB", "mg.hgu133a", "chip mg.hgu133a", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("mg.hgu133a.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(mg.hgu133a_dbconn())
}

