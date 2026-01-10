.onAttach <- function(libname, pkgname) {

    ver <- "0.0.0.9000"

    loadmsg <- paste0("\nLoading the 'BACE' package (version ", ver, "). For an\nintroduction and vignette to the package please see: https://daniel1noble.github.io/BACE/\n
    
    To cite the package, use: citation('BACE')\n")

    packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)

  }
