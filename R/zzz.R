.onAttach <- function(libname, pkgname) {

    ver <- "0.0.0.9000"

    loadmsg <- paste0("\nLoading the 'BACE' package (version ", ver, "). For an\nintroduction and vignette to the package please see: https://daniel1noble.github.io/BACE/")

    packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)
    
    # Print citation information
    packageStartupMessage("\nTo cite BACE in publications, use:\n\n", domain=NULL, appendLF=FALSE)
    cit <- utils::citation("BACE")
    
    # Format each citation with proper spacing
    for (i in seq_along(cit)) {
      packageStartupMessage(format(cit[i], style = "text"), domain=NULL, appendLF=FALSE)
      if (i < length(cit)) {
        packageStartupMessage("\n\n", domain=NULL, appendLF=FALSE)
      }
    }
    packageStartupMessage("", domain=NULL, appendLF=TRUE)

  }
