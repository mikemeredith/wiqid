
# Useful for devel versions.
# Take this file out when submitting to CRAN.

.onAttach <- function(libname, pkgname) {
  version <- try(packageVersion('wiqid'), silent=TRUE)
  if(!inherits(version, "try-error"))
    packageStartupMessage("This is wiqid ", version,
      ". For overview type ?wiqid; for changes do news(p='wiqid').")
}
