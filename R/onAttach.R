

.onAttach <- function(libname, pkgname) {
  version <- try(packageVersion('wiqid'), silent=TRUE)
  if(!inherits(version, "try-error"))
    packageStartupMessage("This is wiqid ", version,
      ". For overview type ?wiqid.")
}
