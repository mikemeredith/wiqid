

.onAttach <- function(libname, pkgname) {
  version <- packageVersion('wiqid')
  packageStartupMessage("This is wiqid ", version,
    ". For overview type ?wiqid.")
}
