# Restart R before testing reverse dependencies

setwd("D:/Github/wiqid_package/reverse_dependencies")

# Install new version of 'wiqid' first
devtools::install_github("mikemeredith/wiqid")
library(wiqid)
packageVersion("wiqid")

# Find which packages are reverse-depends
library(devtools)
( rds <- revdep("wiqid", ignore=c("BEST", "wiqid")) )
# Install them to be sure we have all their dependencies, incl. Suggests
install.packages(rds, dependencies=TRUE)
download.packages(rds, destdir=".", type="source")

( totest <- dir(pattern = ".tar.gz$") )

# Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = FALSE) # not needed if Suggests were installed.
tstcall <- paste("R CMD check", totest)
isok <- rep(NA, length(tstcall))
for(i in seq_along(tstcall)) {
  cat("\n\n***** ", totest[i], "*****\n\n")
  isok[i] <- system(tstcall[i])
}
which(isok != 0)
totest[isok != 0]
sessionInfo()
