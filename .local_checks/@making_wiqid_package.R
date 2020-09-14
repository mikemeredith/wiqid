
# Set directory
# =============
# Choose one:
setwd("C:/Github/wiqid_package") # laptop
setwd("D:/Github/wiqid_package") # desktop
setwd("../..")
dir()  # check...

# Dependencies
# ============
install.packages(c("truncnorm", "coda", "plotrix", "mcmcOutput"))  # 'strong' dependencies
install.packages(c("secr", "shiny", "rjags"))  # Suggests

# Spell check
# ===========
library(spelling)
( out <- spell_check_package(pkg = "wiqid") )
update_wordlist(pkg = "wiqid", confirm = TRUE)

# Development
# ===========
devtools::load_all("C:/GitHub/wiqid_package/wiqid")
system("R CMD INSTALL wiqid") # Use this for a "dev" install.
# test_dir("D:/GitHub/wiqid_package/wiqid/inst/tests/testthat", reporter=ProgressReporter)

# Create the wiqid package
# ========================
unlink(list.files(pattern="Rplots.pdf", recursive=TRUE))
system("R CMD build wiqid")  # Produces the .tar.gz file
pkg <- "wiqid_0.3.0.9000.tar.gz"   # <---- fix version number here

# Check without Suggests packages
# -------------------------------
remove.packages(c("secr", "shiny", "rjags")) # must do before attaching.
# Install Depends and Imports but NOT Suggests
install.packages(c("HDInterval", "truncnorm", "coda"))
Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = FALSE)
system(paste("R CMD check", pkg, "--as-cran"))  # should give NOTEs, no ERRORs.
library(testthat)
test_dir("wiqid/inst/tests/testthat", reporter=ProgressReporter)  # should be no errors

# Now install rjags, secr and shiny and recheck
install.packages(c("rjags", "secr", "shiny"))
Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = TRUE)
system(paste("R CMD check", pkg, "--as-cran"))

# Standard checks
# ---------------
## Pick one to check:
# For laptop (no LaTeX installed)
system(paste("R CMD check", pkg, "--no-manual"))
system(paste("R CMD check", pkg, "--as-cran --no-manual"))

# For desktop
system(paste("R CMD check", pkg))
system(paste("R CMD check", pkg, "--as-cran"))

# Installation
# ------------
## Pick one to install
system(paste("R CMD INSTALL", pkg, "--build"))
system(paste("R CMD INSTALL", pkg))  # use this for R devel

# Test it
# =======
library(testthat)
test_package("wiqid", reporter=ProgressReporter)


# Try it out:
rm(list=ls())
library(wiqid)
?wiqid

data(railSims)
str(railSims) # look at logArea

example(Bsecr0) # takes 12 mins
plotACs(Bout)
plotACs(Bout, 1:10)
plotACs(Bout, 1:10, showLabels=FALSE)
plotACs(Bout, 1:10, howMany=1e6)
plotACs(Bout, NA, howMany=1e6)



# Test showShinyApp (not covered by automated checks)
showShinyApp()
showShinyApp("Beta")
showShinyApp("G")
showShinyApp("Quad")


# Run the examples:
# example("wiqid-package")
example("dgamma2")
example("dbeta2")
example("dt2")
example("closedCapM0")
example("occSS0")
example("BoccSS0")
example("survCJS")
example("weta")
example("richACE")
example("allCombinations")
example("AICc")
example("distBrayCurtis")
example("distShell")
example("distTestData")
example("KillarneyBirds")
example("salamanders")
example("richCurve")
example("richRarefy")
example("Temburong")
example("railSims")
example("occSSrnSite")
example("Bnormal")
example("diagPlot")
example("Bpoisson")
example("Bbinomial")

