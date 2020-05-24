# Individual checks

# setwd("D:/Github/wiqid_package")
setwd("../..")
dir()

# Dependencies
# ============
install.packages(c("truncnorm", "coda", "plotrix", "secr", "shiny", "rjags"))
devtools::install_github("mikemeredith/mcmcOutput")

# Spell check
# ===========
library(devtools)
sIg <- scan("spellcheckIgnore.txt", what='character', comment="#")
tmp <- spell_check("wiqid", ignore=c(hunspell::en_stats, sIg), "en_GB")
length(tmp)  # number of misspellings found
tmp  # error if length == 0

devtools::load_all("C:/GitHub/wiqid_package/wiqid")

# Create the wiqid package
# ========================
unlink(list.files(pattern="Rplots.pdf", recursive=TRUE))
system("R CMD build wiqid")  # Produces the .tar.gz file
# system("R CMD check wiqid_0.2.3.9012.tar.gz")
system("R CMD check --as-cran wiqid_0.2.3.9012.tar.gz")
# system("R CMD check --run-donttest wiqid_0.2.3.9012.tar.gz")
# Sys.setenv(R_ZIPCMD = "C:/Rtools/bin/zip.exe")
system("R CMD INSTALL --build wiqid_0.2.3.9012.tar.gz") # installs and produces the .zip binary
system("R CMD INSTALL wiqid_0.2.3.9012.tar.gz") # installs only

system("R CMD INSTALL wiqid") # Use this for a "dev" install.


# Test it:
library(testthat)
test_package("wiqid", reporter=ProgressReporter)
# test_dir("D:/GitHub/wiqid_package/wiqid/inst/tests/testthat", reporter=ProgressReporter)


# Try it out:
rm(list=ls())
library(wiqid)
?wiqid

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

