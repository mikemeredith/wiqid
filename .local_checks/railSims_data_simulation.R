
# Simulated data for 2-species occupancy

library(wiqid)

# The original version, without covariates:
# =========================================
# The data come from a simulated scenario with the following parameters:
# psiA 	= 0.39 # 	probability of occupancy of species A
# psiBA 	= 0.77 # 	probability of occupancy of B if A is present
# psiBa 	= 0.39 # 	probability of occupancy of B if A is absent
# pA 	= 0.74 # 	probability of detection of species A if B is absent
# pB 	= 0.80 # 	probability of detection of species B if A is absent
# rA 	= 0.60 # 	probability of detection of species A if both are present
# rBA 	= 0.88 # 	probability of detection of species B if both are present and A was detected
# rBa 	= 0.79 # 	probability of detection of species B if both are present and A was not detected
# These values approximate those reported by Richmond et al (2010).

# New version, with covariates
# ============================
# Biological model
nSites <- 160
set.seed(123)
# logArea <- scale(rnorm(nSites))  # returns a 1-column matrix
logArea <- standardize(rnorm(nSites))  # Assume standardised
reeds <- rbinom(nSites, 1, 0.5)==1
beta0 <- 0
beta1 <- 2
psiA <- plogis(beta0 + beta1*logArea)
zA <- rbinom(nSites, 1, psiA)
mean(zA)
psiBa <- 0.77
alpha0 <- -1
alpha1 <- 2
psiBA <- plogis(alpha0 + alpha1*reeds)
psiB <- ifelse(zA, psiBA, psiBa)
zB <- rbinom(nSites, 1, psiB)
mean(zB)

# Detection model
nOccs <- 3
pA <- 0.75 # rA = pA
DHA <- matrix(rbinom(nSites*nOccs, 1, pA*zA), nSites, nOccs)
pB <- 0.80 # rBa = pB,
rBA <- 0.4
rB <- ifelse(DHA, rBA, pB)
DHB <- matrix(rbinom(nSites*nOccs, 1, rB*zB), nSites, nOccs)


# Put together a data frame
railSims <- as.data.frame(cbind(DHA, DHB))
colnames(railSims) <- c("A1","A2","A3", "B1", "B2", "B3")
head(railSims)
railSims$logArea <- round(logArea, 3)
railSims$reeds <- reeds

# compare with old data set
rsnew <- railSims
rm(railSims)
data(railSims)
all.equal(rsnew, railSims, check.attributes=FALSE)
all.equal(rsnew$logArea, c(railSims$logArea))

save("railSims", file="railSims.rda", version=2)
library(tools)
checkRdaFiles("railSims.rda")
resaveRdaFiles("railSims.rda")
checkRdaFiles("railSims.rda")

# export a spreadsheet for PRESENCE
PRESENCE <- railSims
PRESENCE$reeds <- reeds*1
write.csv(PRESENCE, "railSims.csv", row.names=FALSE)
