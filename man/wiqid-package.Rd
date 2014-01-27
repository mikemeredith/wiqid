\name{wiqid-package}
\alias{wiqid-package}
\alias{wiqid}
\docType{package}
\title{
Estimation functions for wildlife parameters
}
\description{
Quick and dirty functions to estimate occupancy, survival, abundance, species richness and diversity, etc. for wildlife populations.
}
\details{

There are a number of sophisticated programs for the analysis of wildlife data, producing estimates of occupancy, survival, abundance, or density. \pkg{wiqid} began as a collection of fast, bare-bones functions which can be run from R suitable for use when you are generating hundreds of simulated data sets.

We now use \pkg{wiqid} in basic wildlife study design and data analysis workshops, where we want to focus on the concepts behind occupancy, mark-recapture, survival, etc. and study design issues. Participants have used \code{lm}, \code{glm} and functions in the \code{secr} and \code{BEST} packages. So  \pkg{wiqid} tries to match the look and feel of these functions.

All functions use standard data frames or matrices for data input. ML estimation function return objects of class \code{wiqid} with parameter estimates on the transformed scale(usually logit functions), variance-covariance matrix, and back-transformed `real' values; there are \code{print} and \code{logLik} methods. Bayesian functions (distinguished by an initial "B") return class \code{Bwiqid} objects, a data frame with a column of MCMC observations for each parameter, plus attributes for \code{sampleSize}; there are \code{print} and \code{plot} methods.

Simulations and bootstraps often generate weird data sets, eg. capture histories with no captures. These functions do not throw errors or give warnings if the data are weird, but return NAs if estimates cannot be calculated. Errors may still occur if the data are impossible, eg. 6 detections in 5 occasions.

The functions are listed by topic below.
}

\section{OCCUPANCY}{

  \bold{Single-season occupancy}
  
  \tabular{ll}{
  \code{\link{occSScov}} \tab general-purpose function; allows site- and survey-specific covariates \cr
  \code{\link{occSS0}} \tab a basic psi(.) p(.) model, faster if this is all you need \cr
  \code{\link{BoccSS0}} \tab a Bayesian implementation of the psi(.) p(.) model \cr
  \code{\link{occSSrn}} \tab Royle-Nichols method \cr
  \code{\link{occSStime}} \tab faster if you have only time effects, also does a plot \cr
  \code{\link{occSScovSite}} \tab faster if you only have site-specific covariates \cr
  \code{\link{occSSautocorA}} \tab autocorrelation among surveys  }

  \bold{Multi-season occupancy}
  
  \tabular{ll}{
  \code{\link{occMScov}} \tab general-purpose function; parameters depend on covariates \cr
  \code{\link{occMS0}} \tab a simple multi-season model with four parameters; faster \cr
  \code{\link{occMSseason}} \tab parameters vary by season; faster  } 
}

\section{DENSITY from spatial capture-recapture data}{
  
We use the \pkg{secr} package for ML estimation of density. For Bayesian estimation, \pkg{wiqid} offers:

  \tabular{ll}{
  \code{\link{Bsecr0}} \tab a Bayesian implementation of the intercept-only model }
}

\section{ABUNDANCE from closed-population capture-recapture data}{
  
Although data for genuinely closed populations are rare, this is an important conceptual stepping-stone from CJS models to robust models.

  \tabular{ll}{
  \code{\link{closedCapM0}} \tab simple model with constant capture probability \cr
  \code{\link{closedCapMt}} \tab capture probability varies with time \cr
  \code{\link{closedCapMtcov}} \tab allows for time-varying covariates \cr
  \code{\link{closedCapMh2}} \tab heterogenity with 2-mixture model \cr
  \code{\link{closedCapMhJK}} \tab jackknife estimator for heterogeneity  }
}

\section{SURVIVAL from capture-recapture data}{

  \bold{Cormack-Jolly-Seber models}
  
  \tabular{ll}{
  \code{\link{survCJS}} \tab model with time-varying covariates \cr
  \code{\link{BsurvCJS}} \tab a Bayesian implementation of the above \cr
  \code{\link{survCJSaj}} \tab allows for different survival for adults and juveniles }
  
Functions for Pollock's robust method are planned.
}

\section{SPECIES RICHNESS from species x sample matrices}{

  \bold{Rarefaction}
  
  \tabular{ll}{
  \code{\link{richRarefy}} \tab Mao's tau estimator for rarefaction \cr
  \code{\link{richCurve}} \tab a shell for plug-in estimators, for example... \cr
  \code{\link{richSobs}} \tab the number of species observed \cr
  \code{\link{richSingle}} \tab the number of singletons observed \cr
  \code{\link{richDouble}} \tab the number of doubletons observed \cr
  \code{\link{richUnique}} \tab the number of uniques observed \cr
  \code{\link{richDuplicate}} \tab the number of duplicates observed }
  
  \bold{Coverage estimators}
  
  \tabular{ll}{
  \code{\link{richACE}} \tab Chao's Abundance-based Coverage Estimator \cr
  \code{\link{richICE}} \tab Chao's Incidence-based Coverage Estimator \cr
  \code{\link{richChao1}} \tab Chao1 estimator \cr
  \code{\link{richChao2}} \tab Chao2 estimator \cr
  \code{\link{richJack1}} \tab first-order jackknife estimator \cr
  \code{\link{richJack2}} \tab second-order jackknife estimator \cr
  \code{\link{richJackA1}} \tab abundance-based first-order jackknife estimator \cr
  \code{\link{richJackA2}} \tab abundance-based second-order jackknife estimator \cr
  \code{\link{richBoot}} \tab bootstrap estimator \cr
  \code{\link{richMM}} \tab Michaelis-Menten estimator \cr
  \code{\link{richRenLau}} \tab Rennolls and Laumonier's estimator  }
}

\section{BIODIVERSITY INDICES}{

  \bold{Alpha diversity}
  
All of these functions express diversity as the number of common species in the assemblage.
  
  \tabular{ll}{
  \code{\link{biodSimpson}} \tab inverse of Simpson's index of dominance \cr
  \code{\link{biodShannon}} \tab exponential form of Shannon's entropy \cr
  \code{\link{biodBerger}} \tab inverse of Berger and Parker's index of dominance \cr
  \code{\link{biodBrillouin}} \tab exponential form of Brillouin's index }

  \bold{Beta diversity / distance}

All of these functions produce distance measures (not similarity) on a scale of 0 to 1. The function \code{\link{distShell}} provides a wrapper to produce a matrix of distance measures across a number of sites.
  
  \tabular{ll}{
  \code{\link{distBrayCurtis}} \tab complement of Bray-Curtis index, aka 'quantitative Sorensen' \cr
  \code{\link{distChaoJaccCorr}} \tab complement of Chao's Jaccard corrected index \cr
  \code{\link{distChaoJaccNaive}} \tab complement of Chao's Jaccard naive index \cr
  \code{\link{distChaoSorCorr}} \tab complement of Chao's Sorensen corrected index \cr
  \code{\link{distChaoSorNaive}} \tab complement of Chao's Sorensen naive index \cr
  \code{\link{distChord}} \tab distance between points on a normalised sphere \cr
  \code{\link{distJaccard}} \tab complement of Jaccard's index of similarity \cr
  \code{\link{distMorisitaHorn}} \tab complement of the Morisita-Horn index of similarity \cr
  \code{\link{distOchiai}} \tab complement of the Ochiai coefficient of similarity \cr
  \code{\link{distPreston}} \tab Preston's coefficient of faunal dissimilarity \cr
  \code{\link{distRogersTanimoto}} \tab complement of the Rogers and Tanimoto's coefficient of similarity \cr
  \code{\link{distSimRatio}} \tab complement of the similarity ratio \cr
  \code{\link{distSorensen}} \tab complement of the Sorensen or Dice index of similarity \cr
  \code{\link{distWhittaker}} \tab Whittaker's index of association }
}

\section{DATA SETS}{

  \tabular{ll}{
  \code{\link{dippers}} \tab Capture-recapture data for European dippers \cr
  \code{\link{distTestData}} \tab artificial data set for distance measures \cr
  \code{\link{KanhaTigers}} \tab camera-trap data for tigers \cr
  \code{\link{KillarneyBirds}} \tab abundance of birds in Irish woodlands \cr
  \code{\link{salamanders}} \tab detection/non-detection data for salamanders \cr
  \code{\link{seedbank}} \tab number of seeds germinating from samples of soil \cr
  \code{\link{Temburong}} \tab counts of tree species in a 1ha plot in Brunei \cr
  \code{\link{TemburongBA}} \tab basal area of tree species in a 1ha plot in Brunei \cr
  \code{\link{weta}} \tab detection/non-detection data and covariates for weta }
}

\section{UTILITY FUNCTIONS}{

  \tabular{ll}{
  \code{\link{AICc}} \tab AIC with small-sample correction \cr
  \code{\link{AICtable}} \tab tabulate AIC for several models \cr
  \code{\link{allCombinations}} \tab model formulae for combinations of covariates \cr
  \code{\link{hdi}} \tab calculate Highest Density Intervals of a distribution \cr
  \code{\link{plotPost}} \tab graphic display of a posterior distribution \cr
  \code{\link{postPriorOverlap}} \tab calculation and display of posterior/prior overlap }
}

\author{
Mike Meredith

Maintainer: Mike Meredith <mmeredith@wcs.org>
}

\keyword{ package }
