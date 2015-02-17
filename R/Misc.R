
# This file contains utilities used in several places in the code
#   and NOT exported:

# signifish : an alternative to signif (added 10-02-2015)
# fixCI : Calculate critical values for CI.
# fixNames : Tidy up the column names in MCMC output: remove [] and ,
# getMARKci : Calculate MARK-style confidence intervals for N
# stdModel : Regularize a list of formulae, ensuring it is a named list of one-sided formulae.
# stddata : Convert a data frame of site and survey data into a list and standardise
# selectCovars : Pull the covars needed for a model matrix into a specific data frame
# AICtable moved to file AICc.R
# ...............................................................................

# A more sensible version of signif
signifish <- function(x, digits=3)
  ifelse(x < 10^digits, signif(x, digits=digits), round(x))

# ...............................................................................

# Deal with confidence interval specification:
fixCI <- function(ci) {
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  return(qnorm(c(alf, 1 - alf)))
}
# .....................................................................

# Tidy up the column names in MCMC output: remove [] and ,
fixNames <- function(x)
  sub(",", "\\.", sub("\\]", "", sub("\\[", "", x)))
# .....................................................................

# Function to calculate the MARK-style confidence intervals for N
# See help for Closed Captures

getMARKci <- function(beta, SE.beta, ci) {
  f0.hat <- exp(beta)
  crit <- qnorm((1 - ci[1]) / 2, lower.tail=FALSE)
  C <- exp(crit * sqrt(log(1 + SE.beta^2))) # See the Burnham et al reference, p212!
  return(c(f0.hat, f0.hat/C, f0.hat*C))
}
# .........................................................................

## Regularize a list of formulae, ensuring it is a named list of one-sided formulae.
## based on Murray Efford's 'stdform' function in 'secr'

# Old version, to be phased out
stdform <- function (flist) {
  warning("stdform is deprecated. Use stdModel instead.")
    LHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==2) '' else trms[2]
    }
    RHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==3) as.formula(paste(trms[c(1,3)])) else form
    }
    lhs <- sapply(flist, LHS)
    temp <- lapply(flist, RHS)
    if (is.null(names(flist))) names(temp) <- lhs
    else names(temp) <- ifelse(names(flist) == '', lhs, names(flist))
    temp
}

# New version:
stdModel <- function (model1, defaultModel) {
  if(is.null(model1))
    return(defaultModel)
  if(inherits(model1, "formula"))
    model1 <- list(model1)
  stopifnot(is.list(model1))
  LHS <- function (form) {
      trms <- as.character (form)
      if (length(trms)==2) '' else trms[2]
  }
  RHS <- function (form) {
      trms <- as.character (form)
      if (length(trms)==3) as.formula(paste(trms[c(1,3)])) else form
  }
  lhs <- sapply(model1, LHS)
  temp <- lapply(model1, RHS)
  if (is.null(names(model1))) {
    names(temp) <- lhs
  } else {
    names(temp) <- ifelse(names(model1) == '', lhs, names(model1))
  }
  newModel <- replace (defaultModel, names(temp), temp)
  return(newModel)
}
# .............................................................................

## Convert a data frame of site and survey data into a list 
# ** Site covars will each have a single column in the data frame,
# ** survey covars will have a column for each survey occasion, and
# column names end with the number of the occasion, eg, temperature
# will be in columns named "temp1", "temp2", etc.

stddata <- function(df, nocc=NULL, scaleBy=0.5)  {
  if (is.null(df))
    return(NULL)
  stopifnot(is.data.frame(df))
  dataList <- as.list(df)
  ## Group variables spread over > 1 column into a single vector
  if (!is.null(nocc)) {
    nocc <- sort(nocc, decreasing=TRUE) # start with biggest
    for (this.nocc in nocc)  {
      # look for names ending with number of occasions
      nam <- names(df)
      clue <- paste0(this.nocc, "$", collapse="")
      clueDo <- grep(clue, nam)
      if(length(clueDo) > 0) {
        for(i in clueDo) {
          # get stem, generate set of names
          stem <- sub(clue, "", nam[i])
          subnames <- paste0(stem, 1:this.nocc)
          subtable <- df[, subnames]
          # check that there's a column for each occasion
          stopifnot(ncol(subtable) == this.nocc)  # do less brutal thing later
          # check that all have same class
          classes <- sapply(subtable, class)
          stopifnot(length(unique(classes)) == 1)  # do less brutal thing later
          # remove original columns from the list:
          dataList <- replace(dataList, subnames, NULL)
          # convert to a matrix, then a vector;
          #   fortunately this also converts factors to character
          tmp <- as.vector(as.matrix(subtable))
          if(is.character(tmp))
            tmp <- as.factor(tmp) # convert to factor AFTER combining columns
          dataList <- c(dataList, list(tmp))
          names(dataList)[length(dataList)] <- stem
        }
      }
    }
  }
  ## Standardize numeric variables to mean = 0 and sd = scaleBy
  if (!is.null(scaleBy)) {
    doScale <- function(x) {
      if (is.numeric(x))
        x <- as.vector(scale(x) * scaleBy)
      return(x)
    }
    dataList <- lapply(dataList, doScale)
  }
  return(dataList)
}
# ...........................................................................

# Pull the covars needed for a model matrix into a specific data frame
selectCovars <- function(formula, dataList, minrows)  {
  wanted <- rownames(attr(terms(formula), "factors"))
  found <- wanted %in% names(dataList)
  wanted <- wanted[found]
  if (length(wanted) > 0)  {
    df <- as.data.frame(dataList[wanted])
    df <- cbind(df, .dummy = rep(NA, minrows))
  } else {
    df <- data.frame(.dummy = rep(NA, minrows))
  }
  stopifnot(nrow(df) %% minrows == 0)
  return(df)
}






