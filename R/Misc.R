
# Tidy up the column names in MCMC output: remove [] and ,
fixNames <- function(x)
  sub(",", "\\.", sub("\\]", "", sub("\\[", "", x)))


# Function to calculate the MARK-style confidence intervals for N
# See help for Closed Captures
# Not exported

getMARKci <- function(beta, SE.beta, ci) {
  f0.hat <- exp(beta)
  crit <- qnorm((1 - ci[1]) / 2, lower.tail=FALSE)
  C <- exp(crit * sqrt(log(1 + SE.beta^2))) # See the Burnham et al reference, p212!
  return(c(f0.hat, f0.hat/C, f0.hat*C))
}

# Creates a table from output of AIC or AICc
# Exported

AICtable <- function(x, digits=3) {
  if(is.vector(x)) {
    name <- deparse(substitute(x))
    IC <- x
    x <- data.frame(x)
    colnames(x) <- name
  } else if (is.data.frame(x) && ncol(x) > 1 )  {
    IC <- x[, 2] 
  } else {
    stop("x must be a vector or a data frame with > 1 column")
  }
  Delta <- IC - min(IC, na.rm=TRUE)
  ModelLik <- exp( - Delta / 2)
  ModelWt <- ModelLik / sum(ModelLik, na.rm=TRUE)
  out <- round(cbind(x, Delta, ModelLik, ModelWt), digits)
  if(!is.null(rownames(out))) { # only sort if rows are named
    ord <- order(IC)
    out <- out[ord, ]
  }
  return(out)
}

## Regularize a list of formulae, ensuring it is a named list
#    of one-sided formulae.
## based on Murray Efford's 'stdform' function in 'secr'

# Not (currently) exported

# Old version, to be phased out
stdform <- function (flist) {
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

## Convert a data frame of site and survey data into a list 
# Site covars will each have a single column in the data frame,
# survey covars will have a column for each survey occasion, and
# column names end with the number of the occasion, eg, temperature
# will be in columns named "temp1", "temp2", etc.

# Not (currently) exported

stddata <- function(df, nocc)  {
  stopifnot(is.data.frame(df))
  dataList <- as.list(df)
  # look for names ending with number of occasions
  nam <- names(df)
  clue <- paste0(nocc, "$", collapse="")
  clueDo <- grep(clue, nam)
  if(length(clueDo) > 0) {
    for(i in clueDo) {
      # get stem, generate set of names
      stem <- sub(clue, "", nam[i])
      subnames <- paste0(stem, 1:nocc)
      subtable <- df[, subnames]
      # check that there's a column for each occasion
      stopifnot(ncol(subtable) == nocc)  # do less brutal thing later
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
  return(dataList)
}
