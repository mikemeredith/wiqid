
# This is intended to become a 'predict' method for 'wiqid' objects

# Additional info to be included in the 'wiqid' fit:
# = model formulae
# = preprocessing information:
#   - for numeric covars, scaling info (mean, SD/2)
#   - for factors, original levels

predict.wiqid <- function(object, newdata, parameter, ci, type=c("link", "response"), ...) {
  
  if(missing(newdata) || !is.data.frame(newdata))
    stop("Please supply a data frame for newdata.")
  if(missing(ci))
    ci <- object$ci
  if(is.null(ci))
    ci <- 0.95
  crit <- fixCI(ci)
  
  # get the model formula
  forms <- object$formulae
  if(is.null(forms))
    stop("No information on models in object.")
  if(missing(parameter))
    parameter <- names(forms)[1]
  parID <- pmatch(parameter, names(forms))  
  if(is.na(parID))
    stop("No submodel found for parameter ", parameter)
  index <- object$index[[parID]]
  formula <- forms[[parID]]
  varsNeeded <- all.vars(formula)
  
  if(length(varsNeeded) == 0) {   # INTERCEPT ONLY MODEL
    stopifnot(length(index) == 1) # check
    cat("This is an intercept-only model, all values identical.\n")    
    intercept <- object$beta[index, 1:2] # est and SE
    intercept <- c(intercept, intercept[1] + intercept[2] * crit)
    lp.mat <- matrix(rep(intercept, each=nrow(newdata)), nrow(newdata))
    rownames(lp.mat) <- rownames(newdata)
    colnames(lp.mat) <- c("est", "SE", "lowCI", "uppCI")

  } else {
    # Get coefficients and varcovar matrix
    coeffs <- coef(object)[index]
    vcv <- object$beta.vcv[index, index]
    # wrangle newdata
    varsMissing <- !varsNeeded %in% names(newdata)
    if(any(varsMissing)) {
      missingText <- paste(varsNeeded[varsMissing], collapse=", ")
      stop("Needed variable(s) missing from new data: ", missingText)
    }
    newdata <- newdata[, names(newdata) %in% varsNeeded, drop=FALSE]
    # scale numeric covars
    newdata <- scaleToMatch(newdata, object$scaling)

    xlev <- object$xlev
    xlev <- xlev[names(xlev)%in% varsNeeded]  # can be empty
    mf <- model.frame(formula, newdata, xlev=xlev)
    modMat <- model.matrix(formula, mf)
    
    
    # Get point estimates and SEs
    lp.mat <- matrix(NA_real_, nrow(modMat), 4)
    rownames(lp.mat) <- rownames(modMat)
    colnames(lp.mat) <- c("est", "SE", "lowCI", "uppCI")
    lp.mat[, 1] <- modMat %*% coeffs
    lp.mat[, 2] <- sqrt(diag(modMat %*% vcv %*% t(modMat)))
    lp.mat[, 3:4] <- sweep(outer(lp.mat[, 2], crit), 1, lp.mat[, 1], "+")
  }
  type <- match.arg(type)
  if(type == "response") {
    link <- object$link
    if(is.null(link))
      link <- "logit"
    plink <- switch(link,
      logit = plogis,
      probit = pnorm,
      stop("Link type ", link, " not recognised."))
    SE <- lp.mat[, 2]
    lp.mat <- plink(lp.mat)
    lp.mat[, 2] <- lp.mat[, 1] * (1 - lp.mat[, 1]) * SE
  }
  attr(lp.mat, "ci") <- ci
  attr(lp.mat, "link") <- object$link
  return(lp.mat)
}