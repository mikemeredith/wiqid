# This is based on AIC in package stats

# There are AICc functions in other packages, which don't all work with
#   occupancy, but AIC does.

AICc <- function(object, ..., nobs) UseMethod("AICc")

AICc.default <- function (object, ..., nobs) 
{
  if (length(list(...))) {
    lls <- lapply(list(object, ...), logLik)
    vals <- sapply(lls, function(el) {
        no <- attr(el, "nobs")
        c(as.numeric(el), attr(el, "df"), if (is.null(no)) NA_integer_ else no)
    })
    #val <- data.frame(df = vals[2L, ], ll = vals[1L, ])
    nos <- na.omit(vals[3L, ])
    if (length(nos) && any(nos != nos[1L])) 
        warning("models are not all fitted to the same number of observations")
    # val <- data.frame(df = val$df, AIC = -2 * val$ll + k * val$df)
    df <- vals[2L, ]
    if(missing(nobs))
      nobs <- vals[3L, ]
    val <- data.frame(df = df,
        AICc = -2 * vals[1L, ] + 2 * df * nobs / pmax(0, nobs - df - 1))
    Call <- match.call()
    Call$nobs <- NULL
    row.names(val) <- as.character(Call[-1L])
  }
  else {
    lls <- logLik(object)
    df <- attr(lls, "df")
    if(missing(nobs))
      nobs <- attr(lls, "nobs")
    if(is.null(nobs))  {
      val <- NA_real_
    } else {
      val <- -2 * as.numeric(lls) + 2 * df * nobs / max(0, nobs - df - 1)
    }
  }
  return(val)
}

# ..........................................................................

# Creates a table from output of AIC or AICc
# Exported

AICtable <- function(x, digits=3, sort) {
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
  if(missing(sort))
    sort <- !is.null(rownames(out)) # sort if rows are named
  if(sort) { 
    ord <- order(IC)
    out <- out[ord, ]
  }
  return(out)
}

