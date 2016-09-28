
# Functions to sanity-check different arguments to 'wiqid' functions.

# Not exported.

# Check detection histories
# -------------------------
verifyDH <- function(DH, allowNA=FALSE) {
  if (!is.matrix(DH) && !is.data.frame(DH))
    stop("Detection history must be a matrix or data frame")
  DH <- as.matrix(DH)
  if(!is.numeric(DH))
    stop("Detection history has non-numeric values.")
  if(!allowNA && any(is.na(DH)))
    stop("Detection history has NA values.")
  if(!any(is.finite(DH)))
    stop("Detection history has no non-missing values.")
  DH <- round(DH)
  range <- range(DH, na.rm=TRUE)
  if(range[1] < 0)
    stop("Detection history has negative values.")
  if(range[2] > 1)
    stop("Detection history has values > 1.")
  return(DH)  
}
