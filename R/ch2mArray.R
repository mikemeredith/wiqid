
# Adapted from Kery & Schaub

# Function to create a m-array based on capture-histories (CH) plus vector of frequencies.

ch2mArray <- function(CH, freq=1){
  if(length(freq) == 1)
    freq <- rep(freq, nrow(CH))
  stopifnot(length(freq) == nrow(CH))
  if(any(freq < 0))
    stop("Sorry, I cannot deal with trap losses (yet).")
  nocc <- ncol(CH)
  ma <- matrix(0, nocc, nocc+1)
    # First column and last row will be removed later
    # Last col is for number never-seen-again
  for(i in 1:nrow(CH)) {
    cht <- which(CH[i, ] != 0) # When was animal caught?
    # Fill in release/recapture data
    for(z in seq_along(cht[-1]))  # Does nothing if length(cht) = 1
      ma[cht[z], cht[z+1]] <- ma[cht[z], cht[z+1]] + freq[i]
  }
  # Marked animals never seen again:
  totCH <- sweep(CH, 1, freq, "*")
  ma[, nocc+1] <- colSums(totCH) - rowSums(ma)
  # Remove 1st col (REcaptures on 1st occasion = all zeros) 
  #  and last row (releases  on last occasion will never be recaptured).
  return(ma[-nocc, -1])
}

