# Print method for class occupancy, ie. output for occSSxxx family of functions.

print.occupancy <- function(x, digits=4, ...)  {
  cat("Call: ")
  print(x$call)
  cat("\nReal values (duplicates ommitted):\n")
  print(unique(x$real), digits=digits, ...)
  cat("\nAIC:", x$AIC, "\n")
}

