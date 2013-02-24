
# From Kery & Schaub

# Function to create a m-array based on capture-histories (CH)

ch2mArray <- function(CH){
   nind <- dim(CH)[1]
   n.occasions <- dim(CH)[2]
   m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
   # Calculate the number of released individuals at each time period
   for (t in 1:n.occasions){
      m.array[t,1] <- sum(CH[,t])
      }
   for (i in 1:nind){
      pos <- which(CH[i,]!=0)
      g <- length(pos)
      for (z in 1:(g-1)){
         m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
         } # z
      } # i
   # Calculate the number of individuals that is never recaptured
   for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
      }
   out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
   return(out)
   }

