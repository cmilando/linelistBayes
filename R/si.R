si <- function(ndays, alpha, beta) {
  prob <- numeric(ndays) # creates a numeric vector
  for (i in 1:ndays){
    prob[i] <- pgamma(    i, shape = alpha, rate = beta) -
      pgamma(i - 1, shape = alpha, rate = beta)
  }
  result <- prob/sum(prob) # normalizes the whole vector
  return(result)
}