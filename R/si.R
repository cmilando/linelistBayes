#' Title
#'
#'  There are different types of serial interval used here
#'  https://www.sciencedirect.com/science/article/pii/S1201971220306111
#'  p[j] is the probability of a serial interval of j days, exactly, and
#'  not the sum of all things before it, so this gives the PDF, not the CDF
#'  right because pgamma() gives the CDF, not the PDF
#'  so this just creates the PDF
#'  so a vector of length 14 with alpha and beta as defined in the doc, for COVID-19
#'
#' @param ndays 
#' @param alpha 
#' @param beta 
#'
#' @return
#' @export
#'
#' @examples
si <- function(ndays, alpha, beta) {
  
  prob <- numeric(ndays) # creates a numeric vector
  
  for (i in 1:ndays){
    prob[i] <- pgamma(i, shape = alpha, rate = beta) -
      pgamma(i - 1, shape = alpha, rate = beta)
  }
  
  result <- prob/sum(prob) # normalizes the whole vector
  
  return(result)
}