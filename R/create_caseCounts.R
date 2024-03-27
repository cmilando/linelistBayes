#' Title
#'
#' @param date_vec 
#' @param location_vec 
#' @param cases_vec 
#' @import lubridate
#' @import magrittr
#' @import dplyr
#' @return
#' @export
#'
#' @examples
create_caseCounts <- function(date_vec, location_vec, cases_vec) {
  
  # check types
  stopifnot(all(is.Date(date_vec)))
  stopifnot(all(is.character(location_vec)))
  stopifnot(all(is.numeric(cases_vec)))
  
  # no negative case numbers
  stopifnot(all(cases_vec >= 0))
  
  # assumes this is just for one location
  stopifnot(length(unique(location_vec)) == 1)
  
  # should all have the same length
  stopifnot(length(date_vec) == length(location_vec))
  stopifnot(length(date_vec) == length(cases_vec))
  
  # create caseCounts
  caseCounts <- data.frame(date = date_vec, 
                           cases = cases_vec, 
                           location = location_vec)
  
  # no is.na
  necessary_columns <- c('date', 'cases', 'location')
  stopifnot(all(! is.na(caseCounts[, necessary_columns])))
  
  # add class attribute
  class(caseCounts) <- c(class(caseCounts), "caseCounts")
  
  return(caseCounts)
  
}