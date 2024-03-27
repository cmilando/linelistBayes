#' Title
#'
#' @param caseCounts 
#' @param reportF 
#' @param reportF_args 
#' @param reportF_missP 
#'
#' @return
#'
#' @examples
#' @import magrittr
#' @import dplyr
#' @import lubridate
#' @export
create_linelist <- function(caseCounts, 
                            reportF = NULL, 
                            reportF_args = NULL,
                            reportF_missP = NULL) {
  
  # ---------------------------------------
  # validate caseCounts
  necessary_columns <- c('date', 'cases', 'location')
  stopifnot(colnames(caseCounts) %in% necessary_columns)
  
  # check types
  stopifnot(all(is.Date(caseCounts$date))) # requires(lubridate)
  stopifnot(all(is.character(caseCounts$location)))
  stopifnot(all(is.numeric(caseCounts$cases)))
  
  # no negative case numbers
  stopifnot(all(caseCounts$cases >= 0))
  
  # assumes this is just for one location
  stopifnot(length(unique(caseCounts$location)) == 1)
  
  # no is.na
  stopifnot(all(! is.na(caseCounts[, necessary_columns])))
  
  # message if reportF is not specific
  if(is.null(reportF)) {
    message("reportF is NULL, setting to runif(n, min = 4, max = 6)")
    reportF = runif
    reportF_args = list(min = 4, max = 6)
  }
  
  # class
  stopifnot("caseCounts" %in% class(caseCounts))
  
  # ---------------------------------------
  # Expand case dates into individuals
  report_date_l <- lapply(1:nrow(caseCounts), function(i) 
    caseCounts$date[rep(i, caseCounts$cases[i])])
  
  report_date_vec <- do.call(c, report_date_l)
  
  stopifnot(identical(length(report_date_vec), sum(caseCounts$cases)))
  
  # ---------------------------------------
  ## now add other columns
  d <- data.frame(report_date = sort(report_date_vec))
  
  ## need to assume some generation distribution here
  tryCatch({
    reportF_args$n <- length(report_date_vec)
    d$delay_int <- round(do.call(reportF, reportF_args))
  }, error = function(cond) {
    stop('Error in arguments to onset function: `n` or user specified')
  })
  
  # checking onset
  stopifnot(all(d$delay_int > 0))
  
  # remove missingness
  if(!is.null(reportF_missP)) {
    stopifnot(is.numeric(reportF_missP))
    stopifnot(reportF_missP >= 0 & reportF_missP < 1)
    remove_delays <- sample(1:nrow(d), size = nrow(d) * reportF_missP)
    d$delay_int[remove_delays] <- NA
  }
  
  d$onset_date = d$report_date - d$delay_int
  
  # add extra columns
  # requires dplyr
  d <- d %>% 
    mutate(
      is_weekend = ifelse(wday(report_date) %in% c(0, 7), 1, 0))
  
  # offset and covert to numbers
  minday = min(d$report_date)
  
  # requires dplyr
  d <- d %>% 
    mutate(report_int = as.vector(report_date - minday + 1), 
           week_int = ceiling(report_int / 7))
  
  # type convert
  d$delay_int <- as.integer(d$delay_int)
  d$is_weekend <- as.integer(d$is_weekend)
  d$report_int <- as.integer(d$report_int)
  d$week_int <- as.integer(d$week_int)
  
  class(d) <- c(class(d), "lineList")

  return(d)
}
