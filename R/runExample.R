#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", package = "linelistBayes")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `linelistBayes`.", call. = FALSE)
  }
  
  rmarkdown::run(paste0(appDir, "/", "v0.Rmd"))
}

#' @export
runBackNowExample <- function() {

  data("halfyear")
  dat_tbl <- halfyear$dat_tbl

  #, so you really just have the black line
  # and only up to a certain point

  # (1)
  # so you want to recreate the red line (onset)

  # expand.grid each person
  dat_tbl <- dat_tbl %>%
    mutate(d1_int = type.convert(d1, as.is = T))
  dim(dat_tbl)
  head(dat_tbl)

  dat_l <- lapply(1:nrow(dat_tbl), function(i)
    dat_tbl$d1_int[rep(i, dat_tbl$Freq[i])])
  dat1_vec <- do.call(c, dat_l)
  stopifnot(identical(length(dat1_vec), sum(dat_tbl$Freq)))

  ## now add other columns
  d <- data.frame(onset = NA, report = dat1_vec)
  head(d)
  min(d$report)

  d <- d %>% arrange(report)
  d <- d %>% mutate(weekend = ifelse(report %% 7 == 0 | report %% 7 == 6, 1, 0))
  d <- d %>% mutate(delay = report - onset,
                    minday = min(report))
  d <- d %>% mutate(report = report - minday + 1, week = ceiling(report / 7))
  head(d)

  # epic <- out$epic
  # repc <- out$repc
  minday <- unique(d$minday)

  # ---------
  ### ok now need to load RCPP


  # so a vector of length 14 with alpha and beta as defined in the doc
  sip <- si(14, 4.29, 1.18)
  sip

  dim(d)
  ## UPDATED `si` to `sip` just to make the function work
  # reduced iter to 100 for testing
  # Notes: only the line-list data are passed here
  # these are evaluated against EPIC and REPC later one
  MAX_ITER <- 2000
  out1 <- backnow(outcome = d$delay, days = d$report, week = d$week,
                  weekend = d$weekend, iter = MAX_ITER, sigma = 0.2,
                  maxdelay = 20, si = sip, size = 6)
  
  print("done")
  
}