#' Title
#'
#' @param input 
#' @param MAX_ITER 
#' @param norm_sigma 
#' @param sip 
#' @param NB_maxdelay 
#' @param NB_size 
#' @param workerID 
#' @param printProgress 
#' @param ... 
#'
#' @return

#' @import coda
#' @import magrittr
#' @import dplyr
#' @export
run_backnow <- function(input, MAX_ITER = 2000, 
                        norm_sigma = NULL, sip = NULL, 
                        NB_maxdelay = NULL, NB_size = NULL, 
                        workerID = NULL, 
                        printProgress = 0, ...) {
  
  # ---------------------------------------------------------
  # input checking
  if(all(c("caseCounts", "lineList") %in% class(input)))    stop()
  if(all(!(c("caseCounts", "lineList") %in% class(input)))) stop()
  
  input_type <- 'NA'
  if("caseCounts" %in% class(input)) {
    timestamp(suffix = "> inputType = caseCounts")
    input_type = 'caseCounts'
    # get lineList, with random each time, from ...
    caseCounts_line <- create_linelist(input, ...)
  } 
  
  if("lineList" %in% class(input)) {
    timestamp(suffix = "> inputType = lineList")
    input_type = 'lineList'
    # this is input.
    caseCounts_line <- input
  } 
  stopifnot(!is.na(input_type))
  
  # maxiter
  stopifnot(is.integer(MAX_ITER))
  stopifnot(MAX_ITER >= 2000)
  if(MAX_ITER >= 20000) warning('`MAX_ITER` >= 20,000 will lead to long run times')
  
  # norm sigma for bayes parameters
  stopifnot(is.numeric(norm_sigma))
  stopifnot(norm_sigma > 0)
  
  # NB_maxdelay, the maximum of the right truncated distribution
  stopifnot(is.integer(NB_maxdelay))
  stopifnot(NB_maxdelay >= 0)
  stopifnot(NB_maxdelay < 50)
  
  # vector for SIP
  stopifnot(is.vector(sip))
  stopifnot(all(is.numeric(sip)))
  stopifnot(all(sip >= 0))
  
  # size of the NB distribut
  stopifnot(is.integer(NB_size))
  stopifnot(NB_size >= 0)
  
  #
  if(is.null(workerID)) workerID <- as.integer(0)
  stopifnot(is.numeric(workerID))
  workerID <- as.integer(workerID)
  stopifnot(workerID >= 0)
  
  #
  stopifnot(printProgress %in% c(0, 1))
  if(!dir.exists(file.path(".", "tmp")) & printProgress == 1) {
    warning("`tmp` dir does not exist, setting printProgress to `0`")
    printProgress <- as.integer(0)
  }
  printProgress <- as.integer(printProgress)
  
  # ---------------------------------------------------------
  # n backnow
  # no progress printing 
  out_list <- 
    backnow_cm(outcome       = caseCounts_line$delay_int, 
               days          = caseCounts_line$report_int, 
               week          = caseCounts_line$week_int,
               weekend       = caseCounts_line$is_weekend, 
               workerID      = workerID,
               printProgress = printProgress,
               iter          = MAX_ITER, 
               sigma         = norm_sigma,
               maxdelay      = NB_maxdelay, 
               si            = sip, 
               size          = NB_size)

  # ---------------------------------------------------------
  # process back-calcution and r(t) across chains
  # after 1000 burn-in
  # also every 2 ...
  N_BURN_IN <- 1000
  probs_to_export <- c(0.025, 0.5, 0.975)
  
  back1  <- out_list$Back[seq(N_BURN_IN + 1, nrow(out_list$Back), by = 2), ]
  
  # depending onthe onset distribution, you may have some NA cols at the tails
  # remove those, but throw an error if ANY NAS remain
  # which columns are ALL NA
  remove_cols <- c()
  for(j in 1:ncol(back1)) {
    if(all(is.na(back1[, j]))) remove_cols <- c(j, remove_cols)
  }
  if(length(remove_cols) > 0) {
    back1 <- back1[, -remove_cols]
  }
  
  for(j in 1:ncol(back1)) {
    if(any(is.na(back1[, j]))) stop("some NA values in `back`")
  }
  
  est_back <- apply(back1, 2, function(x) quantile(x, probs = probs_to_export))
  
  # gewke diagnostics
  gback1 <- geweke.diag(back1)$z
  gback1[is.nan(gback1)] <- 0
  gb1  <- sum(abs(gback1) > 1.96) / length(gback1)
  
  # ---------------------------------------------------------
  # process the r(t)
  r1 <- out_list$R[seq(N_BURN_IN + 1, nrow(out_list$R), by = 2), ]
  
  # which columns are ALL NA
  remove_cols <- c()
  for(j in 1:ncol(r1)) {
    if(all(is.na(r1[, j]))) remove_cols <- c(j, remove_cols)
  }
  if(length(remove_cols) > 0) {
    r1 <- r1[, -remove_cols]
  }
  
  for(j in 1:ncol(r1)) {
    if(any(is.na(r1[, j]))) stop("some NA remains in `r1`")
  }
  
  ##
  est_r <- apply(r1, 2, function(x) quantile(x, probs = probs_to_export))
  
  # gewke diagnostics
  gr1 <- geweke.diag(r1)$z
  gr1[is.nan(gr1)] <- 0
  gr1 <- sum(abs(gr1) > 1.96) / length(gr1)
  
  # ---------------------------------------------------------
  # 
  output_table <- caseCounts_line %>% 
    group_by(report_date) %>% 
    summarize(.groups = 'keep', 
              n = n())
  
  ## DIM OF EST_BACK = length(report_date) + NB_maxdelay
  pre <-  seq.Date(min(output_table$report_date) - NB_maxdelay, 
                   min(output_table$report_date) - 1, by = '1 day')
  
  x_est <- c(pre, output_table$report_date)
  
  ## DIM OF RT is = length(report_date) + NB_maxdelay - NB_size - 1
  # but its the last dates
  rt_size  <- length(x_est) - NB_size - 1
  rt_first <- length(x_est) - rt_size + 1
  rt_last  <- length(x_est)

  
  return(structure(class = "backnow", list(est_back      = est_back, 
                                    est_back_date = x_est, 
                                    est_rt        = est_r, 
                                    est_rt_date   = x_est[rt_first:rt_last], 
                                    geweke_back   = gb1, 
                                    geweke_rt     = gr1, 
                                    report_date   = output_table$report_date,
                                    report_cases  = output_table$n, 
                                    MAX_ITER      = MAX_ITER, 
                                    norm_sigma    = norm_sigma,
                                    NB_maxdelay   = NB_maxdelay, 
                                    si            = sip, 
                                    NB_size       = NB_size)))

}