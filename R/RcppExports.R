# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

findmiss <- function(x) {
    .Call(`_linelistBayes_findmiss`, x)
}

get_mu_vec <- function(x12, beta) {
    .Call(`_linelistBayes_get_mu_vec`, x12, beta)
}

dummy <- function(week, weekend) {
    .Call(`_linelistBayes_dummy`, week, weekend)
}

logLikNB <- function(delay_vec, x12, disp, betaplus, maxdelay) {
    .Call(`_linelistBayes_logLikNB`, delay_vec, x12, disp, betaplus, maxdelay)
}

lambda <- function(curve, si) {
    .Call(`_linelistBayes_lambda`, curve, si)
}

getr <- function(curve, si, size) {
    .Call(`_linelistBayes_getr`, curve, si, size)
}

backnow_cm <- function(outcome, days, week, weekend, iter, sigma, maxdelay, si, size, workerID, printProgress, cd = NULL) {
    .Call(`_linelistBayes_backnow_cm`, outcome, days, week, weekend, iter, sigma, maxdelay, si, size, workerID, printProgress, cd)
}

