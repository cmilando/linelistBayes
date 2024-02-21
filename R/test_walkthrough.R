junk <- function() {
# # ok so starting with this
# library(tidyverse)
# setwd("../../simulation/cm_fcns/")
# load("halfyear.Rdata")
# 
# dat_tbl %>% 
#   mutate(d1_int = type.convert(d1, as.is = T)) %>%
#   ggplot(.) + geom_line(aes(x = d1_int, y = Freq)) +
#   geom_line(aes(x = d1_int, y = Freq), data = o1_tbl %>% 
#               mutate(d1_int = type.convert(d1, as.is = T)),
#             color = 'red') +
#   geom_line(aes(x = d1_int, y = Freq), data = i1_tbl %>% 
#               mutate(d1_int = type.convert(d1, as.is = T)),
#             color = 'blue')
# 
# 
# #, so you really just have the black line
# # and only up to a certain point
# 
# # (1)
# # so you want to recreate the red line (onset)
# 
# # expand.grid each person
# dat_tbl <- dat_tbl %>% 
#   mutate(d1_int = type.convert(d1, as.is = T))
# dim(dat_tbl)
# head(dat_tbl)
# 
# dat_l <- lapply(1:nrow(dat_tbl), function(i) 
#   dat_tbl$d1_int[rep(i, dat_tbl$Freq[i])])
# dat1_vec <- do.call(c, dat_l)
# stopifnot(identical(length(dat1_vec), sum(dat_tbl$Freq)))
# 
# ## now add other columns
# d <- data.frame(onset = NA, report = dat1_vec)
# head(d)
# min(d$report)
# 
# d <- d %>% arrange(report)
# d <- d %>% mutate(weekend = ifelse(report %% 7 == 0 | report %% 7 == 6, 1, 0))
# d <- d %>% mutate(delay = report - onset, 
#                   minday = min(report))
# d <- d %>% mutate(report = report - minday + 1, week = ceiling(report / 7))
# head(d)
# 
# # epic <- out$epic
# # repc <- out$repc
# minday <- unique(d$minday)
# minday
# 
# 
# # ---------
# ### ok now need to load RCPP
# library(RcppArmadillo)
# library(Rcpp)
# library(coda)
# Rcpp::sourceCpp("../../backnow_cm.cpp", verbose = T)
# 
# si <- function(ndays, alpha, beta) {
#   prob <- numeric(ndays) # creates a numeric vector
#   for (i in 1:ndays){
#     prob[i] <- pgamma(    i, shape = alpha, rate = beta) - 
#       pgamma(i - 1, shape = alpha, rate = beta)
#   }
#   result <- prob/sum(prob) # normalizes the whole vector
#   return(result)
# }
# 
# # so a vector of length 14 with alpha and beta as defined in the doc
# sip <- si(14, 4.29, 1.18) 
# sip
# 
# dim(d)
# ## UPDATED `si` to `sip` just to make the function work
# # reduced iter to 100 for testing
# # Notes: only the line-list data are passed here
# # these are evaluated against EPIC and REPC later one
# MAX_ITER <- 2000
# out1 <- backnow(outcome = d$delay, days = d$report, week = d$week, 
#                 weekend = d$weekend, iter = MAX_ITER, sigma = 0.2, 
#                 maxdelay = 20, si = sip, size = 6)
#  # hmm doesn't converge if you don't know 
# # hmmmmmmm but it does ....... interesting ........
# ## i wonder what will happen
# # so i guess i'll try this on the SCC and see what happens
# # it takes a long time to converge beause you don't have onset data
# 
# back1  <- out1$Back[seq(1001, MAX_ITER, by = 2), ] # after 1000 burn-in
# # wait, if you are going for convergence, why are you taking the distribution?
# gback1 <- geweke.diag(back1)$z
# gback1[is.nan(gback1)] <- 0
# gb1  <- sum(abs(gback1) > 1.96) / length(gback1)
# est1 <- apply(back1, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
# dim(est1)
# 
# 
# ## POST PROCESS
# N_REP <- 1
# N_DAY <- 80
# 
# 
# ## Summarizing curve estimates
# ## Collect the lower, median and upper bounds
# cl11 <- mat.or.vec(N_REP, N_DAY)
# cm11 <- mat.or.vec(N_REP, N_DAY)
# cu11 <- mat.or.vec(N_REP, N_DAY)
# 
# # ## Extract known epidemic curve
# # c1 <- mat.or.vec(N_REP, N_DAY)
# # 
# # ## Extract the report curve
# # rc1 <- mat.or.vec(N_REP, N_DAY)
# 
# for (i in 1:N_REP) {
#   md1 <- minday
# 
#   cl11[i, md1:N_DAY] <- est1[1, ]
#   #cl12[i, md1:N_DAY] <- sce1[[i]]$est2[1, ]
#   cm11[i, md1:N_DAY] <- est1[2, ]
#   #cm12[i, md1:N_DAY] <- sce1[[i]]$est2[2, ]
#   cu11[i, md1:N_DAY] <- est1[3, ]
#   
#   # days1 <- as.numeric(names(sce1[[i]]$epic)) + 20
#   # c1[i, days1] <- sce1[[i]]$epic
#   # 
#   # days1 <- as.numeric(names(sce1[[i]]$repc)) + 20
#   # rc1[i, days1] <- sce1[[i]]$repc
# 
# }
# 
# 
# ######################################################
# ## Additional processing for coverage rate and RMSE ##
# ######################################################
# 
# ccov11 <- mat.or.vec(N_REP, N_DAY)
# 
# ce11 <- mat.or.vec(N_REP, N_DAY)
# 
# 
# for (i in 1:N_REP) {
#   ccov11[i, ] <- as.numeric((c1[i, ] >= cl11[i, ]) & (c1[i, ] <= cu11[i, ]))
# 
#   ce11[i, ] <- (cm11[i, ] - c1[i, ])^2
# 
# }
# 
# ccov11 <- apply(ccov11, 2, sum) / N_REP
# 
# ce11 <- sqrt(apply(ce11, 2, mean))
# 
# 
# date <- rep(1:N_DAY, times = 24)
# result <- c(
#   ccov11, #ccov12, #ccov21, ccov22, 
#   #ccov31, ccov32, ccov41, ccov42, 
#   #ccov51, ccov52, ccov61, ccov62,
#   ce11#, ce12#, #ce21, ce22, 
#   #ce31, ce32, ce41, ce42, 
#   #ce51, ce52, ce61, ce62
# )
# 
# ## UPDATE
# N_SCE <- 1 # max = 6
# N_MOD <- 1 # max = 2
# 
# scenario <- rep(1:N_SCE, each = N_DAY * N_MOD, times = N_MOD)
# model <- rep(1:N_MOD, each = N_DAY, times = N_SCE * N_MOD)
# type <- rep(c("coverage", "RMSE"), each = N_DAY * N_SCE * N_MOD)
# cout1 <- data.frame(cbind(date, result, scenario, model, type))
# 
# ################################
# ## Prepare output for figures ##
# ################################
# 
# cl11 <- apply(cl11, 2, mean)
# #cl12 <- apply(cl12, 2, mean)
# cm11 <- apply(cm11, 2, mean)
# # cm12 <- apply(cm12, 2, mean)
# cu11 <- apply(cu11, 2, mean)
# 
# c1 <- apply(c1, 2, mean)
# 
# rc1 <- apply(rc1, 2, mean)
# 
# 
# scenarios <- c(
#   "Single Delay Dist. & Correct Maximum Delay"#, 
#   #"Single Delay Dist. & Incorrect Maximum Delay",
#   #"Two Delay Dist. & Correct Maximum Delay", 
#   #"Two Delay Dist. & Incorrect Maximum Delay",
#   #"Mutiple Delay Dist. & Correct Maximum Delay", 
#   #"Mutiple Delay Dist. & Incorrect Maximum Delay"
# )
# type <- c("1 Dispersion") #, "2 Dispersions")
# low <- c(cl11)#, cl12) #, cl21, cl22, cl31, cl32, cl41, cl42, cl51, cl52, cl61, cl62)
# med <- c(cm11)#, cm12) #, cm21, cm22, cm31, cm32, cm41, cm42, cm51, cm52, cm61, cm62)
# upp <- c(cu11)#, cu12) #, cu21, cu22, cu31, cu32, cu41, cu42, cu51, cu52, cu61, cu62)
# true <- c(c1)#, c1)#, c2, c2, c3, c3, c4, c4, c5, c5, c6, c6)
# report <- c(rc1)#, rc1)#, rc2, rc2, rc3, rc3, rc4, rc4, rc5, rc5, rc6, rc6)
# 
# outc <- data.frame(rep(1:N_DAY, times = N_SCE * N_MOD), 
#                    low, med, upp, true, report, 
#                    rep(type, each = N_DAY, times = N_SCE), 
#                    rep(scenarios, each = N_DAY * N_MOD * N_SCE), 
#                    rep(1:N_SCE, each = N_DAY * N_MOD * N_SCE))
# 
# colnames(outc) <- c("date", "lower", "median", "upper", 
#                     "epic", "repc", "model", "scenario", "sn")
# 
# head(outc)
# 
# 
# library(tidyverse)
# library(ggpubr)
# 
# outc <- outc %>%
#   filter(date > 20) %>%
#   mutate(date = date - 20)
# 
# ## Assuming the first day is 2020-02-01 and the last day is 2020-03-31
# outc <- outc %>% mutate(date = as.Date("2020-01-31") + date)
# 
# outr <- outr %>%
#   filter(date > 20) %>%
#   mutate(date = date - 20)
# 
# outr <- outr %>% mutate(date = as.Date("2020-01-31") + date)
# 
# outc[outc$scenario == "Single Delay Dist. & Correct Maximum Delay", ]$scenario <- "No Delay Improvement,\nCorrect Maximum Delay"
# # outc[outc$scenario == "Two Delay Dist. & Correct Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Correct Maximum Delay"
# # outc[outc$scenario == "Mutiple Delay Dist. & Correct Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Correct Maximum Delay"
# # outc[outc$scenario == "Single Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "No Delay Improvement & Incorrect Maximum Delay"
# # outc[outc$scenario == "Two Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Incorrect Maximum Delay"
# # outc[outc$scenario == "Mutiple Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Incorrect Maximum Delay"
# 
# outr[outr$scenario == "Single Delay Dist. & Correct Maximum Delay", ]$scenario <- "No Delay Improvement,\nCorrect Maximum Delay"
# # outr[outr$scenario == "Two Delay Dist. & Correct Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Correct Maximum Delay"
# # outr[outr$scenario == "Mutiple Delay Dist. & Correct Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Correct Maximum Delay"
# # outr[outr$scenario == "Single Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "No Delay Improvement & Incorrect Maximum Delay"
# # outr[outr$scenario == "Two Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Incorrect Maximum Delay"
# # outr[outr$scenario == "Mutiple Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Incorrect Maximum Delay"
# 
# 
# out1 <- outc %>% filter(sn == 1 | sn == 2)
# 
# head(out1)
# 
# out1 %>% 
#   pivot_longer(cols = c(repc, epic)) %>%
#   ggplot() +
#   annotate("rect", xmin = as.Date("2020-03-11"), 
#            xmax = as.Date("2020-03-31"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
#   geom_line(aes(x = date, y = value, linetype = name), size = 1.2) +
#   geom_smooth(aes(x = date, y = median, ymax = upper, ymin = lower, 
#                   color = model), 
#               size = 1.2, stat = "identity") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(size = 12),
#     strip.text.x = element_text(size = 12),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 11)
#   ) +
#   labs(x = "Date", y = "Daily Counts") +
#   facet_wrap(~scenario)
# 
# 
# 
# ##########################
# ## Read the processed output
# outc <- readRDS("count.rds")
# outr <- readRDS("rest.rds")
# 
# library(tidyverse)
# library(ggpubr)
# 
# outc <- outc %>%
#   filter(date > 20) %>%
#   mutate(date = date - 20)
# 
# ## Assuming the first day is 2020-02-01 and the last day is 2020-03-31
# outc <- outc %>% mutate(date = as.Date("2020-01-31") + date)
# 
# outr <- outr %>%
#   filter(date > 20) %>%
#   mutate(date = date - 20)
# 
# outr <- outr %>% mutate(date = as.Date("2020-01-31") + date)
# 
# outc[outc$scenario == "Single Delay Dist. & Correct Maximum Delay", ]$scenario <- "No Delay Improvement,\nCorrect Maximum Delay"
# # outc[outc$scenario == "Two Delay Dist. & Correct Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Correct Maximum Delay"
# # outc[outc$scenario == "Mutiple Delay Dist. & Correct Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Correct Maximum Delay"
# # outc[outc$scenario == "Single Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "No Delay Improvement & Incorrect Maximum Delay"
# # outc[outc$scenario == "Two Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Incorrect Maximum Delay"
# # outc[outc$scenario == "Mutiple Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Incorrect Maximum Delay"
# 
# outr[outr$scenario == "Single Delay Dist. & Correct Maximum Delay", ]$scenario <- "No Delay Improvement,\nCorrect Maximum Delay"
# # outr[outr$scenario == "Two Delay Dist. & Correct Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Correct Maximum Delay"
# # outr[outr$scenario == "Mutiple Delay Dist. & Correct Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Correct Maximum Delay"
# # outr[outr$scenario == "Single Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "No Delay Improvement & Incorrect Maximum Delay"
# # outr[outr$scenario == "Two Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Sharp Delay Improvement & Incorrect Maximum Delay"
# # outr[outr$scenario == "Mutiple Delay Dist. & Incorrect Maximum Delay", ]$scenario <- "Gradual Delay Improvement & Incorrect Maximum Delay"
# 
# 
# out1 <- outc %>% filter(sn == 1 | sn == 2)
# 
# head(out1)
# 
# out1 %>% 
#   pivot_longer(cols = c(repc, epic)) %>%
#   ggplot() +
#   annotate("rect", xmin = as.Date("2020-03-11"), 
#            xmax = as.Date("2020-03-31"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
#   geom_line(aes(x = date, y = value, linetype = name), size = 1.2) +
#   geom_smooth(aes(x = date, y = median, ymax = upper, ymin = lower, 
#                   color = model), 
#                size = 1.2, stat = "identity") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(size = 12),
#     strip.text.x = element_text(size = 12),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 11)
#   ) +
#   labs(x = "Date", y = "Daily Counts") +
#   facet_wrap(~scenario)
# 
# 
# 
# ### hmm errors out with the following
# # Error: probabilities must be finite and non-negative .....
# 
# 
# ## how do you check convergence?
# 
# 
# 
# 
# 
# 
# 
# # and then predict the future curve? but you'd have to know what r(t) is
# # which .... doesn't make sense

}
