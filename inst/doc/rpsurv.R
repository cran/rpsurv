## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.align = "center",
  fig.width = 6, fig.height = 4.2, out.width = "80%"
)
suppressMessages({
  library(survival)
  library(rpsurv)
})

## ----fit-basic----------------------------------------------------------------
data(brcancer)
brcancer$hormon <- as.numeric(brcancer$hormon)

fit <- rpsurv(Surv(rectime, censrec) ~ hormon, data = brcancer, df = 4, scale = "hazard")
summary(fit)

## ----predict-plot, fig.cap = "Predicted survival by hormonal therapy status."----
plot(fit, newdata = data.frame(hormon = c(0, 1)), col = c("steelblue", "firebrick"),
     main = "rpsurv: predicted survival")

## ----tve-fit------------------------------------------------------------------
fit_tve <- rpsurv(Surv(rectime, censrec) ~ hormon, data = brcancer,
                   df = 4, tve = "hormon", tve.df = 2)
coef(fit_tve)

## ----tve-lrt------------------------------------------------------------------
lrt_stat <- 2 * (fit_tve$loglik - fit$loglik)
lrt_df <- fit_tve$df - fit$df
pchisq(lrt_stat, lrt_df, lower.tail = FALSE)

## ----tvc-example--------------------------------------------------------------
d <- brcancer
d$id <- seq_len(nrow(d))
half <- d$rectime / 2
first  <- data.frame(id = d$id, start = 0, stop = half, status = 0,
                      hormon = d$hormon, x1 = d$x1)
second <- data.frame(id = d$id, start = half, stop = d$rectime, status = d$censrec,
                      hormon = d$hormon, x1 = d$x1 + 1)  # covariate value changes
long <- rbind(first, second)
long <- long[long$start < long$stop, ]

fit_cp <- rpsurv(Surv(start, stop, status) ~ hormon + x1, data = long, df = 4)
fit_cp$counting
coef(fit_cp)

## ----diagnostics, fig.cap = "Cox-Snell residual check: points should lie near the diagonal."----
plot(fit, type = "coxsnell")

## ----deviance-resid-----------------------------------------------------------
dev_resid <- residuals(fit, type = "deviance")
summary(dev_resid)

## ----km-compare, fig.cap = "Model vs. Kaplan-Meier, by hormonal therapy status."----
plot(fit, type = "km", by = "hormon")

## ----validation, eval = requireNamespace("rstpm2", quietly = TRUE)------------
suppressPackageStartupMessages(library(rstpm2))
ref <- stpm2(Surv(rectime, censrec) ~ hormon, data = brcancer, df = 4)
data.frame(
  term = c("hormon", "logLik"),
  rpsurv = c(coef(fit)["hormon"], fit$loglik),
  stpm2  = c(coef(ref)["hormon"], -ref@min)
)

## ----benchmark-table, echo = FALSE--------------------------------------------
res <- readRDS(system.file("extdata", "benchmark_results.rds", package = "rpsurv"))
knitr::kable(res, digits = 3, col.names = c("n", "rpsurv (s)", "stpm2 (s)", "flexsurvspline (s)"))

## ----benchmark-fig, echo = FALSE, out.width = "90%", fig.cap = "Fit time vs. sample size (log-log axes)."----
knitr::include_graphics("figures/benchmark_speed.png")

## ----benchmark-extended, echo = FALSE-----------------------------------------
ext <- readRDS(system.file("extdata", "benchmark_extended.rds", package = "rpsurv"))
ext$speedup <- round(ext$rstpm2 / ext$rpsurv, 1)
knitr::kable(ext, digits = 3,
             col.names = c("model", "n", "rpsurv (s)", "stpm2 (s)", "speedup"))

## ----benchmark-code, eval = FALSE---------------------------------------------
# simulate_data <- function(n, seed = 1) {
#   set.seed(seed)
#   x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#   lp <- 0.5 * x1 - 0.3 * x2
#   u <- runif(n)
#   event_time <- (-log(u) / (0.01 * exp(lp)))^(1 / 1.2)
#   cens_time <- rexp(n, 0.008)
#   data.frame(time = pmin(event_time, cens_time),
#              status = as.numeric(event_time <= cens_time), x1 = x1, x2 = x2)
# }
# d <- simulate_data(200000)
# system.time(rpsurv(Surv(time, status) ~ x1 + x2, data = d, df = 4))
# system.time(stpm2(Surv(time, status) ~ x1 + x2, data = d, df = 4))
# system.time(flexsurvspline(Surv(time, status) ~ x1 + x2, data = d, k = 3))

## ----session-info-------------------------------------------------------------
sessionInfo()

