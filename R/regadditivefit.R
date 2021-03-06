# regadditivefit function
regadditivefit <- function(survtime, cause, comp = FALSE, treatment = NULL, covariates = NULL)
{
  
  Z <- cbind(treatment, covariates)
  if(!comp) fit <- regsurvadditivefit(survtime, cause, Z)
  else fit <- regcompadditivefit(survtime, cause, Z)
  fit <- betavarest(fit)
  fit
}

# regsuradditivefit function
regsurvadditivefit <- function(survtime, cause, Z = NULL)
{
  N <- length(survtime)
  pZ <- ncol(Z)
  comp <- FALSE
  s_zero <- szero(survtime, cause, comp)
  s_one <- sone(survtime, cause, Z, comp)
  Z_bar <- Zbar(s_zero, s_one)
  Z_int <- Zint(survtime, Z_bar)
  score_process <- scoreprocess(Z, Z_bar, cause)
  omega_inv <- omegainv(N, survtime, Z, Z_int, cause)
  coef_est <- coefest(N, omega_inv, score_process)
  baseline_est <- baselineest(cause, s_zero, Z_int, coef_est)
  byprod <- list(N = N,
                 pZ = pZ,
                 survtime = survtime,
                 Z = Z,
                 cause = cause,
                 s_zero = s_zero,
                 s_one = s_one,
                 Z_bar = Z_bar,
                 Z_int = Z_int,
                 score_process = score_process,
                 omega_inv = omega_inv,
                 useIV = FALSE,
                 comp = FALSE)
  res <- list(coef = coef_est,
              baseline = baseline_est,
              byprod = byprod)
  res
}

# regcompadditivefit function
regcompadditivefit <- function(survtime, cause, Z = NULL)
{
  N <- length(survtime)
  pZ <- ncol(Z)
  comp <- TRUE
  #estimate the censoring distribution by kaplan-meier estimator
  km <- survfit(Surv(survtime, cause == 0) ~ 1,  type = "kaplan-meier", conf.type = "log")
  survest <- stepfun(km$time, c(1, km$surv), right = TRUE)
  censorsurv <- survest(survtime)
  s_zero <- szero(survtime, cause, comp, censorsurv)
  s_one <- sone(survtime, cause, Z, comp, censorsurv)
  Z_bar <- Zbar(s_zero, s_one)
  Z_int <- Zint(survtime, Z_bar)
  G_int <- Gint(survtime, censorsurv)
  GZ_int <- GZint(survtime, Z_bar, censorsurv)
  score_process <- scoreprocess(Z, Z_bar, cause)
  omega_inv <- omegainv(N, survtime, Z, Z_int, cause, comp, censorsurv, G_int, GZ_int)
  coef_est <- coefest(N, omega_inv, score_process)
  baseline_est <- baselineest(cause, s_zero, Z_int, coef_est)
  byprod <- list(N = N,
                 pZ = pZ,
                 survtime = survtime,
                 Z = Z,
                 cause = cause,
                 s_zero = s_zero,
                 s_one = s_one,
                 Z_bar = Z_bar,
                 Z_int = Z_int,
                 G_int = G_int,
                 GZ_int = GZ_int,
                 score_process = score_process,
                 omega_inv = omega_inv,
                 censorsurv = censorsurv,
                 useIV = FALSE,
                 comp = TRUE)
  res <- list(coef = coef_est,
              baseline = baseline_est,
              byprod = byprod)
  res
}