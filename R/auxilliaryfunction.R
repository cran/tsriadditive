# sone function
sone <- function(survtime, cause = NULL, Z, comp = FALSE, censorsurv = NULL)
{
  #survtime in increasing order
  s_one <- apply(apply(apply(Z, 2, rev), 2, cumsum), 2, rev)
  if(comp){
    censorest <- apply(Z * (cause >= 2) / censorsurv, 2, cumsum) * censorsurv
    censorest <- rbind(rep(0, ncol(censorest)), censorest[1:(nrow(censorest) - 1), ])
    s_one <- s_one + censorest
  }
  s_one
}

# szero function
szero <- function(survtime, cause = NULL, comp = FALSE, censorsurv = NULL)
{
  #survtime in increasing order
  s_zero <- seq(length(survtime), 1)
  if(comp){
    censorest <- cumsum((cause >= 2) / censorsurv) * censorsurv
    censorest <- c(0, censorest[1:(length(censorest) - 1)])
    s_zero <- s_zero + censorest
  }
  s_zero
}

# Zbar function
Zbar <- function(s_zero, s_one)
{
  Z_bar <- s_one/s_zero
  Z_bar
}

# Zint function
Zint <- function(survtime, Z_bar)
{
  timediff <- c(survtime[1], diff(survtime))
  Z_int <- apply(Z_bar * timediff, 2, cumsum)
  Z_int
}

# Gint function
Gint <- function(survtime, censorsurv)
{
  timediff <- c(survtime[1], diff(survtime))
  temp <- cumsum(rev(timediff * censorsurv))
  temp <- c(0, temp[1:(length(temp) - 1)])
  G_int <- rev(temp)
  G_int
}

# GZ_int function
GZint <- function(survtime, Z_bar, censorsurv)
{
  #This is the integration of our censorsurv * Z_bar from somewhere to infinity
  timediff <- c(survtime[1], diff(survtime))
  temp <- apply(apply(timediff * censorsurv * Z_bar, 2, rev), 2, cumsum)
  temp <- rbind(rep(0, ncol(Z_bar)), temp[1:(nrow(temp) - 1), ])
  GZ_int <- apply(temp, 2, rev)
  GZ_int
}

# A scoreprocess function
scoreprocess <- function(Z, Z_bar, cause)
{
  score_process <- (Z - Z_bar) * (cause == 1)
  score_process
}

# Omega inverse function
# This is for inverse of omega, which is part of the sandwich estimator
omegainv <- function(N, survtime, Z, Z_int, cause, comp = FALSE, censorsurv = NULL, G_int = NULL, GZ_int = NULL)
{
  if(!comp) omega <- t(Z) %*% (Z * survtime - Z_int)
  else omega <- t(Z) %*% (Z * survtime - Z_int + (cause >= 2) * (Z * G_int - GZ_int) / censorsurv)
  omega_inv <- solve(omega / N)
  omega_inv
}

# The finite dimensional coefficients estimator
coefest <- function(N, omega_inv, score_process)
{
  coef_est <- omega_inv %*% apply(score_process / N, 2, sum)
  coef_est
}

# The baseline hazards function
baselineest <- function(cause, s_zero, Z_int, coef_est)
{
  baseline_est <- cumsum((cause == 1) / s_zero) - Z_int %*% coef_est
  #baseline_est <- cummax(pmax(baseline_est, 0))
  baseline_est
}

# A hazard prediction function
hazardpred <- function(fit, newtreatment = NULL, newIV = NULL, newcovariates = NULL)
{
  coef_est <- fit$coef
  baseline_est <- fit$baseline
  useIV <- fit$byprod$useIV
  firstfit <- fit$byprod$firstfit
  survtime <- fit$byprod$survtime
  binary <- fit$byprod$binary
  if(!useIV){
    newobsz <- c(newtreatment, newcovariates) 
    score_pred <- sum(newobsz * coef_est)
    hazard_pred = baseline_est + sum(newobsz * coef_est) * survtime
  }
  else{
    newdata = data.frame(1, newIV, t(newcovariates))
    colnames(newdata) <- names(firstfit$coefficients)
    if(binary) newresidual = newtreatment - predict(firstfit, newdata = newdata, type = "response")
    else newresidual = newtreatment - predict(firstfit, newdata)
    newobsz = c(newtreatment, newresidual, newcovariates) 
    score_pred <- sum(newobsz * coef_est)
    hazard_pred <- baseline_est + score_pred * survtime
  }
  hazard_pred <- cummax(pmax(hazard_pred, 0))
  res <- list(score_pred = score_pred,
              hazard_pred = hazard_pred,
              newobsz = newobsz)
  res
}

# A survival prediction function
survivalpred <- function(fit, newtreatment = NULL, newIV = NULL, newcovariates = NULL)
{
  hazard_pred <- hazardpred(fit, newtreatment, newIV, newcovariates)
  survival_pred <- exp(-hazard_pred$hazard_pred)
  res <- append(list(survival_pred = survival_pred), hazard_pred)
  res
}

