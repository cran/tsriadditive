
betaconfint <- function(coef_est, vcov, alpha)
{
  upper <- coef_est + qnorm(1 - alpha/2) * sqrt(diag(vcov))
  lower <- coef_est - qnorm(1 - alpha/2) * sqrt(diag(vcov))
  res <- list(upper = upper,
              lower = lower)
  res
}


upperconfint <- function(hazard_pred, hazardpredvar_est, newobsz, alpha)
{
  upperconf_int <- pmin(1, exp(-(hazard_pred - qnorm(1 - alpha/2) * sqrt(hazardpredvar_est))))
  upperconf_int
}


lowerconfint <- function(hazard_pred, hazardpredvar_est, newobsz, alpha)
{
  lowerconf_int <- exp(-(hazard_pred + qnorm(1 - alpha/2) * sqrt(hazardpredvar_est)))
  lowerconf_int
}


survprobconfint <- function(hazard_pred, newobsz, fit = NULL, alpha)
{
  fit <- hazardpredvarest(newobsz, fit)
  hazardpredvar_est <- fit$hazardpredvar_est
  upper <- upperconfint(hazard_pred, hazardpredvar_est, newobsz, alpha)
  lower <- lowerconfint(hazard_pred, hazardpredvar_est, newobsz, alpha)
  res <- list(upper = upper,
              lower = lower)
  fit$byprod <- append(fit$byprod, res)
  fit
}
