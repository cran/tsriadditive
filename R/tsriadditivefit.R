# fit an additive hazard using two stage residual inclusion method
tsriadditivefit <- function(survtime, cause, comp = FALSE, treatment = NULL, IV = NULL, covariates = NULL)
{
  #first stage
  binary <- all(treatment == 0 | treatment == 1)
  #if binary is TRUE we treatment as binary treatment o/w continuous
  if(binary){
    firstfit <- glm(treatment ~. , data = data.frame(cbind(treatment, IV, covariates)), family = "binomial")
    residual <- residuals(firstfit, type = "response")
  }
  else{
    firstfit <- lm(treatment ~. , data = data.frame(cbind(treatment, IV, covariates)))
    residual <- residuals(firstfit, type = "response")
  }
  X <- cbind(1, IV, covariates)
  Z <- cbind(treatment, residual, covariates)
  if(!comp) fit <- tsrisurvadditivefit(survtime, cause, Z)
  else fit <- tsricompadditivefit(survtime, cause, Z)
  firststage <- list(X = X, firstfit = firstfit)
  fit$byprod <- append(fit$byprod, firststage)
  fit$byprod$binary <- binary
  fit <- betavarest(fit)
  fit
}

# fit an additive hazard using IV method under survival settings
tsrisurvadditivefit <- function(survtime, cause = NULL, Z)
{
  fit <- regsurvadditivefit(survtime, cause, Z)
  fit$byprod$useIV = TRUE
  fit
}

# fit an additive hazard without using IV method under competing risks settings
tsricompadditivefit <- function(survtime, cause = NULL, Z)
{
  fit <- regcompadditivefit(survtime, cause, Z)
  fit$byprod$useIV <- TRUE
  fit
}
