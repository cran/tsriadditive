#' @title Fitting Additive Hazards Models with Two Stage Residual Inclusion Method
#' @description tsriadditive is used to fit additive hazards models with two stage residual inclusion method.
#' @param survtime the event time
#' @param cause the indicator records the cause. Default to all one. Zero means right censoring. Greater than
#' or equal to two means other cause.
#' @param treatment the treatment variable, can be null
#' @param IV the instrumental variable
#' @param covariates all the observed confounders
#' @return tsriadditive returns an object of class "tsriadditive".
#' An object of class "tsriadditive" is a list containing the following components:
#' \item{coef}{an estimate of the coefficients}
#' \item{baseline}{an estimate of the baseline hazards function}
#' \item{vcov}{an estimate of the variance covariance matrix of coef}
#' \item{byprod}{a byproduct, that will used by other functions}
#' @export tsriadditive
#' @importFrom stats fitted glm lm predict qnorm residuals vcov
#' @importFrom survival survfit
#' @references Ying, A., Xu, R. and Murphy, J. Two-Stage Residual Inclusion for Survival Data and Competing Risks - An
#' Instrumental Variable Approach with Application to SEER- Medicare Linked Data. Statistics in Medicine, 38(10): 1775-1801, 2019.
#' @examples
#' survtime <- rexp(100)
#' cause <- rbinom(100, 1, 0.7)
#' treatment <- rbinom(100, 1, 0.5)
#' IV <- rnorm(100)
#' covariates <- rnorm(100)
#' fit <- tsriadditive(survtime, cause, treatment, IV, covariates)
tsriadditive <- function(survtime, cause = NULL, treatment = NULL, IV = NULL, covariates = NULL)
{
  #if there are ties in surtime, in this function we will just break it down randomly
  #cause: 0 as censored, 1 as event at time, >=2 as other cause event
  #treatment can either be binary or continuous, it can also be None, in which case
  #on IV additive model will be fit even if you pass in IV
  #IV can be None, if so, the non IV additive model will be fit

  #error message leave for later

  #if cause is not provided, we simply assume no censored
  if(is.null(cause)) cause <- rep(1, length(survtime))

  #randomly break ties, order survtime in increasing order
  #reorder all the other variables and save the order into timeorder
  timeorder <- order(survtime, sample(length(survtime)))
  survtime <- survtime[timeorder]
  cause <- cause[timeorder]
  treatment <- treatment[timeorder]
  IV <- IV[timeorder]
  #convert the covariates into matrix object just for later use if it is only dimensional one
  if(is.vector(covariates)) covariates <- matrix(covariates, ncol = 1)
  covariates <- covariates[timeorder, ]
  #comp is the indicator for whether there are competing risks
  comp = !all(cause <= 1)

  #if treatment or IV is None, that means we want to perform regular fit
  if(is.null(treatment) || is.null(IV)) fit <- regadditivefit(survtime, cause, comp, treatment, covariates)
  #so else we have both treatment and IV, we will use 2sri
  else fit <- tsriadditivefit(survtime, cause, comp = comp, treatment, IV, covariates)

  fit$byprod$timeorder <- timeorder
  result <- fit
  class(result) <- "tsriadditive"
  result
}



#' @title Summarizing Additive Hazards Model with Two Stage Residual Inclusion Method Fits
#' @param object an object of class "tsriadditive", usually, a result of a call to tsriadditive.
#' @param x an object of class "summary.tsriadditive", usually, a result of a call to summary.tsriadditive.
#' @param ... further arguments passed to or from other methods.
#' @description summary method for class "tsriadditive".
#' @export tsriadditive
#' @export summary.tsriadditive
#' @importFrom stats pnorm
#' @method summary tsriadditive
#' @S3method summary tsriadditive
#' @examples
#' survtime <- rexp(100)
#' cause <- rbinom(100, 1, 0.7)
#' treatment <- rbinom(100, 1, 0.5)
#' IV <- rnorm(100)
#' covariates <- rnorm(100)
#' fit <- tsriadditive(survtime, cause, treatment, IV, covariates)
#' summary(fit)
#' @return print.summary.lm tries to be smart about formatting coefficients, an estimated variance covariance matrix of
#' the coeffieients, Z-values and the corresponding P-values
summary.tsriadditive <- function(object, ...)
{
  fit = object
  result <- list(Estimate = fit$coef,
                 Std.Error = sqrt(diag(fit$vcov)),
                 zvalue = fit$coef/sqrt(diag(fit$vcov)),
                 Pvalue = 2 * (1 - pnorm(abs(fit$coef)/sqrt(diag(fit$vcov)))))
  class(result) <- "summary.tsriadditive"
  result
}

#' @rdname summary.tsriadditive
#' @S3method print summary.tsriadditive
print.summary.tsriadditive <- function(x, ...){
  result = data.frame(x$Estimate, x$Std.Error, x$zvalue, x$Pvalue)
  colnames(result) <- c("Estimate", "Std. Error", "z value", "Pr(>|t|)")
  print(result)
}

#' @title Predict method for Additive Hazards Model with Two Stage Residual Inclusion Method Fits
#' @description Predicted values based on tsriadditive object.
#' @param object an object of class "tsriadditive", usually, a result of a call to tsriadditive.
#' @param newtreatment a new treatment value.
#' @param newIV a new instrumental variable value.
#' @param newcovariates a new observed covariates.
#' @param ... further arguments passed to or from other methods.
#' @export tsriadditive
#' @export predict.tsriadditive
#' @method predict tsriadditive
#' @S3method predict tsriadditive
#' @return predict.tsriadditive produces a venctor of predictions based on new values.
#' A list with the following components is returned:
#' \item{newobsz}{the vector grouping newtreatment, new IV and newcovariates}
#' \item{score_pred}{the predicted scores}
#' \item{hazard_pred}{the predicted baseline hazards function}
#' \item{surival_pred}{the predicted surival function}
#' @examples
#' survtime <- rexp(100)
#' cause <- rbinom(100, 1, 0.7)
#' treatment <- rbinom(100, 1, 0.5)
#' IV <- rnorm(100)
#' covariates <- rnorm(100)
#' fit <- tsriadditive(survtime, cause, treatment, IV, covariates)
#' predict(fit, 1, 0, 0)
predict.tsriadditive <- function(object, newtreatment = NULL, newIV = NULL, newcovariates = NULL, ...)
{
  fit = object
  res = survivalpred(fit, newtreatment, newIV, newcovariates)
  res
}



#' @title Plotting Predicted Survival Function or Cumulative Incidence Function with Pointwise Confidence Intervals
#' @param x the fitting object after fitting our model
#' @param newtreatment a new treatment value
#' @param newIV a new instrumental variable value
#' @param newcovariates a new observed covariates
#' @param alpha the confidence level 1 - alpha for confidence interval
#' @param unit the time unit we focus
#' @param ... the other arguments you want to put in the built-in plot function
#' @description The function will plot the predicted survival function when fitting a survival model and
#' the predicted cumulative incidence function when fitting a competing risks model.
#' Corresponding pointwise confidence intervals at level alpha are also included.
#' @export tsriadditive
#' @export plot.tsriadditive
#' @importFrom stats stepfun
#' @importFrom graphics plot
#' @method plot tsriadditive
#' @S3method plot tsriadditive
#' @return No return value, called for side effects
#' @examples
#' survtime <- rexp(100)
#' cause <- rbinom(100, 1, 0.7)
#' treatment <- rbinom(100, 1, 0.5)
#' IV <- rnorm(100)
#' covariates <- rnorm(100)
#' fit <- tsriadditive(survtime, cause, treatment, IV, covariates)
#' plot(fit, 1, 0, 0)
plot.tsriadditive <- function(x, newtreatment = NULL, newIV = NULL, newcovariates = NULL, alpha = 0.05, unit = "", ...)
{
  fit <- x
  comp <- fit$byprod$comp
  pred <- survivalpred(fit, newtreatment, newIV, newcovariates)
  survival_pred <- pred$survival_pred
  hazard_pred <- pred$hazard_pred
  newobsz <- pred$newobsz
  survtime <- fit$byprod$survtime
  if(!comp){
    survcurve <- stepfun(fit$byprod$survtime, c(1, survival_pred), right = TRUE)
    plot(survcurve, verticals = TRUE, do.points = FALSE, xlim = c(0, max(survtime)),
         ylim = c(0, 1), lwd = 1.2, xaxs = 'i', yaxs = 'i',
         xlab = paste0("time (", unit, ")"), ylab = "Survival probability", main = NULL, ...)
    fit <- survprobconfint(hazard_pred, newobsz, fit, alpha)
    survuppercurve <- stepfun(fit$byprod$survtime, c(1, fit$byprod$upper), right = TRUE)
    plot(survuppercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3,...)
    survlowercurve <- stepfun(fit$byprod$survtime, c(1, fit$byprod$lower), right = TRUE)
    plot(survlowercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3, ...)
  }
  else{
    distribution_pred <- 1 - survival_pred
    survcurve <- stepfun(fit$byprod$survtime, c(0, distribution_pred), right = TRUE)
    fit <- survprobconfint(hazard_pred, newobsz, fit, alpha)
    survuppercurve <- stepfun(fit$byprod$survtime, c(0, 1 - fit$byprod$upper), right = TRUE)
    survlowercurve <- stepfun(fit$byprod$survtime, c(0, 1 - fit$byprod$lower), right = TRUE)
    plot(survcurve, verticals = TRUE, do.points = FALSE, xlim = c(0, max(survtime)),
         ylim = c(0, max(1 - fit$byprod$lower)), lwd = 1.2, xaxs = 'i', yaxs = 'i',
         xlab = paste0("time (", unit, ")"), ylab = "Cumulative incidence function", main = NULL, ...)
    plot(survuppercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3, ...)
    plot(survlowercurve, verticals = TRUE, do.points = FALSE, add = TRUE, lty = 3, ...)
    #fit$byprod <- append(fit$byprod, list(survcurve = survcurve,
          #                                survuppercurve = survuppercurve,
          #                                survlowercurve = survlowercurve))
  }
  invisible()
}

