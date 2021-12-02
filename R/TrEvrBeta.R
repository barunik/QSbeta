#' @include QsBeta.R
#' @import quantspec
#' @import pbivnorm
NULL

################################################################################
#' Compute Quantile Spectral Betas
#'
#' This function is high-level implementation of \code{getBeta} function and
#' computes general quantile spectral (QS) betas between factor
#' return and asset return for frequencies in (0,pi)
#' for tau1 and tau2 quantile thresholds.
#'
#' @name compQsBeta
#'
#' @param factor time-series of returns of some factor
#' @param asset time-series of returns of some asset
#' @param tau.1 quantile threshold level for the factor's return
#' @param tau.2 quantile threshold level for the asset's return
#'
#' @return Returns vector of QS betas for frequencies between 0 and pi.
#'
#' @export
################################################################################

compQsBeta <- function(factor, asset, tau1, tau2) {
  
  # Estimate smoothed quantile periodogram
  Y <- matrix(c(factor, asset), ncol = 2)
  qqs <- c(tau1, tau2)
  w <- kernelWeight(W = W1, b = 0.5 * dim(Y)[1]^(-1/4))
  sPG <- smoothedPG(Y, levels.1 = qqs, weight = w)
  Beta <- getBeta(sPG)
  freqs <- 1:round(dim(Y)/2)[1]
  
  # Return QS betas
  Re(Beta[freqs,1,1,2,2,1])
}

################################################################################
#' Compute Tail Risk Betas
#'
#' This function applies the function \code{getQsBeta} to compute
#' tail risk (TR) quantile spectral (QS) betas
#' between asset return and factor return for long and short horizon
#' for tau1 and tau2 quantile thresholds. If tau2 is not provided,
#' then it is calculated so that the factor's threshold value is the same
#' as the asset's threshold value.
#'
#' @name compTrBeta
#'
#' @param factor time-series of returns of some factor
#' @param asset time-series of returns of some asset
#' @param tau.1 quantile threshold level for the factor's return
#' @param tau.2 quantile threshold level for the asset's return
#' @param samp.freq sampling frequency of the data in terms
#'                  of number of observations during a year    
#' @param long.horizon length of the long horizon cycle in years
#' @param reference time series from which the threshold
#'                         value for the asset should be computed
#'
#' @return Returns vector of TR QS betas for long and short horizon.
#'
#' @export
################################################################################
compTrBeta <- function(factor, asset, tau1, tau2 = NA,
                       samp.freq = 12, long.horizon = 3,
                       reference = factor) {
  
  # If tau2 is not provided, compute it from the reference time series
  if (missing(tau2)) {
    asset.threshold <- quantile(reference, probs = tau1)
    P <- ecdf(asset)
    tau2 <- P(asset.threshold)
  }
  
  # Compute QS betas for all frequencies between 0 and pi
  betas <- compQsBeta(factor, asset, tau1, tau2)
  
  # Average betas over long and short horizon frequencies
  T <- round(length(factor)/2) # half-length of the time series
  cutoff <- floor(1/(long.horizon * samp.freq) * T) # long and short horizon
  betalong <- mean(betas[1:cutoff])
  betashort <- mean(betas[(cutoff + 1):length(betas)])
  
  list(beta.long = betalong, beta.short = betashort)
  
}

################################################################################
#' Compute Relative Tail Risk Betas
#'
#' This function applies the function \code{getQsBeta} to compute
#' tail risk (TR) quantile spectral (QS) betas
#' between asset return and factor return for long and short horizon
#' for tau1 and tau2 quantile thresholds. If tau2 is not provided,
#' then it is calculated so that the factor's threshold value is the same
#' as the asset's threshold value.
#'
#' @name compRelTrBeta
#'
#' @param factor time-series of returns of some factor
#' @param asset time-series of returns of some asset
#' @param tau.1 quantile threshold level for the factor's return
#' @param tau.2 quantile threshold level for the asset's return
#' @param samp.freq sampling frequency of the data in terms
#'                  of number of observations during a year    
#' @param long.horizon length of the long horizon cycle in years
#' @param reference time series from which the threshold
#'                         value for the asset should be computed
#'
#' @return Returns vector of relative TR QS betas for long and short horizon.
#'
#' @export
################################################################################
compRelTrBeta <- function(factor, asset, tau1, tau2 = NA,
                          samp.freq = 12, long.horizon = 3,
                          reference = factor) {
  
  # If tau2 is not provided, compute it from the reference time series
  if (missing(tau2)) {
    asset.threshold <- quantile(reference, probs = tau1)
    P <- ecdf(asset)
    tau2 <- P(asset.threshold)
  }
  
  # Compute long and short TR betas
  trBetas <- compTrBeta(factor, asset, tau1 = tau1, tau2 = tau2, samp.freq = samp.freq,
                        long.horizon = long.horizon, reference = reference)
  
  # Compute QS betas under Normality assumption
  rho <- cor(factor, asset)
  copula <- pbivnorm(qnorm(tau1,0,1), qnorm(tau2,0,1), rho)
  beta_gauss <- (copula - tau1*tau2)/(tau1*(1-tau1))
  
  betalong <- trBetas$beta.long - beta_gauss
  betashort <- trBetas$beta.short - beta_gauss
  
  list(beta.long = betalong, beta.short = betashort)
  
}

################################################################################
#' Compute Extreme Volatility Risk Betas
#'
#' This function applies the function \code{getQsBeta} to compute
#' extreme volatility risk (EVR) quantile spectral (QS) betas
#' between asset return and negative differences of factor variance
#' for long and short horizon for tau1 and tau2 quantile thresholds.
#' If tau2 is not provided, then it is calculated so that the asset's
#' threshold is tau1-quantile of factor's (reference's) return
#'
#' @name compEvrBeta
#'
#' @param factor time-series of returns of some factor
#'               for which we compute the variance 
#' @param asset time-series of returns of some asset
#' @param tau.1 quantile threshold level for the factor's return
#' @param tau.2 quantile threshold level for the assets's return
#' @param samp.freq sampling frequency of the data in terms
#'                  of number of observations during a year    
#' @param long.horizon length of the long horizon cycle in years
#' @param reference from which time series should the threshold
#'                         value for the asset be computed
#' @param garch.order order of the GARCH(p, q) model which should
#'                    be estimated on the factor series
#' @param mean.mod order of the ARMA(p, q) model which should
#'                    be estimated on the factor series
#' @param distr.mod distribution model of the error term in GARCH,
#'                  for more info, consult the rugarch package
#'
#' @return Returns vector of EVR QS betas for long and short horizon.
#'
#' @export
################################################################################
compEvrBeta <- function(factor, asset, tau1, tau2 = NA,
                        samp.freq = 12, long.horizon = 3,
                        reference = factor, garch.order = c(1, 1),
                        mean.mod = c(1, 0), distr.mod = "norm") {
  
  # Estimate the variance model
  model <- ugarchspec(variance.model = list(garchOrder = garch.order),
                      mean.model = list(armaOrder = mean.mod, include.mean = T),
                      distribution.model = distr.mod)
  fit <- ugarchfit(model, data = factor)
  var.diff.neg <- -diff(fit@fit$var)
  
  # Make both series of the same length
  asset <- asset[-1]
  
  # If tau2 is not provided, compute it from the reference time series
  if (missing(tau2)) {
    asset.threshold <- quantile(reference, probs = tau1)
    P <- ecdf(asset)
    tau2 <- P(asset.threshold)
  }
  
  # Compute QS betas
  betas <- compTrBeta(var.diff.neg, asset, tau1 = tau1,
                      tau2 = tau2, samp.freq = samp.freq,
                      long.horizon = long.horizon)
  
  list(beta.long = betas[[1]], beta.short = betas[[2]])
  
}

