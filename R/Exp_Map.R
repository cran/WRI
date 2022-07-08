#' Numerical implementation of the pointwise Exponential map
#'
#' @description This function implements the Exponential map to calculate \eqn{\hat{F}_t(u)} for a fixed \eqn{u}
#'
#' @param  u The value where the cdf \eqn{\hat{F}} is evaluated
#' @param  forecast A vector that contains the WAR(p) model forecast result
#' @param  cdf The quantile grid used in forecasting
#'
#' @return A numeric value \eqn{\hat{F}_t(u)}
#'
#' @references
#' \cite{Wasserstein Autoregressive Models for Density Time Series, Chao Zhang, Piotr Kokoszka, Alexander Petersen, 2022}
#'
#' @keywords internal

Exp_Map_Barycenter_Method_pw <- function(u, forecast, cdf){

  cutoff <- (forecast <= u)

  temp.1 <- which(cutoff == 1)
  if(is.integer(temp.1) & length(temp.1) == 0L) {
    result <- 0
  }
  else {
    temp.2 <- diff(temp.1)
    if(is.integer(temp.2) & length(temp.2) == 0L) {
      if(temp.1 == 1) {
        result <- cdf[1]
      }
      else {
        result <- cdf[temp.1] - cdf[temp.1-1]
      }
    }
    else {
      temp.3 <- which(temp.2 > 1)
      if(is.integer(temp.3) & length(temp.3) == 0L){
        end <- temp.1[length(temp.1)]
        end.cdf <- cdf[end]
        if(temp.1[1] == 1){
          start.cdf <- 0
        }
        else {
          start.cdf <- cdf[(temp.1[1] - 1)]
        }
      }
      else {
        end <- temp.1[c(temp.3, length(temp.1))]
        start <- temp.1[c(1, temp.3 + 1)] - 1
        if(start[1] == 0){
          start[1] <- 1
          start.cdf <- cdf[start]
          start.cdf[1] <- 0
        }
        else {
          start.cdf <- cdf[start]
        }
        end.cdf <- cdf[end]
      }
      result <- sum(end.cdf - start.cdf)
    }
  }

  return(result)
}


#' Numerical implementation of the Exponential map
#'
#' @description This function implements the Exponential map to calculate \eqn{\hat{F}_t(u)} for all \eqn{u} in the cdf/pdf support
#'
#' @param  density.grid The values where the cdf \eqn{\hat{F}} is evaluated
#' @param  forecast A vector that contains the WAR(p) model forecast result
#' @param  cdf The quantile grid used in forecasting
#'
#' @return A numeric vector that contains \eqn{\hat{F}_t(u)} evaluated over \code{density.grid}
#'
#' @references
#' \cite{Wasserstein Autoregressive Models for Density Time Series, Chao Zhang, Piotr Kokoszka, Alexander Petersen, 2022}
#'
#' @keywords internal

Exp_Map_Barycenter_Method <- function(density.grid, forecast, cdf){
  result<- sapply(density.grid, function(xx) Exp_Map_Barycenter_Method_pw(xx, forecast, cdf))
  return(result)
}
