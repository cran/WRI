#' Function for calculating sample autocovariance
#'
#' @description Calculating sample autocovariance of a specified lag for a centered time series
#'
#' @param  ts A numeric vector consisting of sequentially collected data
#' @param  lag A positive integer indicates the lag of the autocovariance function
#'
#' @return Sample autocovariance
#'
#' @keywords internal

Sample_ACV <- function(ts, lag) {
  l <- length(ts)
  sub.ts.1 <- ts[1:(l-lag)]
  sub.ts.2 <- ts[(1+lag):l]
  result <- sum(sub.ts.1 * sub.ts.2)/l
  return(result)
}


#' Function for calculating sample Wasserstein autocovariance functions
#'
#' @description This function uses a time series of quantile functions to calculate the sample Wasserstein autocovariance functions from order \eqn{0} to \eqn{p} with a specified training window
#'
#' @param  end.day A positive integer, the last index of the training window.
#' @param  training.size A positive integer, the size of the training widnows.
#' @param  quantile A matrix containing all the available quantile functions. Columns represent time indices and rows represent evaluation grid.
#' @param  quantile.grid A numeric vector, the grid over which quantile functions are evaluated.
#' @param  p A positive integer, the maximum order of the sample Wasserstein autocovariance functions.
#'
#' @return A list with
#' \itemize{
#'   \item acvfs - The sample Wasserstein autocovariance functions from order \eqn{0} to \eqn{p}
#'   \item barycenter - The sample average of the quantile functions in the training window
#'   \item quantile.pred - The quantile functions from \eqn{end.day - p + 1} to \eqn{end.day}
#' }
#'
#' @keywords internal

WARp_acvfs <- function(end.day, training.size, quantile, quantile.grid, p){
  start.day <- max(c(0,end.day - training.size + 1))

  # get a subset of quantile functions for training
  quantile.training <- quantile[, start.day:end.day]
  barycenter <- apply(quantile.training, 1, mean)
  centered.quantile <- quantile.training - barycenter

  # sample covariance functions over the grid
  p.w.cov <- lapply(c(0:p),
                    function(xx) apply(centered.quantile,
                                       1,
                                       function(xxx) Sample_ACV(xxx, xx)))

  # interpolate
  n <- length(quantile.grid)
  scale <- quantile.grid[n]
  f.cov <- lapply(p.w.cov, function(xx) splinefun(quantile.grid, xx, method="natural"))

  # take integral
  lambda.hat <- lapply(f.cov, function(xx) integrate(xx, 0, scale, stop.on.error=FALSE))

  # Wasserstein acvfs
  acvfs <- vector(mode="numeric", length=p+1)
  for(i in 1:(p+1)) {
    acvfs[i] <- lambda.hat[[i]]$value
  }
  quantile.pred <- quantile[,(end.day - p + 1):end.day]
  return(list(acvfs, barycenter, quantile.pred))
}

#' Forecast using WAR(p) models
#'
#' @description This is a simplified function that produces the one-step ahead forecasts from WAR(p) models in the tangent space
#'
#' @param  p A positive integer, the order of the WAR(p) model.
#' @param  acvfs A list that is the output of \code{WARp_acvfs()}
#'
#' @return A numeric vector, the one-step ahead forecast produced by the WAR(p) model in the tangent space.
#'
#' @references
#' \cite{Wasserstein Autoregressive Models for Density Time Series, Chao Zhang, Piotr Kokoszka, Alexander Petersen, 2022}
#'
#' @keywords internal

WARp_forecast_tangent <- function(p, acvfs){

  acvf.value <- acvfs[[1]]
  barycenter <- acvfs[[2]]
  last.days <- as.matrix(acvfs[[3]])

  gamma <- matrix(acvf.value[2:(p+1)], nrow = p, ncol = 1, byrow = TRUE)

  # generate elements to populate the upper and lower triangles
  if(p == 1){
    Gamma <- matrix(acvf.value[1],1,1)
  }
  else{
    members <- vector()
    temp.acvf.value <- acvf.value[2:(p+1)]
    for(j in c(1:(p-1))) {
      members <- c(members, temp.acvf.value[1:(p-j)])
    }

    Gamma <- matrix(rep(0,p^2), nrow = p, ncol = p)
    Gamma[lower.tri(Gamma)] <- members
    Gamma <- t(Gamma)
    Gamma[lower.tri(Gamma)] <- members

    diag(Gamma) <- rep(acvf.value[1], times = p)
  }

  phi.hat.matrix <- solve(Gamma, gamma)

  N <- dim(last.days)[2]

  forecast <- (last.days[,N:(N-p+1)] - barycenter) %*% phi.hat.matrix + barycenter

  return(list(forecast,barycenter))
}


#' Calculating innovations in WAR(p) models
#'
#' @description This function calculates innovations in WAR(p) models
#'
#' @param  quantile A matrix containing all the available quantile functions. Columns represent time indices and rows represent evaluation grid.
#' @param  quantile.grid A numeric vector, the grid over which quantile functions are evaluated.
#' @param  p A positive integer, the order of the fitted WAR(p) model.
#'
#' @return A list with
#' \itemize{
#'   \item innovation - The tangent space innovations evaluated over the quantile grid.  Fitting a WAR(p) model for \eqn{n} observations will produce \eqn{n-2p} innovations.
#'   \item cov.surf - The covariance surface
#' }
#'
#' @references
#' \cite{Wasserstein Autoregressive Models for Density Time Series, Chao Zhang, Piotr Kokoszka, Alexander Petersen, 2022}
#'
#' @keywords internal

getInnovation <- function(quantile, quantile.grid, p){

  sample.size <- ncol(quantile)

  innovation.mat <- matrix(0, length(quantile.grid), sample.size-2*p)

  for (i in (2*p):(sample.size-1)) {
    sub.acvf <- WARp_acvfs(i, i, quantile, quantile.grid, p)

    sub.forecast <- WARp_forecast_tangent(p, sub.acvf)

    innovation <- quantile[, i+1] - sub.forecast[[1]]

    innovation.mat[, i-2*p+1] <- innovation
  }

  # Calculate covariance surface
  cov.surf <- fdapace::GetCovSurface(lapply(1:(sample.size-2*p),
                                            function(x) innovation.mat[, x]),
                                     lapply(1:(sample.size-2*p),
                                            function(x) quantile.grid))
  return(list("innovations"=innovation.mat,
              "cov"=cov.surf))
}



#' WAR(p) models: estimation and forecast
#'
#' @description this function produces an object of the WARp class which includes WAR(p) model parameter estimates and relevant quantities (see output list)
#'
#' @details This function takes in a density time series in the form of the corresponding quantile functions as the main input. If the quantile series is not readily available, a general practice is to estimate density functions from samples, then use \code{dens2quantile} from the \code{fdadensity} package to convert density time series to quantile series.
#'
#' @param  quantile A matrix containing all the sample quantile functions. Columns represent time indices and rows represent evaluation grid.
#' @param  quantile.grid A numeric vector, the grid over which quantile functions are evaluated.
#' @param  p A positive integer, the order of the fitted WAR(p) model.
#'
#' @return A \code{WARp} object of:
#' \item{coef}{estimated AR parameters of the fitted WAR(p) model}
#' \item{coef.cov}{covariance matrix of \code{coef}}
#' \item{acvf}{Wasserstein autocovariance function values}
#' \item{Wass.mean}{Wasserstein mean quantile function}
#' \item{quantile}{a matrix containing all the sample quantile functions (columns represent time indices and rows represent evaluation grid)}
#' \item{quantile.grid}{quantile function grid that is utilized in calculation}
#' \item{order}{a positive integer, the order of the fitted WAR(p) model}
#'
#' @references
#' \cite{Wasserstein Autoregressive Models for Density Time Series, Chao Zhang, Piotr Kokoszka, Alexander Petersen, 2022}
#'
#'
#' @examples
#' # Simulate a density time series represented in quantile functions
#' # warSimData$sample.ts: A sample TS of quantile functions of length 100, taken from
#' #            the simulation experiments in Section 4 of Zhang et al. 2022.
#'
#' # warSimData$quantile.grid: The grid over which quantile functions in sample.ts are evaluated.
#'
#' warSimData <- warSim()
#'
#' p <- 3
#' dSup <- seq(-2, 2, 0.02)
#' expSup <- seq(-2, 2, 0.1)
#'
#' # Estimation: fit a WAR(3) model
#' WARp_obj <- WARp(warSimData$sample.ts, warSimData$quantile.grid, p)
#'
#' # Forecast: one-step-ahead forecast
#' forecast_1 <- predict(WARp_obj)               # dSup and expSup are chosen automatically
#' forecast_2 <- predict(WARp_obj, dSup, expSup) # dSup and expSup are chosen by user
#'
#' # Plots
#' par(mfrow=c(1,2))
#'
#' plot(forecast_1$dSup, forecast_1$pred.cdf, type="l", xlab="dSup", ylab="cdf")
#' plot(forecast_1$dSup, forecast_1$pred.pdf, type="l", xlab="dSup", ylab="pdf")
#'
#' plot(forecast_2$dSup, forecast_2$pred.cdf, type="l", xlab="dSup", ylab="cdf")
#' plot(forecast_2$dSup, forecast_2$pred.pdf, type="l", xlab="dSup", ylab="pdf")
#'
#'
#' @export

WARp <- function(quantile, quantile.grid, p){
  sample.size <- ncol(quantile)
  acvf_list <- WARp_acvfs(sample.size, sample.size, quantile, quantile.grid, p)

  # Wasserstein auto covariance functions
  acvf.value <- acvf_list[[1]]
  # Wasserstein mean
  barycenter <- acvf_list[[2]]
  # The last p days of the times series that will be used in prediction
  last.days <- as.matrix(acvf_list[[3]])

  # Prepare components for Yule-Walker estimation
  gamma <- matrix(acvf.value[2:(p+1)], nrow=p, ncol=1, byrow=TRUE)

  # Populate the LHS of Yule-Walker estimation equation
  if( p == 1){
    Gamma <- matrix(acvf.value[1],1,1)
  }
  else{
    members <- vector()
    temp.acvf.value <- acvf.value[2:(p+1)]
    for(j in c(1:(p-1))) {
      members <- c(members, temp.acvf.value[1:(p-j)])
    }

    Gamma <- matrix(rep(0,p^2), nrow = p, ncol = p)
    Gamma[lower.tri(Gamma)] <- members
    Gamma <- t(Gamma)
    Gamma[lower.tri(Gamma)] <- members

    diag(Gamma) <- rep(acvf.value[1], times = p)
  }

  # Solve for autoregressive coefficients
  phi.hat.matrix <- solve(Gamma, gamma)

  # Get autoregressive coefficient covariance matrices
  # Root(s) of autoregressive polynomial
  roots <- Re(polyroot(c(1, -phi.hat.matrix[, 1])))

  coef <- 1/roots

  max.power <- 20
  coef.power <- 0:max.power

  psi.poly <- polynom::polynomial(c(1))

  for (i in 1:p) {
    coef.poly <- polynom::polynomial(coef[i]^coef.power)

    psi.poly <- psi.poly * coef.poly

  }

  psi <- coefficients(psi.poly)
  size.psi <- length(psi)

  psi.sum <- sapply(0:(p-1),
                    function(x){
                      sum(psi[1:(size.psi-x)] * psi[(1+x):size.psi])
                    })

  check.small.roots <- sapply(psi.sum, is.infinite)

  if (sum(check.small.roots) > 0) {
    print("Unable to calculate Wasserstein autoregressive coefficient(s) due
          to small roots.")
  }

  psi.mat <- matrix(0, p, p)

  for(i in 1:p) {
    psi.mat[i, i:p] <- psi.sum[1:(p-i+1)]
  }

  psi.mat <- psi.mat + t(psi.mat) - diag(psi.mat)*diag(1,p)

  # Calculate sigma^2_epsilon in (3.13) Zhang et al. (2022)

  cov.list <- getInnovation(quantile, quantile.grid, p)

  cov <- cov.list[[2]]$cov

  base <- diff(quantile.grid)

  sigma.epsilon.sqrd.num <-  sum(cov[-1, -1]^2 * outer(base, base))

  sigma.epsilon.sqrd.denom <- (sum(diag(cov)[-1] * base))^2

  sigma.epsilon.sqrd <- sigma.epsilon.sqrd.num/sigma.epsilon.sqrd.denom

  coef.cov <- sigma.epsilon.sqrd * solve(psi.mat)

  # Labeling outputs
  names(coef) <- paste0("beta_", 1:p)
  colnames(coef.cov) <- paste0("beta_", 1:p)
  rownames(coef.cov) <- paste0("beta_", 1:p)

  names(acvf.value) <- paste0("lag_", 0:p)

  return(structure(list(
    "coef"=phi.hat.matrix[, 1],
    "coef.cov"=coef.cov,
    "acvf"=acvf.value,
    "Wass.mean"=barycenter,
    "quantile"=quantile,
    "quantile.grid"=quantile.grid,
    "order"=p),
    class="WARp"))

}
