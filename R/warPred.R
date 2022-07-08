#' Prediction by WAR(p) models
#'
#' @description a method of the WARp class which produces a one-step ahead prediction by WAR(p) models
#'
#' @param  object A WARp object, the output of \code{WARp()}.
#' @param  ... Further arguments passed to or from other methods.
#' @param  dSup Optional, a numeric vector, the grid over which forecasted cdf/pdf is evaluated. Should be supplied/ignored with \code{expSup} together.
#' @param  expSup Optional, a numeric vector, the grid over the Exponential map is applied, \code{dSup} should cover and be denser than \code{expSup}.   Should be supplied/ignored with \code{dSup} together.
#'
#' @return A list of:
#'\item{pred.cdf}{predicted cdf}
#'\item{pred.pdf}{predicted pdf}
#'\item{dSup}{support of the predicted cdf/pdf}
#'
#' @references
#' \cite{Wasserstein Autoregressive Models for Density Time Series, Chao Zhang, Piotr Kokoszka, Alexander Petersen, 2022}
#'
#' @seealso \code{\link{WARp}}
#'
#'
#' @export

predict.WARp <- function(object, dSup, expSup, ...){

  quantile <- object$quantile
  quantile.grid <- object$quantile.grid
  p <- object$order

  if(missing(dSup) & missing(expSup)) {
    dSup <- seq(range(quantile)[1],
                range(quantile)[2],
                length.out=5*length(quantile.grid))

    expSup <- seq(range(quantile)[1],
                  range(quantile)[2],
                  length.out=length(quantile.grid))
  }

  else if ((!missing(dSup) & missing(expSup)) | (missing(dSup) & !missing(expSup))) {
    stop("'dSup' and 'expSup' should be supplied or ignored together!")
  }

  if((range(dSup)[1] > range(expSup)[1]) | (range(dSup)[2] < range(expSup)[2]) | length(dSup) < length(expSup)) {
    stop("'dSup' should cover and be denser than 'expSup'.")
  }

  ar.coef.mat <- matrix(object$coef, p, 1)
  sample.size <- ncol(quantile)

  pred.tgt <- (quantile[, sample.size:(sample.size-p+1)] - object$Wass.mean) %*% ar.coef.mat + object$Wass.mean

  pred.cdf <- Exp_Map_Barycenter_Method(expSup, pred.tgt, object$quantile.grid)

  pred.cdf.on.dSup <- approx(expSup, pred.cdf, xout=dSup, rule=2, ties=mean)

  pred.pdf <- diff(pred.cdf.on.dSup$y)/diff(dSup)

  # Remove the initial point due to the effect of numerical differentiation
  result.list <- list("pred.cdf"=pred.cdf.on.dSup$y[-1],
                      "pred.pdf"=pred.pdf,
                      "dSup"=dSup[-1])
}
