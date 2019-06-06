#' @export
ModelMRHLP <- setRefClass(
  "ModelMRHLP",
  fields = list(
    paramMRHLP = "ParamMRHLP",
    statMRHLP = "StatMRHLP"
  ),
  methods = list(
    plot = function() {

      oldpar <- par()[c("mfrow", "mai", "mgp")]
      on.exit(par(oldpar), add = TRUE)

      yaxislim <- c(min(paramMRHLP$fData$Y) - 2 * mean(sqrt(apply(paramMRHLP$fData$Y, 2, var))), max(paramMRHLP$fData$Y) + 2 * mean(sqrt(apply(paramMRHLP$fData$Y, 2, var))))

      # Data, regressors, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      matplot(paramMRHLP$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y", col = gray.colors(paramMRHLP$fData$m), lty = 1)
      title(main = "Time series, MRHLP regimes, and process probabilites")
      colorsvec <- rainbow(paramMRHLP$K)
      for (k in 1:paramMRHLP$K) {
        index <- (statMRHLP$klas == k)
        for (d in 1:paramMRHLP$fData$m) {
          polynomials <- statMRHLP$polynomials[index, d, k]
          lines(statMRHLP$polynomials[, d, k], col = colorsvec[k], lty = "dotted", lwd = 1)
          lines(seq(1:paramMRHLP$fData$n)[index], polynomials, col = colorsvec[k], lwd = 1.5)
        }
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statMRHLP$piik[, 1], type = "l", xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)), col = colorsvec[1], lwd = 1.5)
      if (paramMRHLP$K > 1) {
        for (k in 2:paramMRHLP$K) {
          lines(statMRHLP$piik[, k], col = colorsvec[k], lwd = 1.5)
        }
      }

      # Data, regression model, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      matplot(paramMRHLP$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y", col = gray.colors(paramMRHLP$fData$m), lty = 1)
      title(main = "Time series, estimated MRHLP model, and segmentation")

      # Transition time points
      tk = which(diff(statMRHLP$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], col = "red", lty = "dotted", lwd = 1.5)
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statMRHLP$klas, type = "l", xlab = "", ylab = "Estimated class labels", col = "red", lwd = 2)

      # # Model log-likelihood during EM
      # par(mfrow = c(1, 1))
      # plot.default(unlist(statMRHLP$stored_loglik), type = "l", xlab = "EM iteration number", ylab = "log-lokelihodd", col = "blue")
    }
  )
)

