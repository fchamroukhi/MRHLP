FittedMRHLP <- setRefClass(
  "FittedMRHLP",
  fields = list(
    modelMRHLP = "ModelMRHLP",
    paramMRHLP = "ParamMRHLP",
    statMRHLP = "StatMRHLP"
  ),
  methods = list(
    plot = function() {

      yaxislim <- c(min(modelMRHLP$Y) - 2 * mean(sqrt(apply(modelMRHLP$Y, 2, var))), max(modelMRHLP$Y) + 2 * mean(sqrt(apply(modelMRHLP$Y, 2, var))))

      # Data, regressors, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      matplot(modelMRHLP$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y", col = gray.colors(modelMRHLP$m), lty = 1)
      title(main = "Time series, MRHLP regimes, and process probabilites")
      colorsvec <- rainbow(modelMRHLP$K)
      for (k in 1:modelMRHLP$K) {
        index <- (statMRHLP$klas == k)
        for (d in 1:modelMRHLP$m) {
          polynomials <- statMRHLP$polynomials[index, d, k]
          lines(statMRHLP$polynomials[, d, k], col = colorsvec[k], lty = "dotted", lwd = 1)
          lines(seq(1:modelMRHLP$n)[index], polynomials, col = colorsvec[k], lwd = 1.5)
        }
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statMRHLP$piik[, 1], type = "l", xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)), col = colorsvec[1], lwd = 1.5)
      if (modelMRHLP$K > 1) {
        for (k in 2:modelMRHLP$K) {
          lines(statMRHLP$piik[, k], col = colorsvec[k], lwd = 1.5)
        }
      }

      # Data, regression model, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      matplot(modelMRHLP$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y", col = gray.colors(modelMRHLP$m), lty = 1)
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

FittedMRHLP <- function(modelMRHLP, paramMRHLP, statMRHLP) {
  new("FittedMRHLP", modelMRHLP = modelMRHLP, paramMRHLP = paramMRHLP, statMRHLP = statMRHLP)
}
