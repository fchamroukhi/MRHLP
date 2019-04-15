FittedMRHLP <- setRefClass(
  "FittedMRHLP",
  fields = list(
    modelMRHLP = "ModelMRHLP",
    paramMRHLP = "ParamMRHLP",
    statMRHLP = "StatMRHLP"
  ),
  methods = list(
    plot = function() {
      yaxislim <- c(min(modelMRHLP$Y) - 2 * mean(sqrt(apply(modelMRHLP$Y, 2, var))), min(modelMRHLP$Y) + 2 * mean(sqrt(apply(modelMRHLP$Y, 2, var))))

      # Data, regressors, and segmentation
      par(mfrow = c(2, 1))
      matplot(modelMRHLP$Y, type = "l", lty = 1, col = gray.colors(modelMRHLP$m), ylab = "y", xlab = "")
      title(main = "Time series, MRHLP regimes, and process probabilites")
      colors = rainbow(modelMRHLP$K)
      for (k in 1:modelMRHLP$K) {
        index <- (statMRHLP$klas == k)
        for (d in 1:modelMRHLP$m) {
          polynomials <- statMRHLP$polynomials[(statMRHLP$klas == k), d, k]
          lines(statMRHLP$polynomials[, d, k], lty = "dotted", lwd = 2, col = colors[k])
          lines(seq(1:modelMRHLP$n)[index], polynomials, lwd = 2, col = colors[k])
        }
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statMRHLP$piik[, 1], type = "l", lwd = 2, col = colors[1], xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)))
      if (modelMRHLP$K > 1) {
        for (k in 2:modelMRHLP$K) {
          lines(statMRHLP$piik[, k], type = "l", lwd = 2, col = colors[k])
        }
      }

      # Data, regression model, and segmentation
      par(mfrow = c(2, 1))
      matplot(modelMRHLP$Y, type = "l", lty = 1, col = gray.colors(modelMRHLP$m), ylab = "y", xlab = "")
      title(main = "Time series, estimated MRHLP model, and segmentation")

      # Transition time points
      tk = which(diff(statMRHLP$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], lty = "dotted", lwd = 2, col = "red")
      }
      # Probablities of the hidden process (segmentation)
      plot.default(statMRHLP$klas, type = "l", lwd = 2, col = "red", xlab = "", ylab = "Estimated class labels")

      # Model log-likelihood during EM
      par(mfrow = c(1, 1))
      plot.default(unlist(statMRHLP$stored_loglik)[5:35], type = "l", col = "blue", xlab = "EM iteration number", ylab = "log-lokelihodd")
    }
  )
)

FittedMRHLP <- function(modelMRHLP, paramMRHLP, statMRHLP) {
  new("FittedMRHLP", modelMRHLP = modelMRHLP, paramMRHLP = paramMRHLP, statMRHLP = statMRHLP)
}
