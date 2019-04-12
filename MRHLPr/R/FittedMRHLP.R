FittedMRHLP <- setRefClass(
  "FittedMRHLP",
  fields = list(
    modelMRHLP = "ModelMRHLP",
    paramMRHLP = "ParamMRHLP",
    statMRHLP = "StatMRHLP"
  ),
  methods = list(
    plot = function() {
      par(mfrow = c(2, 1))
      plot.default(modelMRHLP$Y, type = "l", ylab = "y", xlab = "")
      title(main = "Time series, MRHLP regimes and process probabilities")
      colors = rainbow(modelMRHLP$K)
      for (k in 1:modelMRHLP$K) {
        index <- (StatMRHLP$klas == k)
        polynomials <- StatMRHLP$polynomials[(StatMRHLP$klas == k), k]
        lines(StatMRHLP$polynomials[, k], lty = "dotted", lwd = 2, col = colors[k])
        lines(seq(1:modelMRHLP$m)[index], polynomials, lwd = 2, col = colors[k])
      }

      plot.default(StatMRHLP$piik[, 1], type = "l", lwd = 2, col = colors[1], xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)))
      if (K > 1) {
        for (k in 2:modelMRHLP$K) {
          lines(StatMRHLP$piik[, k], type = "l", lwd = 2, col = colors[k])
        }
      }

      par(mfrow = c(2, 1))
      plot.default(modelMRHLP$Y, type = "l", ylab = "y", xlab = "")
      title(main = "Time series, estimated RHLP model, and segmentation")

      tk = which(diff(StatMRHLP$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], lty = "dotted", lwd = 2, col = "red")
      }
      plot.default(StatMRHLP$klas, type = "l", lwd = 2, col = "red", xlab = "", ylab = "Estimated class labels")
    }
  )
)

FittedMRHLP <- function(modelMRHLP, paramMRHLP, statMRHLP) {
  new("FittedMRHLP", modelMRHLP = modelMRHLP, paramMRHLP = paramMRHLP, statMRHLP = statMRHLP)
}
