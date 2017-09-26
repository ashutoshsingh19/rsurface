#' @title Response surface regression analysis of CCD data
#'
#' @description For the simultaneous production of response surface analysis output by
#'    \pkg{rsm} used in combination with \pkg{graphics} in the second-order
#'    polynomial approach. The output includes regression model fitting and plot of the
#'    fitted response surface.
#'
#' @param x the matrix of experimental data that contains columns with the uncoded levels
#'    for each experimental factor and the observed values for the response variable in the
#'    rightmost column.
#'
#' @return
#'
#' The user will be prompted to enter \dQuote{1} for a 3-D plot of the response surface,
#' or \dQuote{2} to plot the contour of the predicted variance of the response
#'
#' \dQuote{Data.For.Analysis}, includes the data set and the coding coefficients
#'     for the transformation of the independent factors
#'
#' \dQuote{Response.Surface.Summary}, includes the response surface for variable, hypothesis
#'    tests for linear, quadratic, and crossproduct terms, lack of fit test, parameter
#'    estimates, the factor ANOVA table, canonical analysis, and eigenvectors
#'
#' @references
#' Mead, R., Gilmour, S. G., and Mead, A. 2012. Statistical Principles for the Design of
#'    Experiments: Applications to Real Experiments. Cambridge University Press, Cambridge.
#'
#' Panneton, B., Philion, H., Dutilleul, P., Theriault, R., and Khelifi, M. 1999. Full
#'    factorial design versus central composite design: Statistical comparison and
#'    experimental implications for spray droplet deposition. Transactions of the American
#'    Society of Agricultural Engineers 42:877-883.
#'
#' @examples
#' if(interactive()){
#'  ccd_analysis(ExampleData)
#'  }
#'
#' @export

ccd_analysis <- function(x) {
  Coefficients <- rsm::coded.data(x)
  colnames(x)[ncol(x)] <- "Response"
  for (i in 1:(ncol(x) - 1)) {
    names(x[, i]) <- paste("Factor", i, sep = "")
  }
  rsmodel <- c(NA)
  for (i in 1:(ncol(x) - 1)) {
    rsmodel <- c(rsmodel, paste("Factor", i, sep = ""))
  }
  rsmodel <- rsmodel[-1]
  rsmodel <- paste(rsmodel, collapse = ',')
  rsmodel1 <-
    stats::as.formula(paste("Response ~ SO(", rsmodel, ")", sep = ""))
  ResponseSurface <- rsm::rsm(rsmodel1, data = x)
  rsmodel2 <- gsub(",", "+", rsmodel)
  rsmodel2 <- stats::as.formula(paste("~", rsmodel2, sep = ""))
  NormalizedFactors <- scale(x[-ncol(x)])
  plotprompt <-
    paste(
      "For a 3-D plot of the response surface press 1, to plot the contour of the predicted variance of the response, press any other key: "
    )
  plotval <- readline(plotprompt)
  if (plotval == 1) {
    graphics::persp(
      ResponseSurface,
      rsmodel2,
      at = rsm::canonical(ResponseSurface),
      col = grDevices::rainbow(50),
      contours = "colors"
    )
  } else {
    X <- cbind(matrix(1, nrow(x), 1), as.matrix(x[1:nrow(x), 1:2]))
    X <- cbind(X, matrix(NA, nrow(X), 3))
    for (i in 1:nrow(X)) {
      X[i, 4] <- X[i, 2]^2
      X[i, 5] <- X[i, 3]^2
      X[i, 6] <- X[i, 2] * X[i, 3]
    }
    BIGX <- matrix(NA, 1, 6)
    for (i in seq(min(X[1:nrow(X), 2]), max(X[1:nrow(X), 2]), 0.1)) {
      newrow <- matrix(NA, 1, 6)
      for (o in seq(min(X[1:nrow(X), 3]), max(X[1:nrow(X), 3]), 0.1)) {
        newrow[1, 1] <- 1
        newrow[1, 2] <- i
        newrow[1, 3] <- o
        newrow[1, 4] <- i * i
        newrow[1, 5] <- o * o
        newrow[1, 6] <- i * o
        BIGX <- rbind(BIGX, newrow)
      }
    }
    BIGX <- BIGX[2:nrow(BIGX), 1:6]
    Factor1 <- BIGX[1:nrow(BIGX), 2]
    Factor2 <- BIGX[1:nrow(BIGX), 3]
    A <- solve(crossprod(X))
    Predicted_Variance <- matrix(NA, nrow(BIGX), 1)
    for (i in 1:nrow(BIGX)) {
      b <- BIGX[i, 1:6]
      Predicted_Variance[i, 1] <- t(b) %*% A %*% b
    }
      ForPlot <- as.data.frame(cbind(Factor1, Factor2, Predicted_Variance))
      colnames(ForPlot) <- c("Factor1", "Factor2", "Predicted_Variance")
      p <- plotly::plot_ly(ForPlot,
                           x = ~Factor1,
                           y = ~Factor2,
                           z = ~Predicted_Variance,
                           type = "contour",
                           colorscale = "Portland",
                           contours = list(
                             showlabels = TRUE),
                           line = list(smoothing = 0)
      )
      print(p)
  }
  list(
    Data.For.Analysis = Coefficients,
    Response.Surface.Summary = summary(ResponseSurface)
  )
}
