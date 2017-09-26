#' @title Defines a CCD for k quantitative factors
#'
#' @description Defines a rotatable central composite design (CCD) for experimentation with
#'    k quantitative factors and j replicates of the central point
#'    (which corresponds to the \sQuote{average treatment combination}).
#'
#' @param j the number of replicates of the central point
#'
#' @param k the number of quantitative factors used or studied in the experimentation.
#'
#' @return Three output tables containing the level of replication (number of replicates) and the experimental uncoded
#'    values of the quantitative factors to be used for experimentation and one plot showing
#'    the corresponding variance of the predicted response.
#'
#' \dQuote{Factorial.Points}, the first table which contains the treatment combinations for a 2^\emph{k} factorial design
#'    (which, in coded form, corresponds to the vertices of a square, a cube, or a hyper-cube when \emph{k} = 2, 3 and
#'    more).
#'
#' \dQuote{Axial.Points}, the second table which contains 2k axial or \dQuote{star} points.
#'
#' \dQuote{Central.Point}, the third table which contains the number of replicates for the central point, coded (0, 0, 0).
#'
#' @examples
#' #Enter the function as shown below. The user will be prompted to input smallest
#' #and greatest values for each factor that will be used for experimentation.
#'
#' if(interactive()){
#'  design_ccd(5, 3)
#'  }
#'
#' @references
#' Mead, R., Gilmour, S. G., and Mead, A. 2012. Statistical Principles for the Design of Experiments:
#'    Applications to Real Experiments. Cambridge University Press, Cambridge.
#'
#' Panneton, B., Philion, H., Dutilleul, P., Theriault, R., and Khelifi, M. 1999. Full factorial design
#'    versus central composite design: Statistical comparison and experimental implications for spray
#'    droplet deposition. Transactions of the American Society of Agricultural Engineers 42:877-883.
#'
#' @export

design_ccd <- function(j, k) {
  MinStarPoints <- matrix(NA, nrow = 1, ncol = k)
  MaxStarPoints <- matrix(NA, nrow = 1, ncol = k)
  CentrePoints <- matrix(NA, nrow = 1, ncol = k)
  AlphaPoints <- matrix(NA, nrow = 1, ncol = k)
  MinusLambdaPoints <- matrix(NA, nrow = 1, ncol = k)
  PlusLambdaPoints <- matrix(NA, nrow = 1, ncol = k)
  for (i in 1:k) {
    prompt1 <- paste("What is the smallest value of factor ", i, "? ", sep = "")
    prompt2 <- paste("What is the greatest value of factor ", i, "? ", sep = "")
    minval <- readline(prompt1)
    maxval <- readline(prompt2)
    MinStarPoints[1, i] <- as.numeric(minval)
    MaxStarPoints[1, i] <- as.numeric(maxval)
    CentrePoints[1, i] <- (MinStarPoints[1, i] + MaxStarPoints[1, i]) / 2
    AlphaPoints[1, i] <- abs(MaxStarPoints[1, i] - MinStarPoints[1, i]) / (2 * 2^(k / 4))
    PlusLambdaPoints[1, i] <- CentrePoints[1, i] + AlphaPoints[1, i]
    MinusLambdaPoints[1, i] <- CentrePoints[1, i] - AlphaPoints[1, i]
  }
  FactorialDesign <- MinusLambdaPoints
    for (i in 1:k) {
    FactorialDesign2 <- FactorialDesign
    FactorialDesign <- rbind(FactorialDesign, FactorialDesign2)
    FactorialDesign[(2^(i - 1) + 1): 2^i, i] <- PlusLambdaPoints[1, i]
    }
  headings <- c("Replication")
  for (i in 1: k) {
   headings <- c(headings, paste("Factor", i, sep = ""))
  }
  ones <- matrix(1, nrow = 2^k, ncol = 1)
  ones <- as.data.frame(ones)
  FactorialDesign <- cbind(ones, FactorialDesign)
  colnames(FactorialDesign) <- headings
  AxialDesign <- matrix(NA, nrow = (2 * k), ncol = k)
  for (i in 1:k) {
    AxialDesign1 <- CentrePoints
    AxialDesign1[1, i] <- MinStarPoints[1, i]
    AxialDesign2 <- CentrePoints
    AxialDesign2[1, i] <- MaxStarPoints[1, i]
    AxialDesign[(2 * i - 1), 1: k] <- AxialDesign1
    AxialDesign[(2 * i), 1: k] <- AxialDesign2
  }
  AxialOnes <- matrix(1, ncol = 1, nrow = 2*k)
  AxialOnes <- as.data.frame(AxialOnes)
  AxialDesign <- cbind(AxialOnes, AxialDesign)
  colnames(AxialDesign) <- headings
  CentreDesign <- as.data.frame(cbind(j, CentrePoints))
  colnames(CentreDesign) <- headings
  CentreRows <- CentrePoints
  for (i in 1:j) {
    row2 <- CentrePoints
    CentreRows <- rbind(CentreRows, row2)
  }
  CentreRows <- CentreRows[-nrow(CentreRows),]
  CentreRows<- as.data.frame(cbind(1, CentreRows))
  colnames(CentreRows) <- headings

  X <- as.matrix(rbind(FactorialDesign[1:nrow(FactorialDesign), 1:3], AxialDesign[1:nrow(AxialDesign), 1:3], CentreRows[1:nrow(CentreRows), 1:3]))
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
if (k == '2') {
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
list(Factorial.Points = FactorialDesign,
       Axial.Points = AxialDesign,
       Central.Point = CentreDesign)
}
