#' Example data for the \pkg{rsurface} package
#'
#' This example uses experimental data published in Czitrom and Spagon (1997),
#'     \emph{Statistical Case Studies for Industrial Process Improvement} that
#'     describes a semiconductor wafer processing experiment. A goal of this experiment
#'     was to fit response surface models to the deposition layer stress
#'     as a function of two particular controllable factors of the chemical vapor deposition
#'     (CVD) reactor process. These factors were pressure (measured in torr)
#'     and the ratio of the gaseous reactants hydrogen gas and tungsten(VI) fluoride.
#'
#' @format A data frame with three columns and ten rows of values
#' \describe{
#'   \item{Factor1}{Pressure measured in torr}
#'   \item{Factor2}{The ratio of gaseous reactants.
#'       The smallest and greatest values for the ratios of hydrogen gas to tungsten(VI) fluoride were chosen to be 2 and 10.}
#'   \item{Response}{Deposition layer stress}
#'   }
#'
#' @references
#' Czitrom, V., and Spagon, P. D., (1997), Statistical Case Studies for Industrial Process Improvement, Philadelphia, PA, ASA-SIAM Series on Statistics and Applied Probability.
"ExampleData"
