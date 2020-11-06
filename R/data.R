#' Random Rain
#'
#' Randomly simulated synthetic rain fields with Matern covariances
#'
#' @docType data
#'
#' @usage data(rrain)
#'
#' @format A \code{4x10x128x128} matrix
#'
#' @keywords datasets
#'
#' @details These fields were used in Buschow et al. (2019) <doi:10.5194/gmd-12-3401-2019>. The first array corresponds to the four model configurations from that paper (different roughness nu and scale sc), the second dimension contains ten realizations for each model.
#' @source simulated using the 'RandomFields' package, code available at <10.5281/zenodo.3257511>
#'
#' @examples
#' data(rrain)
"rrain"  
