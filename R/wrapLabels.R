#' @title Wrap Labels
#'
#' @description Takes in a string and adds line breaks at specified breakpoints
#'
#' @param x string of interest
#' @param width character length for breakpoint
#'
#' @return input string with new line at every width characters
#'
#' @export
wrapLabels <- function(x, width = 60) {
  sapply(x, function(str) {
    if (nchar(str) > width) {
      paste(strwrap(str, width = width), collapse = "\n")
    } else {
      str
    }
  }, USE.NAMES = FALSE)
}
