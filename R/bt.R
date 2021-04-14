#' Perform a binomial test
#'
#' @param a A number
#' @param b A number
#' @param p A number
#' @return Returns the binomial test pvalue for an event that occurred \code{a} times out of a total of \code{b} with a null probability assumption of \code{p}
#' @examples
#' bt(2,10,.5)
#' @export
#'

bt <- function(a, b, p = 0.5) {binom.test(a, b, p, alternative=
                                            c("greater"), conf.level = 0.95)$p.value}
