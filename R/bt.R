#' @export


bt <- function(a, b, p = 0.5) {binom.test(a, b, p, alternative=
                                            c("greater"), conf.level = 0.95)$p.value}
