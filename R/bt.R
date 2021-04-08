bt <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=
                                            c("two.sided"), conf.level = 0.95)$p.value}
