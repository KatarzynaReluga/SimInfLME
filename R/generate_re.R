#' Generate random elements
#'
#' Function to generate centered and scaled random effects or errors from chi-squared,
#' Student's t or normal distributions
#'
#' @param n Number of observations
#' @param dg Degrees of freedom
#' @param scaling_factor Scaling factor
#' @param type Distribution of random element, \code{chisquare}: chi-squared distribution,
#' \code{student_t}: Student's t-distribution, \code{normal}: normal distribution
#' @param start_seed Seed to reproduce simulations
#'
#' @importFrom stats rnorm rt rchisq
#'
#' @return \item{random_elements}{Vector of random elements}
#'
#' @export
#'
#' @examples
#'
#' # Generate scaled and centered Student's t-distributed random elements
#'
#' t_student_re <- generate_re(n = 100, dg = 6, scaling_factor = 1,
#' type = "student_t")
#'
#'


generate_re <- function(n,
                        dg = 5,
                        scaling_factor,
                        type = c("chisquared", "student_t", "normal"),
                        start_seed = 1) {

  set.seed(start_seed)

  if (type == "chisquare") {
    random_elements <-
      ((rchisq(n, dg) - dg) / sqrt(2 * dg)) * scaling_factor
  } else if (type == "student_t") {
    random_elements <-
      (rt(n, dg) / sqrt(dg / (dg - 2))) * scaling_factor
  } else {
    random_elements <-
      rnorm(n, 0, scaling_factor)
  }

  return(random_elements)

}
