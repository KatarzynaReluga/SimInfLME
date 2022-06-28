#' Generate centered and scaled random elements
#'
#' @param n Number of observations
#' @param dg Degrees of freedom
#' @param scaling_factor Scaling factor
#' @param type Type of random element, \code{chisquare}: chi-squared distributed,
#' \code{student_t}: Student's t-distributed
#' @param start_seed Seed to reprodcue simulations
#'
#' @importFrom stats rt rchisq
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
                        dg,
                        scaling_factor,
                        type = c("chisquare", "student_t"),
                        start_seed = 1) {

  set.seed(start_seed)

  if (type == "chisquare") {
    random_elements <-
      ((rchisq(n, dg) - dg) / sqrt(2 * dg)) * scaling_factor
  } else {
    random_elements <-
      (rt(n, dg) / sqrt(dg / (dg - 2))) * scaling_factor
  }

  return(random_elements)

}
