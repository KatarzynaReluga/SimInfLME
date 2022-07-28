#' Generate NERM samples
#'
#' Function to generate data samples with outcomes from nested error regression model (NERM)
#' in a parametric way
#'
#' @param generate_u List to generate random effects, that is
#' \itemize{
#'  \item scaling_factor -- scaling factor
#'  \item type_dist -- distribution of random element, \code{chisquare}: chi-squared distribution,
#'  \code{student_t}: Student's t-distribution, \code{normal}: normal distribution
#'  \item dg -- degrees of freedom.
#' }
#' @param generate_e List to generate random effects, that is
#' \itemize{
#'  \item scaling_factor -- scaling factor
#'  \item type_dist -- distribution of random element, \code{chisquare}: chi-squared distribution,
#'  \code{student_t}: Student's t-distribution, \code{normal}: normal distribution
#'  \item dg -- degrees of freedom.
#' }
#' @param beta Vector with fixed parameters
#' @param X Matrix with covariates
#' @param id_cluster Vector with cluster labels
#' @param no_sim Number of samples
#' @param start_seed Seed to reproduce simulations
#' @param cluster_means Cluster-level covariates for fixed parameters
#' @param return_u Return random effects? Default: return_u = FALSE
#'
#' @return
#' \item{samples}{Data frame or list with data frames with samples from NERM}
#'
#' @export
#'
#' @examples
#'
#' # Basic setup -------------------------------------------------------------------
#' set.seed(11992)
#' m = 25
#' n_j = 5
#' n = m * n_j
#' scaling_factor_e = 1
#' scaling_factor_u = sqrt(0.5)
#' x_ij = runif(n, 0, 1)
#' X = cbind(X0 = rep(1, n), X1 = x_ij)
#' t_X  = t(X)
#' id_cluster = rep(1:m, each = n_j)
#' beta = c(1, 1)
#' cluster_means <- as.matrix(aggregate(X,
#' list(id_cluster),
#' FUN = mean)[, -1])
#'
#' # Generate NERM sample ---------------------------------------------------------
#' generate_samples <- generate_NERM(generate_u = list(type_dist = "chisquared",
#'                                                     scaling_factor = 1,
#'                                                     dg = 6),
#'                                   generate_e = list(type_dist = "chisquared",
#'                                                     scaling_factor = 1,
#'                                                     dg = 6),
#'                                   beta = beta,
#'                                   X = X,
#'                                   id_cluster = id_cluster,
#'                                   start_seed = 1,
#'                                   no_sim = 10,
#'                                   cluster_means = cluster_means,
#'                                   return_u = FALSE)
#'


generate_NERM <- function(generate_u = list(type_dist = "chisquared",
                                            scaling_factor = 1,
                                            dg = 6),
                          generate_e = list(type_dist = "chisquared",
                                            scaling_factor = 1,
                                            dg = 6),
                          beta,
                          X,
                          cluster_means,
                          id_cluster,
                          start_seed = 1,
                          no_sim = 10,
                          return_u = FALSE) {
  generate_sample_wrapper <- function(start_seed, return_u) {
    set.seed(start_seed)

    re_u <- generate_re(
      n = length(as.numeric(table(id_cluster))),
      dg = generate_u$dg,
      scaling_factor = generate_u$scaling_factor,
      type_dist = generate_u$type_dist,
      start_seed = start_seed
    )
    re_u_repeat <- rep(re_u, as.numeric(table(id_cluster)))

    re_e <- generate_re(
      n = dim(X)[1],
      dg = generate_e$dg,
      scaling_factor = generate_e$scaling_factor,
      type_dist = generate_e$type_dist,
      start_seed = start_seed
    )

    y  = X %*% beta + re_u_repeat + re_e
    mu = crossprod(t(cluster_means), beta) + re_u
    mu_repeat <- rep(mu, as.numeric(table(id_cluster)))



    if (return_u) {
      data_sample <- data.frame(y, X, id_cluster, mu = mu_repeat,
                                u = re_u_repeat)
    } else {
      data_sample <- data.frame(y, X, id_cluster, mu = mu_repeat)
    }
    return(data_sample)
  }



  if (no_sim == 1) {
    samples <- lapply(start_seed, generate_sample_wrapper, return_u)
    samples <- data.frame(samples)
  } else {
    sim_seed <- list()
    for (j in 1:no_sim)
      sim_seed[[j]] <- start_seed * j

    samples <- lapply(sim_seed, generate_sample_wrapper, return_u)

  }
  return(samples)

}
