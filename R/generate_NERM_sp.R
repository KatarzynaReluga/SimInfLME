#' Generate NERM samples
#'
#' Function to generate data samples with outcomes from nested error regression model (NERM)
#' in a semiparametric way
#'
#' @inheritParams generate_NERM
#' @param empirical_u Estimated random effects
#' @param var_u_est Estimated variance of random effects
#' @param empirical_e Estimated errors
#' @param var_e_est Estimated errors of random effects
#
#' @return
#' \item{samples}{Data frame or list with data frames with samples}
#'
#' @export
#'
#' @examples
#'
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
#'
#' cluster_means <- as.matrix(aggregate(X,
#' list(id_cluster),
#' FUN = mean)[, -1])
#'
#' #Formula to construct intervals
#' formula_y <-  y ~ 1 + X1 + (1| id_cluster)
#'
#' data_sample <- generate_NERM(generate_u = list(type = "chisquared",
#'                                                     scaling_factor = 1,
#'                                                     dg = 6),
#'                                   generate_e = list(type = "chisquared",
#'                                                     scaling_factor = 1,
#'                                                     dg = 6),
#'                                   beta = beta,
#'                                   X = X,
#'                                   cluster_means,
#'                                   id_cluster = id_cluster,
#'                                   start_seed = 1,
#'                                   no_sim = 1)
#'
#'
#' fitted_NERM <- fit_NERM(formula_y, data_sample,
#'                         id_cluster, cluster_means)
#'
#'
#' generate_NERM_sp <- function(empirical_u = fitted_NERM$u_hat,
#'                              var_u_est = fitted_NERM$var_u,
#'                              empirical_e = fitted_NERM$e_hat,
#'                              var_e_est = fitted_NERM$var_e,
#'                              beta = fitted_NERM$beta_hat, X = X,
#'                              cluster_means,
#'                              id_cluster = id_cluster,
#'                              start_seed = 1, no_sim = 10)
#'
#'


generate_NERM_sp <- function(empirical_u, var_u_est,
                             empirical_e, var_e_est,
                             beta,
                             X,
                             cluster_means,
                             id_cluster,
                             start_seed = 1,
                             no_sim = 10,
                             return_u = FALSE) {


  generate_sample_wrapper <- function(start_seed, return_mu = FALSE) {

    set.seed(start_seed)

    u_scale = empirical_u * sqrt(var_u_est) / sqrt(mean(empirical_u ^ 2))
    u_scale_center = u_scale - mean(u_scale)
    u_star = sample(u_scale_center, length(unique(id_cluster)),
                         replace = TRUE)
    re_u_repeat <- rep(u_star, as.numeric(table(id_cluster)))

    e_scale = empirical_e * sqrt(var_e_est) / sqrt(mean(empirical_e ^ 2))
    e_scale_center = e_scale - mean(e_scale)
    re_e = sample(e_scale_center,  length(id_cluster), replace = TRUE)


    y  = X %*% beta + re_u_repeat + re_e
    mu = crossprod(t(cluster_means), beta) + u_star
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
    samples <- lapply(start_seed, generate_sample_wrapper)
    samples <- data.frame(samples)
  } else {
    sim_seed <- list()
    for (j in 1:no_sim)
      sim_seed[[j]] <- start_seed * j

    samples <- lapply(sim_seed, generate_sample_wrapper)

  }
}


