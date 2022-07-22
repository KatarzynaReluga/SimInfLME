#' Bootstrap NERM
#'
#' Generic function to generate bootstrap samples
#' in a parametric or semiparametric way
#'
#' @inheritParams generate_NERM
#' @inheritParams generate_NERM_sp
#' @inheritParams construct_intervals
#' @param boot_object Object defining bootstrap type
#' @param beta_est Estimated beta parameter
#' @param boot_seed Seed to generate bootstrap samples
#' @param n_boot Number of bootstrap samples
#' @param formula_y Formula to fit NERM
#' @param data_sample Data frame
#' @param ... Additional parameters
#'
#' @return List with following parameters
#' \item{boot_NERMS}{Bootstrap NERM samples}
#' \item{fit_boot_NERMS}{Fitted NERM models to bootstrap samples}
#' \item{double_boot_NERMS}{Double bootstrap NERM samples}
#' \item{fit_double_boot_NERMS}{Fitted NERM models to double-bootstrap samples}
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
#' formula_y <-  y ~ -1 + X0 + X1 + (1| id_cluster)
#'
#' data_sample <- generate_NERM(generate_u = list(type = "chisquared",
#'                                                scaling_factor = 1,
#'                                                dg = 6),
#'                              generate_e = list(type = "chisquared",
#'                                                scaling_factor = 1,
#'                                                dg = 6),
#'                              beta = beta,
#'                              X = X,
#'                              id_cluster = id_cluster,
#'                              start_seed = 1,
#'                              no_sim = 1,
#'                              cluster_means = cluster_means)
#'
#'
#'
#' fitted_NERM <- fit_NERM(formula_y, data_sample,
#'                         id_cluster, cluster_means)
#'
#' type_method = c("parametric")
#' class(type_method) <- type_method
#'
#' bootstrap_NERM <- bootstrap_NERM(type_method,
#'                                  var_u_est = fitted_NERM$var_u,
#'                                  var_e_est = fitted_NERM$var_e,
#'                                  beta_est = fitted_NERM$beta_hat,
#'                                  boot_seed = 2,
#'                                  n_boot = 100, formula_y,
#'                                  data_sample = data_sample,
#'                                  id_cluster = id_cluster,
#'                                  cluster_means = cluster_means,
#'                                  type_var_estimator = "var")
#'
#' type_method = c("semiparametric")
#' class(type_method) <- type_method
#'
#' bootstrap_NERM <- bootstrap_NERM(type_method,
#'                                  var_u_est = fitted_NERM$var_u,
#'                                  var_e_est = fitted_NERM$var_e,
#'                                  empirical_u = fitted_NERM$u_hat,
#'                                  empirical_e = fitted_NERM$e_hat,
#'                                  beta_est = fitted_NERM$beta_hat,
#'                                  boot_seed = 2,
#'                                  n_boot = 100, formula_y,
#'                                  data_sample = data_sample,
#'                                  id_cluster = id_cluster,
#'                                  cluster_means = cluster_means,
#'                                  type_var_estimator = "var")
#'
#'


bootstrap_NERM <- function(...)
  UseMethod("bootstrap_NERM")



#'
#' @describeIn bootstrap_NERM Generate bootstrap samples in a parametric way
#' @export
#'

bootstrap_NERM.parametric <- function(boot_object,
                                      var_u_est,
                                      var_e_est,
                                      beta_est,
                                      boot_seed = 10,
                                      n_boot = 100,
                                      formula_y,
                                      data_sample,
                                      id_cluster,
                                      type_var_estimator,
                                      cluster_means,
                                      ...) {
  X <- format_data_matrix(data = data_sample,
                          name_col = "X",
                          select_col = "x|X")

  generate_samples <-
    generate_NERM(
      generate_u = list(type = "normal",
                        scaling_factor = sqrt(var_u_est)),
      generate_e = list(type = "normal",
                        scaling_factor = sqrt(var_e_est)),
      beta = beta_est,
      X = X,
      id_cluster = id_cluster,
      start_seed = boot_seed,
      no_sim = n_boot,
      cluster_means = cluster_means,
      return_u = TRUE
    )


  NERMS_wrapper <- function(data_sample) {
    boot_estimates <- fit_NERM(formula_y, data_sample,
                               id_cluster, cluster_means = cluster_means)
    return(boot_estimates)
  }


  fit_boot_NERMS <- lapply(generate_samples, NERMS_wrapper)

  output <- list(boot_NERMS = generate_samples,
                 fit_boot_NERMS = fit_boot_NERMS)

  # Estimate NERM for each sample
  if (type_var_estimator == "mse_bc") {
    double_boot_seeds = boot_seed * 1:length(fit_boot_NERMS)

    for (i in 1:length(fit_boot_NERMS)) {
      fit_boot_NERMS[[i]]$double_boot_seed <- double_boot_seeds[i]
    }

    db_wrapper <- function(estimated_NERM_db) {
      generate_samples <- generate_NERM(
        generate_u = list(
          type = "normal",
          scaling_factor = sqrt(estimated_NERM_db$var_u)
        ),
        generate_e = list(
          type = "normal",
          scaling_factor = sqrt(estimated_NERM_db$var_e)
        ),
        beta = estimated_NERM_db$beta_hat,
        X = X,
        cluster_means = cluster_means,
        id_cluster = id_cluster,
        start_seed = estimated_NERM_db$double_boot_seed,
        no_sim = 1,
        return_u = TRUE
      )

    }
    generate_db_samples <- lapply(fit_boot_NERMS, db_wrapper)

    fit_double_boot_NERMS <-
      lapply(generate_db_samples, NERMS_wrapper)

    output <- list(
      boot_NERMS = generate_samples,
      fit_boot_NERMS = fit_boot_NERMS,
      double_boot_NERMS = generate_db_samples,
      fit_double_boot_NERMS = fit_double_boot_NERMS
    )

  }
  return(output)
}




#'
#' @describeIn bootstrap_NERM Generate bootstrap samples in a semiparametric way
#' @export
#'

bootstrap_NERM.semiparametric <- function(boot_object,
                                          var_u_est,
                                          var_e_est,
                                          empirical_u = NULL,
                                          empirical_e = NULL,
                                          beta_est,
                                          boot_seed = 10,
                                          n_boot = 100,
                                          formula_y,
                                          data_sample,
                                          id_cluster,
                                          type_var_estimator,
                                          cluster_means,
                                          ...) {
  X <- format_data_matrix(data = data_sample,
                          name_col = "X",
                          select_col = "x|X")

  # Fix beta_est and beta_hat, decide on one of them.
  generate_samples <- generate_NERM_sp(
    empirical_u = empirical_u,
    var_u_est = var_u_est,
    empirical_e = empirical_e,
    var_e_est = var_e_est,
    beta = beta_est,
    X = X,
    cluster_means = cluster_means,
    id_cluster = id_cluster,
    start_seed = boot_seed,
    no_sim = n_boot,
    return_u = TRUE
  )

  NERMS_wrapper <- function(data_sample) {
    boot_estimates <- fit_NERM(formula_y, data_sample,
                               id_cluster, cluster_means)
    return(boot_estimates)
  }


  fit_boot_NERMS <- lapply(generate_samples, NERMS_wrapper)

  output <- list(boot_NERMS = generate_samples,
                 fit_boot_NERMS = fit_boot_NERMS)

  # Estimate NERM for each sample
  if (type_var_estimator == "mse_bc") {
    double_boot_seeds = boot_seed * 1:length(fit_boot_NERMS)

    for (i in 1:length(fit_boot_NERMS)) {
      fit_boot_NERMS[[i]]$double_boot_seed <- double_boot_seeds[i]
    }

    db_wrapper <- function(estimated_NERM_db) {
      generate_samples <-
        generate_NERM_sp(
          empirical_u = estimated_NERM_db$u_hat,
          var_u_est = estimated_NERM_db$var_u,
          empirical_e = estimated_NERM_db$e_hat,
          var_e_est = estimated_NERM_db$var_e,
          beta = estimated_NERM_db$beta_hat,
          X = X,
          cluster_means = cluster_means,
          id_cluster = id_cluster,
          start_seed = estimated_NERM_db$double_boot_seed,
          no_sim = 1,
          return_u = TRUE
        )

    }
    generate_db_samples <- lapply(fit_boot_NERMS, db_wrapper)

    fit_double_boot_NERMS <-
      lapply(generate_db_samples, NERMS_wrapper)

    output <- list(
      boot_NERMS = generate_samples,
      fit_boot_NERMS = fit_boot_NERMS,
      double_boot_NERMS = generate_db_samples,
      fit_double_boot_NERMS = fit_double_boot_NERMS
    )

  }
  return(output)



}
