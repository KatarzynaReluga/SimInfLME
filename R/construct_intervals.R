#' Intervals for linear mixed effects
#'
#' Function construct_intLME is a wrapper function to
#' construct simultaneous and individual intervals for
#' linear mixed effects. For the moment, only NERM is implemented,
#' but the function can be easily extended to other model types
#'
#' @inheritParams generate_NERM
#' @param type_method Method to construct intervals
#' @param type_var_estimator Type of variance estimator
#' @param model_type Type of model, for now only NERM is supported
#' @param formula_y Formula for fitting NERM
#' @param data_sample Data frame
#' @param alpha Level alpha to construct intervals
#' @param n_boot Number of bootstrap samples, default n_boot = 1000
#' @param boot_seed Vector with cluster labels
#' @param ... Additional parameters
#'
#' @return List with following parameters:
#' \item{int_up}{Upper boundary of individual interval}
#' \item{int_down}{Lower boundary of individual interval}
#' \item{length_int}{Length of individual interval}
#' \item{int_sim_up}{Upper boundary of simultaneous interval}
#' \item{int_sim_down}{Lower boundary of simultaneous interval}
#' \item{length_sim}{Length of simulaneous interval}
#'
#'
#'
#' @export
#'
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
#' intervals_3t <- construct_intervals(type_method = c("parametric"),
#'                                     type_var_estimator = c("mse_3t"),
#'                                     model_type = c("NERM"),
#'                                     formula_y, data_sample,
#'                                     cluster_means = cluster_means,
#'                                     id_cluster = id_cluster,
#'                                     alpha = 0.05,
#'                                     n_boot = 100,
#'                                     boot_seed = 1)
#'
#'
#' intervals_bc <- construct_intervals(type_method = c("parametric"),
#'                                     type_var_estimator = c("mse_bc"),
#'                                     model_type = c("NERM"),
#'                                     formula_y, data_sample,
#'                                     cluster_means = cluster_means,
#'                                     id_cluster = id_cluster,
#'                                     alpha = 0.05,
#'                                     n_boot = 100,
#'                                     boot_seed = 1)
#'
#' intervals_b <- construct_intervals(type_method = c("parametric"),
#'                                    type_var_estimator = c("mse_b"),
#'                                    model_type = c("NERM"),
#'                                    formula_y, data_sample,
#'                                    cluster_means = cluster_means,
#'                                    id_cluster = id_cluster,
#'                                    alpha = 0.05,
#'                                    n_boot = 100,
#'                                    boot_seed = 1)
#'
#' intervals_spa <- construct_intervals(type_method = c("parametric"),
#'                                      type_var_estimator = c("mse_spa"),
#'                                      model_type = c("NERM"),
#'                                      formula_y, data_sample,
#'                                      cluster_means = cluster_means,
#'                                      id_cluster = id_cluster,
#'                                      alpha = 0.05,
#'                                      n_boot = 100,
#'                                      boot_seed = 1)
#'

construct_intervals <- function(type_method = c("parametric",
                                                "semiparametric",
                                                "analytical"),
                                type_var_estimator = c("g1", "mse_a",
                                                       "mse_b", "mse_3t",
                                                       "mse_spa", "mse_bc"),
                                model_type = c("NERM"),
                                formula_y,
                                data_sample,
                                cluster_means,
                                id_cluster,
                                alpha = 0.05,
                                n_boot = 1000,
                                boot_seed = 1,
                                ...) {
  # Define the type of model
  model_type <- match.arg(model_type)
  class(model_type) <- model_type

  # Define the type of method
  type_method <- match.arg(type_method)
  class(type_method) <- type_method

  intervals <- construct_intLME(
    model_obj = model_type,
    type_method = type_method,
    type_var_estimator = type_var_estimator,
    formula_y,
    data_sample,
    cluster_means = cluster_means,
    id_cluster = id_cluster,
    alpha = 0.05,
    n_boot = 100,
    boot_seed = 1
  )

  return(intervals)

}

#' Intervals for linear mixed models under NERM
#'
#' Generic function to construct individual and simultaneous intervals
#' for linear mixed models
#'
#' @inheritParams construct_intervals
#' @param model_obj Object of class model
#'
#' @importFrom postcAIC compute_corrected_mse
#' @importFrom stats qnorm quantile
#'
#' @return List with following parameters:
#' \item{int_up}{Upper boundary of individual interval}
#' \item{int_down}{Lower boundary of individual interval}
#' \item{length_int}{Length of individual interval}
#' \item{int_sim_up}{Upper boundary of simultaneous interval}
#' \item{int_sim_down}{Lower boundary of simultaneous interval}
#' \item{length_sim}{Length of simulaneous interval}
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
#' model_obj = c("NERM")
#' # For now, only NERM is supported, but other model choices can be easily programmed
#' class(model_obj) <- model_obj
#'
#' construct_LMEintervals <- construct_intLME(model_obj,
#'                                            type_method = c("parametric"),
#'                                            type_var_estimator = c("mse_bc"),
#'                                            formula_y, data_sample,
#'                                            cluster_means = cluster_means, id_cluster,
#'                                            alpha = 0.05, n_boot = 100,
#'                                            boot_seed = 1)
#'
#'



construct_intLME <- function(...)
  UseMethod("construct_intLME")



#'
#' @describeIn construct_intLME Method to construct individual and simultanous intervvals for NERM
#' @export
#'



construct_intLME.NERM <- function(model_obj,
                                  type_method = c("parametric",
                                                  "semiparametric",
                                                  "analytical"),
                                  type_var_estimator = c("var", "mse_a",
                                                         "mse_b", "mse_3t",
                                                         "mse_spa", "mse_bc"),
                                  formula_y,
                                  data_sample,
                                  cluster_means,
                                  id_cluster,
                                  alpha,
                                  n_boot,
                                  boot_seed = 1,
                                  ...) {
  # Method to construct intervals
  type_method <- match.arg(type_method)
  class(type_method) <- type_method

  # Variance method
  type_var_estimator <- match.arg(type_var_estimator)
  class(type_var_estimator) <- type_var_estimator

  fitted_NERM <-
    fit_NERM(formula_y, data_sample, id_cluster, cluster_means)
  m = length(fitted_NERM$mu_hat)
  var_u_est  = fitted_NERM$var_u
  var_e_est  = fitted_NERM$var_e

  if (type_method == "analytical") {
    q_ind_alpha = qnorm(1 - alpha / 2)
    int_up <-
      fitted_NERM$mu_hat + sqrt(fitted_NERM$mse) *  q_ind_alpha
    int_down <-
      fitted_NERM$mu_hat - sqrt(fitted_NERM$mse) *  q_ind_alpha

    q_sim_alpha = qnorm(1 - alpha / (2 * m))
    int_sim_up <-
      fitted_NERM$mu_hat + sqrt(fitted_NERM$mse) *  q_sim_alpha
    int_sim_down <-
      fitted_NERM$mu_hat - sqrt(fitted_NERM$mse) *  q_sim_alpha


  } else {
    boot_samples <- bootstrap_NERM(
      type_method,
      var_u_est = var_u_est,
      var_e_est = var_e_est,
      beta_est = fitted_NERM$beta_hat,
      boot_seed = boot_seed,
      n_boot = n_boot,
      formula_y = formula_y,
      data_sample = data_sample,
      id_cluster = id_cluster,
      cluster_means = cluster_means,
      type_var_estimator = type_var_estimator
    )

    # Retrieve mu_boot and mu_hat_boot
    mu_matrix <-
      matrix(
        0,
        ncol = length(boot_samples$fit_boot_NERMS[[1]]$mu_hat),
        nrow = length(boot_samples$fit_boot_NERMS)
      )

    mu_hat_matrix <-
      matrix(
        0,
        ncol = length(boot_samples$fit_boot_NERMS[[1]]$mu_hat),
        nrow = length(boot_samples$fit_boot_NERMS)
      )

    sd_mixed <-
      matrix(
        0,
        ncol = length(boot_samples$fit_boot_NERMS[[1]]$mu_hat),
        nrow = length(boot_samples$fit_boot_NERMS)
      )

    for (i in 1:length(boot_samples$boot_NERMS)) {
      mu_matrix[i, ] <- unique(boot_samples$boot_NERMS[[i]]$mu)
      mu_hat_matrix[i, ] <- boot_samples$fit_boot_NERMS[[i]]$mu_hat
      sd_mixed[i, ] <- sqrt(boot_samples$fit_boot_NERMS[[i]]$g1)
    }

    # Bootstrap differences
    dif <- mu_hat_matrix - mu_matrix

    rmse_mixed <- boot_MSE(type_var_estimator,
                           boot_samples,
                           cluster_means,
                           var_u_est,
                           var_e_est)$rmse_mixed

    dif_M <- matrix(0, ncol = length(rmse_mixed), nrow = n_boot)
    dif_t <- matrix(0, ncol = length(rmse_mixed), nrow = n_boot)


    if (type_var_estimator == "g1") {
      for (i in 1:n_boot) {
        dif_M[i, ] = abs(dif[i, ]) / sd_mixed[i, ]
        dif_t[i, ] = dif[i, ] / sd_mixed[i, ]
      }
    } else {
      for (i in 1:m) {
        dif_M[, i] = abs(dif[, i]) / rmse_mixed[i]
        dif_t[, i] = dif[, i] / rmse_mixed[i]
      }

    }

    q_sim_alpha = quantile(apply(dif_M, 1L, max) ,
                           prob = 1 - alpha,
                           type = 8)
    q_ind_alpha = apply(dif_t, 2L, quantile, probs = c(alpha / 2, 1 - alpha / 2))


    if (type_var_estimator == "g1") {
      int_up <-
        fitted_NERM$mu_hat + sqrt(fitted_NERM$g1) *  q_ind_alpha[2]
      int_down <-
        fitted_NERM$mu_hat + sqrt(fitted_NERM$g1) *  q_ind_alpha[1]
      length_int <- int_up - int_down

      int_sim_up <-
        fitted_NERM$mu_hat + sqrt(fitted_NERM$g1) * q_sim_alpha
      int_sim_down <-
        fitted_NERM$mu_hat - sqrt(fitted_NERM$g1) * q_sim_alpha
      length_sim <- int_sim_up - int_sim_down

    } else {
      int_up <- fitted_NERM$mu_hat + rmse_mixed *  q_ind_alpha[2]
      int_down <- fitted_NERM$mu_hat + rmse_mixed *  q_ind_alpha[1]
      length_int <- int_up - int_down

      int_sim_up <- fitted_NERM$mu_hat + rmse_mixed * q_sim_alpha
      int_sim_down <- fitted_NERM$mu_hat - rmse_mixed * q_sim_alpha
      length_sim <- int_sim_up - int_sim_down
    }

  }

  output <-
    list(
      int_up = int_up,
      int_down = int_down,
      length_int = length_int,
      int_sim_up = int_sim_up,
      int_sim_down = int_sim_down,
      length_sim = length_sim
    )
  return(output)

}
