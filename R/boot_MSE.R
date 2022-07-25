#' Bootstrap MSE
#'
#' Generic function to obtain bootstrap MSE estimates
#'
#' @inheritParams generate_NERM
#' @param obj_var Variance object
#' @param boot_samples List of bootstrap samples
#' @param var_u_est Variance of random effects
#' @param var_e_est Variance of errors
#' @param ... Additional parameters
#'
#' @return List with following parameters:
#' \item{mse_mixed}{Mean squared error}
#' \item{rmse_mixed}{Root mean squared error}
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
#' fitted_NERM <- fit_NERM(formula_y, data_sample,
#'                         id_cluster, cluster_means)
#'
#' type_method = c("parametric")
#' class(type_method) <- type_method
#'
#' # Bootstrap MSE
#'
#' boot_samples <- bootstrap_NERM(type_method,
#'                                var_u_est = fitted_NERM$var_u,
#'                                var_e_est = fitted_NERM$var_e,
#'                                beta_est = fitted_NERM$beta_hat,
#'                                boot_seed = 2,
#'                                n_boot = 100, formula_y,
#'                                data_sample = data_sample,
#'                                id_cluster = id_cluster,
#'                                cluster_means = cluster_means,
#'                                type_var_estimator = "mse_b")
#'
#' type_var_estimator = "mse_b"
#' class(type_var_estimator) <- "mse_b"
#'
#' boot_mixed <- boot_MSE(type_var_estimator, boot_samples,
#'                        cluster_means, var_u_est, var_e_est)
#'
#' # Bias-corrected bootstrap MSE
#'
#' boot_samples <- bootstrap_NERM(type_method,
#'                                var_u_est = fitted_NERM$var_u,
#'                                var_e_est = fitted_NERM$var_e,
#'                                beta_est = fitted_NERM$beta_hat,
#'                                boot_seed = 2,
#'                                n_boot = 100, formula_y,
#'                                data_sample = data_sample,
#'                                id_cluster = id_cluster,
#'                                cluster_means = cluster_means,
#'                                type_var_estimator = "mse_bc")
#'
#' type_var_estimator = "mse_bc"
#' class(type_var_estimator) <- "mse_bc"
#'
#' boot_mixed <- boot_MSE(type_var_estimator, boot_samples,
#'                        cluster_means, var_u_est, var_e_est)
#'


boot_MSE <- function(...)
  UseMethod("boot_MSE")


#'
#' @describeIn boot_MSE Function to compute bootstrap MSE estimates
#' @export
#'


boot_MSE.mse_b <- function(obj_var, boot_samples, ...) {
  boot_NERMS  = boot_samples$boot_NERMS
  fit_boot_NERMS  = boot_samples$fit_boot_NERMS
  mu_matrix <- matrix(0,
                      ncol = length(fit_boot_NERMS[[1]]$mu_hat),
                      nrow = length(fit_boot_NERMS))
  mu_hat_matrix <-
    matrix(0,
           ncol = length(fit_boot_NERMS[[1]]$mu_hat),
           nrow = length(fit_boot_NERMS))

  for (i in 1:length(boot_samples$boot_NERMS)) {
    mu_matrix[i,] <- unique(boot_NERMS[[i]]$mu)
    mu_hat_matrix[i,] <- fit_boot_NERMS[[i]]$mu_hat
  }

  mse_mixed <- colMeans((mu_hat_matrix - mu_matrix) ^ 2)
  rmse_mixed = sqrt(mse_mixed)

  output <- list(mse_mixed = mse_mixed,
                 rmse_mixed = rmse_mixed)
  return(output)

}


#'
#' @describeIn boot_MSE Function to compute 3-term bootstrap MSE estimates
#' @export
#'


boot_MSE.mse_3t <- function(obj_var,
                            boot_samples,
                            cluster_means,
                            var_u_est,
                            var_e_est,
                            ...) {
  id_cluster <- boot_samples$boot_NERMS[[1]]$id_cluster
  fit_boot_NERMS <- boot_samples$fit_boot_NERMS
  data_sample <- boot_samples$boot_NERMS[[1]]
  X <- format_data_matrix(data = data_sample,
                          name_col = "X",
                          select_col = "x|X")


  for (i in 1:length(boot_samples$boot_NERMS)) {
    fit_boot_NERMS[[i]]$y_star <- boot_samples$boot_NERMS[[i]]$y
    fit_boot_NERMS[[i]]$id_cluster <-
      boot_samples$boot_NERMS[[i]]$id_cluster
  }


  # Compute bootstrap parameters
  boot_params <-
    lapply(fit_boot_NERMS,
           compute_boot_params,
           var_e_est,
           var_u_est,
           cluster_means,
           X)

  # Compute mse 3 terms
  mu_matrix <- matrix(0,
                      ncol = length(fit_boot_NERMS[[1]]$mu_hat),
                      nrow = length(fit_boot_NERMS))

  mu_hat_boot_mixed_matrix <-
    matrix(0,
           ncol = length(fit_boot_NERMS[[1]]$mu_hat),
           nrow = length(fit_boot_NERMS))

  mu_hat_boot_matrix <-
    matrix(0,
           ncol = length(fit_boot_NERMS[[1]]$mu_hat),
           nrow = length(fit_boot_NERMS))

  for (i in 1:length(boot_samples$boot_NERMS)) {
    mu_matrix[i,] <- unique(boot_samples$boot_NERMS[[i]]$mu)
    mu_hat_boot_mixed_matrix[i,] <-
      boot_params[[i]]$mu_hat_boot_mixed
    mu_hat_boot_matrix[i,] <-  boot_params[[i]]$mu_hat_boot
  }

  mse_mixed = colMeans((mu_hat_boot_mixed_matrix - mu_matrix) ^ 2) +
    colMeans((mu_hat_boot_mixed_matrix - mu_hat_boot_matrix) ^ 2) +
    2 * colMeans((mu_hat_boot_mixed_matrix - mu_matrix) * (mu_hat_boot_mixed_matrix - mu_hat_boot_matrix)
    )

  rmse_mixed = sqrt(mse_mixed)

  output <- list(mse_mixed = mse_mixed,
                 rmse_mixed = rmse_mixed)
  return(output)


}



#'
#' @describeIn boot_MSE Function to compute semiparametric bootstrap estimates
#' @export
#'


boot_MSE.mse_spa <- function(obj_var,
                             boot_samples,
                             cluster_means,
                             var_u_est,
                             var_e_est,
                             ...) {
  id_cluster <- boot_samples$boot_NERMS[[1]]$id_cluster
  fit_boot_NERMS <- boot_samples$fit_boot_NERMS

  for (i in 1:length(boot_samples$boot_NERMS)) {
    fit_boot_NERMS[[i]]$y_star <- boot_samples$boot_NERMS[[i]]$y
    fit_boot_NERMS[[i]]$id_cluster <-
      boot_samples$boot_NERMS[[i]]$id_cluster
  }
  data_sample <- boot_samples$boot_NERMS[[1]]
  X <- format_data_matrix(data = data_sample,
                          name_col = "X",
                          select_col = "x|X")

  # Compute bootstrap parameters
  boot_params <-
    lapply(fit_boot_NERMS,
           compute_boot_params,
           var_e_est,
           var_u_est,
           cluster_means,
           X)

  var_mixed <- boot_params[[1]]$var_mixed
  af_var_mixed <- boot_params[[1]]$af_var_mixed


  # Compute mse 3 terms
  mu_matrix <- matrix(0,
                      ncol = length(fit_boot_NERMS[[1]]$mu_hat),
                      nrow = length(fit_boot_NERMS))

  mu_hat_boot_mixed_matrix <-
    matrix(0,
           ncol = length(fit_boot_NERMS[[1]]$mu_hat),
           nrow = length(fit_boot_NERMS))

  mu_hat_boot_matrix <-
    matrix(0,
           ncol = length(fit_boot_NERMS[[1]]$mu_hat),
           nrow = length(fit_boot_NERMS))

  var_mixed_boot <- matrix(0,
                    ncol = length(fit_boot_NERMS[[1]]$mu_hat),
                    nrow = length(fit_boot_NERMS))

  af_var_mixed_boot <- matrix(0,
                    ncol = length(fit_boot_NERMS[[1]]$mu_hat),
                    nrow = length(fit_boot_NERMS))

  for (i in 1:length(boot_samples$boot_NERMS)) {
    mu_matrix[i,] <- unique(boot_samples$boot_NERMS[[i]]$mu)

    mu_hat_boot_mixed_matrix[i,] <-
      boot_params[[i]]$mu_hat_boot_mixed

    mu_hat_boot_matrix[i,] <-  boot_params[[i]]$mu_hat_boot

    var_mixed_boot[i,] <- boot_params[[i]]$var_mixed_boot

    af_var_mixed_boot[i,] <- boot_params[[i]]$af_var_mixed_boot
  }

  var_mixed2_boot <-  colMeans(var_mixed_boot + af_var_mixed_boot)

  mse_mixed = 2 * (var_mixed + af_var_mixed) - var_mixed2_boot +
    colMeans((mu_hat_boot_mixed_matrix - mu_hat_boot_matrix) ^ 2) +
    2 * colMeans((mu_hat_boot_mixed_matrix - mu_matrix) * (mu_hat_boot_mixed_matrix - mu_hat_boot_matrix)
    )

  rmse_mixed = sqrt(mse_mixed)

  output <- list(mse_mixed = mse_mixed,
                 rmse_mixed = rmse_mixed)
  return(output)

}


#'
#' @describeIn boot_MSE Function to compute bias-corrected
#' bootstrap MSE estimate
#' @export
#'


boot_MSE.mse_bc <- function(obj_var, boot_samples,
                            ...) {
  # First order MSE
  var_obj_init <- "mse_b"
  class(var_obj_init) <- var_obj_init

  mse_b <- boot_MSE(obj_var = var_obj_init, boot_samples)$mse_mixed

  # Second order MSE
  double_boot_NERMS  = boot_samples$double_boot_NERMS
  fit_double_boot_NERMS  = boot_samples$fit_double_boot_NERMS

  # Matrices to collect mu
  mu_db_matrix <-
    matrix(
      0,
      ncol = length(fit_double_boot_NERMS[[1]]$mu_hat),
      nrow = length(double_boot_NERMS)
    )
  mu_db_hat_matrix <-
    matrix(
      0,
      ncol = length(fit_double_boot_NERMS[[1]]$mu_hat),
      nrow = length(double_boot_NERMS)
    )

  for (i in 1:length(boot_samples$double_boot_NERMS)) {
    mu_db_matrix[i,] <- unique(double_boot_NERMS[[i]]$mu)
    mu_db_hat_matrix[i,] <- fit_double_boot_NERMS[[i]]$mu_hat
  }

  mse_bd <- colMeans((mu_db_hat_matrix - mu_db_matrix) ^ 2)
  mse_bc0 <- 2 * mse_b - mse_bd

  mse_mixed = numeric(length(fit_double_boot_NERMS[[1]]$mu_hat))

  for (i in 1:length(mse_bc0)) {
    if (mse_b[i] >= mse_bc0[i]) {
      mse_mixed[i] = mse_bc0[i]
    } else {
      mse_mixed[i] = mse_b[i] * exp((mse_b[i] - mse_bc0[i]) / mse_bc0[i])
    }
  }

  rmse_mixed = sqrt(mse_mixed)

  output <- list(mse_mixed = mse_mixed,
                 rmse_mixed = rmse_mixed)
  return(output)

}


#' Compute bootstrap parameters
#'
#' Function to compute parameters from bootstrap samples
#'
#' @inheritParams generate_NERM
#' @param fit_boot_NERMS List of parameters with fitted values from NERM
#' @param var_u_est Variance of random effects
#' @param var_e_est Variance of errors
#'
#' @return List with following parameters
#' \item{mu_hat_boot_mixed}{Parameter mu computed with elements from an original sample and a bootstrap sample}
#' \item{mu_hat_boot}{Parameter mu computed with elements from a bootstrap sample}
#' \item{var_mixed_boot}{Variance of bootstrap mixed effect}
#' \item{af_var_mixed_boot}{Adjustemnt factor for variance of bootstrap mixed effect}
#' \item{var_mixed}{Variance of mixed effect}
#' \item{af_var_mixed}{Adjustemnt factor for variance of mixed effect}
#'
#'
#'
#' @importFrom postcAIC create_Z
#'
#'
#'

compute_boot_params <- function(fit_boot_NERMS,
                                var_e_est,
                                var_u_est,
                                cluster_means,
                                X) {
  # Clusters, Z matrix, n_d and m
  id_cluster  = fit_boot_NERMS$id_cluster
  Z = create_Z("NERM", id_cluster)
  n_d  = as.numeric(table(id_cluster))
  m = length(n_d)
  n_total = length(id_cluster)

  # Bootstrap quantities
  var_u_boot  = fit_boot_NERMS$var_u
  var_e_boot  = fit_boot_NERMS$var_e
  beta_boot = fit_boot_NERMS$beta_hat

  gamma_boot = var_u_boot / (var_u_boot + var_e_boot / n_d)
  var_mixed_boot  = (1 - gamma_boot) * var_u_boot

  R_boot = var_e_boot * diag(n_total)
  G_boot = var_u_boot * diag(m)

  V_boot <- R_boot + kronecker(G_boot,  matrix(1, n_d, n_d))
  solve_V_boot <-
    1 / var_e_boot * diag(n_total) - kronecker(diag(m) * gamma_boot * 1 / (n_d * var_e_boot), matrix(1, n_d, n_d))

  y_star <- fit_boot_NERMS$y_star

  temp1_af_var_mixed_boot = cluster_means - gamma_boot * cluster_means
  temp2_af_var_mixed_boot = solve(crossprod(X, crossprod(solve_V_boot, X)))
  af_var_mixed_boot <- numeric(m)

  for (i in 1:m) {
    af_var_mixed_boot[i] = t(temp1_af_var_mixed_boot[i,]) %*% temp2_af_var_mixed_boot %*% (temp1_af_var_mixed_boot[i,])
  }

  # NERM quantities
  gamma_est <- var_u_est / (var_u_est + var_e_est / n_d)
  solve_V = 1 / var_e_est * diag(n_total) - kronecker(diag(m) * gamma_est * 1 / (n_d * var_e_est), matrix(1, n_d, n_d))
  G_est = var_u_est * diag(m)

  var_mixed = (1 - gamma_est) * var_u_est

  temp1_var_mixed = cluster_means - gamma_est * cluster_means
  temp2_var_mixed = solve(crossprod(X, crossprod(solve_V, X)))
  af_var_mixed <- numeric(m)

  for (i in 1:m) {
    af_var_mixed[i] = t(temp1_var_mixed[i,]) %*% temp2_var_mixed %*% (temp1_var_mixed[i,])
  }


  beta_mixed <- crossprod(t(tcrossprod(solve(crossprod(
    t(crossprod(X, solve_V)), X
  )), X)), crossprod(solve_V, y_star))


  mu_hat_boot_mixed <- crossprod(t(cluster_means), beta_mixed) +
    crossprod(t(tcrossprod(G_est, Z)), crossprod(solve_V, (y_star - crossprod(t(
      X
    ), beta_mixed))))

  mu_hat_boot <- crossprod(t(cluster_means), beta_boot) +
    crossprod(t(tcrossprod(G_boot, Z)), crossprod(solve_V_boot, (y_star - crossprod(t(
      X
    ), beta_boot))))

  output <- list(
    mu_hat_boot_mixed = mu_hat_boot_mixed,
    mu_hat_boot = mu_hat_boot,
    var_mixed_boot = var_mixed_boot,
    af_var_mixed_boot = af_var_mixed_boot,
    var_mixed = var_mixed,
    af_var_mixed = af_var_mixed
  )
  return(output)

}
