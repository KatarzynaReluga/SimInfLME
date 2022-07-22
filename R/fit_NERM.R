#' Fit NERM
#'
#' Function to fit nested error regression model
#'
#' @inheritParams construct_intervals
#' @inheritParams generate_NERM
#'
#' @importFrom postcAIC compute_corrected_mse
#' @importFrom lme4 lmer VarCorr fixef ranef
#' @importFrom stats predict
#'
#' @return List with following parameters:
#' \item{var_u}{Estimated variance of random effects}
#' \item{var_e}{Estimated variance of errors}
#' \item{beta_hat}{Estimated regression parameters}
#' \item{u_hat}{Estimated random effects}
#' \item{mu_hat}{Estimated mixed effects}
#' \item{e_hat}{Estimated errors}
#' \item{mse}{Analytical mse}
#' \item{g1}{Analytical variance of mixed effect}
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
#'                              cluster_means)
#'
#'
#' fitted_NERM <- fit_NERM(formula_y, data_sample,
#'                         id_cluster, cluster_means)
#'


fit_NERM <- function(formula_y,
                     data_sample,
                     id_cluster,
                     cluster_means) {
  m = length(as.numeric(table(id_cluster)))
  id_cluster = data_sample$id_cluster
  X <- format_data_matrix(data = data_sample,
                          name_col = "X",
                          select_col = "x|X")

  fit_y = lmer(formula_y, data = data_sample)
  var_u = ((as.data.frame(VarCorr(fit_y)))$vcov)[1]
  var_e = ((as.data.frame(VarCorr(fit_y)))$vcov)[2]
  beta_hat = matrix(t(fixef(fit_y)), ncol = 1, nrow = 2)
  u_hat = as.vector(unlist(ranef(fit_y)))
  mu_hat = crossprod(t(cluster_means), beta_hat) + u_hat
  e_hat = data_sample$y - (predict(fit_y) + rep(u_hat, as.numeric(table(id_cluster))))

  C_cluster <- cbind(cluster_means, diag(m))

  #Compute analytical MSE
  mse = compute_corrected_mse(
    C_cluster = C_cluster,
    X = X,
    sig_u = var_u,
    sig_e = var_e,
    clusterID = id_cluster
  )$mse

  temp = var_u / (var_u + var_e / as.numeric(table(id_cluster)))

  # Compute g1
  g1 = (1 - temp) * var_u

  output <- list(
    var_u = var_u,
    var_e = var_e,
    beta_hat = beta_hat,
    u_hat = u_hat,
    mu_hat = mu_hat,
    e_hat = e_hat,
    mse = mse,
    g1 = g1
  )
  return(output)
}
