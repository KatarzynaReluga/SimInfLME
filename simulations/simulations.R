###################################################
## Simulations to reproduce results in the paper ##
###################################################
# Basic setup
set.seed(11992)
m = 25
n_j = 5
n = m * n_j
scaling_factor_e = 1
scaling_factor_u = sqrt(0.5)
x_ij = runif(n, 0, 1)
X = cbind(X0 = rep(1, n), X1 = x_ij)
t_X  = t(X)
id_cluster = rep(1:m, each = n_j)
beta = c(1, 1)
no_sim = 1000

cluster_means <- as.matrix(aggregate(X,
                           list(id_cluster), FUN = mean)[, -1])

#Formula to construct intervals
formula_y <-  y ~ -1 + X0 + X1 + (1| id_cluster)

# Normally distributed errors and random effects

data_samples <- generate_NERM(generate_u = list(type = "normal",
                                               scaling_factor = scaling_factor_u),
                              generate_e = list(type = "normal",
                                                scaling_factor = scaling_factor_e),
                              beta = beta,
                              X = X,
                              id_cluster = id_cluster,
                              start_seed = 1,
                              no_sim = no_sim,
                              cluster_means = cluster_means)

# Analytical intervals

analytical_interval <- function(data_sample) {

  interval <- construct_intervals(type_method = "analytical", type_var_estimator = "mse_a",
                                             model_type = c("NERM"),
                                             formula_y,
                                             data_sample = data_sample,
                                             cluster_means = cluster_means, id_cluster,
                                             alpha = 0.05, n_boot = 1000,
                                             boot_seed = 1)
  interval
}

analytical_intervals <- lapply(data_samples, analytical_interval)

# Add true mixed parameter

for (i in 1:length(analytical_intervals)) {
  analytical_intervals[[i]]$mu <- unique(data_samples[[i]]$mu)
}
















