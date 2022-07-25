##############################################
## Code to reproduce simulation results in  ##
## in Realuga and Sperlich (2022)           ##
##############################################
# Basic setup --------------------------------------
set.seed(10)
m = 25
n_j = 5
n = m * n_j
scaling_factor_e = sqrt(0.5)
scaling_factor_u = 1
x_ij = runif(n, 0, 1)
X = cbind(X0 = rep(1, n), X1 = x_ij)
t_X  = t(X)
id_cluster = rep(1:m, each = n_j)
beta = c(1, 1)
no_sim = 10
alpha = 0.05
n_boot = 1000
boot_seed = 10

cluster_means <- as.matrix(aggregate(X,
                                     list(id_cluster), FUN = mean)[, -1])

#Formula to construct intervals ---------------------

formula_y <-  y ~ -1 + X0 + X1 + (1| id_cluster)

##################################################
# Normally distributed errors and random effects #
##################################################

data_samples <- generate_NERM(generate_u = list(type = "student_t",
                                                scaling_factor = scaling_factor_u, 
                                                dg = 6),
                              generate_e = list(type = "student_t",
                                                scaling_factor = scaling_factor_e, 
                                                dg = 6),
                              beta = beta,
                              X = X,
                              id_cluster = id_cluster,
                              start_seed = 1,
                              no_sim = no_sim,
                              cluster_means = cluster_means)


######################################
# Semiparametric bootstrap intervals #
######################################

## Compute intervals for mixed effects ------------------------------

## 1. Variability of mixed parameter: var_mixed
spb_int_var_mixed <- lapply(data_samples, construct_intervals,
                            type_method = "semiparametric",
                            type_var_estimator = "var_mixed",
                            model_type = c("NERM"),
                            formula_y = formula_y,
                            cluster_means = cluster_means,
                            id_cluster = id_cluster,
                            alpha = alpha,
                            n_boot = n_boot,
                            boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(spb_int_var_mixed)) {
  spb_int_var_mixed[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_var_mixed <- compute_coverage_length(spb_int_var_mixed)

########################################################################

## 2. Variability of mixed parameter: msa_a
spb_int_mse_a <- lapply(data_samples, construct_intervals,
                        type_method = "semiparametric",
                        type_var_estimator = "mse_a",
                        model_type = c("NERM"),
                        formula_y = formula_y,
                        cluster_means = cluster_means,
                        id_cluster = id_cluster,
                        alpha = alpha,
                        n_boot = n_boot,
                        boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(spb_int_mse_a)) {
  spb_int_mse_a[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_a <- compute_coverage_length(spb_int_mse_a)

########################################################################

## 3. Variability of mixed parameter: msa_b
spb_int_mse_b <- lapply(data_samples, construct_intervals,
                        type_method = "semiparametric",
                        type_var_estimator = "mse_b",
                        model_type = c("NERM"),
                        formula_y = formula_y,
                        cluster_means = cluster_means,
                        id_cluster = id_cluster,
                        alpha = alpha,
                        n_boot = n_boot,
                        boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(spb_int_mse_b)) {
  spb_int_mse_b[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_b <- compute_coverage_length(spb_int_mse_b)

########################################################################

## 4. Variability of mixed parameter: msa_bc
spb_int_mse_bc <- lapply(data_samples, construct_intervals,
                         type_method = "semiparametric",
                         type_var_estimator = "mse_bc",
                         model_type = c("NERM"),
                         formula_y = formula_y,
                         cluster_means = cluster_means,
                         id_cluster = id_cluster,
                         alpha = alpha,
                         n_boot = n_boot,
                         boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(spb_int_mse_bc)) {
  spb_int_mse_bc[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_bc <- compute_coverage_length(spb_int_mse_bc)

########################################################################

## 5. Variability of mixed parameter: msa_3t
spb_int_mse_3t <- lapply(data_samples, construct_intervals,
                         type_method = "semiparametric",
                         type_var_estimator = "mse_3t",
                         model_type = c("NERM"),
                         formula_y = formula_y,
                         cluster_means = cluster_means,
                         id_cluster = id_cluster,
                         alpha = alpha,
                         n_boot = n_boot,
                         boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(spb_int_mse_3t)) {
  spb_int_mse_3t[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_3t <- compute_coverage_length(spb_int_mse_3t)

########################################################################

## 6. Variability of mixed parameter: msa_3t
spb_int_mse_spa <- lapply(data_samples, construct_intervals,
                          type_method = "parametric",
                          type_var_estimator = "mse_spa",
                          model_type = c("NERM"),
                          formula_y = formula_y,
                          cluster_means = cluster_means,
                          id_cluster = id_cluster,
                          alpha = alpha,
                          n_boot = n_boot,
                          boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(spb_int_mse_spa)) {
  spb_int_mse_spa[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_spa <- compute_coverage_length(spb_int_mse_spa)


##################################
# Parametric bootstrap intervals #
##################################

## Compute intervals for mixed effects ------------------------------

## 1. Variability of mixed parameter: var_mixed
pb_int_var_mixed <- lapply(data_samples, construct_intervals,
                           type_method = "parametric",
                           type_var_estimator = "var_mixed",
                           model_type = c("NERM"),
                           formula_y = formula_y,
                           cluster_means = cluster_means,
                           id_cluster = id_cluster,
                           alpha = alpha,
                           n_boot = n_boot,
                           boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(pb_int_var_mixed)) {
  pb_int_var_mixed[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_var_mixed <- compute_coverage_length(pb_int_var_mixed)

########################################################################

## 2. Variability of mixed parameter: msa_a
pb_int_mse_a <- lapply(data_samples, construct_intervals,
                       type_method = "parametric",
                       type_var_estimator = "mse_a",
                       model_type = c("NERM"),
                       formula_y = formula_y,
                       cluster_means = cluster_means,
                       id_cluster = id_cluster,
                       alpha = alpha,
                       n_boot = n_boot,
                       boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(pb_int_mse_a)) {
  pb_int_mse_a[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_a <- compute_coverage_length(pb_int_mse_a)

########################################################################

## 3. Variability of mixed parameter: msa_b
pb_int_mse_b <- lapply(data_samples, construct_intervals,
                       type_method = "parametric",
                       type_var_estimator = "mse_b",
                       model_type = c("NERM"),
                       formula_y = formula_y,
                       cluster_means = cluster_means,
                       id_cluster = id_cluster,
                       alpha = alpha,
                       n_boot = n_boot,
                       boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(pb_int_mse_b)) {
  pb_int_mse_b[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_b <- compute_coverage_length(pb_int_mse_b)

########################################################################

## 4. Variability of mixed parameter: msa_bc
pb_int_mse_bc <- lapply(data_samples, construct_intervals,
                        type_method = "parametric",
                        type_var_estimator = "mse_bc",
                        model_type = c("NERM"),
                        formula_y = formula_y,
                        cluster_means = cluster_means,
                        id_cluster = id_cluster,
                        alpha = alpha,
                        n_boot = n_boot,
                        boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(pb_int_mse_bc)) {
  pb_int_mse_bc[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_bc <- compute_coverage_length(pb_int_mse_bc)

########################################################################

## 5. Variability of mixed parameter: msa_3t
pb_int_mse_3t <- lapply(data_samples, construct_intervals,
                        type_method = "parametric",
                        type_var_estimator = "mse_3t",
                        model_type = c("NERM"),
                        formula_y = formula_y,
                        cluster_means = cluster_means,
                        id_cluster = id_cluster,
                        alpha = alpha,
                        n_boot = n_boot,
                        boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(pb_int_mse_3t)) {
  pb_int_mse_3t[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_3t <- compute_coverage_length(pb_int_mse_3t)

########################################################################

## 6. Variability of mixed parameter: msa_3t
pb_int_mse_spa <- lapply(data_samples, construct_intervals,
                         type_method = "parametric",
                         type_var_estimator = "mse_spa",
                         model_type = c("NERM"),
                         formula_y = formula_y,
                         cluster_means = cluster_means,
                         id_cluster = id_cluster,
                         alpha = alpha,
                         n_boot = n_boot,
                         boot_seed = boot_seed)

## Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(pb_int_mse_spa)) {
  pb_int_mse_spa[[i]]$mu <- unique(data_samples[[i]]$mu)
}

## Compute coverage and lengths ----------------------------------------
cov_len_mse_spa <- compute_coverage_length(pb_int_mse_spa)

###########################
## Analytical intervals  ##
###########################
analytical_interval <- function(data_sample) {
  
  interval <- construct_intervals(type_method = "analytical",
                                  type_var_estimator = "mse_a",
                                  model_type = c("NERM"),
                                  formula_y,
                                  data_sample = data_sample,
                                  cluster_means = cluster_means, id_cluster,
                                  alpha = 0.05, n_boot = 100,
                                  boot_seed = 1)
  interval
}

analytical_intervals <- lapply(data_samples, analytical_interval)


# Add true mixed parameter to the list with intervals -----------------

for (i in 1:length(analytical_intervals)) {
  analytical_intervals[[i]]$mu <- unique(data_samples[[i]]$mu)
}

# Compute coverage and lengths ----------------------------------------

coverage_length <- compute_coverage_length(analytical_intervals)
