% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/construct_intervals.R
\name{construct_intLME}
\alias{construct_intLME}
\alias{construct_intLME.NERM}
\title{Intervals for linear mixed models under NERM}
\usage{
construct_intLME(...)

\method{construct_intLME}{NERM}(
  model_obj,
  type_method = c("parametric", "semiparametric", "analytical"),
  type_var_estimator = c("var_mixed", "mse_a", "mse_b", "mse_3t", "mse_spa", "mse_bc"),
  formula_y,
  data_sample,
  cluster_means,
  id_cluster,
  alpha,
  n_boot,
  boot_seed = 1,
  ...
)
}
\arguments{
\item{...}{Additional parameters}

\item{model_obj}{Object of class model}

\item{type_method}{Method to construct intervals}

\item{type_var_estimator}{Type of variance estimator}

\item{formula_y}{Formula for fitting NERM}

\item{data_sample}{Data frame}

\item{cluster_means}{Cluster-level covariates for fixed parameters}

\item{id_cluster}{Vector with cluster labels}

\item{alpha}{Level alpha to construct intervals}

\item{n_boot}{Number of bootstrap samples, default n_boot = 1000}

\item{boot_seed}{Vector with cluster labels}
}
\value{
List with following parameters:
\item{mu_hat}{Estimated mixed effects}
\item{int_up}{Upper boundary of individual interval}
\item{int_down}{Lower boundary of individual interval}
\item{length_ind}{Length of individual interval}
\item{int_sim_up}{Upper boundary of simultaneous interval}
\item{int_sim_down}{Lower boundary of simultaneous interval}
\item{length_sim}{Length of simulaneous interval}
}
\description{
Generic function to construct individual and simultaneous intervals
for linear mixed models
}
\section{Methods (by class)}{
\itemize{
\item \code{NERM}: Method to construct individual and simultanous intervvals for NERM
}}

\examples{

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

cluster_means <- as.matrix(aggregate(X,
list(id_cluster),
FUN = mean)[, -1])

#Formula to construct intervals
formula_y <-  y ~ -1 + X0 + X1 + (1| id_cluster)

data_sample <- generate_NERM(generate_u = list(type_dist = "chisquared",
                                               scaling_factor = 1,
                                               dg = 6),
                             generate_e = list(type_dist = "chisquared",
                                               scaling_factor = 1,
                                               dg = 6),
                             beta = beta,
                             X = X,
                             id_cluster = id_cluster,
                             start_seed = 1,
                             no_sim = 1,
                             cluster_means = cluster_means)

model_obj = c("NERM")
# For now, only NERM is supported, but other model choices can be easily programmed
class(model_obj) <- model_obj

construct_LMEintervals_pb <- construct_intLME(model_obj,
                                           type_method = c("parametric"),
                                           type_var_estimator = c("mse_bc"),
                                           formula_y, data_sample,
                                           cluster_means = cluster_means, id_cluster,
                                           alpha = 0.05, n_boot = 100,
                                           boot_seed = 1)

construct_LMEintervals_spb <- construct_intLME(model_obj,
                                           type_method = c("semiparametric"),
                                           type_var_estimator = c("mse_bc"),
                                           formula_y, data_sample,
                                           cluster_means = cluster_means, id_cluster,
                                           alpha = 0.05, n_boot = 100,
                                           boot_seed = 1)

}
