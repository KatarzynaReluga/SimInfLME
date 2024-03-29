% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_MSE.R
\name{compute_boot_params}
\alias{compute_boot_params}
\title{Compute bootstrap parameters}
\usage{
compute_boot_params(fit_boot_NERMS, var_e_est, var_u_est, cluster_means, X)
}
\arguments{
\item{fit_boot_NERMS}{List of parameters with fitted values from NERM}

\item{var_e_est}{Variance of errors}

\item{var_u_est}{Variance of random effects}

\item{cluster_means}{Cluster-level covariates for fixed parameters}

\item{X}{Matrix with covariates}
}
\value{
List with following parameters
\item{mu_hat_boot_mixed}{Parameter mu computed with elements from an original sample and a bootstrap sample}
\item{mu_hat_boot}{Parameter mu computed with elements from a bootstrap sample}
\item{var_mixed_boot}{Variance of bootstrap mixed effect}
\item{af_var_mixed_boot}{Adjustemnt factor for variance of bootstrap mixed effect}
\item{var_mixed}{Variance of mixed effect}
\item{af_var_mixed}{Adjustemnt factor for variance of mixed effect}
}
\description{
Function to compute parameters from bootstrap samples
}
