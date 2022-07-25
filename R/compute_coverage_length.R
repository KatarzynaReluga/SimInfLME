#' Compute coverage and length
#'
#' Function to calculate the empirical coverage, the average length and the variance
#' of lengths of individual and simultaneous intervals
#' over simulation runs
#'
#' @param intervals List of individual and simultaneous intervals from simulations
#'
#' @importFrom stats var
#'
#' @return List with following parameters:
#' \item{cov_ind}{Empircial coverage of individual intervals}
#' \item{average_length_ind}{Average length of individual intervals}
#' \item{var_length_ind}{Variance of length of individual intervals}
#' \item{cov_sim}{Empircial coverage of simultaneous intervals}
#' \item{average_length_sim}{Average length of simultaneous intervals}
#' \item{var_length_sim}{Variance of length of simultaneous intervals}
#'
#'
#' @export
#'


compute_coverage_length <- function(intervals) {


  coverage_list <- lapply(intervals, empirical_coverage)

  int_ind_matrix <- matrix(0, nrow = length(intervals),
                           ncol = length(intervals[[1]]$mu))

  int_sim_matrix <- matrix(0, nrow = length(intervals),
                           ncol = length(intervals[[1]]$mu))

  tests_ind <- matrix(0, nrow = length(intervals),
                    ncol = length(intervals[[1]]$mu))

  tests_sim <- numeric(length(intervals))

  for (i in 1:length(intervals)) {
    int_ind_matrix[i, ] <- intervals[[i]]$length_ind
    int_sim_matrix[i, ] <- intervals[[i]]$length_sim
    tests_ind[i, ] <- coverage_list[[i]]$test_ind
    tests_sim[i] <- coverage_list[[i]]$test_sim
  }

  cov_ind <- mean(colMeans(tests_ind))
  average_length_ind = mean(colMeans(int_ind_matrix))
  var_length_ind = mean(apply(int_ind_matrix, 2, var))

  cov_sim <- mean(tests_sim)
  average_length_sim = mean(int_sim_matrix)
  var_length_sim = mean(apply(int_sim_matrix, 2, var))

  output <- list(cov_ind = cov_ind,
                 average_length_ind = average_length_ind,
                 var_length_ind = var_length_ind,
                 cov_sim = cov_sim,
                 average_length_sim = average_length_sim,
                 var_length_sim = var_length_sim)
  output

}

#' Compute empirical coverage
#'
#' Function to calculate the empirical coverage of individual and simultaneous intervals
#' over simulation runs
#'
#' @param intervals_one_sample List of individual and simultaneous intervals from
#' one simulation simulations
#'
#' @return List with following parameters:
#' \item{test_ind}{Empircial coverage of individual intervals}
#' \item{test_sim}{Empircial coverage of simultaneous intervals}
#'
#'


# call it test
empirical_coverage <- function(intervals_one_sample) {

  mu = intervals_one_sample$mu
  int_down = intervals_one_sample$int_down
  int_up = intervals_one_sample$int_up

  int_sim_down = intervals_one_sample$int_sim_down
  int_sim_up = intervals_one_sample$int_sim_up

  test_ind <- (int_down <= mu) & (mu <= int_up)
  test_sim = sum(((int_sim_down <= mu) & (mu <= int_sim_up))) == length(mu)

  output <- list(test_ind = test_ind,
                 test_sim = test_sim)
  output
}
