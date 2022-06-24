#### Dataset documentation: est_log_ors ####

#' Effect size estimates of 10 SNPs
#'
#' A vector of length 10 containing the estimates of 10 SNP effect sizes
#' obtained from a multi-SNP logistic regression model.
#'
#' @source {Data is simulated using Hapgen2 and a multi-SNP logistic regression
#'   model is fitted. These are the estimated log-odds ratios from the
#'   regression output.}
#' @examples
#' data(est_log_ors)
"est_log_ors"

#### Dataset documentation: est_log_ors ####

#' covariance matrix of 10 SNPs
#'
#' A 10 by 10 covariance matrix containing the pairwise covariances of the
#' estimates of 10 SNP effect sizes obtained from a multi-SNP logistic
#' regression model.
#'
#' @source {Data is simulated using Hapgen2 and a multi-SNP logistic regression
#'   model is fitted. This is the covariance matrix of the 10 estimated log-odds
#'   ratios obtained from the regression output.}
#' @examples
#' data(est_cov_mat)
"est_cov_mat"
