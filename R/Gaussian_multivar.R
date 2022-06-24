#' Calculate PPI using the Gaussian prior
#'
#' @inheritParams laplace_mv
#' @param w The Gaussian effect size prior variance.
#' @return A list consisting of (in order) the posterior inclusion probability
#'   by SNP, the posterior probability of each model, the prior_probabilities
#'   used (these are the same as those specified in the function arguments).
#'   Additionally if \code{allow_no_causal} is \code{TRUE} then the posterior
#'   probability of the null model (no causal SNPs) is also returned as the
#'   fourth list element.
#' @export
#' @examples
#' out <- gauss_mv(num_snps = 10, beta_hats = est_log_ors,
#' V_inv = solve(est_cov_mat), maxk = 2, w = 0.04,
#' prior_probs = c(0.5, 0.5), allow_no_causal = FALSE)
#' out[[1]] # PPIs of the 10 SNPs
#' out[[2]] # Posterior probabilities of the 1-SNP and 2-SNP models
#' out[[3]] # prior probabilities of the number of causal SNPs

gauss_mv <- function(num_snps, beta_hats, V_inv, maxk, w, prior_probs,
                     allow_no_causal){
  # check consistency in prior probs and 'allow_no_causal'
  # define a new prior vector without allowing no causal SNPs
  if(allow_no_causal == T){
    prior_probs_no_zero <- prior_probs[-1]
    if(length(prior_probs) != (maxk + 1)){
      stop("number of prior probs not consistent with values of other arguments")
    }
  } else{
    prior_probs_no_zero <- prior_probs
    if(length(prior_probs) != maxk){
      stop("number of prior probs not consistent with values of other arguments")
    }
  }

  # set up marginal likelihood list for all the non-null models
  weighted_ml <- vector("list", length = maxk)

  for(k in 1:maxk){
    #cat(k,"causal snps in the model" ,"\n")
    W = rep(w, k)
    det_W <- prod(W)
    W_inv <- matrix(1/W, nrow = k, ncol= k)

    # determine all possible combinations of causal SNPs
    snp_combs <- combine_snps(num_snps, k)
    num_combs <- dim(snp_combs)[1]

    #set up marginal likelihood matrix for given number of causal SNPs in the model
    weighted_ml[[k]] <- cbind(snp_combs, NA)

    # loop through all SNP combinations
    for(i in 1:num_combs){
      causal_snps_in_model <- snp_combs[i, ]
      sigma_c <- V_inv[causal_snps_in_model, causal_snps_in_model, drop = F]
      Omega <- W_inv + sigma_c
      V_inv_star <- V_inv[, causal_snps_in_model, drop = F]
      Q <- V_inv_star %*% solve(Omega) %*% t(V_inv_star)
      R <- t(beta_hats) %*% (V_inv - Q) %*% beta_hats
      ml_gauss <- (det(Omega) * det_W) ^ (-0.5) * exp(-R / 2)
      # prior * marginal likelihood
      weighted_ml[[k]][i, (k + 1)] <- ml_gauss *
        prior_probs_no_zero[k] / num_combs
    } # end of i loop
  } # end of k loop

  # initialise post_model_prob to be the same as weighted_ml
  post_model_prob <- weighted_ml

  #P(D) [denominator]
  if(allow_no_causal == F){
    tot_weighted_ml <- as.vector(sum(sapply(weighted_ml, sum_last_col)))

  } else {
    ml_null_model <- exp(-0.5 * t(beta_hats) %*% V_inv %*% beta_hats)
    weighted_ml_null <- as.vector(ml_null_model * prior_probs[1])
    tot_weighted_ml <- as.vector(weighted_ml_null +
                                 sum(sapply(weighted_ml, sum_last_col)))
    post_null_prob <- weighted_ml_null / tot_weighted_ml
  }

  for(k in 1 :maxk){
    post_model_prob[[k]][, k + 1] <- weighted_ml[[k]][, k + 1] / tot_weighted_ml
  }

  # check that sum = 1
  # cat("sum of posterior models probabilities = ",
  #    as.vector(post_null_prob + sum(sapply(post_model_prob, sum_last_col))), "\n")

  # calculate posterior probability for each snp
  pp_by_snp <- rep(NA, num_snps) # vector of snp posterior probs
  comps <- sapply(1:num_snps, function(y)
    sapply(weighted_ml, sum_last_col_for_snp, snp = y))
  if(maxk == 1) {pp_by_snp <- comps / tot_weighted_ml} else {
    pp_by_snp <- apply(comps, 2, sum) / tot_weighted_ml
  }

  if(allow_no_causal == F){
    return(list(pp_by_snp, post_model_prob, prior_probs))
  } else {
    return(list(pp_by_snp, post_model_prob, prior_probs, post_null_prob))
  }
} # end of function loop



