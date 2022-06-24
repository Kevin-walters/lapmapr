#' Calculate PPI using the Laplace prior
#'
#' @param num_snps The number of SNPs in your dataset.
#' @param beta_hats A vector of estimates of the SNP effect sizes (log-odds
#'   ratios).
#' @param V_inv The inverse of the covariance matrix (of dim \code{num_snps} by
#'   \code{num_snps}) of the estimated SNP effect sizes.
#' @param maxk The maximum number of causal SNPs allowed in the model.
#' @param lambda The rate parameter of the Laplace prior.
#' @param prior_probs A vector of prior probabilities for the number of causal
#'   SNPs. For example if you want to allow no more than two causal SNPs with
#'   either one and two causal SNPs equally likely a priori then use
#'   \code{c(0.5, 0.5)}. You can allow for the possibility that there are no
#'   causal SNPs but the flag \code{allow_no_causal} must be \code{TRUE}. For
#'   example to specify that the prior probability of no causal SNPs is 0.5 and
#'   that the prior probability of either one or two causal SNPs is 0.25 use
#'   \code{c(0.5, 0.25, 0.25)}.
#' @param allow_no_causal Are you allowing for the possibility that there are no
#'   causal SNPs? I.e. is the first element of \code{prior_probs} the
#'   probability that there are no causal SNPs a priori?
#' @return A list consisting of (in order) the posterior inclusion probability
#'   by SNP, the posterior probability of each model, the prior_probabilities
#'   used (these are the same as those specified in the function arguments).
#'   Additionally if \code{allow_no_causal} is \code{TRUE} then the posterior
#'   probability of the null model (no causal SNPs) is also returned as the
#'   fourth list element.
#' @export
#' @examples
#' out <- laplace_mv(num_snps = 10, beta_hats = est_log_ors,
#' V_inv = solve(est_cov_mat), maxk = 2, lambda = 5,
#' prior_probs = c(0.5, 0.5), allow_no_causal = FALSE)
#' out[[1]] # PPIs of the 10 SNPs
#' out[[2]] # Posterior probabilities of the 1-SNP and 2-SNP models
#' out[[3]] # prior probabilities of the number of causal SNPs


laplace_mv <- function(num_snps, beta_hats, V_inv, maxk, lambda, prior_probs,
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

  #loop through the number of causal snps in the model
  for(k in 1:maxk){
    #cat(k,"causal snps in the model" ,"\n")

    # determine all possible combinations of causal SNPs
    snp_combs <- combine_snps(num_snps, k)
    num_combs <- dim(snp_combs)[1]

    #set up marginal likelihood matrix for given number of causal SNPs
    # in the model
    weighted_ml[[k]] <- cbind(snp_combs, NA)

    #determine all combinations of A
    A_comb <- gtools::permutations(2, k, c(-1,1), set = F, repeats.allowed = T)
    num_A_combs <- dim(A_comb)[1]

    #set up matrix for lower and upper limit for each combination of A to
    # be used in obtaining cdf
    lower_limit <- A_comb * 0
    lower_limit <- replace(lower_limit, A_comb == -1, -Inf)
    upper_limit <- A_comb * 0
    upper_limit <- replace(upper_limit, A_comb == 1, Inf)

    #loop through all snp combinations
    for(i in 1:num_combs){
      #cat("i= ", i, "\t")
      #cat("snp combs = ", snp_combs[i,], "\n")
      causal_snps_in_model <- snp_combs[i, ]
      sigma_c <- V_inv[causal_snps_in_model, causal_snps_in_model, drop = F]
      nu_inv <- solve(sigma_c)
      V_inv_star <- V_inv[, causal_snps_in_model, drop = F]
      U <- V_inv_star %*%  nu_inv
      T_1 <- t(beta_hats) %*% (V_inv - (U %*% t(V_inv_star))) %*% beta_hats
      gdata::lowerTriangle(nu_inv) <- gdata::upperTriangle(nu_inv, byrow = T)


      #loop through combinations of A
      A_combs <- vector("numeric", length = num_A_combs)
      for(j in 1:num_A_combs){
        #cat("j= ", j, " ")
        A <- t(A_comb[j,, drop = F])
        mu <- nu_inv %*% (t(V_inv_star) %*% beta_hats - (lambda * A))
        T_2 <- 2 * lambda * t(beta_hats) %*% U %*% A -
          (lambda ^ 2) * t(A) %*% nu_inv %*% A
        t <- T_1 + T_2

        #calculates CDF
        cdf_mvnorm <- mvtnorm::pmvnorm(lower = as.vector(lower_limit[j,]),
                                       upper=as.vector(upper_limit[j,]),
                                       mean = as.vector(mu),
                                       sigma = nu_inv)
        #sigma = solve(sigma_c))
        #if(is.nan(log(cdf_mvnorm))) cat(i,"\t " ,j, "\t", cdf_mvnorm, "\n")
        A_combs[j] <- ifelse(cdf_mvnorm> 0, exp(-0.5 * t + log(cdf_mvnorm)), 0)
      } #end of A combs loop
      #cat("\n")
      ml_lap <- (lambda/2) ^ k * (2*pi) ^ (k/2) * det(sigma_c) ^ (-0.5) *
        sum(A_combs)
      weighted_ml[[k]][i, (k + 1)] <- ml_lap * prior_probs_no_zero[k] / num_combs
    } #end of snp combs loop
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
  #cat("sum of post model probs = ",
  #   as.vector(post_null_prob + sum(sapply(post_model_prob, sum_last_col))), "\n")

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
}
