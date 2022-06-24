#' @importFrom utils combn

# obtain snps combination for each model
combine_snps <- function(p, k) t((combn(seq(1:p), k)))

# adds together the elements in the final column of a matrix
sum_last_col <-  function(x) sum(x[, dim(x)[2]])

# sums the final column of x for rows that contain the number 'snp'
sum_last_col_for_snp <-  function(x, snp) {
  log_vec <- apply(x, 1, function(z) any(z == snp))
  final_col <- dim(x)[2]
  sum(x[log_vec, final_col])
}
