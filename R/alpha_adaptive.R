# Adaptive alpha adjustment for tree-structured hypothesis testing
#
# Implements Algorithm 1 from Appendix D. The key idea: at each tree level,
# adjust the significance threshold to account for both the number of tests
# (multiplicity) and the estimated power at that level. When cumulative
# power is low, the probability of reaching a node by chance is already
# small, so natural gating may provide adequate FWER control in practice
# (but the formal guarantee requires tau = 0; see Remark in Appendix D).


#' Compute Adaptive Alpha Levels by Tree Depth
#'
#' Implements Algorithm 1 from the adaptive alpha appendix. Computes
#' adjusted significance levels for each tree level based on estimated
#' power decay through the tree.
#'
#' @param k Branching factor. Either a scalar (constant k at all levels)
#'   or an integer vector of length \code{max_depth} where \code{k[ell]}
#'   is the branching factor at level \code{ell}.
#' @param delta_hat Estimated standardized effect size (e.g., Cohen's d).
#'   Conservative (larger) values produce more stringent adjustment,
#'   which preserves the FWER guarantee. Use an upper bound on the
#'   true effect size.
#' @param N_total Total sample size at the root level.
#' @param tau Cumulative power threshold (default 0.1). When cumulative
#'   power drops below \code{tau}, natural gating is deemed sufficient
#'   and nominal alpha is used.
#' @param max_depth Maximum depth to compute (default 20).
#' @param thealpha Nominal significance level (default 0.05).
#'
#' @return Named numeric vector of adjusted alpha levels, one per depth
#'   (1 through \code{max_depth}). Names are depth levels as characters.
#'
#' @details
#' The formula at level \eqn{\ell} is:
#' \deqn{\alpha_\ell^{adj} = \alpha / (k^{(\ell-1)} \cdot \prod_{j=1}^{\ell-1} \hat\theta_j)}
#' when cumulative power exceeds \code{tau}, and \eqn{\alpha_\ell^{adj} = \alpha}
#' otherwise.
#'
#' The FWER guarantee (Theorem in Appendix D) requires that power
#' is not underestimated (i.e., \eqn{\hat\theta_j \geq \theta_j}).
#' In practice this means using a conservatively large \code{delta_hat}.
#'
#' @examples
#' # Alpha schedule for a 4-ary tree with moderate effect
#' compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000)
#'
#' # Binary tree with very high power — approaches Bonferroni
#' compute_adaptive_alphas(k = 2, delta_hat = 0.5, N_total = 1e6,
#'                         tau = 0, max_depth = 5)
#'
#' @importFrom stats pnorm qnorm
#' @export
compute_adaptive_alphas <- function(k, delta_hat, N_total, tau = 0.1,
                                     max_depth = 20L, thealpha = 0.05) {
  stopifnot(length(delta_hat) == 1, delta_hat > 0,
            length(N_total) == 1, N_total > 0,
            length(tau) == 1, tau >= 0, tau <= 1,
            length(thealpha) == 1, thealpha > 0, thealpha < 1)

  # Handle scalar or vector k
  if (length(k) == 1L) {
    stopifnot(k >= 2)
    k_vec <- rep(as.integer(k), max_depth)
  } else {
    stopifnot(length(k) >= max_depth, all(k >= 1))
    k_vec <- as.integer(k[seq_len(max_depth)])
  }

  z_crit <- qnorm(1 - thealpha / 2)
  alphas <- numeric(max_depth)

  # Running products tracking the path from root to current level:
  # path_power = product of theta_j for j = 1..(ell-1)
  # cum_n_divisor = product of k_j for j = 1..(ell-1), used for sample size
  # num_tests = product of k_j for j = 1..(ell-1), the multiplicity at level ell
  # All doubles to avoid integer overflow at large depth
  path_power <- 1.0
  cum_n_divisor <- 1.0
  num_tests <- 1.0

  for (ell in seq_len(max_depth)) {
    # Sample size at this level under equal splitting
    n_ell <- N_total / cum_n_divisor
    # Estimated power at this level
    theta_ell <- pnorm(delta_hat * sqrt(n_ell) - z_crit)
    theta_ell <- max(min(theta_ell, 1), 0)

    if (ell == 1L) {
      # Root: single test, no adjustment
      alphas[ell] <- thealpha
    } else {
      cum_power_through_ell <- path_power * theta_ell

      if (cum_power_through_ell > tau) {
        # Adjust for multiplicity and cumulative power.
        # The min() guard ensures we never exceed nominal alpha.
        alphas[ell] <- min(thealpha, thealpha / (num_tests * path_power))
      } else {
        # Natural gating suffices: low cumulative power means reaching
        # this node by chance is already unlikely enough
        alphas[ell] <- thealpha
      }
    }

    # Update running products for next level
    path_power <- path_power * theta_ell
    cum_n_divisor <- cum_n_divisor * k_vec[ell]
    if (ell < max_depth) {
      num_tests <- num_tests * k_vec[ell]
    }
  }

  names(alphas) <- as.character(seq_len(max_depth))
  return(alphas)
}


#' Adaptive Alpha Adjustment Based on Power Decay
#'
#' Factory function that creates an alpha adjustment function for use
#' with \code{\link{find_blocks}}. The returned function adjusts
#' significance levels at each tree depth based on estimated power
#' decay (Algorithm 1 from the adaptive alpha appendix).
#'
#' @inheritParams compute_adaptive_alphas
#'
#' @return A function with signature
#'   \code{function(pval, batch, nodesize, thealpha, thew0, depth)}
#'   conforming to the \code{alphafn} interface used by
#'   \code{\link{find_blocks}}.
#'
#' @details
#' The returned function uses the \code{depth} parameter (passed by
#' \code{find_blocks}) to look up the pre-computed alpha for each
#' node's tree depth. The \code{pval}, \code{batch}, \code{nodesize},
#' and \code{thew0} parameters are accepted for interface compatibility
#' but are not used — unlike online FDR methods, the adaptive alpha
#' depends only on tree structure, not on observed p-values.
#'
#' Results are cached internally: the vector of adjusted alphas is
#' computed once per unique value of \code{thealpha} and reused on
#' subsequent calls.
#'
#' @examples
#' # Create an adaptive alpha function for a 4-ary tree
#' my_alpha <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)
#'
#' # Use with find_blocks
#' # find_blocks(idat, bdat, ..., alphafn = my_alpha)
#'
#' # Inspect the alpha schedule it will use
#' compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000)
#'
#' @export
alpha_adaptive <- function(k, delta_hat, N_total, tau = 0.1,
                           max_depth = 20L) {
  # Validate at factory time so errors surface early, not deep
  # inside a find_blocks run
  stopifnot(delta_hat > 0, N_total > 0, tau >= 0, tau <= 1)
  if (length(k) == 1L) {
    stopifnot(k >= 2)
  } else {
    stopifnot(all(k >= 1))
  }

  # Cache for computed alphas, keyed by thealpha
  cache_env <- new.env(parent = emptyenv())
  cache_env$alphas <- NULL
  cache_env$thealpha <- NULL

  # Return closure conforming to the alphafn interface
  function(pval, batch, nodesize, thealpha = 0.05, thew0 = 0.05 - 0.001,
           depth = NULL) {
    # Lazy compute and cache the alpha schedule
    if (is.null(cache_env$alphas) || !identical(thealpha, cache_env$thealpha)) {
      cache_env$alphas <- compute_adaptive_alphas(
        k = k, delta_hat = delta_hat, N_total = N_total,
        tau = tau, max_depth = max_depth, thealpha = thealpha
      )
      cache_env$thealpha <- thealpha
    }
    alpha_by_depth <- cache_env$alphas

    if (is.null(depth)) {
      warning("alpha_adaptive requires depth; returning nominal alpha")
      return(rep(thealpha, length(pval)))
    }

    # Look up pre-computed alpha for each node's depth
    depths <- as.integer(depth)
    # Clamp to valid range
    depths <- pmax(1L, pmin(depths, length(alpha_by_depth)))

    return(alpha_by_depth[depths])
  }
}
