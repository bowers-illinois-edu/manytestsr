# Adaptive alpha adjustment for tree-structured hypothesis testing
#
# Implements the error load framework from Appendix B of the supplement.
# The key idea: at each tree level, the "error load" G_ell counts the
# expected number of all-null sibling groups that the procedure tests.
# When the total error load sum(G_ell) <= 1, the unadjusted procedure
# controls FWER at level alpha ("natural gating"). When it exceeds 1,
# we adjust alpha at each level to compensate (Algorithm 1).


#' Compute Error Load for Natural Gating Assessment
#'
#' Computes the error load at each tree level, which measures the expected
#' number of all-null sibling groups tested. When the total error load is
#' at most 1, the unadjusted procedure (fixed alpha at every node) controls
#' FWER at level alpha --- no adjustment needed. When it exceeds 1, adaptive
#' alpha adjustment is required.
#'
#' Two interfaces are provided. The \strong{parametric} interface
#' (\code{k}, \code{delta_hat}, \code{N_total}) assumes a regular k-ary
#' tree with equal splitting: sample size at level ell is
#' \eqn{N / k^{ell-1}}. The \strong{tree} interface
#' (\code{node_dat}) accepts per-node sample sizes from an actual
#' (possibly irregular) tree, as returned by \code{\link{find_blocks}}.
#'
#' @param k Branching factor. Either a scalar (constant k at all levels)
#'   or an integer vector of length \code{max_depth} where \code{k[ell]}
#'   is the branching factor at level \code{ell}. Used in parametric mode;
#'   ignored when \code{node_dat} is provided.
#' @param delta_hat Estimated standardized effect size (e.g., Cohen's d).
#'   Used to compute power at each level via the normal approximation.
#' @param N_total Total sample size at the root level. Used in parametric
#'   mode; ignored when \code{node_dat} is provided.
#' @param node_dat Optional data.frame or data.table with columns
#'   \code{nodenum}, \code{parent}, \code{depth}, and \code{nodesize}.
#'   When provided, the function computes per-node power from
#'   \code{nodesize} and aggregates error load by depth. This supports
#'   irregular trees (e.g., DPP design with unequal splits).
#' @param max_depth Maximum tree depth to compute. In parametric mode,
#'   defaults to 20. In tree mode, inferred from \code{node_dat}.
#' @param thealpha Nominal significance level (default 0.05).
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{G}}{Numeric vector of error load at each level
#'     (length \code{max_depth}).}
#'   \item{\code{sum_G}}{Total error load: \code{sum(G)}.}
#'   \item{\code{needs_adjustment}}{Logical: \code{TRUE} when
#'     \code{sum_G > 1}.}
#'   \item{\code{thetas}}{Estimated power (conditional rejection
#'     probability) at each level.}
#'   \item{\code{critical_level}}{First level where theta < 1/k
#'     (where natural gating begins to dominate), or \code{NA} if
#'     theta never drops below 1/k.}
#'   \item{\code{n_by_level}}{Sample size at each level.}
#' }
#'
#' @details
#' In \strong{parametric mode}, the error load at level ell is:
#' \deqn{G_\ell = k^{\ell-1} \prod_{j=0}^{\ell-1} \theta_j}
#' where \eqn{\theta_j = \Phi(\delta \sqrt{n_j} - z_{1-\alpha/2})}
#' and \eqn{n_j = N / k^{j}}.
#'
#' In \strong{tree mode}, the function computes per-node power from
#' the actual sample sizes. For each node at depth d, it computes
#' \eqn{\theta = \Phi(\delta \sqrt{n_{\text{node}}} - z_{1-\alpha/2})}.
#' The error load at depth d is the sum over all nodes at depth d of the
#' product of theta values along the path from root to that node's parent,
#' which generalizes the regular-tree formula to irregular branching.
#'
#' @examples
#' # Parametric mode: regular 3-ary tree, moderate effect
#' compute_error_load(k = 3, delta_hat = 0.3, N_total = 500)
#'
#' # High power, wide tree: needs adjustment
#' res <- compute_error_load(k = 10, delta_hat = 0.5, N_total = 5000,
#'                           max_depth = 3)
#' res$sum_G     # likely > 1
#' res$needs_adjustment
#'
#' @importFrom stats pnorm qnorm
#' @export
compute_error_load <- function(k = NULL, delta_hat, N_total = NULL,
                               node_dat = NULL,
                               max_depth = 20L, thealpha = 0.05) {

  stopifnot(length(delta_hat) == 1, delta_hat > 0,
            length(thealpha) == 1, thealpha > 0, thealpha < 1)

  z_crit <- qnorm(1 - thealpha / 2)

  if (!is.null(node_dat)) {
    # --- Tree mode: use actual per-node sample sizes ---
    return(.error_load_from_tree(node_dat, delta_hat, z_crit, thealpha))
  }

  # --- Parametric mode: regular k-ary tree with equal splits ---
  stopifnot(!is.null(k), !is.null(N_total),
            length(N_total) == 1, N_total > 0)

  # Handle scalar or vector k
  if (length(k) == 1L) {
    stopifnot(k >= 2)
    k_vec <- rep(as.integer(k), max_depth)
  } else {
    stopifnot(length(k) >= max_depth, all(k >= 1))
    k_vec <- as.integer(k[seq_len(max_depth)])
  }

  G <- numeric(max_depth)
  thetas <- numeric(max_depth)
  n_by_level <- numeric(max_depth)
  cum_n_divisor <- 1.0
  # cum_theta tracks prod_{j=0}^{ell-1} theta_j (product through current level)
  cum_theta <- 1.0
  # cum_k tracks k^{ell-1} = prod of k_vec[1..(ell-1)]
  cum_k <- 1.0
  critical_level <- NA_integer_

  for (ell in seq_len(max_depth)) {
    n_ell <- N_total / cum_n_divisor
    n_by_level[ell] <- n_ell
    theta_ell <- pnorm(delta_hat * sqrt(n_ell) - z_crit)
    theta_ell <- max(min(theta_ell, 1), 0)
    thetas[ell] <- theta_ell

    # G_ell = k^{ell-1} * prod_{j=0}^{ell-1} theta_j
    # In 1-indexed terms: cum_k * cum_theta * theta_ell
    # (cum_k = k^{ell-2} before update, cum_theta = prod thetas through ell-1)
    # Actually: we update cum_k and cum_theta AFTER using them, so:
    # G[ell] = cum_k * cum_theta * theta_ell
    G[ell] <- cum_k * cum_theta * theta_ell

    # Track critical level: theta < 1/k means error load shrinks
    if (is.na(critical_level) && ell > 1L && theta_ell < 1.0 / k_vec[ell - 1L]) {
      critical_level <- ell
    }

    # Update running products for next level
    cum_theta <- cum_theta * theta_ell
    cum_n_divisor <- cum_n_divisor * k_vec[ell]
    cum_k <- cum_k * k_vec[ell]
  }

  names(G) <- names(thetas) <- names(n_by_level) <- as.character(seq_len(max_depth))

  return(list(
    G = G,
    sum_G = sum(G),
    needs_adjustment = sum(G) > 1,
    thetas = thetas,
    critical_level = critical_level,
    n_by_level = n_by_level
  ))
}


#' Compute error load from an actual tree structure (internal)
#'
#' @param node_dat data.frame with nodenum, parent, depth, nodesize
#' @param delta_hat Standardized effect size
#' @param z_crit Critical z-value
#' @param thealpha Nominal alpha
#' @return Same structure as compute_error_load
#' @keywords internal
.error_load_from_tree <- function(node_dat, delta_hat, z_crit, thealpha) {

  nd <- as.data.frame(node_dat)
  required_cols <- c("nodenum", "parent", "depth", "nodesize")
  missing_cols <- setdiff(required_cols, names(nd))
  if (length(missing_cols) > 0L) {
    stop("node_dat must have columns: ", paste(missing_cols, collapse = ", "))
  }

  # Compute per-node power from nodesize
  nd$theta <- pnorm(delta_hat * sqrt(nd$nodesize) - z_crit)
  nd$theta <- pmax(pmin(nd$theta, 1), 0)

  # Compute path_power for each node: product of theta values along
  # the path from root to this node's parent.
  # For the root, path_power = 1 (no ancestors).
  # For a child at depth 2, path_power = theta_root.
  # We traverse by depth order.

  nd$path_power <- NA_real_
  depths <- sort(unique(nd$depth))

  # Index by nodenum for fast lookup
  node_idx <- stats::setNames(seq_len(nrow(nd)), nd$nodenum)

  for (d in depths) {
    rows_at_d <- which(nd$depth == d)
    for (r in rows_at_d) {
      if (nd$parent[r] == 0) {
        # Root node: no ancestors
        nd$path_power[r] <- 1.0
      } else {
        parent_row <- node_idx[as.character(nd$parent[r])]
        nd$path_power[r] <- nd$path_power[parent_row] * nd$theta[parent_row]
      }
    }
  }

  # Error load at depth d: sum over all nodes at depth d of
  # path_power * theta (ancestors all rejected AND this node rejects).
  # This matches G_ell = k^{ell-1} * prod_{j=0}^{ell-1} theta_j
  # from the supplement definition.
  max_depth <- max(nd$depth)
  G <- numeric(max_depth)
  thetas_by_level <- numeric(max_depth)
  n_by_level <- numeric(max_depth)

  for (d in seq_len(max_depth)) {
    rows_at_d <- which(nd$depth == d)
    G[d] <- sum(nd$path_power[rows_at_d] * nd$theta[rows_at_d])
    # Report average theta and n at this level for diagnostics
    thetas_by_level[d] <- mean(nd$theta[rows_at_d])
    n_by_level[d] <- mean(nd$nodesize[rows_at_d])
  }

  # Critical level: first depth where average theta < 1/(avg branching factor)
  # For irregular trees, use the ratio G_{d+1}/G_d as the effective k*theta
  critical_level <- NA_integer_
  for (d in seq_len(max_depth - 1L)) {
    if (G[d] > 0 && G[d + 1L] / G[d] < 1) {
      critical_level <- d + 1L
      break
    }
  }

  names(G) <- names(thetas_by_level) <- names(n_by_level) <- as.character(seq_len(max_depth))

  return(list(
    G = G,
    sum_G = sum(G),
    needs_adjustment = sum(G) > 1,
    thetas = thetas_by_level,
    critical_level = critical_level,
    n_by_level = n_by_level,
    node_detail = nd[, c("nodenum", "parent", "depth", "nodesize", "theta", "path_power")]
  ))
}


#' Compute Adaptive Alpha Levels by Tree Depth
#'
#' Implements Algorithm 1 from Appendix B of the supplement. First checks
#' whether natural gating suffices (total error load \eqn{\le 1}); if so,
#' returns nominal alpha at every level. Otherwise, computes adjusted
#' significance levels that compensate for the error load at each depth.
#'
#' @param k Branching factor. Either a scalar (constant k at all levels)
#'   or an integer vector of length \code{max_depth} where \code{k[ell]}
#'   is the branching factor at level \code{ell}.
#' @param delta_hat Estimated standardized effect size (e.g., Cohen's d).
#'   Conservative (larger) values produce more stringent adjustment,
#'   which preserves the FWER guarantee. Use an upper bound on the
#'   true effect size.
#' @param N_total Total sample size at the root level.
#' @param max_depth Maximum depth to compute (default 20).
#' @param thealpha Nominal significance level (default 0.05).
#'
#' @return Named numeric vector of adjusted alpha levels, one per depth
#'   (1 through \code{max_depth}). Names are depth levels as characters.
#'   Has attribute \code{"error_load"} containing the
#'   \code{\link{compute_error_load}} result, so the caller can inspect
#'   whether adjustment was needed.
#'
#' @details
#' The function first calls \code{\link{compute_error_load}} to assess
#' whether natural gating suffices. When \eqn{\sum G_\ell \le 1}, no
#' adjustment is needed and nominal \code{thealpha} is returned at every
#' level.
#'
#' When adjustment is needed, the formula at level \eqn{\ell} is:
#' \deqn{\alpha_\ell^{adj} = \min\left\{\alpha,\;
#'   \frac{\alpha}{k^{(\ell-1)} \cdot \prod_{j=1}^{\ell-1} \hat\theta_j}
#'   \right\}}
#'
#' The FWER guarantee (Theorem in the supplement) requires that power
#' is not underestimated (i.e., \eqn{\hat\theta_j \geq \theta_j}).
#' In practice this means using a conservatively large \code{delta_hat}.
#'
#' @examples
#' # Natural gating sufficient: all alphas = 0.05
#' compute_adaptive_alphas(k = 3, delta_hat = 0.2, N_total = 100,
#'                         max_depth = 4)
#'
#' # Needs adjustment: alphas shrink at deeper levels
#' compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000,
#'                         max_depth = 5)
#'
#' @importFrom stats pnorm qnorm
#' @export
compute_adaptive_alphas <- function(k, delta_hat, N_total,
                                     max_depth = 20L, thealpha = 0.05) {
  stopifnot(length(delta_hat) == 1, delta_hat > 0,
            length(N_total) == 1, N_total > 0,
            length(thealpha) == 1, thealpha > 0, thealpha < 1)

  # Check error load first: if natural gating suffices, no adjustment needed
  el <- compute_error_load(k = k, delta_hat = delta_hat, N_total = N_total,
                           max_depth = max_depth, thealpha = thealpha)

  if (!el$needs_adjustment) {
    alphas <- rep(thealpha, max_depth)
    names(alphas) <- as.character(seq_len(max_depth))
    attr(alphas, "error_load") <- el
    return(alphas)
  }

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
      # Adjust for multiplicity and cumulative power.
      # The min() guard ensures we never exceed nominal alpha.
      alphas[ell] <- min(thealpha, thealpha / (num_tests * path_power))
    }

    # Update running products for next level
    path_power <- path_power * theta_ell
    cum_n_divisor <- cum_n_divisor * k_vec[ell]
    if (ell < max_depth) {
      num_tests <- num_tests * k_vec[ell]
    }
  }

  names(alphas) <- as.character(seq_len(max_depth))
  attr(alphas, "error_load") <- el
  return(alphas)
}


#' Compute Adaptive Alpha Levels from an Actual Tree
#'
#' Tree-mode counterpart to \code{\link{compute_adaptive_alphas}}. Instead
#' of assuming a regular k-ary tree with equal splits, this function takes
#' the actual tree structure (with per-node sample sizes) and computes
#' per-depth adjusted significance levels. This handles irregular trees
#' where branching factor and sample sizes vary across nodes.
#'
#' @param node_dat A data.frame or data.table with columns
#'   \code{nodenum}, \code{parent}, \code{depth}, and \code{nodesize}.
#'   Typically extracted from a \code{\link{find_blocks}} result.
#'   The root node must have \code{parent = 0} and \code{depth = 1}.
#' @param delta_hat Estimated standardized effect size (e.g., Cohen's d).
#'   Conservative (larger) values produce more stringent adjustment,
#'   which preserves the FWER guarantee.
#' @param max_depth Maximum depth to compute. Defaults to the maximum
#'   depth present in \code{node_dat}.
#' @param thealpha Nominal significance level (default 0.05).
#'
#' @return Named numeric vector of adjusted alpha levels, one per depth
#'   (1 through \code{max_depth}). Names are depth levels as characters.
#'   Has attribute \code{"error_load"} containing the
#'   \code{\link{compute_error_load}} result from the tree.
#'
#' @details
#' The algorithm mirrors \code{\link{compute_adaptive_alphas}} but uses
#' actual per-node power instead of the parametric assumption:
#'
#' \enumerate{
#'   \item Compute per-node power \eqn{\theta} and path power
#'     (product of ancestor thetas) via \code{compute_error_load}.
#'   \item If total error load \eqn{\sum G_\ell \le 1}: natural gating
#'     suffices, return nominal \code{thealpha} at every depth.
#'   \item Otherwise, at depth \eqn{d}:
#'     \deqn{\alpha_d = \min\left\{\alpha,\;
#'       \frac{\alpha}{\sum_{\text{nodes at depth } d}
#'         \text{path\_power}(\text{node})} \right\}}
#' }
#'
#' The denominator at depth \eqn{d} is the sum of path powers over all
#' nodes at that depth --- i.e., the expected number of tests the
#' procedure conducts at depth \eqn{d}. For a regular k-ary tree with
#' equal sample sizes, this reduces to
#' \eqn{k^{d-1} \prod_{j=1}^{d-1} \theta_j}, matching the parametric
#' formula in \code{\link{compute_adaptive_alphas}}.
#'
#' @examples
#' # Build a small irregular tree
#' nd <- data.frame(
#'   nodenum  = 1:7,
#'   parent   = c(0, 1, 1, 2, 2, 3, 3),
#'   depth    = c(1, 2, 2, 3, 3, 3, 3),
#'   nodesize = c(500, 250, 250, 125, 125, 100, 150)
#' )
#' compute_adaptive_alphas_tree(node_dat = nd, delta_hat = 0.5)
#'
#' @param budget_weights Controls how the error budget is allocated across
#'   depths (Theorems B.3 and B.5 in the supplement). Options:
#'   \describe{
#'     \item{\code{NULL}}{(default) No budget weights --- use the
#'       telescoping formula \eqn{\alpha_\ell = \alpha / G_\ell}. This
#'       is the tightest bound for regular k-ary trees (Theorem B.2).}
#'     \item{\code{"equal"}}{Equal weights: \eqn{w_\ell = 1/(L-1)} for
#'       all depths \eqn{\ell \ge 2}.}
#'     \item{\code{"proportional"}}{Error-load-proportional:
#'       \eqn{w_\ell = G_\ell / \sum G}, giving a uniform adjusted
#'       alpha of \eqn{\alpha / \sum G} at every depth.}
#'     \item{numeric vector}{Custom weights for depths 2 through L.
#'       Must have length \eqn{L - 1} and sum to at most 1.}
#'   }
#'
#' @importFrom stats pnorm qnorm
#' @export
compute_adaptive_alphas_tree <- function(node_dat, delta_hat,
                                         max_depth = NULL,
                                         thealpha = 0.05,
                                         budget_weights = NULL) {
  stopifnot(length(delta_hat) == 1, delta_hat > 0,
            length(thealpha) == 1, thealpha > 0, thealpha < 1)

  z_crit <- qnorm(1 - thealpha / 2)

  # Compute per-node power and path_power via the existing tree helper
  el <- .error_load_from_tree(node_dat, delta_hat, z_crit, thealpha)

  if (is.null(max_depth)) {
    max_depth <- max(el$node_detail$depth)
  }

  # If natural gating suffices, return nominal alpha everywhere
  # (budget weights are irrelevant when the tree protects itself)
  if (!el$needs_adjustment) {
    alphas <- rep(thealpha, max_depth)
    names(alphas) <- as.character(seq_len(max_depth))
    attr(alphas, "error_load") <- el
    return(alphas)
  }

  # Compute per-depth error loads G_ell for weight resolution
  nd <- el$node_detail
  G_by_depth <- numeric(max_depth)
  for (d in seq_len(max_depth)) {
    rows_at_d <- which(nd$depth == d)
    G_by_depth[d] <- sum(nd$path_power[rows_at_d])
  }

  # Resolve budget weights into a numeric vector (NULL = telescoping)
  weights <- .resolve_budget_weights(budget_weights, G_by_depth, max_depth)

  # Compute per-depth alpha.
  # With weights: alpha_ell = min(alpha, w_ell * alpha / G_ell)
  # Without (NULL): alpha_ell = min(alpha, alpha / G_ell) --- telescoping
  alphas <- numeric(max_depth)
  for (d in seq_len(max_depth)) {
    if (d == 1L || G_by_depth[d] == 0) {
      # Root always gets nominal alpha; zero error load needs no shrinkage
      alphas[d] <- thealpha
    } else if (is.null(weights)) {
      alphas[d] <- min(thealpha, thealpha / G_by_depth[d])
    } else {
      alphas[d] <- min(thealpha, weights[d] * thealpha / G_by_depth[d])
    }
  }

  names(alphas) <- as.character(seq_len(max_depth))
  attr(alphas, "error_load") <- el
  return(alphas)
}


# Resolve user-facing budget_weights argument into a numeric vector.
# Returns NULL when budget_weights is NULL (caller uses telescoping).
# Otherwise returns a numeric vector of length max_depth where
# weights[1] = 0 (root gets nominal alpha) and sum(weights[2:L]) <= 1.
.resolve_budget_weights <- function(budget_weights, G_by_depth, max_depth) {
  if (is.null(budget_weights)) {
    return(NULL)
  }

  L <- max_depth
  n_depths <- L - 1L  # depths 2 through L

  if (is.character(budget_weights)) {
    budget_weights <- match.arg(budget_weights, c("equal", "proportional"))

    if (budget_weights == "equal") {
      w <- rep(1 / n_depths, n_depths)
    } else {
      # Proportional: w_ell = G_ell / sum(G), so alpha_ell = alpha / sum(G)
      G_sub <- G_by_depth[2:L]
      total_G <- sum(G_sub)
      if (total_G <= 0) {
        w <- rep(1 / n_depths, n_depths)
      } else {
        w <- G_sub / total_G
      }
    }
    return(c(0, w))
  }

  if (is.numeric(budget_weights)) {
    if (length(budget_weights) != n_depths) {
      stop("budget_weights must have length ", n_depths,
           " (one weight per depth from 2 to ", L, ")")
    }
    if (any(budget_weights < 0)) {
      stop("budget_weights must be non-negative")
    }
    if (sum(budget_weights) > 1 + 1e-10) {
      stop("budget_weights must sum to at most 1 (got ",
           round(sum(budget_weights), 4), ")")
    }
    return(c(0, budget_weights))
  }

  stop("budget_weights must be NULL, 'equal', 'proportional', or a numeric vector")
}


#' Adaptive Alpha Adjustment Based on Power Decay
#' with \code{\link{find_blocks}}. The returned function adjusts
#' significance levels at each tree depth based on estimated power
#' decay (Algorithm 1 from Appendix B of the supplement).
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
#' subsequent calls. When the error load is at most 1 (natural gating
#' suffices), nominal alpha is returned at every level.
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
alpha_adaptive <- function(k, delta_hat, N_total, max_depth = 20L) {
  # Validate at factory time so errors surface early, not deep
  # inside a find_blocks run
  stopifnot(delta_hat > 0, N_total > 0)
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
        max_depth = max_depth, thealpha = thealpha
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


#' Adaptive Alpha from an Actual Tree Structure
#'
#' Factory function counterpart to \code{\link{alpha_adaptive}} for
#' irregular trees. Creates an alpha adjustment function for use with
#' \code{\link{find_blocks}}, using actual per-node sample sizes to
#' compute the alpha schedule instead of assuming a regular k-ary tree.
#'
#' @inheritParams compute_adaptive_alphas_tree
#'
#' @return A function with signature
#'   \code{function(pval, batch, nodesize, thealpha, thew0, depth)}
#'   conforming to the \code{alphafn} interface used by
#'   \code{\link{find_blocks}}.
#'
#' @details
#' The returned function behaves identically to the closure from
#' \code{\link{alpha_adaptive}}: it uses the \code{depth} parameter
#' to look up pre-computed alphas, ignoring \code{pval}, \code{batch},
#' and other arguments. The difference is that the alpha schedule comes
#' from \code{\link{compute_adaptive_alphas_tree}} rather than the
#' parametric formula, so it handles irregular branching correctly.
#'
#' Results are cached internally: the alpha schedule is computed once
#' per unique value of \code{thealpha} and reused on subsequent calls.
#'
#' @examples
#' # Build a tree template from the DPP design
#' nd <- data.frame(
#'   nodenum  = 1:7,
#'   parent   = c(0, 1, 1, 2, 2, 3, 3),
#'   depth    = c(1, 2, 2, 3, 3, 3, 3),
#'   nodesize = c(500, 250, 250, 125, 125, 100, 150)
#' )
#' my_alpha <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)
#'
#' # Use with find_blocks
#' # find_blocks(idat, bdat, ..., alphafn = my_alpha)
#'
#' @inheritParams compute_adaptive_alphas_tree
#' @export
alpha_adaptive_tree <- function(node_dat, delta_hat, max_depth = NULL,
                                budget_weights = NULL) {
  # Validate at factory time
  stopifnot(delta_hat > 0)
  nd <- as.data.frame(node_dat)
  required_cols <- c("nodenum", "parent", "depth", "nodesize")
  missing_cols <- setdiff(required_cols, names(nd))
  if (length(missing_cols) > 0L) {
    stop("node_dat must have columns: ", paste(missing_cols, collapse = ", "))
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
      cache_env$alphas <- compute_adaptive_alphas_tree(
        node_dat = nd, delta_hat = delta_hat,
        max_depth = max_depth, thealpha = thealpha,
        budget_weights = budget_weights
      )
      cache_env$thealpha <- thealpha
    }
    alpha_by_depth <- cache_env$alphas

    if (is.null(depth)) {
      warning("alpha_adaptive_tree requires depth; returning nominal alpha")
      return(rep(thealpha, length(pval)))
    }

    # Look up pre-computed alpha for each node's depth
    depths <- as.integer(depth)
    # Clamp to valid range
    depths <- pmax(1L, pmin(depths, length(alpha_by_depth)))

    return(alpha_by_depth[depths])
  }
}


#' Find All Descendants of Given Nodes in a Tree
#'
#' BFS traversal using the parent column to collect every node
#' descended from the given starting nodes.
#'
#' @param node_dat data.frame with at least \code{nodenum} and
#'   \code{parent} columns.
#' @param nodenums Integer vector of node IDs whose descendants
#'   to find.
#' @return Integer vector of descendant node IDs (not including
#'   the starting nodes themselves). Empty integer vector if no
#'   descendants exist.
#' @keywords internal
.get_all_descendants <- function(node_dat, nodenums) {
  nd <- as.data.frame(node_dat)
  # BFS: start with children of the given nodes, then expand
  descendants <- integer(0)
  frontier <- nd$nodenum[nd$parent %in% nodenums]

  while (length(frontier) > 0L) {
    descendants <- c(descendants, frontier)
    frontier <- nd$nodenum[nd$parent %in% frontier]
  }

  return(descendants)
}


#' Adaptive Alpha with Branch Pruning
#'
#' Factory function that creates an adaptive alpha adjustment system
#' supporting branch pruning. Unlike \code{\link{alpha_adaptive_tree}},
#' which pre-computes a fixed schedule, this version can recompute the
#' schedule on a pruned subtree after each depth --- giving more alpha
#' to surviving branches when dead branches are removed.
#'
#' @inheritParams alpha_adaptive_tree
#' @param budget_weights Controls depth-wise budget allocation. Accepts
#'   the same values as \code{\link{compute_adaptive_alphas_tree}}, plus
#'   \code{"remaining"}: a sequential spending process where a fraction
#'   \code{spending_fraction} of the remaining budget is spent at each
#'   depth. See Details.
#' @param budget_total Initial error budget (default 1.0). The constraint
#'   \eqn{\sum w_\ell \le} \code{budget_total} guarantees FWER control.
#' @param spending_fraction Fraction of remaining budget to spend at each
#'   depth when \code{budget_weights = "remaining"} (default 0.5). At
#'   depth \eqn{\ell}, the weight is
#'   \eqn{w_\ell = f \times B_\ell} where \eqn{f} is the spending
#'   fraction and \eqn{B_\ell} is the remaining budget.
#' @param switching Logical (default \code{FALSE}). When \code{TRUE},
#'   implements the switching corollary (Corollary B.1): after each
#'   update, if the remaining pruned error load fits within the
#'   remaining budget, all deeper depths revert to nominal alpha.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{alphafn}}{A function with the standard alphafn
#'     interface: \code{function(pval, batch, nodesize, thealpha,
#'     thew0, depth)}. Looks up the current alpha schedule by depth.}
#'   \item{\code{update}}{A function
#'     \code{function(pruned_node_dat, thealpha)} that recomputes
#'     the alpha schedule on the given (pruned) tree. Called by
#'     \code{find_blocks} after each depth's testable decisions.}
#'   \item{\code{reset}}{A function \code{function(thealpha)} that
#'     restores the alpha schedule to the full (unpruned) tree.
#'     Called by \code{find_blocks} at the start of each run to
#'     ensure independence across simulation iterations.}
#' }
#'
#' @details
#' The FWER guarantee follows from Theorem B.5 in the supplement:
#' predictable budget weights with data-dependent denominators. The
#' weights are "predictable" because \eqn{w_\ell} depends only on the
#' testing history through depth \eqn{\ell - 1}, not on depth-\eqn{\ell}
#' outcomes. The union bound across depths gives FWER \eqn{\le \alpha}
#' whenever \eqn{\sum w_\ell \le 1}.
#'
#' The \strong{switching corollary} (when \code{switching = TRUE}):
#' after pruning narrows the surviving tree, if the remaining error load
#' \eqn{\sum_{\ell \ge s} D_\ell \le B_s} (remaining budget), then
#' \eqn{\alpha_\ell = \alpha} for all \eqn{\ell \ge s}. This works
#' by setting \eqn{w_\ell = D_\ell}, so the \eqn{D_\ell} in numerator
#' and denominator cancel, leaving nominal alpha.
#'
#' When \code{find_blocks} detects a list-valued \code{alphafn},
#' it extracts these three components and calls \code{reset} at
#' the start of each run and \code{update} after each depth.
#'
#' @examples
#' nd <- data.frame(
#'   nodenum  = 1:7,
#'   parent   = c(0, 1, 1, 2, 2, 3, 3),
#'   depth    = c(1, 2, 2, 3, 3, 3, 3),
#'   nodesize = c(500, 250, 250, 125, 125, 100, 150)
#' )
#' obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)
#' # obj$alphafn -- pass to find_blocks
#' # obj$update(pruned_nd, 0.05) -- recompute on surviving tree
#' # obj$reset(0.05) -- restore full-tree schedule
#'
#' @export
alpha_adaptive_tree_pruned <- function(node_dat, delta_hat,
                                       max_depth = NULL,
                                       budget_weights = NULL,
                                       budget_total = 1.0,
                                       spending_fraction = 0.5,
                                       switching = FALSE) {
  # Validate at factory time
  stopifnot(delta_hat > 0, budget_total > 0, budget_total <= 1,
            spending_fraction > 0, spending_fraction <= 1,
            is.logical(switching))
  full_nd <- as.data.frame(node_dat)
  required_cols <- c("nodenum", "parent", "depth", "nodesize")
  missing_cols <- setdiff(required_cols, names(full_nd))
  if (length(missing_cols) > 0L) {
    stop("node_dat must have columns: ", paste(missing_cols, collapse = ", "))
  }

  # Determine whether to use the remaining-budget sequential process
  use_remaining <- identical(budget_weights, "remaining")
  # For non-"remaining" budget_weights, pass through to compute function
  bw_for_compute <- if (use_remaining) NULL else budget_weights

  # Shared environment: alphafn reads, update/reset write.
  state <- new.env(parent = emptyenv())
  state$alphas <- NULL
  state$thealpha <- NULL
  state$budget_remaining <- budget_total

  # Compute schedule on a given tree and store in shared state.
  # When use_remaining = TRUE, the budget-weight for this depth is
  # spending_fraction * budget_remaining, and budget_remaining is
  # decremented. When switching = TRUE, check if the remaining error
  # load fits within the remaining budget.
  .compute_and_store <- function(nd, thealpha, is_update = FALSE) {
    z_crit <- stats::qnorm(1 - thealpha / 2)
    el <- .error_load_from_tree(nd, delta_hat, z_crit, thealpha)
    cur_max_depth <- if (is.null(max_depth)) max(el$node_detail$depth) else max_depth

    # --- Switching corollary check ---
    # If remaining error load fits within remaining budget,
    # all remaining depths get nominal alpha
    if (switching && is_update && !el$needs_adjustment) {
      # Pruned error load < 1 means natural gating on surviving tree
      alphas <- rep(thealpha, cur_max_depth)
      names(alphas) <- as.character(seq_len(cur_max_depth))
      attr(alphas, "error_load") <- el
      state$alphas <- alphas
      state$thealpha <- thealpha
      return(invisible(NULL))
    }

    if (switching && is_update && el$sum_G <= state$budget_remaining) {
      # Remaining error load fits in remaining budget: switch to nominal
      alphas <- rep(thealpha, cur_max_depth)
      names(alphas) <- as.character(seq_len(cur_max_depth))
      attr(alphas, "error_load") <- el
      state$alphas <- alphas
      state$thealpha <- thealpha
      return(invisible(NULL))
    }

    # --- Remaining-budget sequential spending ---
    if (use_remaining && is_update) {
      # Spend a fraction of the remaining budget at this update
      c_ell <- spending_fraction * state$budget_remaining
      state$budget_remaining <- state$budget_remaining - c_ell

      # Build a single-depth weight vector: c_ell for the current
      # shallowest untested depth, then distribute remaining budget
      # evenly across deeper depths as a conservative default
      nd_detail <- el$node_detail
      G_by_depth <- numeric(cur_max_depth)
      for (d in seq_len(cur_max_depth)) {
        rows_at_d <- which(nd_detail$depth == d)
        G_by_depth[d] <- sum(nd_detail$path_power[rows_at_d])
      }

      # Compute alphas using the current spending amount
      alphas <- numeric(cur_max_depth)
      for (d in seq_len(cur_max_depth)) {
        if (d == 1L || G_by_depth[d] == 0) {
          alphas[d] <- thealpha
        } else {
          # Use c_ell as the weight for all remaining depths
          alphas[d] <- min(thealpha, c_ell * thealpha / G_by_depth[d])
        }
      }
      names(alphas) <- as.character(seq_len(cur_max_depth))
      attr(alphas, "error_load") <- el
      state$alphas <- alphas
      state$thealpha <- thealpha
      return(invisible(NULL))
    }

    # --- Default: pass through to compute_adaptive_alphas_tree ---
    state$alphas <- compute_adaptive_alphas_tree(
      node_dat = nd, delta_hat = delta_hat,
      max_depth = max_depth, thealpha = thealpha,
      budget_weights = bw_for_compute
    )
    state$thealpha <- thealpha
  }

  # --- alphafn: standard closure interface ---
  alphafn <- function(pval, batch, nodesize, thealpha = 0.05,
                      thew0 = 0.05 - 0.001, depth = NULL) {
    # Lazy initialization on first call
    if (is.null(state$alphas) || !identical(thealpha, state$thealpha)) {
      .compute_and_store(full_nd, thealpha, is_update = FALSE)
    }
    alpha_by_depth <- state$alphas

    if (is.null(depth)) {
      warning("alpha_adaptive_tree_pruned requires depth; returning nominal alpha")
      return(rep(thealpha, length(pval)))
    }

    depths <- as.integer(depth)
    depths <- pmax(1L, pmin(depths, length(alpha_by_depth)))
    return(alpha_by_depth[depths])
  }

  # --- update: recompute schedule on pruned tree ---
  update_fn <- function(pruned_node_dat, thealpha = 0.05) {
    .compute_and_store(pruned_node_dat, thealpha, is_update = TRUE)
  }

  # --- reset: restore full-tree schedule and budget ---
  reset_fn <- function(thealpha = 0.05) {
    state$budget_remaining <- budget_total
    .compute_and_store(full_nd, thealpha, is_update = FALSE)
  }

  return(list(
    alphafn = alphafn,
    update = update_fn,
    reset = reset_fn
  ))
}
