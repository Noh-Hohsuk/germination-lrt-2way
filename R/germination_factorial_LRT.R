# germination_factorial_LRT_final.R
#
#
# Release v1
#   Unified two-factor germination analysis with selectable interaction scale:
#     * scale = "probability", method = "mle"
#         constrained binomial likelihood ratio tests on the probability scale
#     * scale = "logit", method = "mle"
#         ordinary binomial logistic regression with likelihood-ratio tests
#     * scale = "logit", method = "firth"
#         Firth logistic regression with penalized likelihood-ratio tests
#
#   Automatic post-hoc workflow:
#     * nonsignificant interaction -> pairwise comparisons for significant main effects
#     * significant interaction    -> simple-effect pairwise comparisons
#     * Holm multiplicity adjustment
#     * compact letter displays
#
#   Additional features:
#     * heuristic separation diagnostics for ordinary logistic MLE
#     * fitted probability matrices
#     * consistent overall-test table across analysis modes
#     * concise print(), detailed summary(), and extractor functions
#
# v4 refinements
#   * cleaner user-facing output
#   * very small p-values displayed in scientific notation rather than 0
#   * explicit distinction between primary logit-scale effect estimates and
#     descriptive probability differences in post-hoc output
#   * compact logistic-separation diagnostics
# ------------------------------------------------------------
# Two-factor likelihood-based analysis for germination data
# supporting probability-scale and logit-scale interaction tests.
#
# Main features
#   * scale = "probability": proposed probability-scale Binomial LRT
#       - constrained additive model:
#             p_ij = mu + alpha_i + beta_j
#             0 <= p_ij <= 1
#       - nloptr::NLOPT_LD_SLSQP with multiple starting values
#       - automatic post-hoc testing, Holm adjustment, and CLD
#   * scale = "logit", method = "mle": ordinary binomial logistic regression
#       - additive null on the logit scale
#       - likelihood-ratio tests for Factor A, Factor B, and interaction
#       - automatic logit-scale post-hoc testing, Holm adjustment, and CLD
#       - heuristic diagnostics for possible separation
#   * scale = "logit", method = "firth": Firth logistic regression
#       - penalized likelihood-ratio tests for the same logit-scale hypotheses
#       - automatic logit-scale post-hoc testing, Holm adjustment, and CLD
#       - main-effect pairwise Firth comparisons use logistf's PLR comparison
#         because merging factor levels changes formula term labels
#       - useful when ordinary logistic fitting is affected by separation
#
# Required packages
#   probability scale: nloptr, multcompView
#   Firth logit:       logistf
#   install.packages(c("nloptr", "multcompView", "logistf"))
#
# Required data columns
#   Fac_A, Fac_B, g_ijk, n_ijk
#
# Factor levels must be coded 1,2,...,I and 1,2,...,J.
# ------------------------------------------------------------

if (!requireNamespace("nloptr", quietly = TRUE)) {
  stop("Package 'nloptr' is required. Install it with install.packages('nloptr').")
}

if (!requireNamespace("multcompView", quietly = TRUE)) {
  stop("Package 'multcompView' is required. Install it with install.packages('multcompView').")
}


# ============================================================
# Utilities
# ============================================================

.make_L <- function(I, J) {
  q <- 1 + (I - 1) + (J - 1)
  L <- matrix(0, nrow = I * J, ncol = q)

  r <- 1L
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      L[r, 1] <- 1
      if (i > 1) L[r, i] <- 1
      if (j > 1) L[r, I + j - 1] <- 1
      r <- r + 1L
    }
  }
  L
}


.prepare_data <- function(df) {
  req <- c("Fac_A", "Fac_B", "g_ijk", "n_ijk")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }

  df <- df[, req]
  df$Fac_A <- as.integer(df$Fac_A)
  df$Fac_B <- as.integer(df$Fac_B)
  df$g_ijk <- as.numeric(df$g_ijk)
  df$n_ijk <- as.numeric(df$n_ijk)

  if (anyNA(df)) stop("Missing values are not allowed.")
  if (any(df$g_ijk < 0) || any(df$n_ijk <= 0) ||
      any(df$g_ijk > df$n_ijk)) {
    stop("Invalid binomial counts.")
  }

  I <- max(df$Fac_A)
  J <- max(df$Fac_B)

  if (!setequal(sort(unique(df$Fac_A)), seq_len(I))) {
    stop("Fac_A must be coded consecutively as 1,...,I.")
  }
  if (!setequal(sort(unique(df$Fac_B)), seq_len(J))) {
    stop("Fac_B must be coded consecutively as 1,...,J.")
  }

  tab <- table(df$Fac_A, df$Fac_B)
  if (any(tab == 0)) {
    stop("Every A x B treatment combination must be represented.")
  }

  df <- df[order(df$Fac_A, df$Fac_B), ]
  rownames(df) <- NULL

  list(df = df, I = I, J = J)
}


.cell_totals <- function(df, I, J) {
  z <- aggregate(
    cbind(g_ijk, n_ijk) ~ Fac_A + Fac_B,
    data = df,
    FUN = sum
  )
  z <- z[order(z$Fac_A, z$Fac_B), ]
  if (nrow(z) != I * J) stop("Unexpected number of treatment cells.")
  z
}


.safe_xlog <- function(x, y) {
  ans <- numeric(length(x))
  use <- x > 0
  ans[use] <- x[use] * log(y[use])
  ans
}


.full_binom_loglik <- function(df, p_cell, J) {
  idx <- (df$Fac_A - 1L) * J + df$Fac_B
  sum(
    dbinom(
      df$g_ijk,
      size = df$n_ijk,
      prob = p_cell[idx],
      log = TRUE
    )
  )
}


.pair_matrix <- function(K) {
  if (K < 2) return(matrix(integer(0), ncol = 2))
  t(combn(seq_len(K), 2))
}


# ============================================================
# Constrained additive model M_{A+B}
# ============================================================

.fit_additive_binomial <- function(
    df,
    maxeval = 10000,
    xtol_rel = 1e-10,
    ftol_rel = 1e-12,
    feasibility_tol = 1e-8,
    gradient_eps = 1e-12,
    print_level = 0
) {
  prep <- .prepare_data(df)
  df <- prep$df
  I <- prep$I
  J <- prep$J

  L <- .make_L(I, J)
  ct <- .cell_totals(df, I, J)

  G <- as.numeric(ct$g_ijk)
  N <- as.numeric(ct$n_ijk)
  F <- N - G

  # Negative log-likelihood without combinatorial constants.
  nll <- function(theta) {
    p <- as.vector(L %*% theta)

    if (any(!is.finite(p)) || any(p < 0) || any(p > 1)) {
      return(1e100)
    }

    if (any(p == 0 & G > 0) || any(p == 1 & F > 0)) {
      return(1e100)
    }

    val <- -sum(.safe_xlog(G, p) + .safe_xlog(F, 1 - p))
    if (!is.finite(val)) 1e100 else val
  }

  # Analytic gradient.  Clipping is used only for stable gradient
  # evaluation near p = 0 or p = 1.
  grad_nll <- function(theta) {
    p <- as.vector(L %*% theta)
    pg <- pmin(pmax(p, gradient_eps), 1 - gradient_eps)

    score_p <- G / pg - F / (1 - pg)
    as.vector(-crossprod(L, score_p))
  }

  # nloptr uses constraints g(theta) <= 0.
  eval_g_ineq <- function(theta) {
    p <- as.vector(L %*% theta)
    c(p - 1, -p)
  }

  eval_jac_g_ineq <- function(theta) {
    rbind(L, -L)
  }

  q <- ncol(L)
  lb <- c(0, rep(-Inf, q - 1))
  ub <- c(1, rep( Inf, q - 1))

  # Multiple starting values
  btot <- aggregate(
    cbind(g_ijk, n_ijk) ~ Fac_B,
    data = df, FUN = sum
  )
  btot <- btot[order(btot$Fac_B), ]
  pB <- btot$g_ijk / btot$n_ijk

  atot <- aggregate(
    cbind(g_ijk, n_ijk) ~ Fac_A,
    data = df, FUN = sum
  )
  atot <- atot[order(atot$Fac_A), ]
  pA <- atot$g_ijk / atot$n_ijk

  p0 <- sum(df$g_ijk) / sum(df$n_ijk)

  # Keep starting values slightly in the interior.
  cap <- function(x) pmin(pmax(x, 1e-6), 1 - 1e-6)

  pB0 <- cap(pB)
  pA0 <- cap(pA)
  p00 <- cap(p0)

  theta_B <- c(
    pB0[1],
    rep(0, I - 1),
    pB0[2:J] - pB0[1]
  )

  theta_A <- c(
    pA0[1],
    pA0[2:I] - pA0[1],
    rep(0, J - 1)
  )

  theta_0 <- c(
    p00,
    rep(0, I - 1),
    rep(0, J - 1)
  )

  starts <- list(
    B_start = theta_B,
    A_start = theta_A,
    common_start = theta_0
  )

  fits <- vector("list", length(starts))
  names(fits) <- names(starts)

  for (sname in names(starts)) {
    ans <- try(
      nloptr::nloptr(
        x0 = starts[[sname]],
        eval_f = nll,
        eval_grad_f = grad_nll,
        lb = lb,
        ub = ub,
        eval_g_ineq = eval_g_ineq,
        eval_jac_g_ineq = eval_jac_g_ineq,
        opts = list(
          algorithm = "NLOPT_LD_SLSQP",
          maxeval = maxeval,
          xtol_rel = xtol_rel,
          ftol_rel = ftol_rel,
          print_level = print_level
        )
      ),
      silent = TRUE
    )

    if (inherits(ans, "try-error")) {
      fits[[sname]] <- list(
        start = sname,
        ok = FALSE,
        feasible = FALSE,
        logLik = -Inf,
        error = as.character(ans)
      )
      next
    }

    theta_hat <- as.numeric(ans$solution)
    p_hat <- as.vector(L %*% theta_hat)

    max_violation <- max(c(0, -p_hat, p_hat - 1))
    feasible <- is.finite(max_violation) &&
                max_violation <= feasibility_tol

    # Numerical cleanup only after feasibility assessment.
    p_eval <- pmin(pmax(p_hat, 0), 1)

    ll <- if (feasible) {
      .full_binom_loglik(df, p_eval, J)
    } else {
      -Inf
    }

    fits[[sname]] <- list(
      start = sname,
      ok = TRUE,
      feasible = feasible,
      status = ans$status,
      message = ans$message,
      iterations = ans$iterations,
      theta = theta_hat,
      p_cell = p_eval,
      logLik = ll,
      max_constraint_violation = max_violation
    )
  }

  llvec <- sapply(fits, function(z) z$logLik)

  if (all(!is.finite(llvec))) {
    stop("The constrained additive model did not yield a feasible solution.")
  }

  best_name <- names(which.max(llvec))
  best <- fits[[best_name]]

  pmat <- matrix(
    best$p_cell,
    nrow = I,
    ncol = J,
    byrow = TRUE,
    dimnames = list(
      paste0("A", seq_len(I)),
      paste0("B", seq_len(J))
    )
  )

  diagnostics <- do.call(
    rbind,
    lapply(fits, function(z) {
      data.frame(
        start = z$start,
        status = if (!is.null(z$status)) z$status else NA_integer_,
        feasible = isTRUE(z$feasible),
        logLik = if (is.finite(z$logLik)) z$logLik else NA_real_,
        max_constraint_violation =
          if (!is.null(z$max_constraint_violation))
            z$max_constraint_violation else NA_real_,
        stringsAsFactors = FALSE
      )
    })
  )

  list(
    theta = best$theta,
    p_cell = best$p_cell,
    p_matrix = pmat,
    logLik = best$logLik,
    start = best$start,
    status = best$status,
    message = best$message,
    iterations = best$iterations,
    max_constraint_violation = best$max_constraint_violation,
    diagnostics = diagnostics,
    L = L
  )
}


# ============================================================
# Overall tests
# ============================================================

.fit_overall_models <- function(
    df,
    maxeval = 10000,
    xtol_rel = 1e-10,
    ftol_rel = 1e-12
) {
  prep <- .prepare_data(df)
  df <- prep$df
  I <- prep$I
  J <- prep$J

  # M_A
  a <- aggregate(
    cbind(g_ijk, n_ijk) ~ Fac_A,
    data = df, FUN = sum
  )
  a <- a[order(a$Fac_A), ]
  pA <- a$g_ijk / a$n_ijk
  llA <- sum(
    dbinom(
      df$g_ijk,
      size = df$n_ijk,
      prob = pA[df$Fac_A],
      log = TRUE
    )
  )

  # M_B
  b <- aggregate(
    cbind(g_ijk, n_ijk) ~ Fac_B,
    data = df, FUN = sum
  )
  b <- b[order(b$Fac_B), ]
  pB <- b$g_ijk / b$n_ijk
  llB <- sum(
    dbinom(
      df$g_ijk,
      size = df$n_ijk,
      prob = pB[df$Fac_B],
      log = TRUE
    )
  )

  # M_{A+B}
  add <- .fit_additive_binomial(
    df,
    maxeval = maxeval,
    xtol_rel = xtol_rel,
    ftol_rel = ftol_rel
  )
  llAdd <- add$logLik

  # Saturated M_{A x B}
  sat <- .cell_totals(df, I, J)
  pSat <- sat$g_ijk / sat$n_ijk
  idx <- (df$Fac_A - 1L) * J + df$Fac_B

  llSat <- sum(
    dbinom(
      df$g_ijk,
      size = df$n_ijk,
      prob = pSat[idx],
      log = TRUE
    )
  )

  LR_A   <- max(0, 2 * (llAdd - llB))
  LR_B   <- max(0, 2 * (llAdd - llA))
  LR_Int <- max(0, 2 * (llSat - llAdd))

  df_A   <- I - 1
  df_B   <- J - 1
  df_Int <- (I - 1) * (J - 1)

  tests <- data.frame(
    Effect = c("Factor A", "Factor B", "Interaction"),
    LR = c(LR_A, LR_B, LR_Int),
    df = c(df_A, df_B, df_Int),
    p_value = c(
      pchisq(LR_A, df_A, lower.tail = FALSE),
      pchisq(LR_B, df_B, lower.tail = FALSE),
      pchisq(LR_Int, df_Int, lower.tail = FALSE)
    ),
    stringsAsFactors = FALSE
  )

  list(
    tests = tests,
    additive_fit = add,
    fitted = list(
      M_A = pA,
      M_B = pB,
      M_AplusB = add$p_matrix,
      M_AxB = matrix(
        pSat,
        nrow = I,
        ncol = J,
        byrow = TRUE,
        dimnames = list(
          paste0("A", seq_len(I)),
          paste0("B", seq_len(J))
        )
      )
    ),
    logLik = c(
      M_A = llA,
      M_B = llB,
      M_AplusB = llAdd,
      M_AxB = llSat
    )
  )
}


# ============================================================
# Main-effect pairwise tests under M_{A+B}
# ============================================================

.fit_pairwise_null_additive <- function(
    df,
    factor = c("A", "B"),
    level1,
    level2,
    maxeval = 10000,
    xtol_rel = 1e-10,
    ftol_rel = 1e-12
) {
  factor <- match.arg(factor)
  d <- df

  if (factor == "A") {
    d$Fac_A[d$Fac_A == level2] <- level1
    vals <- sort(unique(d$Fac_A))
    d$Fac_A <- match(d$Fac_A, vals)
  } else {
    d$Fac_B[d$Fac_B == level2] <- level1
    vals <- sort(unique(d$Fac_B))
    d$Fac_B <- match(d$Fac_B, vals)
  }

  .fit_additive_binomial(
    d,
    maxeval = maxeval,
    xtol_rel = xtol_rel,
    ftol_rel = ftol_rel
  )
}


.posthoc_main_effect <- function(
    df,
    overall,
    factor = c("A", "B"),
    alpha = 0.05,
    adjust = "holm",
    maxeval = 10000,
    xtol_rel = 1e-10,
    ftol_rel = 1e-12
) {
  factor <- match.arg(factor)
  prep <- .prepare_data(df)
  I <- prep$I
  J <- prep$J
  K <- if (factor == "A") I else J

  pairs <- .pair_matrix(K)
  ll_alt <- overall$additive_fit$logLik
  Padd <- overall$additive_fit$p_matrix

  out <- vector("list", nrow(pairs))

  for (k in seq_len(nrow(pairs))) {
    l1 <- pairs[k, 1]
    l2 <- pairs[k, 2]

    nullfit <- .fit_pairwise_null_additive(
      df,
      factor = factor,
      level1 = l1,
      level2 = l2,
      maxeval = maxeval,
      xtol_rel = xtol_rel,
      ftol_rel = ftol_rel
    )

    LR <- max(0, 2 * (ll_alt - nullfit$logLik))
    p <- pchisq(LR, df = 1, lower.tail = FALSE)

    if (factor == "A") {
      est <- mean(Padd[l1, ] - Padd[l2, ])
    } else {
      est <- mean(Padd[, l1] - Padd[, l2])
    }

    out[[k]] <- data.frame(
      Factor = factor,
      Level1 = l1,
      Level2 = l2,
      Contrast = paste0(factor, l1, " - ", factor, l2),
      Estimate = est,
      LR = LR,
      df = 1,
      p_value = p,
      stringsAsFactors = FALSE
    )
  }

  ans <- do.call(rbind, out)
  ans$p_adjusted <- p.adjust(ans$p_value, method = adjust)
  ans$significant <- ans$p_adjusted < alpha
  rownames(ans) <- NULL
  ans
}


# ============================================================
# Simple-effect pairwise tests
# ============================================================

.simple_pairwise_LRT <- function(
    df,
    target = c("A_within_B", "B_within_A"),
    alpha = 0.05,
    adjust = "holm"
) {
  target <- match.arg(target)
  prep <- .prepare_data(df)
  df <- prep$df
  I <- prep$I
  J <- prep$J

  out <- list()
  z <- 0L

  if (target == "A_within_B") {
    for (j in seq_len(J)) {
      pairs <- .pair_matrix(I)

      for (k in seq_len(nrow(pairs))) {
        i1 <- pairs[k, 1]
        i2 <- pairs[k, 2]

        G1 <- sum(df$g_ijk[df$Fac_A == i1 & df$Fac_B == j])
        N1 <- sum(df$n_ijk[df$Fac_A == i1 & df$Fac_B == j])
        G2 <- sum(df$g_ijk[df$Fac_A == i2 & df$Fac_B == j])
        N2 <- sum(df$n_ijk[df$Fac_A == i2 & df$Fac_B == j])

        p1 <- G1 / N1
        p2 <- G2 / N2
        p0 <- (G1 + G2) / (N1 + N2)

        ll_alt <- sum(
          dbinom(
            c(G1, G2),
            size = c(N1, N2),
            prob = c(p1, p2),
            log = TRUE
          )
        )

        ll_null <- sum(
          dbinom(
            c(G1, G2),
            size = c(N1, N2),
            prob = c(p0, p0),
            log = TRUE
          )
        )

        LR <- max(0, 2 * (ll_alt - ll_null))

        z <- z + 1L
        out[[z]] <- data.frame(
          Simple_effect = paste0("A within B", j),
          Fixed_level = j,
          Level1 = i1,
          Level2 = i2,
          Contrast = paste0("A", i1, " - A", i2),
          Estimate = p1 - p2,
          LR = LR,
          df = 1,
          p_value = pchisq(LR, 1, lower.tail = FALSE),
          stringsAsFactors = FALSE
        )
      }
    }

  } else {
    for (i in seq_len(I)) {
      pairs <- .pair_matrix(J)

      for (k in seq_len(nrow(pairs))) {
        j1 <- pairs[k, 1]
        j2 <- pairs[k, 2]

        G1 <- sum(df$g_ijk[df$Fac_A == i & df$Fac_B == j1])
        N1 <- sum(df$n_ijk[df$Fac_A == i & df$Fac_B == j1])
        G2 <- sum(df$g_ijk[df$Fac_A == i & df$Fac_B == j2])
        N2 <- sum(df$n_ijk[df$Fac_A == i & df$Fac_B == j2])

        p1 <- G1 / N1
        p2 <- G2 / N2
        p0 <- (G1 + G2) / (N1 + N2)

        ll_alt <- sum(
          dbinom(
            c(G1, G2),
            size = c(N1, N2),
            prob = c(p1, p2),
            log = TRUE
          )
        )

        ll_null <- sum(
          dbinom(
            c(G1, G2),
            size = c(N1, N2),
            prob = c(p0, p0),
            log = TRUE
          )
        )

        LR <- max(0, 2 * (ll_alt - ll_null))

        z <- z + 1L
        out[[z]] <- data.frame(
          Simple_effect = paste0("B within A", i),
          Fixed_level = i,
          Level1 = j1,
          Level2 = j2,
          Contrast = paste0("B", j1, " - B", j2),
          Estimate = p1 - p2,
          LR = LR,
          df = 1,
          p_value = pchisq(LR, 1, lower.tail = FALSE),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  ans <- do.call(rbind, out)

  # Multiplicity family: all pairwise comparisons within one fixed level
  # of the conditioning factor.
  ans$p_adjusted <- ave(
    ans$p_value,
    ans$Fixed_level,
    FUN = function(x) p.adjust(x, method = adjust)
  )

  ans$significant <- ans$p_adjusted < alpha
  rownames(ans) <- NULL
  ans
}


# ============================================================
# Compact letter displays
# ============================================================

.make_cld <- function(level_labels, pairwise_table, alpha = 0.05) {
  if (length(level_labels) == 1) {
    return(data.frame(Level = level_labels, Group = "a"))
  }

  pvec <- pairwise_table$p_adjusted
  names(pvec) <- paste(
    level_labels[pairwise_table$Level1],
    level_labels[pairwise_table$Level2],
    sep = "-"
  )

  letters <- multcompView::multcompLetters(
    pvec,
    threshold = alpha
  )$Letters

  data.frame(
    Level = names(letters),
    Group = unname(letters),
    stringsAsFactors = FALSE
  )
}


.cld_main_effect <- function(
    overall,
    factor = c("A", "B"),
    pairwise_table,
    alpha = 0.05
) {
  factor <- match.arg(factor)
  P <- overall$additive_fit$p_matrix

  if (factor == "A") {
    labels <- rownames(P)
    avg <- rowMeans(P)

    cld <- .make_cld(labels, pairwise_table, alpha)
    est <- data.frame(
      Level = labels,
      Average_fitted_probability = as.numeric(avg),
      stringsAsFactors = FALSE
    )

  } else {
    labels <- colnames(P)
    avg <- colMeans(P)

    cld <- .make_cld(labels, pairwise_table, alpha)
    est <- data.frame(
      Level = labels,
      Average_fitted_probability = as.numeric(avg),
      stringsAsFactors = FALSE
    )
  }

  ans <- merge(est, cld, by = "Level", sort = FALSE)
  ans <- ans[match(labels, ans$Level), ]
  rownames(ans) <- NULL
  ans
}


.cld_simple_effect <- function(
    df,
    target = c("A_within_B", "B_within_A"),
    pairwise_table,
    alpha = 0.05
) {
  target <- match.arg(target)
  prep <- .prepare_data(df)
  df <- prep$df
  I <- prep$I
  J <- prep$J

  out <- list()
  z <- 0L

  if (target == "A_within_B") {
    for (j in seq_len(J)) {
      tab <- pairwise_table[pairwise_table$Fixed_level == j, ]
      labels <- paste0("A", seq_len(I))
      cld <- .make_cld(labels, tab, alpha)

      probs <- sapply(seq_len(I), function(i) {
        G <- sum(df$g_ijk[df$Fac_A == i & df$Fac_B == j])
        N <- sum(df$n_ijk[df$Fac_A == i & df$Fac_B == j])
        G / N
      })

      tmp <- data.frame(
        Fixed_factor = "B",
        Fixed_level = j,
        Level = labels,
        Observed_probability = probs,
        stringsAsFactors = FALSE
      )
      tmp <- merge(tmp, cld, by = "Level", sort = FALSE)
      tmp <- tmp[match(labels, tmp$Level), ]

      z <- z + 1L
      out[[z]] <- tmp
    }

  } else {
    for (i in seq_len(I)) {
      tab <- pairwise_table[pairwise_table$Fixed_level == i, ]
      labels <- paste0("B", seq_len(J))
      cld <- .make_cld(labels, tab, alpha)

      probs <- sapply(seq_len(J), function(j) {
        G <- sum(df$g_ijk[df$Fac_A == i & df$Fac_B == j])
        N <- sum(df$n_ijk[df$Fac_A == i & df$Fac_B == j])
        G / N
      })

      tmp <- data.frame(
        Fixed_factor = "A",
        Fixed_level = i,
        Level = labels,
        Observed_probability = probs,
        stringsAsFactors = FALSE
      )
      tmp <- merge(tmp, cld, by = "Level", sort = FALSE)
      tmp <- tmp[match(labels, tmp$Level), ]

      z <- z + 1L
      out[[z]] <- tmp
    }
  }

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}



# ============================================================
# Logit-scale models: ordinary MLE and Firth penalized likelihood
# ============================================================

.prepare_logit_data <- function(df) {
  prep <- .prepare_data(df)
  d <- prep$df
  d$A <- factor(d$Fac_A, levels = seq_len(prep$I),
                labels = paste0("A", seq_len(prep$I)))
  d$B <- factor(d$Fac_B, levels = seq_len(prep$J),
                labels = paste0("B", seq_len(prep$J)))
  list(df = d, I = prep$I, J = prep$J)
}

.logit_cell_grid <- function(I, J) {
  expand.grid(
    A = factor(paste0("A", seq_len(I)), levels = paste0("A", seq_len(I))),
    B = factor(paste0("B", seq_len(J)), levels = paste0("B", seq_len(J))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

.as_probability_matrix <- function(p, I, J) {
  # expand.grid varies A fastest; convert explicitly to A x B matrix.
  matrix(
    p,
    nrow = I,
    ncol = J,
    byrow = FALSE,
    dimnames = list(paste0("A", seq_len(I)), paste0("B", seq_len(J)))
  )
}

.fit_overall_logit_mle <- function(df) {
  prep <- .prepare_logit_data(df)
  d <- prep$df
  I <- prep$I
  J <- prep$J

  form_B   <- cbind(g_ijk, n_ijk - g_ijk) ~ B
  form_A   <- cbind(g_ijk, n_ijk - g_ijk) ~ A
  form_Add <- cbind(g_ijk, n_ijk - g_ijk) ~ A + B
  form_Int <- cbind(g_ijk, n_ijk - g_ijk) ~ A * B

  fit_B <- stats::glm(form_B, data = d, family = stats::binomial(link = "logit"))
  fit_A <- stats::glm(form_A, data = d, family = stats::binomial(link = "logit"))
  fit_Add <- stats::glm(form_Add, data = d, family = stats::binomial(link = "logit"))
  fit_Int <- stats::glm(form_Int, data = d, family = stats::binomial(link = "logit"))

  llB <- as.numeric(stats::logLik(fit_B))
  llA <- as.numeric(stats::logLik(fit_A))
  llAdd <- as.numeric(stats::logLik(fit_Add))
  llInt <- as.numeric(stats::logLik(fit_Int))

  LR_A <- max(0, 2 * (llAdd - llB))
  LR_B <- max(0, 2 * (llAdd - llA))
  LR_Int <- max(0, 2 * (llInt - llAdd))

  df_A <- I - 1
  df_B <- J - 1
  df_Int <- (I - 1) * (J - 1)

  tests <- data.frame(
    Effect = c("Factor A", "Factor B", "Interaction"),
    LR = c(LR_A, LR_B, LR_Int),
    df = c(df_A, df_B, df_Int),
    p_value = c(
      stats::pchisq(LR_A, df_A, lower.tail = FALSE),
      stats::pchisq(LR_B, df_B, lower.tail = FALSE),
      stats::pchisq(LR_Int, df_Int, lower.tail = FALSE)
    ),
    stringsAsFactors = FALSE
  )

  grid <- .logit_cell_grid(I, J)
  fitted <- list(
    M_A = .as_probability_matrix(
      stats::predict(fit_A, newdata = grid, type = "response"), I, J
    ),
    M_B = .as_probability_matrix(
      stats::predict(fit_B, newdata = grid, type = "response"), I, J
    ),
    M_AplusB = .as_probability_matrix(
      stats::predict(fit_Add, newdata = grid, type = "response"), I, J
    ),
    M_AxB = .as_probability_matrix(
      stats::predict(fit_Int, newdata = grid, type = "response"), I, J
    )
  )

  model_list <- list(M_B = fit_B, M_A = fit_A, M_AplusB = fit_Add, M_AxB = fit_Int)
  coef_all <- unlist(lapply(model_list, stats::coef))

  se_all <- unlist(lapply(model_list, function(fit) {
    V <- try(stats::vcov(fit), silent = TRUE)
    if (inherits(V, "try-error")) return(NA_real_)
    sqrt(diag(V))
  }))

  fitted_all <- unlist(lapply(model_list, stats::fitted))
  min_fitted <- suppressWarnings(min(fitted_all, na.rm = TRUE))
  max_fitted <- suppressWarnings(max(fitted_all, na.rm = TRUE))
  max_abs_coef <- suppressWarnings(max(abs(coef_all), na.rm = TRUE))
  max_se <- suppressWarnings(max(se_all, na.rm = TRUE))

  # Heuristic diagnostic only.  A TRUE value is not a formal proof of
  # separation, but flags the common numerical symptoms encountered when
  # fitted probabilities are driven very close to 0 or 1.
  possible_separation <-
    any(!is.finite(coef_all)) ||
    is.finite(max_abs_coef) && max_abs_coef > 10 ||
    is.finite(max_se) && max_se > 100 ||
    is.finite(min_fitted) && min_fitted < 1e-8 ||
    is.finite(max_fitted) && max_fitted > 1 - 1e-8

  diagnostics <- list(
    converged = c(
      M_B = isTRUE(fit_B$converged),
      M_A = isTRUE(fit_A$converged),
      M_AplusB = isTRUE(fit_Add$converged),
      M_AxB = isTRUE(fit_Int$converged)
    ),
    all_coefficients_finite = all(is.finite(coef_all)),
    max_abs_coefficient = max_abs_coef,
    max_standard_error = max_se,
    min_fitted_probability = min_fitted,
    max_fitted_probability = max_fitted,
    possible_separation = possible_separation,
    note = paste(
      "possible_separation is a heuristic diagnostic based on large/non-finite",
      "coefficients, very large standard errors, or fitted probabilities extremely",
      "close to 0 or 1. If TRUE, consider scale = 'logit', method = 'firth'."
    )
  )

  list(
    tests = tests,
    fitted = fitted,
    logLik = c(M_B = llB, M_A = llA, M_AplusB = llAdd, M_AxB = llInt),
    models = list(M_B = fit_B, M_A = fit_A, M_AplusB = fit_Add, M_AxB = fit_Int),
    diagnostics = diagnostics
  )
}

.expand_binomial_for_firth <- function(df) {
  prep <- .prepare_logit_data(df)
  d <- prep$df

  rows <- vector("list", 2L * nrow(d))
  z <- 0L

  for (r in seq_len(nrow(d))) {
    g <- d$g_ijk[r]
    f <- d$n_ijk[r] - d$g_ijk[r]

    if (g > 0) {
      z <- z + 1L
      rows[[z]] <- data.frame(
        y = 1,
        A = d$A[r],
        B = d$B[r],
        weight = g
      )
    }

    if (f > 0) {
      z <- z + 1L
      rows[[z]] <- data.frame(
        y = 0,
        A = d$A[r],
        B = d$B[r],
        weight = f
      )
    }
  }

  if (z == 0L) stop("No observations available for Firth logistic regression.")

  out <- do.call(rbind, rows[seq_len(z)])
  out$A <- factor(out$A, levels = paste0("A", seq_len(prep$I)))
  out$B <- factor(out$B, levels = paste0("B", seq_len(prep$J)))

  list(df = out, I = prep$I, J = prep$J)
}

.firth_nested_test <- function(full_fit, reduced_fit) {
  z <- stats::anova(full_fit, reduced_fit, method = "nested")
  chisq <- unname(z$chisq)
  df <- unname(z$df)

  # Recompute the upper-tail probability from the reported PLR statistic.
  # This avoids display values of exactly 0 that can arise from numerical
  # underflow in package output for very strong effects.
  p_value <- stats::pchisq(chisq, df = df, lower.tail = FALSE)

  c(chisq = chisq, df = df, p_value = p_value)
}

.fit_overall_logit_firth <- function(df) {
  if (!requireNamespace("logistf", quietly = TRUE)) {
    stop(
      "Package 'logistf' is required for method = 'firth'. ",
      "Install it with install.packages('logistf')."
    )
  }

  prep <- .expand_binomial_for_firth(df)
  d <- prep$df
  I <- prep$I
  J <- prep$J

  # pl = FALSE avoids computing profile-likelihood intervals for every fit;
  # the factorial tests below use logistf's nested penalized LRT via anova().
  fit_B <- logistf::logistf(y ~ B, data = d, weights = weight, firth = TRUE, pl = FALSE)
  fit_A <- logistf::logistf(y ~ A, data = d, weights = weight, firth = TRUE, pl = FALSE)
  fit_Add <- logistf::logistf(y ~ A + B, data = d, weights = weight, firth = TRUE, pl = FALSE)
  fit_Int <- logistf::logistf(y ~ A * B, data = d, weights = weight, firth = TRUE, pl = FALSE)

  tA <- .firth_nested_test(fit_Add, fit_B)
  tB <- .firth_nested_test(fit_Add, fit_A)
  tI <- .firth_nested_test(fit_Int, fit_Add)

  tests <- data.frame(
    Effect = c("Factor A", "Factor B", "Interaction"),
    LR = c(tA["chisq"], tB["chisq"], tI["chisq"]),
    df = as.integer(c(tA["df"], tB["df"], tI["df"])),
    p_value = c(tA["p_value"], tB["p_value"], tI["p_value"]),
    stringsAsFactors = FALSE
  )
  rownames(tests) <- NULL

  grid <- .logit_cell_grid(I, J)
  fitted <- list(
    M_A = .as_probability_matrix(
      stats::predict(fit_A, newdata = grid, type = "response"), I, J
    ),
    M_B = .as_probability_matrix(
      stats::predict(fit_B, newdata = grid, type = "response"), I, J
    ),
    M_AplusB = .as_probability_matrix(
      stats::predict(fit_Add, newdata = grid, type = "response"), I, J
    ),
    M_AxB = .as_probability_matrix(
      stats::predict(fit_Int, newdata = grid, type = "response"), I, J
    )
  )

  # logistf stores penalized log-likelihoods as c(null, full).
  pll <- function(fit) unname(fit$loglik["full"])

  list(
    tests = tests,
    fitted = fitted,
    penalized_logLik = c(
      M_B = pll(fit_B),
      M_A = pll(fit_A),
      M_AplusB = pll(fit_Add),
      M_AxB = pll(fit_Int)
    ),
    models = list(M_B = fit_B, M_A = fit_A, M_AplusB = fit_Add, M_AxB = fit_Int),
    diagnostics = list(
      note = paste(
        "Factorial tests are nested penalized likelihood-ratio tests from",
        "logistf::anova.logistf(method = 'nested'). Frequency weights are used",
        "to represent the original binomial counts without expanding to one row per seed."
      )
    )
  )
}


# ============================================================
# Logit-scale post-hoc comparisons
# ============================================================

.recode_pairwise_factor <- function(df, factor = c("A", "B"), level1, level2) {
  factor <- match.arg(factor)
  d <- .prepare_data(df)$df

  if (factor == "A") {
    d$Fac_A[d$Fac_A == level2] <- level1
    vals <- sort(unique(d$Fac_A))
    d$Fac_A <- match(d$Fac_A, vals)
  } else {
    d$Fac_B[d$Fac_B == level2] <- level1
    vals <- sort(unique(d$Fac_B))
    d$Fac_B <- match(d$Fac_B, vals)
  }

  d
}


.fit_logit_additive_reduced_mle <- function(
    df,
    factor = c("A", "B"),
    level1,
    level2
) {
  factor <- match.arg(factor)
  d0 <- .recode_pairwise_factor(df, factor, level1, level2)
  prep <- .prepare_logit_data(d0)
  d <- prep$df

  stats::glm(
    cbind(g_ijk, n_ijk - g_ijk) ~ A + B,
    data = d,
    family = stats::binomial(link = "logit")
  )
}


.fit_logit_additive_reduced_firth <- function(
    df,
    factor = c("A", "B"),
    level1,
    level2
) {
  factor <- match.arg(factor)

  if (!requireNamespace("logistf", quietly = TRUE)) {
    stop(
      "Package 'logistf' is required for method = 'firth'. ",
      "Install it with install.packages('logistf')."
    )
  }

  d0 <- .recode_pairwise_factor(df, factor, level1, level2)
  prep <- .expand_binomial_for_firth(d0)
  d <- prep$df

  logistf::logistf(
    y ~ A + B,
    data = d,
    weights = weight,
    firth = TRUE,
    pl = FALSE
  )
}


.posthoc_logit_main_effect <- function(
    df,
    overall,
    factor = c("A", "B"),
    method = c("mle", "firth"),
    alpha = 0.05,
    adjust = "holm"
) {
  factor <- match.arg(factor)
  method <- match.arg(method)

  prep <- .prepare_data(df)
  I <- prep$I
  J <- prep$J
  K <- if (factor == "A") I else J

  pairs <- .pair_matrix(K)
  Padd <- overall$fitted$M_AplusB
  altfit <- overall$models$M_AplusB

  out <- vector("list", nrow(pairs))

  for (k in seq_len(nrow(pairs))) {
    l1 <- pairs[k, 1]
    l2 <- pairs[k, 2]

    if (method == "mle") {
      redfit <- .fit_logit_additive_reduced_mle(
        df, factor = factor, level1 = l1, level2 = l2
      )
      LR <- max(
        0,
        2 * (
          as.numeric(stats::logLik(altfit)) -
          as.numeric(stats::logLik(redfit))
        )
      )
      p <- stats::pchisq(LR, df = 1, lower.tail = FALSE)

    } else {
      redfit <- .fit_logit_additive_reduced_firth(
        df, factor = factor, level1 = l1, level2 = l2
      )

      # The reduced pairwise model is created by merging two factor levels.
      # Although its design space is nested in the additive model, the two
      # fitted logistf objects have different formula term labels.  Therefore
      # anova.logistf(method = "nested") cannot recover the nesting from the
      # formulas and may fail inside terms.formula().
      #
      # Use logistf's PLR comparison for this specific merged-level comparison.
      # The models differ by one degree of freedom, and the test still compares
      # the same logit-scale equality hypothesis for the two factor levels.
      z <- stats::anova(altfit, redfit, method = "PLR")
      LR <- unname(z$chisq)

      # Recompute from the chi-square statistic so extremely small
      # p-values remain representable instead of being returned as 0.
      p <- stats::pchisq(LR, df = 1, lower.tail = FALSE)
    }

    if (factor == "A") {
      log_odds_diff <- mean(stats::qlogis(Padd[l1, ]) - stats::qlogis(Padd[l2, ]))
      prob_diff <- mean(Padd[l1, ] - Padd[l2, ])
    } else {
      log_odds_diff <- mean(stats::qlogis(Padd[, l1]) - stats::qlogis(Padd[, l2]))
      prob_diff <- mean(Padd[, l1] - Padd[, l2])
    }

    out[[k]] <- data.frame(
      Factor = factor,
      Level1 = l1,
      Level2 = l2,
      Contrast = paste0(factor, l1, " - ", factor, l2),
      Log_odds_difference = log_odds_diff,
      Odds_ratio = exp(log_odds_diff),
      Average_probability_difference = prob_diff,
      LR = LR,
      df = 1,
      p_value = p,
      stringsAsFactors = FALSE
    )
  }

  ans <- do.call(rbind, out)
  ans$p_adjusted <- stats::p.adjust(ans$p_value, method = adjust)
  ans$significant <- ans$p_adjusted < alpha
  rownames(ans) <- NULL
  ans
}


.fit_simple_logit_pair_mle <- function(
    df,
    target = c("A_within_B", "B_within_A"),
    fixed_level,
    level1,
    level2
) {
  target <- match.arg(target)
  prep <- .prepare_data(df)
  d <- prep$df

  if (target == "A_within_B") {
    d <- d[
      d$Fac_B == fixed_level & d$Fac_A %in% c(level1, level2),
      , drop = FALSE
    ]
    d$G <- factor(
      d$Fac_A,
      levels = c(level1, level2),
      labels = c("g1", "g2")
    )
  } else {
    d <- d[
      d$Fac_A == fixed_level & d$Fac_B %in% c(level1, level2),
      , drop = FALSE
    ]
    d$G <- factor(
      d$Fac_B,
      levels = c(level1, level2),
      labels = c("g1", "g2")
    )
  }

  fit0 <- stats::glm(
    cbind(g_ijk, n_ijk - g_ijk) ~ 1,
    data = d,
    family = stats::binomial(link = "logit")
  )

  fit1 <- stats::glm(
    cbind(g_ijk, n_ijk - g_ijk) ~ G,
    data = d,
    family = stats::binomial(link = "logit")
  )

  LR <- max(
    0,
    2 * (
      as.numeric(stats::logLik(fit1)) -
      as.numeric(stats::logLik(fit0))
    )
  )

  b <- stats::coef(fit1)
  # Gg2 is level2 - level1; report level1 - level2.
  lod <- -unname(b["Gg2"])

  list(
    LR = LR,
    p_value = stats::pchisq(LR, 1, lower.tail = FALSE),
    log_odds_difference = lod,
    odds_ratio = exp(lod),
    models = list(null = fit0, alternative = fit1)
  )
}


.fit_simple_logit_pair_firth <- function(
    df,
    target = c("A_within_B", "B_within_A"),
    fixed_level,
    level1,
    level2
) {
  target <- match.arg(target)

  if (!requireNamespace("logistf", quietly = TRUE)) {
    stop(
      "Package 'logistf' is required for method = 'firth'. ",
      "Install it with install.packages('logistf')."
    )
  }

  prep <- .prepare_data(df)
  d <- prep$df

  if (target == "A_within_B") {
    d <- d[
      d$Fac_B == fixed_level & d$Fac_A %in% c(level1, level2),
      , drop = FALSE
    ]
    group_num <- ifelse(d$Fac_A == level1, 1L, 2L)
  } else {
    d <- d[
      d$Fac_A == fixed_level & d$Fac_B %in% c(level1, level2),
      , drop = FALSE
    ]
    group_num <- ifelse(d$Fac_B == level1, 1L, 2L)
  }

  rows <- list()
  z <- 0L
  for (r in seq_len(nrow(d))) {
    g <- d$g_ijk[r]
    f <- d$n_ijk[r] - d$g_ijk[r]

    if (g > 0) {
      z <- z + 1L
      rows[[z]] <- data.frame(
        y = 1,
        G = factor(group_num[r], levels = c(1, 2), labels = c("g1", "g2")),
        weight = g
      )
    }

    if (f > 0) {
      z <- z + 1L
      rows[[z]] <- data.frame(
        y = 0,
        G = factor(group_num[r], levels = c(1, 2), labels = c("g1", "g2")),
        weight = f
      )
    }
  }

  ed <- do.call(rbind, rows[seq_len(z)])
  ed$G <- factor(ed$G, levels = c("g1", "g2"))

  fit0 <- logistf::logistf(
    y ~ 1, data = ed, weights = weight, firth = TRUE, pl = FALSE
  )
  fit1 <- logistf::logistf(
    y ~ G, data = ed, weights = weight, firth = TRUE, pl = FALSE
  )

  tt <- .firth_nested_test(fit1, fit0)
  b <- stats::coef(fit1)
  lod <- -unname(b["Gg2"])

  list(
    LR = unname(tt["chisq"]),
    p_value = unname(tt["p_value"]),
    log_odds_difference = lod,
    odds_ratio = exp(lod),
    models = list(null = fit0, alternative = fit1)
  )
}


.posthoc_logit_simple_effect <- function(
    df,
    target = c("A_within_B", "B_within_A"),
    method = c("mle", "firth"),
    alpha = 0.05,
    adjust = "holm"
) {
  target <- match.arg(target)
  method <- match.arg(method)

  prep <- .prepare_data(df)
  d <- prep$df
  I <- prep$I
  J <- prep$J

  out <- list()
  z <- 0L

  if (target == "A_within_B") {
    for (j in seq_len(J)) {
      pairs <- .pair_matrix(I)

      for (k in seq_len(nrow(pairs))) {
        i1 <- pairs[k, 1]
        i2 <- pairs[k, 2]

        fit <- if (method == "mle") {
          .fit_simple_logit_pair_mle(
            d, "A_within_B", fixed_level = j, level1 = i1, level2 = i2
          )
        } else {
          .fit_simple_logit_pair_firth(
            d, "A_within_B", fixed_level = j, level1 = i1, level2 = i2
          )
        }

        G1 <- sum(d$g_ijk[d$Fac_A == i1 & d$Fac_B == j])
        N1 <- sum(d$n_ijk[d$Fac_A == i1 & d$Fac_B == j])
        G2 <- sum(d$g_ijk[d$Fac_A == i2 & d$Fac_B == j])
        N2 <- sum(d$n_ijk[d$Fac_A == i2 & d$Fac_B == j])

        z <- z + 1L
        out[[z]] <- data.frame(
          Simple_effect = paste0("A within B", j),
          Fixed_level = j,
          Level1 = i1,
          Level2 = i2,
          Contrast = paste0("A", i1, " - A", i2),
          Log_odds_difference = fit$log_odds_difference,
          Odds_ratio = fit$odds_ratio,
          Observed_probability_difference = G1 / N1 - G2 / N2,
          LR = fit$LR,
          df = 1,
          p_value = fit$p_value,
          stringsAsFactors = FALSE
        )
      }
    }

  } else {
    for (i in seq_len(I)) {
      pairs <- .pair_matrix(J)

      for (k in seq_len(nrow(pairs))) {
        j1 <- pairs[k, 1]
        j2 <- pairs[k, 2]

        fit <- if (method == "mle") {
          .fit_simple_logit_pair_mle(
            d, "B_within_A", fixed_level = i, level1 = j1, level2 = j2
          )
        } else {
          .fit_simple_logit_pair_firth(
            d, "B_within_A", fixed_level = i, level1 = j1, level2 = j2
          )
        }

        G1 <- sum(d$g_ijk[d$Fac_A == i & d$Fac_B == j1])
        N1 <- sum(d$n_ijk[d$Fac_A == i & d$Fac_B == j1])
        G2 <- sum(d$g_ijk[d$Fac_A == i & d$Fac_B == j2])
        N2 <- sum(d$n_ijk[d$Fac_A == i & d$Fac_B == j2])

        z <- z + 1L
        out[[z]] <- data.frame(
          Simple_effect = paste0("B within A", i),
          Fixed_level = i,
          Level1 = j1,
          Level2 = j2,
          Contrast = paste0("B", j1, " - B", j2),
          Log_odds_difference = fit$log_odds_difference,
          Odds_ratio = fit$odds_ratio,
          Observed_probability_difference = G1 / N1 - G2 / N2,
          LR = fit$LR,
          df = 1,
          p_value = fit$p_value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  ans <- do.call(rbind, out)
  ans$p_adjusted <- ave(
    ans$p_value,
    ans$Fixed_level,
    FUN = function(x) stats::p.adjust(x, method = adjust)
  )
  ans$significant <- ans$p_adjusted < alpha
  rownames(ans) <- NULL
  ans
}


.cld_logit_main_effect <- function(
    overall,
    factor = c("A", "B"),
    pairwise_table,
    alpha = 0.05
) {
  factor <- match.arg(factor)
  P <- overall$fitted$M_AplusB

  if (factor == "A") {
    labels <- rownames(P)
    avg <- rowMeans(P)
  } else {
    labels <- colnames(P)
    avg <- colMeans(P)
  }

  cld <- .make_cld(labels, pairwise_table, alpha)
  est <- data.frame(
    Level = labels,
    Average_fitted_probability = as.numeric(avg),
    stringsAsFactors = FALSE
  )

  ans <- merge(est, cld, by = "Level", sort = FALSE)
  ans <- ans[match(labels, ans$Level), ]
  rownames(ans) <- NULL
  ans
}


.cld_logit_simple_effect <- function(
    overall,
    target = c("A_within_B", "B_within_A"),
    pairwise_table,
    alpha = 0.05
) {
  target <- match.arg(target)
  P <- overall$fitted$M_AxB
  I <- nrow(P)
  J <- ncol(P)

  out <- list()
  z <- 0L

  if (target == "A_within_B") {
    for (j in seq_len(J)) {
      tab <- pairwise_table[pairwise_table$Fixed_level == j, ]
      labels <- rownames(P)
      cld <- .make_cld(labels, tab, alpha)

      tmp <- data.frame(
        Fixed_factor = "B",
        Fixed_level = j,
        Level = labels,
        Fitted_probability = as.numeric(P[, j]),
        stringsAsFactors = FALSE
      )

      tmp <- merge(tmp, cld, by = "Level", sort = FALSE)
      tmp <- tmp[match(labels, tmp$Level), ]
      z <- z + 1L
      out[[z]] <- tmp
    }

  } else {
    for (i in seq_len(I)) {
      tab <- pairwise_table[pairwise_table$Fixed_level == i, ]
      labels <- colnames(P)
      cld <- .make_cld(labels, tab, alpha)

      tmp <- data.frame(
        Fixed_factor = "A",
        Fixed_level = i,
        Level = labels,
        Fitted_probability = as.numeric(P[i, ]),
        stringsAsFactors = FALSE
      )

      tmp <- merge(tmp, cld, by = "Level", sort = FALSE)
      tmp <- tmp[match(labels, tmp$Level), ]
      z <- z + 1L
      out[[z]] <- tmp
    }
  }

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}



# ============================================================
# Main user function
# ============================================================
#
# germination_factorial_LRT(
#   data,
#   scale = c("probability", "logit"),
#   method = c("mle", "firth"),
#   alpha = 0.05,
#   posthoc = TRUE,
#   adjust = "holm",
#   simple_effects = c("both", "A_within_B", "B_within_A"),
#   cld = TRUE,
#   maxeval = 10000,
#   xtol_rel = 1e-10,
#   ftol_rel = 1e-12
# )
#
# Required columns:
#   Fac_A, Fac_B, g_ijk, n_ijk
#
# Interpretation:
#   scale = "probability"
#     Interaction is defined through constancy of absolute differences
#     in germination probability.
#
#   scale = "logit"
#     Interaction is defined through constancy of differences in log-odds.
#
#   method = "firth"
#     Available only with scale = "logit".
#
# The default call
#   germination_factorial_LRT(data)
# reproduces the probability-scale analysis used in the manuscript.
#
# ============================================================
# User-facing analysis function
# ============================================================

germination_factorial_LRT <- function(
    data,
    scale = c("probability", "logit"),
    method = c("mle", "firth"),
    alpha = 0.05,
    posthoc = TRUE,
    adjust = "holm",
    simple_effects = c("both", "A_within_B", "B_within_A"),
    cld = TRUE,
    maxeval = 10000,
    xtol_rel = 1e-10,
    ftol_rel = 1e-12
) {
  scale <- match.arg(scale)
  method <- match.arg(method)
  simple_effects <- match.arg(simple_effects)

  if (scale == "probability" && method != "mle") {
    stop(
      "method = 'firth' is available only for scale = 'logit'. ",
      "For probability-scale analysis, use method = 'mle'."
    )
  }

  if (scale == "probability") {
    overall <- .fit_overall_models(
      data,
      maxeval = maxeval,
      xtol_rel = xtol_rel,
      ftol_rel = ftol_rel
    )
    overall$tests$significant <- overall$tests$p_value < alpha

    pA <- overall$tests$p_value[overall$tests$Effect == "Factor A"]
    pB <- overall$tests$p_value[overall$tests$Effect == "Factor B"]
    pI <- overall$tests$p_value[overall$tests$Effect == "Interaction"]

    ph <- list()
    letters <- list()

    if (!posthoc) {
      rule <- "Post-hoc comparisons were not requested."

    } else if (pI >= alpha) {
      rule <- paste0(
        "Because the interaction was not significant at alpha = ",
        alpha,
        ", pairwise comparisons for significant main effects were ",
        "performed under the additive probability-scale model."
      )

      if (pA < alpha) {
        ph$Factor_A <- .posthoc_main_effect(
          data, overall, factor = "A",
          alpha = alpha, adjust = adjust,
          maxeval = maxeval,
          xtol_rel = xtol_rel,
          ftol_rel = ftol_rel
        )

        if (cld) {
          letters$Factor_A <- .cld_main_effect(
            overall, "A", ph$Factor_A, alpha
          )
        }
      }

      if (pB < alpha) {
        ph$Factor_B <- .posthoc_main_effect(
          data, overall, factor = "B",
          alpha = alpha, adjust = adjust,
          maxeval = maxeval,
          xtol_rel = xtol_rel,
          ftol_rel = ftol_rel
        )

        if (cld) {
          letters$Factor_B <- .cld_main_effect(
            overall, "B", ph$Factor_B, alpha
          )
        }
      }

    } else {
      rule <- paste0(
        "Because the interaction was significant at alpha = ",
        alpha,
        ", post-hoc comparisons were performed as simple-effect tests ",
        "rather than as overall main-effect comparisons."
      )

      if (simple_effects %in% c("both", "A_within_B")) {
        ph$A_within_B <- .simple_pairwise_LRT(
          data,
          target = "A_within_B",
          alpha = alpha,
          adjust = adjust
        )

        if (cld) {
          letters$A_within_B <- .cld_simple_effect(
            data,
            target = "A_within_B",
            pairwise_table = ph$A_within_B,
            alpha = alpha
          )
        }
      }

      if (simple_effects %in% c("both", "B_within_A")) {
        ph$B_within_A <- .simple_pairwise_LRT(
          data,
          target = "B_within_A",
          alpha = alpha,
          adjust = adjust
        )

        if (cld) {
          letters$B_within_A <- .cld_simple_effect(
            data,
            target = "B_within_A",
            pairwise_table = ph$B_within_A,
            alpha = alpha
          )
        }
      }
    }

    result <- list(
      overall_tests = overall$tests,
      fitted_probabilities = overall$fitted,
      logLik = overall$logLik,
      additive_fit = overall$additive_fit,
      posthoc_rule = rule,
      posthoc = ph,
      cld = letters,
      diagnostics = overall$additive_fit$diagnostics,
      scale = scale,
      method = method,
      interaction_null = paste(
        "Absolute differences in germination probability are constant",
        "across levels of the other factor."
      ),
      alpha = alpha,
      p_adjust_method = adjust,
      simple_effects = simple_effects
    )

  } else {
    if (method == "mle") {
      overall <- .fit_overall_logit_mle(data)
      ll_name <- "logLik"
    } else {
      overall <- .fit_overall_logit_firth(data)
      ll_name <- "penalized_logLik"
    }

    overall$tests$significant <- overall$tests$p_value < alpha

    pA <- overall$tests$p_value[overall$tests$Effect == "Factor A"]
    pB <- overall$tests$p_value[overall$tests$Effect == "Factor B"]
    pI <- overall$tests$p_value[overall$tests$Effect == "Interaction"]

    ph <- list()
    letters <- list()

    if (!posthoc) {
      rule <- "Post-hoc comparisons were not requested."

    } else if (pI >= alpha) {
      rule <- paste0(
        "Because the interaction was not significant at alpha = ",
        alpha,
        ", pairwise comparisons for significant main effects were ",
        "performed under the additive logit-scale model."
      )

      if (pA < alpha) {
        ph$Factor_A <- .posthoc_logit_main_effect(
          data,
          overall,
          factor = "A",
          method = method,
          alpha = alpha,
          adjust = adjust
        )

        if (cld) {
          letters$Factor_A <- .cld_logit_main_effect(
            overall, "A", ph$Factor_A, alpha
          )
        }
      }

      if (pB < alpha) {
        ph$Factor_B <- .posthoc_logit_main_effect(
          data,
          overall,
          factor = "B",
          method = method,
          alpha = alpha,
          adjust = adjust
        )

        if (cld) {
          letters$Factor_B <- .cld_logit_main_effect(
            overall, "B", ph$Factor_B, alpha
          )
        }
      }

    } else {
      rule <- paste0(
        "Because the interaction was significant at alpha = ",
        alpha,
        ", post-hoc comparisons were performed as logit-scale simple-effect ",
        "tests rather than as overall main-effect comparisons."
      )

      if (simple_effects %in% c("both", "A_within_B")) {
        ph$A_within_B <- .posthoc_logit_simple_effect(
          data,
          target = "A_within_B",
          method = method,
          alpha = alpha,
          adjust = adjust
        )

        if (cld) {
          letters$A_within_B <- .cld_logit_simple_effect(
            overall,
            target = "A_within_B",
            pairwise_table = ph$A_within_B,
            alpha = alpha
          )
        }
      }

      if (simple_effects %in% c("both", "B_within_A")) {
        ph$B_within_A <- .posthoc_logit_simple_effect(
          data,
          target = "B_within_A",
          method = method,
          alpha = alpha,
          adjust = adjust
        )

        if (cld) {
          letters$B_within_A <- .cld_logit_simple_effect(
            overall,
            target = "B_within_A",
            pairwise_table = ph$B_within_A,
            alpha = alpha
          )
        }
      }
    }

    result <- list(
      overall_tests = overall$tests,
      fitted_probabilities = overall$fitted,
      logLik = if (ll_name == "logLik") overall$logLik else NULL,
      penalized_logLik = if (ll_name == "penalized_logLik") overall$penalized_logLik else NULL,
      models = overall$models,
      posthoc_rule = rule,
      posthoc = ph,
      cld = letters,
      diagnostics = overall$diagnostics,
      scale = scale,
      method = method,
      interaction_null = paste(
        "Differences in log-odds of germination are constant",
        "across levels of the other factor."
      ),
      alpha = alpha,
      p_adjust_method = adjust,
      simple_effects = simple_effects
    )
  }

  structure(result, class = "germination_factorial_LRT")
}


# ============================================================
# Printing helpers
# ============================================================

.format_p_value <- function(p, digits = 4) {
  out <- rep(NA_character_, length(p))

  ok <- !is.na(p)
  if (!any(ok)) return(out)

  x <- p[ok]

  out_ok <- ifelse(
    x < .Machine$double.eps,
    paste0("< ", format(.Machine$double.eps, scientific = TRUE, digits = 3)),
    ifelse(
      x < 1e-4,
      format(x, scientific = TRUE, digits = digits),
      formatC(x, format = "f", digits = digits)
    )
  )

  out[ok] <- out_ok
  out
}


.prepare_table_for_print <- function(tab) {
  out <- tab

  if ("p_value" %in% names(out)) {
    out$p_value <- .format_p_value(out$p_value)
  }

  if ("p_adjusted" %in% names(out)) {
    out$p_adjusted <- .format_p_value(out$p_adjusted)
  }

  out
}


.print_posthoc_estimand_note <- function(x, nm) {
  if (identical(x$scale, "probability")) {
    if (nm %in% c("Factor_A", "Factor_B")) {
      cat(
        "Effect estimate: difference in fitted germination probability ",
        "under the additive probability-scale model.\n",
        sep = ""
      )
    } else {
      cat(
        "Effect estimate: within-stratum difference in germination probability.\n",
        sep = ""
      )
    }

  } else {
    if (nm %in% c("Factor_A", "Factor_B")) {
      cat(
        "Primary logit-scale estimates: log-odds difference and odds ratio.\n",
        "Average probability difference is reported only as a descriptive ",
        "probability-scale summary.\n",
        sep = ""
      )
    } else {
      cat(
        "Primary logit-scale estimates: within-stratum log-odds difference ",
        "and odds ratio.\n",
        "Observed probability difference is reported only as a descriptive ",
        "probability-scale summary.\n",
        sep = ""
      )
    }
  }
}


# ============================================================
# Print method
# ============================================================

print.germination_factorial_LRT <- function(x, ...) {
  cat("\nTwo-factor germination analysis\n")
  cat("================================\n\n")

  if (identical(x$scale, "probability")) {
    method_label <- "Binomial LRT with constrained probability-scale additivity"
  } else if (identical(x$method, "firth")) {
    method_label <- "Firth logistic regression with penalized LRTs"
  } else {
    method_label <- "Ordinary binomial logistic regression with LRTs"
  }

  cat("Interaction scale : ", x$scale, "\n", sep = "")
  cat("Method            : ", method_label, "\n", sep = "")
  cat("Interaction null  : ", x$interaction_null, "\n\n", sep = "")

  cat("Overall tests\n")
  print(
    .prepare_table_for_print(x$overall_tests),
    row.names = FALSE,
    digits = 6
  )

  cat("\nPost-hoc rule\n")
  cat(x$posthoc_rule, "\n")

  if (length(x$posthoc) > 0) {
    cat(
      "\nDetailed post-hoc tables are stored in $posthoc",
      " and compact letter displays in $cld.\n",
      sep = ""
    )
    cat(
      "Use summary(fit, detail = TRUE) for full post-hoc output.\n",
      sep = ""
    )
  }

  if (identical(x$scale, "logit") && identical(x$method, "mle")) {
    d <- x$diagnostics
    cat("\nLogistic-fit diagnostic\n")
    cat(
      "  Possible separation : ",
      d$possible_separation,
      "\n",
      sep = ""
    )
    if (isTRUE(d$possible_separation)) {
      cat(
        "  WARNING: Ordinary logistic MLE may be unstable; ",
        "consider method = 'firth'.\n",
        sep = ""
      )
    }
  }

  invisible(x)
}



# ============================================================
# Summary method
# ============================================================

summary.germination_factorial_LRT <- function(
    object,
    detail = FALSE,
    diagnostics = TRUE,
    ...
) {
  x <- object

  out <- list(
    scale = x$scale,
    method = x$method,
    interaction_null = x$interaction_null,
    overall_tests = x$overall_tests,
    posthoc_rule = x$posthoc_rule,
    posthoc = if (detail) x$posthoc else NULL,
    cld = if (detail) x$cld else NULL,
    diagnostics = if (diagnostics) x$diagnostics else NULL
  )

  class(out) <- "summary.germination_factorial_LRT"
  out
}


print.summary.germination_factorial_LRT <- function(x, ...) {
  cat("\nTwo-factor germination analysis summary\n")
  cat("=======================================\n\n")

  cat("Interaction scale : ", x$scale, "\n", sep = "")
  cat("Method            : ", x$method, "\n", sep = "")
  cat("Interaction null  : ", x$interaction_null, "\n\n", sep = "")

  cat("Overall tests\n")
  print(
    .prepare_table_for_print(x$overall_tests),
    row.names = FALSE,
    digits = 6
  )

  cat("\nPost-hoc rule\n")
  cat(x$posthoc_rule, "\n")

  if (!is.null(x$posthoc) && length(x$posthoc) > 0) {
    for (nm in names(x$posthoc)) {
      cat("\nPairwise comparisons: ", nm, "\n", sep = "")
      tmp_obj <- list(scale = x$scale)
      .print_posthoc_estimand_note(tmp_obj, nm)
      print(
        .prepare_table_for_print(x$posthoc[[nm]]),
        row.names = FALSE,
        digits = 6
      )

      if (!is.null(x$cld) && nm %in% names(x$cld)) {
        cat("\nCompact letter display: ", nm, "\n", sep = "")
        print(x$cld[[nm]], row.names = FALSE, digits = 6)
      }
    }
  }

  if (!is.null(x$diagnostics) && identical(x$scale, "logit") &&
      identical(x$method, "mle")) {
    d <- x$diagnostics
    cat("\nLogistic-fit diagnostic\n")
    cat(
      "  All models converged       : ",
      paste(names(d$converged), d$converged, sep = "=", collapse = ", "),
      "\n",
      sep = ""
    )
    cat("  All coefficients finite   : ", d$all_coefficients_finite, "\n", sep = "")
    cat("  Max |coefficient|         : ", format(d$max_abs_coefficient, digits = 6), "\n", sep = "")
    cat("  Max standard error        : ", format(d$max_standard_error, digits = 6), "\n", sep = "")
    cat("  Min fitted probability    : ", format(d$min_fitted_probability, scientific = TRUE, digits = 6), "\n", sep = "")
    cat("  Max fitted probability    : ", format(d$max_fitted_probability, scientific = TRUE, digits = 6), "\n", sep = "")
    cat("  Possible separation       : ", d$possible_separation, "\n", sep = "")
  }

  invisible(x)
}


# ============================================================
# Convenience extractor
# ============================================================

summary_tables <- function(fit) {
  if (!inherits(fit, "germination_factorial_LRT")) {
    stop("'fit' must be an object returned by germination_factorial_LRT().")
  }

  list(
    overall = fit$overall_tests,
    pairwise = fit$posthoc,
    cld = fit$cld
  )
}


# ============================================================
# EXAMPLES
# ============================================================

# dat <- read.csv("one_dataset.csv")

# 1. Probability-scale analysis (default; manuscript method)
# fit_prob <- germination_factorial_LRT(
#   data = dat,
#   scale = "probability",
#   method = "mle",
#   posthoc = TRUE
# )
# print(fit_prob)
# summary(fit_prob, detail = TRUE)

# 2. Ordinary logistic regression
# fit_logit <- germination_factorial_LRT(
#   data = dat,
#   scale = "logit",
#   method = "mle",
#   posthoc = TRUE
# )
# print(fit_logit)
# summary(fit_logit, detail = TRUE)

# 3. Firth logistic regression
# fit_firth <- germination_factorial_LRT(
#   data = dat,
#   scale = "logit",
#   method = "firth",
#   posthoc = TRUE
# )
# print(fit_firth)
# summary(fit_firth, detail = TRUE)

# 4. Direct access to result components
# fit_prob$overall_tests
# fit_prob$posthoc
# fit_prob$cld
# fit_logit$diagnostics

# 5. Convenience extractor
# summary_tables(fit_prob)
