# factorial_binomial_LRT.R
# GitHub implementation for:
# Assessing treatment interactions on the probability scale
# in factorial germination experiments
#
# Main user function:
#   germination_factorial_LRT()
#
# ------------------------------------------------------------
# Probability-scale Binomial LRT for two-factor germination data
# with automatic post-hoc testing and Holm adjustment.
#
# Main features
#   * Overall LRTs for Factor A, Factor B, and A x B interaction
#   * Constrained additive model:
#         p_ij = mu + alpha_i + beta_j
#         0 <= p_ij <= 1
#   * R implementation uses nloptr::NLOPT_LD_SLSQP
#   * Multiple starting values; best feasible likelihood retained
#   * If interaction is not significant:
#         pairwise main-effect LRTs for significant factors
#   * If interaction is significant:
#         simple-effect pairwise LRTs
#   * Holm multiplicity adjustment
#   * Compact letter display (CLD) via multcompView
#
# Required packages
#   install.packages(c("nloptr", "multcompView"))
#
# Required data columns
#   Fac_A, Fac_B, g_ijk, n_ijk
#
# Fac_A and Fac_B may be numeric, character, or factor variables.
# Original level labels are preserved in the returned output.
# For factor variables, the declared factor-level order is retained.
# For numeric variables, levels are ordered numerically; for character
# variables, the order of first appearance is retained.
# ------------------------------------------------------------

if (!requireNamespace("nloptr", quietly = TRUE)) {
  stop("Package 'nloptr' is required. Install it with install.packages('nloptr').")
}

if (!requireNamespace("multcompView", quietly = TRUE)) {
  stop("Package 'multcompView' is required. Install it with install.packages('multcompView').")
}



# ============================================================
# Original factor-label handling
# ============================================================

.user_levels <- function(x) {
  if (is.factor(x)) {
    return(levels(droplevels(x)))
  }

  if (is.numeric(x) || is.integer(x)) {
    return(sort(unique(x)))
  }

  unique(as.character(x))
}


.encode_user_data <- function(data) {
  req <- c("Fac_A", "Fac_B", "g_ijk", "n_ijk")
  miss <- setdiff(req, names(data))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }

  A_levels <- .user_levels(data$Fac_A)
  B_levels <- .user_levels(data$Fac_B)

  # Use character matching so that numeric labels such as 2016/2017
  # and character labels are handled by the same machinery.
  A_chr <- as.character(data$Fac_A)
  B_chr <- as.character(data$Fac_B)

  A_key <- as.character(A_levels)
  B_key <- as.character(B_levels)

  encoded <- data
  encoded$Fac_A <- match(A_chr, A_key)
  encoded$Fac_B <- match(B_chr, B_key)

  if (anyNA(encoded$Fac_A) || anyNA(encoded$Fac_B)) {
    stop("Could not encode one or more factor levels.")
  }

  list(
    data = encoded,
    A_levels = as.character(A_levels),
    B_levels = as.character(B_levels)
  )
}


.relabel_probability_matrix <- function(mat, A_levels, B_levels) {
  rownames(mat) <- A_levels
  colnames(mat) <- B_levels
  mat
}


.relabel_main_pairwise <- function(tab, factor, A_levels, B_levels) {
  if (is.null(tab) || nrow(tab) == 0) return(tab)

  lev <- if (factor == "A") A_levels else B_levels

  tab$Level1 <- lev[as.integer(tab$Level1)]
  tab$Level2 <- lev[as.integer(tab$Level2)]
  tab$Contrast <- paste(tab$Level1, tab$Level2, sep = " - ")
  tab
}


.relabel_simple_pairwise <- function(tab, target, A_levels, B_levels) {
  if (is.null(tab) || nrow(tab) == 0) return(tab)

  if (target == "A_within_B") {
    tab$Fixed_level <- B_levels[as.integer(tab$Fixed_level)]
    tab$Level1 <- A_levels[as.integer(tab$Level1)]
    tab$Level2 <- A_levels[as.integer(tab$Level2)]
    tab$Contrast <- paste(tab$Level1, tab$Level2, sep = " - ")
    tab$Simple_effect <- paste0("A within B=", tab$Fixed_level)
  } else {
    tab$Fixed_level <- A_levels[as.integer(tab$Fixed_level)]
    tab$Level1 <- B_levels[as.integer(tab$Level1)]
    tab$Level2 <- B_levels[as.integer(tab$Level2)]
    tab$Contrast <- paste(tab$Level1, tab$Level2, sep = " - ")
    tab$Simple_effect <- paste0("B within A=", tab$Fixed_level)
  }

  tab
}


.relabel_main_cld <- function(tab, factor, A_levels, B_levels) {
  if (is.null(tab) || nrow(tab) == 0) return(tab)

  lev <- if (factor == "A") A_levels else B_levels

  # Internal CLD levels are A1,A2,... or B1,B2,...
  idx <- suppressWarnings(as.integer(sub("^[AB]", "", tab$Level)))
  tab$Level <- lev[idx]

  if ("Average_fitted_probability" %in% names(tab)) {
    names(tab)[names(tab) == "Average_fitted_probability"] <-
      "Mean_fitted_probability"
  }

  tab
}


.relabel_simple_cld <- function(tab, target, A_levels, B_levels) {
  if (is.null(tab) || nrow(tab) == 0) return(tab)

  if (target == "A_within_B") {
    idx <- suppressWarnings(as.integer(sub("^A", "", tab$Level)))
    tab$Level <- A_levels[idx]
    tab$Fixed_level <- B_levels[as.integer(tab$Fixed_level)]
  } else {
    idx <- suppressWarnings(as.integer(sub("^B", "", tab$Level)))
    tab$Level <- B_levels[idx]
    tab$Fixed_level <- A_levels[as.integer(tab$Fixed_level)]
  }

  tab
}


.relabel_fit_output <- function(overall, A_levels, B_levels) {
  overall$additive_fit$p_matrix <-
    .relabel_probability_matrix(
      overall$additive_fit$p_matrix,
      A_levels,
      B_levels
    )

  overall$fitted$M_AplusB <-
    .relabel_probability_matrix(
      overall$fitted$M_AplusB,
      A_levels,
      B_levels
    )

  overall$fitted$M_AxB <-
    .relabel_probability_matrix(
      overall$fitted$M_AxB,
      A_levels,
      B_levels
    )

  names(overall$fitted$M_A) <- A_levels
  names(overall$fitted$M_B) <- B_levels

  overall
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
# User-facing analysis function
# ============================================================

germination_factorial_LRT <- function(
    data,
    alpha = 0.05,
    posthoc = TRUE,
    adjust = "holm",
    simple_effects = c("both", "A_within_B", "B_within_A"),
    cld = TRUE,
    maxeval = 10000,
    xtol_rel = 1e-10,
    ftol_rel = 1e-12
) {
  simple_effects <- match.arg(simple_effects)

  # Preserve the user's original factor labels and work internally
  # with consecutive integer codes.
  enc <- .encode_user_data(data)
  dat <- enc$data
  A_levels <- enc$A_levels
  B_levels <- enc$B_levels

  overall <- .fit_overall_models(
    dat,
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
      tmpA <- .posthoc_main_effect(
        dat, overall, factor = "A",
        alpha = alpha, adjust = adjust,
        maxeval = maxeval,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel
      )

      if (cld) {
        letters$Factor_A <- .cld_main_effect(
          overall, "A", tmpA, alpha
        )
      }

      ph$Factor_A <- .relabel_main_pairwise(
        tmpA, "A", A_levels, B_levels
      )

      if (cld) {
        letters$Factor_A <- .relabel_main_cld(
          letters$Factor_A, "A", A_levels, B_levels
        )
      }
    }

    if (pB < alpha) {
      tmpB <- .posthoc_main_effect(
        dat, overall, factor = "B",
        alpha = alpha, adjust = adjust,
        maxeval = maxeval,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel
      )

      if (cld) {
        letters$Factor_B <- .cld_main_effect(
          overall, "B", tmpB, alpha
        )
      }

      ph$Factor_B <- .relabel_main_pairwise(
        tmpB, "B", A_levels, B_levels
      )

      if (cld) {
        letters$Factor_B <- .relabel_main_cld(
          letters$Factor_B, "B", A_levels, B_levels
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
      tmpAB <- .simple_pairwise_LRT(
        dat,
        target = "A_within_B",
        alpha = alpha,
        adjust = adjust
      )

      if (cld) {
        letters$A_within_B <- .cld_simple_effect(
          dat,
          target = "A_within_B",
          pairwise_table = tmpAB,
          alpha = alpha
        )
      }

      ph$A_within_B <- .relabel_simple_pairwise(
        tmpAB, "A_within_B", A_levels, B_levels
      )

      if (cld) {
        letters$A_within_B <- .relabel_simple_cld(
          letters$A_within_B,
          "A_within_B",
          A_levels,
          B_levels
        )
      }
    }

    if (simple_effects %in% c("both", "B_within_A")) {
      tmpBA <- .simple_pairwise_LRT(
        dat,
        target = "B_within_A",
        alpha = alpha,
        adjust = adjust
      )

      if (cld) {
        letters$B_within_A <- .cld_simple_effect(
          dat,
          target = "B_within_A",
          pairwise_table = tmpBA,
          alpha = alpha
        )
      }

      ph$B_within_A <- .relabel_simple_pairwise(
        tmpBA, "B_within_A", A_levels, B_levels
      )

      if (cld) {
        letters$B_within_A <- .relabel_simple_cld(
          letters$B_within_A,
          "B_within_A",
          A_levels,
          B_levels
        )
      }
    }
  }

  # Relabel fitted-probability output after all computations that
  # require the internal A1/A2/... and B1/B2/... naming convention.
  overall <- .relabel_fit_output(overall, A_levels, B_levels)

  structure(
    list(
      overall_tests = overall$tests,
      fitted_probabilities = overall$fitted,
      logLik = overall$logLik,
      additive_fit = overall$additive_fit,
      posthoc_rule = rule,
      posthoc = ph,
      cld = letters,
      alpha = alpha,
      p_adjust_method = adjust,
      simple_effects = simple_effects,
      level_labels = list(
        Factor_A = A_levels,
        Factor_B = B_levels
      )
    ),
    class = "germination_factorial_LRT"
  )
}


# ============================================================
# Print method
# ============================================================

print.germination_factorial_LRT <- function(x, ...) {
  cat("\nProbability-scale Binomial LRT for a two-factor germination experiment\n")
  cat("=====================================================================\n\n")

  cat("Overall tests\n")
  print(x$overall_tests, row.names = FALSE, digits = 6)

  cat("\nPost-hoc rule\n")
  cat(x$posthoc_rule, "\n")

  if (length(x$posthoc) > 0) {
    for (nm in names(x$posthoc)) {
      cat("\nPairwise comparisons: ", nm, "\n", sep = "")
      print(x$posthoc[[nm]], row.names = FALSE, digits = 6)

      if (nm %in% names(x$cld)) {
        cat("\nCompact letter display: ", nm, "\n", sep = "")
        print(x$cld[[nm]], row.names = FALSE, digits = 6)
      }
    }
  }

  cat(
    "\nMultiplicity adjustment: ",
    x$p_adjust_method,
    "\n",
    sep = ""
  )

  if ("A_within_B" %in% names(x$posthoc) ||
      "B_within_A" %in% names(x$posthoc)) {
    cat(
      "For simple effects, multiplicity adjustment is applied separately ",
      "within each fixed level of the conditioning factor.\n",
      sep = ""
    )
  }

  cat(
    "CLD note: levels sharing a letter are not significantly different ",
    "according to the adjusted pairwise tests at the selected alpha level.\n",
    sep = ""
  )

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

# Example 1: default automatic workflow
#
# dat <- read.csv("one_dataset.csv")
# # Fac_A and Fac_B may contain their original labels, e.g.
# # Fac_A = c("N4", "N14", ...), Fac_B = c(2016, 2017).
#
# fit <- germination_factorial_LRT(
#   data = dat,
#   alpha = 0.05,
#   posthoc = TRUE,
#   adjust = "holm",
#   simple_effects = "both",
#   cld = TRUE
# )
#
# print(fit)


# Example 2: if interaction is significant and only A-within-B
# simple effects are scientifically relevant
#
# fit <- germination_factorial_LRT(
#   data = dat,
#   simple_effects = "A_within_B"
# )


# Example 3: extract results for tables or export
#
# fit$overall_tests
# fit$posthoc$Factor_A
# fit$posthoc$Factor_B
# fit$posthoc$A_within_B
# fit$posthoc$B_within_A
# fit$cld$Factor_A
# fit$cld$Factor_B
# fit$cld$A_within_B
# fit$cld$B_within_A
#
# write.csv(fit$overall_tests, "overall_tests.csv", row.names = FALSE)
