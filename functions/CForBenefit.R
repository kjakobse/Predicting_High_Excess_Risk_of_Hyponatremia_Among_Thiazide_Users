# Using model-based c-for-benefit, we can provide tau_hat predicted from model trained on training data and p_0 and p_1
# from model trained on test data to obtain a valid estimate of discriminative performance of the training data model.
# We can further estimate discriminative performance on a validation data set by providing data from such a set, using
# the training set model to obtain predictions. Further, Y_hat and W_hat should be derived from the validation sample.
# p_0 and p_1 in a validation sample is computed based entirely on a model trained on the validation sample, 
# provided using validation_forest argument.
CForBenefit <- function(forest,
                        validation_forest = NULL,
                        match = c("none", "covariates", "CATE", "control_risk", "treated_risk", 
                                  "CATE_control_risk", "CATE_treated_risk", "combined_outcome_risk"),
                        match_method = "nearest",
                        match_distance = "mahalanobis",
                        tau_hat_method = c("risk_diff", "tau_avg", "tau_treated", "tau_control"),
                        CI = c("simple", "bootstrap", "none"),
                        level = 0.95,
                        n_quantiles = 4,
                        n_bootstraps = 999L,
                        bootstrap_honesty = TRUE,
                        n_cores = 1L,
                        time_limit = Inf,
                        time_limit_CI = Inf,
                        mem_limit = NULL,
                        verbose = TRUE,
                        compile = TRUE,
                        ...) {
  # compile c-for-benefit functions
  stopifnot(
    "compile must be boolean (TRUE or FALSE)" =
      isTRUE(compile) || isFALSE(compile)
  )
  if (compile) {
  Rcpp::sourceCpp("functions/mbcb.cpp")
  Rcpp::sourceCpp("functions/rcorr.cpp")
  }
  # ensure causal forest
  stopifnot(
    "forest must be a causal_forest object" =
      inherits(forest, c("causal_forest", "causal_survival_forest"))
  )
  stopifnot(
    "validation_forest must be NULL or a causal_forest object" =
      is.null(validation_forest) || 
      inherits(validation_forest, c("causal_forest", "causal_survival_forest"))
  )
  
  
  # Compute expected potential outcomes in training and validation samples
  p_0 <- as.numeric(forest$Y.hat - forest$W.hat * forest$predictions)
  p_1 <- as.numeric(forest$Y.hat + (1 - forest$W.hat) * forest$predictions)
  if (!is.null(validation_forest)) {
    val_p_0 <- as.numeric(validation_forest$Y.hat - validation_forest$W.hat * validation_forest$predictions)
    val_p_1 <- as.numeric(validation_forest$Y.hat + (1 - validation_forest$W.hat) * validation_forest$predictions)
  }
  
  subclass <- NULL
  
  # ensure correct data types
  stopifnot(
    "level must be a scalar between 0 and 1" = 
      is.numeric(level) && length(level) == 1 && 0 < level && level < 1
  )
  stopifnot(
    "n_bootstraps must be a positive integer" = 
      is.numeric(n_bootstraps) && length(n_bootstraps) == 1 && 1 <= n_bootstraps
  )
  stopifnot(
    "n_cores must be a positive integer" = 
      is.numeric(n_cores) && length(n_cores) == 1 && 1 <= n_cores
  )
  stopifnot("match_method must be a character" = is.character(match_method))
  stopifnot("match_distance must be a character" = is.character(match_distance))
  stopifnot(
    "verbose must be boolean (TRUE or FALSE)" =
      isTRUE(verbose) || isFALSE(verbose)
  )
  match <- match.arg(match)
  tau_hat_method <- match.arg(tau_hat_method)
  CI <- match.arg(CI)
  if (match == "control_risk" && tau_hat_method != "tau_treated" ||
      match != "control_risk" && tau_hat_method == "tau_treated") {
    warning(paste0(
      "To obtain valid inference of C-for-benefit, matching on control outcome risk (match = 'control_risk')\n", 
      "and concordance evaluation using estimated effect on the treated (tau_hat_method = 'tau_treated')\n",
      "should be used together."
    ), immediate. = TRUE)
  }
  if (match == "treated_risk" && tau_hat_method != "tau_control" ||
      match != "treated_risk" && tau_hat_method == "tau_control") {
    warning(paste0(
      "To obtain valid inference of C-for-benefit, matching on treated outcome risk (match = 'treated_risk')\n", 
      "and concordance evaluation using estimated effect on the control (tau_hat_method = 'tau_control')\n",
      "should be used together."
    ), immediate. = TRUE)
  }
  
  # combine data into a tibble
  if (is.null(validation_forest)) {
    data_tbl <- dplyr::tibble(
      match_id = seq_along(forest$Y.orig),
      X = forest$X,
      p_0 = p_0,
      p_1 = p_1,
      tau_hat = forest$predictions[, 1]
    )
    data_tbl <- dplyr::mutate(
      data_tbl,
      W = forest$W.orig,
      Y = forest$Y.orig
    )
    if ("sample.weights" %in% names(forest)) {
      data_tbl <- dplyr::mutate(data_tbl, sample_weights = forest$sample.weights)
    } else {
      data_tbl <- dplyr::mutate(data_tbl, sample_weights = 1)
    }
  } else {
    data_tbl <- dplyr::tibble(
      match_id = seq_along(validation_forest$Y.orig),
      X = validation_forest$X,
      p_0 = val_p_0,
      p_1 = val_p_1,
      tau_hat = validation_forest$predictions[, 1]
    )
    data_tbl <- dplyr::mutate(
      data_tbl,
      W = validation_forest$W.orig,
      Y = validation_forest$Y.orig
    )
    if ("sample.weights" %in% names(validation_forest)) {
      data_tbl <- dplyr::mutate(data_tbl, sample_weights = validation_forest$sample.weights)
    } else {
      data_tbl <- dplyr::mutate(data_tbl, sample_weights = 1)
    }
  }
  
  # create diagnostic tables
  # 2x2 counts tables for W and Y within quantiles of tau_hat
  quantile_counts <- data_tbl |>
    dplyr::mutate(
      quantile = cut(
        tau_hat, 
        breaks = quantile(tau_hat, seq(0, 1, length.out = n_quantiles + 1)), 
        labels = seq_len(n_quantiles),
        include.lowest = TRUE
      )
    ) |>
    dplyr::group_by(quantile) |> 
    dplyr::count(W, Y) |>
    dplyr::ungroup()
  
  # probability of each level of observed benefit in each quantile
  qbf <- quantile_benefit_probabilities <- quantile_counts |>
    dplyr::group_by(quantile) |>
    dplyr::summarise(
      p_obs_harm = sum((1-W)*(1-Y)*n)/(sum((1-W)*(1-Y)*n)+sum((1-W)*Y*n)) * sum(W*Y*n)/(sum(W*Y*n)+sum(W*(1-Y)*n)),
      p_obs_no_effect = sum((1-W)*(1-Y)*n)/(sum((1-W)*(1-Y)*n)+sum((1-W)*Y*n)) * sum(W*(1-Y)*n)/(sum(W*Y*n)+sum(W*(1-Y)*n)) +
        sum((1-W)*Y*n)/(sum((1-W)*(1-Y)*n)+sum((1-W)*Y*n)) * sum(W*Y*n)/(sum(W*Y*n)+sum(W*(1-Y)*n)),
      p_obs_benefit = 1 - p_obs_harm - p_obs_no_effect,
      #sum((1-W)*Y*n)/(sum((1-W)*Y*n)+sum((1-W)*(1-Y)*n)) * sum(W*(1-Y)*n)/(sum(W*Y*n)+sum(W*(1-Y)*n)),
      diff_harm_benefit = p_obs_harm - p_obs_benefit,
      #sum(W*Y*n)/(sum(W*Y*n)+sum(W*(1-Y)*n)) - sum((1-W)*Y*n)/(sum((1-W)*Y*n)+sum((1-W)*(1-Y)*n)),
      ratio_harm_benefit = p_obs_harm / p_obs_benefit
      #(sum((1-W)*Y*n)+sum((1-W)*(1-Y)*n))*sum(W*Y*n) / (sum((1-W)*Y*n)*(sum(W*Y*n)+sum(W*(1-Y)*n)))
    )
  
  # Probability of observing concordance among pairs of pairs from different quantiles
  nrow_qbp <- nrow(qbf)
  quantile_cfb <- dplyr::tibble(
    pair = vector("character", (nrow_qbp - 1) * nrow_qbp / 2),
    cfb = vector("numeric", (nrow_qbp - 1) * nrow_qbp / 2)
  )
  k <- 1
  for (i in seq(nrow(qbf), 2, -1)) {
    for (j in seq(i - 1, 1, -1)) {
      quantile_cfb[k, "pair"] <- paste("quantile", j, "+", i)
      quantile_cfb[k, "cfb"] <- (
        qbf$p_obs_harm[i] * qbf$p_obs_no_effect[j] +
          qbf$p_obs_harm[i] * qbf$p_obs_benefit[j] + 
          qbf$p_obs_no_effect[i]*qbf$p_obs_benefit[j]
      ) / (
        1 - (
          qbf$p_obs_harm[i] * qbf$p_obs_harm[j] +
            qbf$p_obs_no_effect[i] * qbf$p_obs_no_effect[j] +
            qbf$p_obs_benefit[i] * qbf$p_obs_benefit[j]
        )
      )
      k <- k + 1
    }
  }
  
  
  # set time limit on calculating C-for-benefit
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  })
  setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)
  
  # model based c-for-benefit (no matching)
  if (match == "none") {
    # c-for-benefit
    if (is.null(validation_forest)) {
      cindex <- mbcb(p_0, p_1, forest$predictions, forest$sample.weights)
    } else {
      tau_hat <- predict(forest, newdata = validation_forest$X, num.threads = n_cores)
      cindex <- mbcb(val_p_0, val_p_1, tau_hat, validation_forest$sample.weights)
    }
    c_for_benefit <- cindex$`C Index`
    c_for_benefit_se <- cindex$S.D. / 2
    
    if (CI == "simple") { # Simple confidence interval
      lower_CI <- c_for_benefit + qnorm(0.5 - level / 2) * c_for_benefit_se
      upper_CI <- c_for_benefit + qnorm(0.5 + level / 2) * c_for_benefit_se
    } else if (CI == "bootstrap") { # bootstrap confidence intervals
      CB_for_CI <- c()
      B <- 1
      setTimeLimit(cpu = time_limit_CI, elapsed = time_limit_CI, transient = TRUE)
      if (n_cores == 1L) {
        if (verbose) cli::cli_progress_bar("Bootstrapping", total = n_bootstraps)
        while (B <= n_bootstraps) {
          tryCatch(
            {
              # sample bootstrap indices with replacement
              boot_index <- sample(length(forest$Y.orig), length(forest$Y.orig), replace = TRUE)
              
              # fit cf model to bootstrap sample
              boot_cf <- grf::causal_forest(
                X = dplyr::slice(forest$X.orig, boot_index),
                Y = forest$Y.orig[boot_index],
                W = forest$W.orig[boot_index],
                Y.hat = forest$Y.hat[boot_index],
                W.hat = forest$W.hat[boot_index],
                num.trees = forest$`_num_trees`,
                sample.weights = forest$sample.weights[boot_index],
                clusters = forest$clusters,
                equalize.cluster.weights = forest$equalize.cluster.weights,
                sample.fraction = forest$tunable.params$sample.fraction,
                mtry = forest$tunable.params$mtry,
                min.node.size = forest$tunable.params$min.node.size,
                honesty = bootstrap_honesty,
                honesty.fraction = forest$tunable.params$honesty.fraction,
                honesty.prune.leaves = forest$tunable.params$honesty.prune.leaves,
                alpha = forest$tunable.params$alpha,
                imbalance.penalty = forest$tunable.params$imbalance.penalty,
                ci.group.size = forest$`_ci_group_size`,
                num.threads = n_cores
              )
              
              # extract predicted treatment effect and expected potential outcomes
              if (is.null(validation_forest)) {
                if (exists(forest$sample.weights)) {
                  sample_weights_boot <- forest$sample.weights
                } else {
                  sample_weights_boot <- rep(1, length(forest$Y.orig))
                }
                tau_hat_boot <- as.numeric(boot_cf$predictions)
                p_0_boot <- as.numeric(forest$Y.hat[boot_index] - forest$W.hat[boot_index] * tau_hat_boot)
                p_1_boot <- as.numeric(forest$Y.hat[boot_index] + (1 - forest$W.hat[boot_index]) * tau_hat_boot)
              } else {
                if (exists(validation_forest$sample.weights)) {
                  sample_weights_boot <- validation_forest$sample.weights
                } else {
                  sample_weights_boot <- rep(1, length(validation_forest$Y.orig))
                }
                tau_hat_boot <- predict(boot_cf, newdata = validation_forest$X, num.threads = n_cores)
                p_0_boot <- as.numeric(validation_forest$Y.hat - validation_forest$W.hat * validation_forest$predictions)
                p_1_boot <- as.numeric(validation_forest$Y.hat + (1 - validation_forest$W.hat) * validation_forest$predictions)
              }
              
              # C-for-benefit in bootstrap sample
              CB_for_CI[B] <- mbcb(p_0_boot, p_1_boot, tau_hat_boot, sample_weights_boot)$`C Index`
              
              # increment bootstrap index
              B <- B + 1
              if (verbose) cli::cli_progress_update()
            },
            error = function(e) {
              if (
                grepl(
                  "reached elapsed time limit|reached CPU time limit",
                  e$message
                )
              ) {
                input <- readline(
                  "Time limit reached. Do you want execution to continue (y/n) "
                )
                if (input %in% c("y", "n")) {
                  if (input == "n") stop("Time limit reached, execution stopped.")
                } else {
                  input <- readline(
                    "Please input either 'y' to continue or 'n' to stop execution. "
                  )
                  if (input %in% c("y", "n")) {
                    if (input == "n") stop("Time limit reached, execution stopped.")
                  } else {
                    stop("Answer not 'y' or 'n'. Execution stopped.")
                  }
                }
              } else {
                stop(e)
              }
            }
          )
        }
      } else {
        future::plan(multisession, workers = n_cores)
        CB_for_CI <- furrr::future_map_dbl(
          .x = seq_len(n_bootstraps),
          .options = furrr_options(
            packages = c("dplyr", "grf", "Rcpp"),
            globals = c("forest", "validation_forest")
          ),
          .progress = TRUE,
          .f = \(B) {
            Rcpp::sourceCpp("functions/mbcb.cpp")
            
            # sample bootstrap indices with replacement
            boot_index <- sample(length(forest$Y.orig), length(forest$Y.orig), replace = TRUE)
            
            # fit cf model to bootstrap sample
            boot_cf <- grf::causal_forest(
              X = dplyr::slice(forest$X.orig, boot_index),
              Y = forest$Y.orig[boot_index],
              W = forest$W.orig[boot_index],
              Y.hat = forest$Y.hat[boot_index],
              W.hat = forest$W.hat[boot_index],
              num.trees = forest$`_num_trees`,
              sample.weights = forest$sample.weights[boot_index],
              clusters = forest$clusters,
              equalize.cluster.weights = forest$equalize.cluster.weights,
              sample.fraction = forest$tunable.params$sample.fraction,
              mtry = forest$tunable.params$mtry,
              min.node.size = forest$tunable.params$min.node.size,
              honesty = bootstrap_honesty,
              honesty.fraction = forest$tunable.params$honesty.fraction,
              honesty.prune.leaves = forest$tunable.params$honesty.prune.leaves,
              alpha = forest$tunable.params$alpha,
              imbalance.penalty = forest$tunable.params$imbalance.penalty,
              ci.group.size = forest$`_ci_group_size`,
              num.threads = 1L
            )
            
            # extract predicted treatment effect and expected potential outcomes
            if (is.null(validation_forest)) {
              if (exists(forest$sample.weights)) {
                sample_weights_boot <- forest$sample.weights
              } else {
                sample_weights_boot <- rep(1, length(forest$Y.orig))
              }
              tau_hat_boot <- as.numeric(boot_cf$predictions)
              p_0_boot <- as.numeric(forest$Y.hat[boot_index] - forest$W.hat[boot_index] * tau_hat_boot)
              p_1_boot <- as.numeric(forest$Y.hat[boot_index] + (1 - forest$W.hat[boot_index]) * tau_hat_boot)
            } else {
              if (exists(validation_forest$sample.weights)) {
                sample_weights_boot <- validation_forest$sample.weights
              } else {
                sample_weights_boot <- rep(1, length(validation_forest$Y.orig))
              }
              tau_hat_boot <- predict(boot_cf, newdata = validation_forest$X, num.threads = n_cores)
              p_0_boot <- as.numeric(validation_forest$Y.hat - validation_forest$W.hat * validation_forest$predictions)
              p_1_boot <- as.numeric(validation_forest$Y.hat + (1 - validation_forest$W.hat) * validation_forest$predictions)
            }
            
            # C-for-benefit in bootstrap sample
            return(mbcb(p_0_boot, p_1_boot, tau_hat_boot, sample_weights_boot)$`C Index`)
          }
        )
        future::plan(sequential)
      }
      lower_CI <- as.numeric(quantile(CB_for_CI, 0.5 - level / 2))
      upper_CI <- as.numeric(quantile(CB_for_CI, 0.5 + level / 2))
    } else {
      lower_CI <- NA
      upper_CI <- NA
    } 
    
    # return results
    return(list(
      c_for_benefit = c_for_benefit,
      c_for_benefit_se = c_for_benefit_se,
      lower_CI = lower_CI,
      upper_CI = upper_CI,
      cfb_tbl = dplyr::tibble(
        "ci_{0.5 - level / 2}" := lower_CI,
        "estimate" = c_for_benefit,
        "ci_{0.5 + level / 2}" := upper_CI
      ),
      diagnostic_quantiles = list(
        quantile_counts = quantile_counts,
        quantile_benefit_probabilities = quantile_benefit_probabilities,
        quantile_c_for_benefit = quantile_cfb
      )
    ))
  }
  
  # each patient can only match to one patient from the other treatment arm
  if (!("estimand" %in% ...names())) {
    # ATT: all treated patients are matched to a control patient
    # ATC: all control patients are matched to a treated patient
    if (sum(data_tbl$W == 1) <= sum(data_tbl$W == 0)) estimand_method <- "ATT" else estimand_method <- "ATC"
  } else {
    estimand_method <- dots_list(...)$estimand
  }
  
  # match on covariates, predicted effect, control risk, or treated risk
  if (match == "covariates") {
    if ("estimand" %in% ...names()) {
      matched <- MatchIt::matchit(
        DF2formula(bind_cols(tibble(W = data_tbl$W), as_tibble(data_tbl$X))),
        data = unnest(data_tbl, "X"),
        method = match_method,
        distance = match_distance,
        ...
      )
    } else {
      matched <- MatchIt::matchit(
        DF2formula(bind_cols(tibble(W = data_tbl$W), as_tibble(data_tbl$X))),
        data = unnest(data_tbl, "X"),
        method = match_method,
        distance = match_distance,
        estimand = estimand_method,
        ...
      )
    }
  } else if (match == "CATE") {
    if ("estimand" %in% ...names()) {
      matched <- MatchIt::matchit(
        W ~ tau_hat,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        ...
      )
    } else {
      matched <- MatchIt::matchit(
        W ~ tau_hat,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        estimand = estimand_method,
        ...
      )
    }
  } else if (match == "control_risk") {
    if ("estimand" %in% ...names()) {
      matched <- MatchIt::matchit(
        W ~ p_0,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        ...
      )
    } else {
      matched <- MatchIt::matchit(
        W ~ p_0,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        estimand = estimand_method,
        ...
      )
    }
  } else if (match == "treated_risk") {
    if ("estimand" %in% ...names()) {
      matched <- MatchIt::matchit(
        W ~ p_1,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        ...
      )
    } else {
      matched <- MatchIt::matchit(
        W ~ p_1,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        estimand = estimand_method,
        ...
      )
    }
  } else if (match == "CATE_control_risk") {
      if ("estimand" %in% ...names()) {
        matched <- MatchIt::matchit(
          W ~ p_0 + tau_hat,
          data = data_tbl,
          method = match_method,
          distance = match_distance,
          ...
        )
      } else {
        matched <- MatchIt::matchit(
          W ~ p_0 + tau_hat,
          data = data_tbl,
          method = match_method,
          distance = match_distance,
          estimand = estimand_method,
          ...
        )
      }
  } else if (match == "CATE_treated_risk") {
    if ("estimand" %in% ...names()) {
      matched <- MatchIt::matchit(
        W ~ p_1 + tau_hat,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        ...
      )
    } else {
      matched <- MatchIt::matchit(
        W ~ p_1 + tau_hat,
        data = data_tbl,
        method = match_method,
        distance = match_distance,
        estimand = estimand_method,
        ...
      )
    }
  }
  if (match != "combined_outcome_risk") {
    matched_patients <- MatchIt::match.data(matched)
    matched_patients$subclass <- as.numeric(matched_patients$subclass)
    matched_patients <- dplyr::as_tibble(matched_patients)
  }
  
  if(match == "combined_outcome_risk") {
    args <- rlang::dots_list(...)
    if(!is.null(args$estimand)) args$estimand <- NULL
    matched_trt <- rlang::inject(MatchIt::matchit(
      W ~ p_0,
      data = data_tbl,
      method = match_method,
      distance = match_distance,
      estimand = "ATT",
      !!!args
    ))
    matched_ctr <- rlang::inject(MatchIt::matchit(
      W ~ p_1,
      data = data_tbl,
      method = match_method,
      distance = match_distance,
      estimand = "ATC",
      !!!args
    ))
    matched_patients_trt <- data_tbl |>
      filter(matched_trt$treat == 1) |>
      mutate(subclass = seq_len(n())) |>
      bind_rows(
        data_tbl |>
          slice(as.integer(matched_trt$match.matrix[, 1])) |>
          mutate(subclass = seq_len(n()))
      ) |>
      arrange(subclass, W) |>
      mutate(
        sample_weights = sum((W == 1) * sample_weights),
        .by = subclass
        ) |>
      mutate(
        Y = Y - p_0,
        weights = 1
      )
    matched_patients_ctr <- data_tbl |>
      filter(matched_ctr$treat == 0) |>
      mutate(subclass = seq_len(n()) + sum(data_tbl$W == 1)) |>
      bind_rows(
        data_tbl |>
          slice(as.integer(matched_ctr$match.matrix[, 1])) |>
          mutate(subclass = seq_len(n()) + sum(data_tbl$W == 1))
      ) |>
      arrange(subclass, W) |>
      mutate(
        sample_weights = sum((W == 0) * sample_weights),
        .by = subclass
      ) |>
      mutate(
        Y = Y - p_1,
        weights = 1
      )
    matched_patients <- bind_rows(
      matched_patients_trt,
      matched_patients_ctr
    )
  }

  # sort on subclass and W
  matched_patients <- matched_patients |>
    dplyr::arrange(subclass, W)
  
  # matched observed treatment effect
  observed_te <- matched_patients |>
    dplyr::select(subclass, Y, W, weights) |>
    dplyr::summarise(
      Y = sum(weights * Y) / sum(weights),
      .by = c(subclass, W)
    ) |>
    dplyr::summarise(
      matched_tau_obs = sum((W == 1) * Y - (W == 0) * Y),
      .by = subclass
    )
  matched_patients <- matched_patients |>
    inner_join(observed_te, by = "subclass")
  
  # matched predicted treatment effect
  if (tau_hat_method == "risk_diff") {
    # matched p_0 = P(Y = 1| W = 0)
    # matched p_1 = P(Y = 1| W = 1)
    matched_patients <- matched_patients |>
      dplyr::mutate(
        matched_p_0 = sum((W == 0) * weights * p_0) / sum((W == 0) * weights),
        matched_p_1 = sum((W == 1) * weights * p_1) / sum((W == 1) * weights),
        matched_tau_hat = matched_p_1 - matched_p_0,
        .by = "subclass"
      )
  } else if (tau_hat_method == "tau_avg") {
    matched_patients <- matched_patients |>
      dplyr::mutate(
        matched_tau_hat = (
          sum((W == 0) * weights * tau_hat) / sum((W == 0) * weights) + 
            sum((W == 1) * weights * tau_hat) / sum((W == 1) * weights)
          ) / 2,
        .by = "subclass"
      )
  } else if (tau_hat_method == "tau_treated") {
    matched_patients <- matched_patients |>
      dplyr::mutate(
        matched_tau_hat = sum((W == 1) * weights * tau_hat) / sum((W == 1) * weights),
        .by = "subclass"
      )
  } else if (tau_hat_method == "tau_control") {
    matched_patients <- matched_patients |>
      dplyr::mutate(
        matched_tau_hat = sum((W == 0) * weights * tau_hat) / sum((W == 0) * weights),
        .by = "subclass"
      )
  }
  # wait if memory usage exceeds limit
  if (!is.null(mem_limit)) {
    time <- 0
    while(mem_used_lgl(mem_limit$memory_limit_Gb, mem_limit$session_id)) {
      if (time > mem_limit$timeout) {
        stop("time limit reached without memory becoming available")
      }
      Sys.sleep(mem_limit$sleep_time)
      time <- time + sleep_time
    }
  }
  # C-for-benefit
  if (match == "combined_outcome_risk") {
    matched_obs_pred <- matched_patients |> 
      slice_head(n = 1, by = subclass) |> 
      select("matched_tau_hat", "matched_tau_obs", "sample_weights")
    cindex <- rcorr(
      matched_obs_pred |> pull("matched_tau_hat"),
      matched_obs_pred |> pull("matched_tau_obs"),
      matched_obs_pred |> pull("sample_weights")
    )
  } else {
    matched_obs_pred <- matched_patients |> 
      slice_head(n = 1, by = subclass) |> 
      select("matched_tau_hat", "matched_tau_obs")
    cindex <- rcorr(
      matched_obs_pred |> pull("matched_tau_hat"),
      matched_obs_pred |> pull("matched_tau_obs"),
      matched_patients |> 
        summarise(
          sample_weights = ifelse(estimand_method == "ATT", sum((W == 1) * sample_weights), sum((W == 0) * sample_weights)),
          .by = subclass
        ) |> 
        pull("sample_weights")
    )
  }
  c_for_benefit <- cindex["C Index"][[1]]
  c_for_benefit_se <- cindex["S.D."][[1]] / 2
  
  # Confidence interval
  if (CI == "simple") {
    lower_CI <- c_for_benefit + qnorm(0.5 - level / 2) * c_for_benefit_se
    upper_CI <- c_for_benefit + qnorm(0.5 + level / 2) * c_for_benefit_se
  } else if (CI == "bootstrap") {
    CB_for_CI <- c()
    B <- 0
    if (verbose) cli::cli_progress_bar("Bootstrapping", total = n_bootstraps)
    setTimeLimit(cpu = time_limit_CI, elapsed = time_limit_CI, transient = TRUE)
    while (B < n_bootstraps) {
      tryCatch(
        {
          # bootstrap matched patient pairs
          subclass_IDs <- unique(matched_patients$subclass)
          sample_subclass <- sample(
            subclass_IDs,
            length(subclass_IDs),
            replace = TRUE
          )
          # create bootstrap sample from pairs with sampled IDs
          bootstrap_matched_patients <- matched_patients |>
            dplyr::slice(2 * sample_subclass)
          # C-for-benefit in bootstrap sample
          bootstrap_cindex <- rcorr(
            bootstrap_matched_patients$matched_tau_hat,
            bootstrap_matched_patients$matched_tau_obs,
            bootstrap_matched_patients$sample_weights
          )
          CB_for_CI <- c(CB_for_CI, bootstrap_cindex["C Index"][[1]])
          B <- B + 1
          if (verbose) cli::cli_progress_update()
        },
        error = function(e) {
          if (
            grepl(
              "reached elapsed time limit|reached CPU time limit",
              e$message
            )
          ) {
            input <- readline(
              "Time limit reached. Do you want execution to continue (y/n) "
            )
            if (input %in% c("y", "n")) {
              if (input == "n") stop("Time limit reached, execution stopped.")
            } else {
              input <- readline(
                "Please input either 'y' to continue or 'n' to stop execution. "
              )
              if (input %in% c("y", "n")) {
                if (input == "n") stop("Time limit reached, execution stopped.")
              } else {
                stop("Answer not 'y' or 'n'. Execution stopped.")
              }
            }
          } else {
            stop(e)
          }
        }
      )
    }
    lower_CI <- as.numeric(quantile(CB_for_CI, 0.5 - level / 2))
    upper_CI <- as.numeric(quantile(CB_for_CI, 0.5 + level / 2))
  } else if (CI == "none") {
    lower_CI <- NA
    upper_CI <- NA
  }
  
  # scatterplot of predicted treatment effect within matched groups
  limits <- range(matched_patients$tau_hat)
  x_ann <- limits[1] + diff(limits) * 0.1
  y_ann <- limits[1] + diff(limits) * 0.9
  plot_data <- dplyr::tibble(
    x = matched_patients$tau_hat[seq(1, nrow(matched_patients), 2)], 
    y = matched_patients$tau_hat[seq(2, nrow(matched_patients), 2)]
  )
  plot_fct <- \(plot_data) {
    ggplot2::ggplot(plot_data) + 
      ggplot2::geom_point(ggplot2::aes(x = x, y = y)) + 
      ggplot2::geom_abline(intercept = 0, slope = 1) + 
      ggplot2::xlab("Predicted harm in control group") +
      ggplot2::ylab("Predicted harm in treated group") +
      ggplot2::scale_x_continuous(limits = limits) +
      ggplot2::scale_y_continuous(limits = limits) +
      ggplot2::annotate(
        "label", 
        x = x_ann, 
        y = y_ann, 
        label = glue::glue(
          "Correlation: {sprintf('%.2f', cor(plot_data$x, plot_data$y))}"
        )
      )
  }
  plot_fct <- as.list(plot_fct)
  plot_fct[[2]][[2]][[3]][[3]] <- x_ann
  plot_fct[[2]][[2]][[3]][[4]] <- y_ann
  plot_fct[[2]][[2]][[2]][[2]][[3]][[2]] <- limits
  plot_fct[[2]][[2]][[2]][[3]][[2]] <- limits
  plot_fct <- as.function(plot_fct)
  rlang::fn_env(plot_fct) <- .GlobalEnv
  
  return(
    list(
      matched_patients = matched_patients,
      c_for_benefit = c_for_benefit,
      c_for_benefit_se = c_for_benefit_se,
      lower_ci = lower_CI,
      upper_ci = upper_CI,
      cfb_tbl = dplyr::tibble(
        "ci_{0.5 - level / 2}" := lower_CI,
        "estimate" = c_for_benefit,
        "ci_{0.5 + level / 2}" := upper_CI
      ),
      matched_tau_hat_corr_plot_data = plot_data,
      matched_tau_hat_corr_plot = plot_fct,
      diagnostic_quantiles = list(
        quantile_counts = quantile_counts,
        quantile_benefit_probabilities = quantile_benefit_probabilities,
        quantile_c_for_benefit = quantile_cfb
      )
    )
  )
}
