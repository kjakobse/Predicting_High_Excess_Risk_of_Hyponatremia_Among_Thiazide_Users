CFRiskStrategy <- function(cf, 
                           cate_all,
                           W_all,
                           type = c("train", "test"),
                           level = 0.95,
                           excess_risk_cutoff = 0,
                           excess_risk_quantiles = NULL,
                           tol = 1e-5) {
  if (min(cf$W.hat) < tol) {
    cf$W.hat[which(cf$W.hat < tol)] <- tol
    warning("lowest W.hat was closer to 0 than the tolerance of ", tol, ".\n",
            "All W.hat lower than tol are set to tol. COnsider if the results are trustworthy!")
  }
  if (max(cf$W.hat) > 1-tol) {
    cf$W.hat[which(cf$W.hat > 1-tol)] <- 1-tol
    warning("highest W.hat was closer to 1 than the tolerance of ", tol, ".\n",
            "All W.hat lower than tol are set to tol. COnsider if the results are trustworthy!")
  }
  
  # input checks 
  type <- match.arg(type)
  
  
  # lower and upper confidence limits
  lower_lvl <- 0.5 - level / 2
  upper_lvl <- 0.5 + level / 2
  risk_strategy <- dplyr::tibble()
  
  # Expected population excess risk among the treated
  pop_excess_risk_simple <- tryCatch(
    dplyr::tibble(
      name = "population_excess_risk_in_treated_simple",
      "lower_{lower_lvl}" := NA_real_,
      estimate = weighted.mean(cf$predictions[cf$W.orig == 1], cf$sample.weights[cf$W.orig == 1]),
      "upper_{upper_lvl}" := NA_real_,
      std_err = NA_real_
    ),
    error = \(e) {
      dplyr::tibble(
        name = "population_excess_risk_in_treated_simple",
        "lower_{lower_lvl}" := NA_real_,
        estimate = NA_real_,
        "upper_{upper_lvl}" := NA_real_,
        std_err = NA_real_,
        error = e$message
      )
    }
  )
  pop_excess_risk_robust <- tryCatch(
    ifelse(
      type == "train",
      grf::average_treatment_effect(cf, target.sample = "treated"),
      with(
        cf,
        ATEAll(
          Y.orig = Y.orig, 
          Y.hat = Y.hat,
          W.orig = W.orig,
          W.hat = W.hat,
          tau.hat.pointwise = predictions,
          sample.weights = sample.weights,
          target.sample = "treated"
        )
      )
    ) |>
      (\(x) {
        dplyr::tibble(
          name = "population_excess_risk_in_treated_robust",
          "lower_{lower_lvl}" := x["estimate"] + qnorm(lower_lvl) * x["std.err"],
          estimate = x["estimate"],
          "upper_{upper_lvl}" := x["estimate"] + qnorm(upper_lvl) * x["std.err"],
          std_err = x["std.err"]
        )
      })(),
    error = \(e) {
      dplyr::tibble(
        name = "population_excess_risk_in_treated_robust",
        "lower_{lower_lvl}" := NA_real_,
        estimate = NA_real_,
        "upper_{upper_lvl}" := NA_real_,
        std_err = NA_real_,
        error = e$message
      )
    }
  ) 
  risk_strategy <- dplyr::bind_rows(risk_strategy, pop_excess_risk_simple)
  risk_strategy <- dplyr::bind_rows(risk_strategy, pop_excess_risk_robust)
  
  # Units to treat with non-thiazide instead of thiazide to avoid 1 case of hyponatremia
  # This indicates how many must receive the non-thiazide treatment instead of thiazide treatment
  # if we move fully to non-thiazide in the population to avoid one case of hypo-natremia.
  unt_pop <- tryCatch(
    ifelse(
      type == "train",
      grf::average_treatment_effect(cf, target.sample = "treated"),
      with(
        cf,
        ATEAll(
          Y.orig = Y.orig, 
          Y.hat = Y.hat,
          W.orig = W.orig,
          W.hat = W.hat,
          tau.hat.pointwise = predictions,
          sample.weights = sample.weights,
          target.sample = "treated"
        )
      )
    ) |>
      (\(x) {
        dplyr::tibble(
          name = "units_needed_to_treat",
          "lower_{lower_lvl}" := 1 / (x["estimate"] + qnorm(lower_lvl) * x["std.err"]),
          estimate = 1 / x["estimate"],
          "upper_{upper_lvl}" := 1 / (x["estimate"] + qnorm(upper_lvl) * x["std.err"]),
          std_err = NA_real_
        )
      })(),
    error = \(e) {
      dplyr::tibble(
        name = "units_needed_to_treat",
        "lower_{lower_lvl}" := NA_real_,
        estimate = NA_real_,
        "upper_{upper_lvl}" := NA_real_,
        std_err = NA_real_,
        error = e$message
      )
    }
  )
  risk_strategy <- dplyr::bind_rows(risk_strategy, unt_pop)
  
  if (!is.null(excess_risk_quantiles)) {
    excess_risk_quantiles <- sort(excess_risk_quantiles)
    excess_risk_cutoff <- dplyr::tibble(
      predictions = cf$predictions[, 1],
      W = cf$W.orig,
      sample_weight = cf$sample.weights
    ) |>
      dplyr::filter(W == 1) |>
      dplyr::arrange(predictions) |>
      dplyr::mutate(cum_weight = cumsum(sample_weight)) |>
      (\(x) dplyr::slice(x, purrr::map_dbl(excess_risk_quantiles, \(quan) sum(x$cum_weight <= quan * sum(x$sample_weight)))))() |>
      dplyr::pull(predictions) 
    if (any(excess_risk_quantiles == 0)) {
      excess_risk_cutoff <- c(min(cf$predictions[, 1]), excess_risk_cutoff)
    }
  } else {
    # weighted ecdf
    excess_risk_quantiles <- purrr::map_dbl(
      excess_risk_cutoff,
      \(cutoff) sum(cf$sample.weights[cf$predictions <= cutoff]) / sum(cf$sample.weights)
    )
  }
  
  for (i in seq_along(excess_risk_cutoff)) {
    quantile <- excess_risk_quantiles[i]
    cutoff <- excess_risk_cutoff[i]
    # cutoff for treatment strategy
    cutoff_tbl <- dplyr::tibble(
      name = glue::glue("excess_risk_cutoff_>{sprintf('%.0f', 100*quantile)}%"),
      "lower_{lower_lvl}" := NA_real_,
      estimate = cutoff,
      "upper_{upper_lvl}" := NA_real_,
      std_err = NA_real_
    )
    
    # expected excess risk among those treated above x% CATE
    # This provides the average excess risk among the subpopulation whom we would advice to take the 
    # alterative (non-thiazide) treatment, using a strategy of avoiding thiazides when excess risk is above
    # a pre-determined cut-off (here we consider 5%, 1%, and 0%)
    excess_risk_simple <- tryCatch(
      dplyr::tibble(
        name = glue::glue("excess_risk_simple_>{sprintf('%.0f', 100*quantile)}%"),
        "lower_{lower_lvl}" := NA_real_,
        estimate = ifelse(
          any(cf$predictions > cutoff & cf$W.orig == 1),
          weighted.mean(
            cf$predictions[cf$predictions >= cutoff & cf$W.orig == 1],
            cf$sample.weights[cf$predictions >= cutoff & cf$W.orig == 1]
          ),
          NA
        ),
        "upper_{upper_lvl}" := NA_real_,
        std_err = NA_real_
      ),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("excess_risk_tau_>{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_,
          error = e$message
        )
      }
    )
    
    excess_risk_robust <- tryCatch(
      ifelse(
        type == "train",
        grf::average_treatment_effect(cf, target.sample = "treated", subset = cf$predictions >= cutoff),
        with(
          cf,
          ATEAll(
            Y.orig = Y.orig[predictions >= cutoff], 
            Y.hat = Y.hat[predictions >= cutoff],
            W.orig = W.orig[predictions >= cutoff],
            W.hat = W.hat[predictions >= cutoff],
            tau.hat.pointwise = predictions[predictions >= cutoff],
            sample.weights = sample.weights[predictions >= cutoff],
            target.sample = "treated"
          )
        )
      ) |>
        (\(x) {
          dplyr::tibble(
            name = glue::glue("excess_risk_robust_>{sprintf('%.0f', 100*quantile)}%"),
            "lower_{lower_lvl}" := x["estimate"] + qnorm(lower_lvl) * x["std.err"],
            estimate = x["estimate"],
            "upper_{upper_lvl}" := x["estimate"] + qnorm(upper_lvl) * x["std.err"],
            std_err = x["std.err"]
          )
        })(),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("excess_risk_robust_>{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_,
          error = e$message
        )
      }
    )
    
    # Number of treated persons with an excess risk above x%
    num_patients <- sum(cf$predictions >= cutoff & cf$W.orig == 1) |>
      (\(x) {
        dplyr::tibble(
          name = glue::glue("num_patients_>{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = as.double(x),
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_
        )
      })()
    num_patients_ipcw <- sum(cf$sample.weights * (cf$predictions >= cutoff & cf$W.orig == 1)) |>
      (\(x) {
        dplyr::tibble(
          name = glue::glue("num_patients_ipcw_>{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = as.double(x),
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_
        )
      })()
    num_patients_all <- sum(cate_all >= cutoff & W_all == 1) |>
      (\(x) {
        dplyr::tibble(
          name = glue::glue("num_patients_all_>{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = as.double(x),
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_
        )
      })()
    
    # Difference in average excess risk between observed treatment strategy, and a treatment strategy 
    # where patients are treated only below a cut-off excess risk of x%. This provides the difference in
    # average excess risk between the two strategies. This imagines one exposure (non-thiazide treatment)
    # as having a "null" effect, while the main exposure (thiazide treatment) has some effect.
    # Note that the updated strategy will always treat fewer individuals. The average risk for the treated 
    # below and above the cut-off is also useful to get some context for the differences in average excess risk. 
    ate_treated_cutoff_simple <- tryCatch(
      dplyr::tibble(
        name = glue::glue("excess_risk_simple_<{sprintf('%.0f', 100*quantile)}%"),
        "lower_{lower_lvl}" := NA_real_,
        estimate = weighted.mean(
          (\(cf) {
            pred <- cf$predictions
            pred[pred >= cutoff] <- 0
            return(pred[cf$W.orig == 1])
          })(cf = cf),
          cf$sample.weights[cf$W.orig == 1]
        ),
        "upper_{upper_lvl}" := NA_real_,
        std_err = NA_real_
      ),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("excess_risk_simple_<{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_,
          error = e$message
        )
      }
    )
    
    ate_treated_cutoff_robust <- tryCatch(
      ifelse(
        type == "train",
        grf::average_treatment_effect(cf, target.sample = "treated", subset = cf$predictions < cutoff),
        with(
          cf,
          ATEAll(
            Y.orig = Y.orig[predictions < cutoff], 
            Y.hat = Y.hat[predictions < cutoff],
            W.orig = W.orig[predictions < cutoff],
            W.hat = W.hat[predictions < cutoff],
            tau.hat.pointwise = predictions[predictions < cutoff],
            sample.weights = sample.weights[predictions < cutoff],
            target.sample = "treated"
          )
        )
      ) |>
        (\(x) {
          dplyr::tibble(
            name = glue::glue("excess_risk_robust_<{sprintf('%.0f', 100*quantile)}%"),
            "lower_{lower_lvl}" := x["estimate"] + qnorm(0.025) * x["std.err"],
            estimate = x["estimate"],
            "upper_{upper_lvl}" := x["estimate"] + qnorm(0.975) * x["std.err"],
            std_err = x["std.err"]
          )
        })(),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("excess_risk_robust_<{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_,
          error = e$message
        )
      }
    ) 
    ate_treated <- ifelse(
        type == "train",
        grf::average_treatment_effect(cf, target.sample = "treated"),
        with(
          cf,
          ATEAll(
            Y.orig = Y.orig, 
            Y.hat = Y.hat,
            W.orig = W.orig,
            W.hat = W.hat,
            tau.hat.pointwise = predictions,
            sample.weights = sample.weights,
            target.sample = "treated"
          )
        )
      ) |>
      (\(x) {
        dplyr::tibble(
          name = "excess_risk_robust",
          "lower_{lower_lvl}" := x["estimate"] + qnorm(0.025) * x["std.err"],
          estimate = x["estimate"],
          "upper_{upper_lvl}" := x["estimate"] + qnorm(0.975) * x["std.err"],
          std_err = x["std.err"]
        )
      })()

    # The avoided risk is the difference between the excess risk in the treated population and the 
    # excess risk in the population where the top % of excess risk are changed to non-thiazide treatment.
    # Since this changes the excess risk to zero, the left-over from the difference are the top % of excess 
    # risks. 
    risk_avoid <- tryCatch(
      dplyr::tibble(
        name = glue::glue("risk_avoid_{sprintf('%.0f', 100*quantile)}%"),
        estimate = (1 - quantile) * excess_risk_robust$estimate,
        std_err = (1 - quantile) * excess_risk_robust$std_err,
        "lower_{lower_lvl}" := estimate + qnorm(lower_lvl) * std_err,
        "upper_{upper_lvl}" := estimate + qnorm(upper_lvl) * std_err
      ),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("risk_avoid_{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_,
          error = e$message
        )
      }
    ) |>
      dplyr::select("name", paste0("lower_", lower_lvl), "estimate", paste0("upper_", upper_lvl), "std_err")
    
    # Expected total decrease in hyponatremia cases by not treating patients with excess risk >x%
    # This is again based on a strategy of treating "as usual", except if excess risk is >x%, in which case
    # we treat with non-thiazide. This value gives an idea of the total impact of such a strategy with 
    # different cut-offs, as it gives us the number of people who we expect would have experienced hypo-
    # natremia with the current strategy but avoided it had they been treated using the new strategy. 
    # Effectively (since only these receive a different treatment) this measures how many of those >x% risk 
    # we expect would have avoided experiencing hypo-natremia had they received a non-thiazide treatment. 
    cases_avoid <- tryCatch(
      ifelse(
        type == "train",
        grf::average_treatment_effect(cf, target.sample = "treated", subset = cf$predictions >= cutoff),
        with(
          cf,
          ATEAll(
            Y.orig = Y.orig[predictions >= cutoff], 
            Y.hat = Y.hat[predictions >= cutoff],
            W.orig = W.orig[predictions >= cutoff],
            W.hat = W.hat[predictions >= cutoff],
            tau.hat.pointwise = predictions[predictions >= cutoff],
            sample.weights = sample.weights[predictions >= cutoff],
            target.sample = "treated"
          )
        )
      ) |>
        (\(x) {
          dplyr::tibble(
            name = glue::glue("cases_avoid_{sprintf('%.0f', 100*quantile)}%"),
            # place number of patients with risk above cutoff in lower_{lower_level} for access when aggregating
            "lower_{lower_lvl}" := sum(cf$predictions >= cutoff & cf$W.orig == 1),
            estimate = x["estimate"],
            "upper_{upper_lvl}" := sum(cf$predictions >= cutoff & cf$W.orig == 1),
            std_err = x["std.err"]
          )
        })(),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("cases_avoid_{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_
        )
      }
    )
    
    # Units to treat with non-thiazide instead of thiazide to avoid 1 case of hyponatremia
    # This indicates how many must receive the non-thiazide treatment instead of thiazide treatment from the
    # sub-population with a excess risk >x% to on average avoid one case of hypo-natremia.
    units_nt <- tryCatch(
      ifelse(
        type == "train",
        grf::average_treatment_effect(cf, target.sample = "treated", subset = cf$predictions >= cutoff),
        with(
          cf,
          ATEAll(
            Y.orig = Y.orig[predictions >= cutoff], 
            Y.hat = Y.hat[predictions >= cutoff],
            W.orig = W.orig[predictions >= cutoff],
            W.hat = W.hat[predictions >= cutoff],
            tau.hat.pointwise = predictions[predictions >= cutoff],
            sample.weights = sample.weights[predictions >= cutoff],
            target.sample = "treated"
          )
        )
      ) |>
        (\(x) {
          dplyr::tibble(
            name = glue::glue("units_nt_{sprintf('%.0f', 100*quantile)}%"),
            "lower_{lower_lvl}" := 1 / (x["estimate"] + qnorm(lower_lvl) * x["std.err"]),
            estimate = 1 / x["estimate"],
            "upper_{upper_lvl}" := 1 / (x["estimate"] + qnorm(upper_lvl) * x["std.err"])
          )
        })(),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("units_nt_{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_
        )
      }
    )
    
    # population level benefit of treatment strategy where treatment if below a cutoff
    # Comparison of cf treatment strategy vs. inverse strategy
    # This looks at a treatment strategy where we treat only based on the predicted excess risk, i.e. 
    # everyone below a cutoff receives thiazides and everyone above the cutoff receives non-thiazide. This is 
    # a more intrusive strategy, since we completely change how treatment is chosen, instead of recommending
    # against a particular treatment option for a select subpopulation. 
    # Here we compare the treatment strategy of thiazide to all below a cut-off and non-thiazide to those 
    # above the cut-off with the most extreme inverse strategy, i.e. flipping the treatment assignment on 
    # its head. This gives an idea of the average gain in the population from  the chosen strategy compared
    # to the alternative strategy for each individual. 
    w_treat <- sum(cf$sample.weights[cf$predictions < cutoff]) / sum(cf$sample.weights)
    w_notreat <- sum(cf$sample.weights[cf$predictions >= cutoff]) / sum(cf$sample.weights)
    ate_treat <- tryCatch(
      ifelse(
        type == "train",
        grf::average_treatment_effect(cf, subset = cf$predictions < cutoff),
        with(
          cf,
          ATEAll(
            Y.orig = Y.orig[predictions < cutoff], 
            Y.hat = Y.hat[predictions < cutoff],
            W.orig = W.orig[predictions < cutoff],
            W.hat = W.hat[predictions < cutoff],
            tau.hat.pointwise = predictions[predictions < cutoff],
            sample.weights = sample.weights[predictions < cutoff],
            target.sample = "all"
          )
        )
      ), 
      error = \(e) NA_real_
      )
    ate_notreat <- tryCatch(
      ifelse(
        type == "train",
        grf::average_treatment_effect(cf, subset = cf$predictions >= cutoff),
        with(
          cf,
          ATEAll(
            Y.orig = Y.orig[predictions >= cutoff], 
            Y.hat = Y.hat[predictions >= cutoff],
            W.orig = W.orig[predictions >= cutoff],
            W.hat = W.hat[predictions >= cutoff],
            tau.hat.pointwise = predictions[predictions >= cutoff],
            sample.weights = sample.weights[predictions >= cutoff],
            target.sample = "all"
          )
        )
      ), 
      error = \(e) NA_real_
      )
    PB_hat <- tryCatch(
      tibble(
        name = glue::glue("PB_hat_{sprintf('%.0f', 100*quantile)}%"),
        "lower_{lower_lvl}" := (w_treat * ate_treat["estimate"] - w_notreat * ate_notreat["estimate"]) + 
          qnorm(lower_lvl) * sqrt((w_treat * ate_treat["std.err"])^2 + (w_notreat * ate_notreat["std.err"])^2),
        estimate = w_treat * ate_treat["estimate"] - w_notreat * ate_notreat["estimate"],
        "upper_{upper_lvl}" := (ate_treat["estimate"] - ate_notreat["estimate"]) + 
          qnorm(upper_lvl) * sqrt((w_treat * ate_treat["std.err"])^2 + (w_notreat * ate_notreat["std.err"])^2),
        std_err = sqrt((w_treat * ate_treat["std.err"])^2 + (w_notreat * ate_notreat["std.err"])^2)
      ),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("PB_hat_{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_
        )
      }
    )
    
    # cf treatment strategy vs. non treated (contains zeroes in ATE computation)
    # Here, instead of comparing the treatment strategy against the inverse strategy, we compare it against 
    # the null-strategy of not giving thiazides to anyone. This corresponds to the excess risk among those
    # receiving thiazides averages over the full population.
    PB0 <- tryCatch(
      dplyr::tibble(
        name = glue::glue("PB0_{sprintf('%.0f', 100*quantile)}%"),
        "lower_{lower_lvl}" := (w_treat * ate_treat["estimate"]) + 
          qnorm(lower_lvl) * w_treat * ate_treat["std.err"],
        estimate = w_treat * ate_treat["estimate"],
        "upper_{upper_lvl}" := (ate_treat["estimate"] - ate_notreat["estimate"]) + 
          qnorm(upper_lvl) * w_treat * ate_treat["std.err"],
        std_err = w_treat * ate_treat["std.err"]
      ),
      error = \(e) {
        dplyr::tibble(
          name = glue::glue("PB0_{sprintf('%.0f', 100*quantile)}%"),
          "lower_{lower_lvl}" := NA_real_,
          estimate = NA_real_,
          "upper_{upper_lvl}" := NA_real_,
          std_err = NA_real_
        )
      }
    )
    
    tmp <- bind_rows(
      cutoff_tbl,
      excess_risk_tau,
      excess_risk_robust,
      num_patients,
      num_patients_ipcw,
      num_patients_all,
      ate_treated_cutoff_simple,
      ate_treated_cutoff,
      risk_avoid,
      cases_avoid,
      units_nt,
      PB_hat,
      PB0
    )
    
    risk_strategy <- bind_rows(risk_strategy, tmp)
  }
  return(risk_strategy)
}