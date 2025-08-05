TIH_DynamicSubgroups <- function(temp_path = "E:/workdata/708445/KIJA/causal_forest/results/data/temp/",
                                 remove_temp = c("all", "dynamic_subgroups", "none"),
                                 cf_append = "",
                                 num.trees = 2000,
                                 model_name,
                                 setup_name,
                                 cohort,
                                 cluster_name,
                                 n_folds,
                                 n_rankings,
                                 continuous,
                                 discrete,
                                 tune_parameters,
                                 data_index,
                                 num_threads,
                                 num_clusters,
                                 num_fold_par,
                                 horizon,
                                 thiazide_limit) {
  # Sample fold indices to use with all imputed datasets -------------------------------------------------------------
  remove_temp <- match.arg(remove_temp)
  # set option used to retrive available cores for futuremice
  options("parallelly.availableCores.methods" = "system")
  
  if (length(list.files(temp_path, pattern = paste0("^indices_dyn_", setup_name, ".rds$"))) == 0 && data_index == 1) {
    indices_dyn <- (\(cohort) {
      clusters <- pull(cohort, cluster_name)
      folds <- rep(0, length(clusters))
      for (i in unique(clusters)) {
        folds[clusters == i] <- sample(
          sort(seq_along(folds[clusters == i]) %% n_folds) + 1
        )
      }
      return(split(seq_along(folds), folds))
    })(
      cohort = cohort
    )
    saveRDS(indices_dyn, paste0(temp_path ,"indices_dyn_", setup_name, ".rds"))
  } else {
    while (length(list.files(temp_path, pattern = paste0("^indices_dyn_", setup_name, ".rds$"))) == 0) {
      Sys.sleep(5)
    }
    indices_dyn <- readRDS(paste0(temp_path, "indices_dyn_", setup_name, ".rds"))
  }
  
  # Sample weights for dynamic subgroups -------------------------------------------------------------------------------
  # compute sample weights in each fold using the remaining folds to train survival forest
  if (num_fold_par > 1) {
    plan(multisession, workers = min(num_fold_par, length(indices_dyn)))
  }
  sample_weights_dyn <- furrr::future_map2(
    .x = indices_dyn,
    .y = seq_along(indices_dyn),
    .options = furrr_options(
      packages = c("dplyr", "grf", "future", "furrr"),
      globals = c(
        "num_threads", "num_clusters", "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "temp_path", 
        "horizon", "thiazide_limit", "num.trees", "cluster_name"
        ),
      seed = TRUE
    ),
    .f = function(idx, fold, data_index, setup_name, cohort, continuous, discrete) {
      if (length(list.files(temp_path, pattern = paste0("^sample_weights_", fold, "_", data_index, "_", setup_name, ".rds$"))) > 0) {
        return(NULL)
      }
      continuous <- continuous$sample_weight
      discrete <- discrete$sample_weight
      forest_C <- GRFAnalysisWrapper(
        (cohort |> mutate(D_C = 1 - !!sym(paste0("D_", thiazide_limit))))[-idx,],
        type = "survival",
        continuous = c(!!continuous, "W"),
        discrete = !!discrete,
        Y = paste0("Y_", thiazide_limit), D = "D_C",
        clusters = cohort[[cluster_name]][-idx],
        num.trees = num.trees, alpha = 0.0002,
        seed = sample(-999999999:999999999, 1),
        num.threads = max(floor(num_threads / num_clusters / num_fold_par), 1)
      )
      newdata <- cohort |> 
        slice(idx) |>
        filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
        select(!!continuous, "W") |>
        bind_cols(
          cohort |>
            slice(idx) |>
            filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
            select(!!discrete) |>
            DiscreteCovariatesToOneHot()
        )
      Y_obs <- cohort |> 
        slice(idx) |>
        filter(!is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon)))) |>
        pull(paste0("Y_", thiazide_limit))
      predict_C <- predict(
        forest_C, 
        newdata = newdata,
        failure.times = pmin(horizon, Y_obs), 
        prediction.times = "time"
      )
      censoring_prob <- predict_C$predictions[, 1]
      sample_weights <- 1 / censoring_prob
      # naming: sample_weights_(fold number)_(imputation number)_(setup)
      saveRDS(
        sample_weights, 
        paste0(temp_path, "sample_weights_", fold, "_", data_index, "_", setup_name, ".rds")
      )
      return(NULL)
    },
    cohort = cohort,
    continuous = continuous,
    discrete = discrete,
    setup_name = setup_name,
    data_index = data_index
  )
  plan(sequential)
  
  for (
    j in seq_along(list.files(temp_path, pattern = paste0("^sample_weights_\\d{1,}_", data_index, "_", setup_name, ".rds$")))
  ) {
    sample_weights_dyn[[j]] <- readRDS(paste0(temp_path, "sample_weights_", j, "_", data_index, "_", setup_name, ".rds"))
  }
  
  # Y_hat and W_hat for dynamic subgroups ------------------------------------------------------------------------------ 
  # compute Y_hat and W_hat in each fold using the remaining folds to train regression forests
  sample_weights_dyn_inv <- list()
  for(i in seq_along(sample_weights_dyn)) sample_weights_dyn_inv[[i]] <- unlist(sample_weights_dyn[-i])
  
  if (num_fold_par > 1) {
    plan(multisession, workers = min(num_fold_par, length(indices_dyn)))
  }
  Y_hat_dyn <- furrr::future_pmap(
    .l = list(
      idx = indices_dyn,
      fold = seq_along(indices_dyn),
      sample_weights = sample_weights_dyn_inv # need all sample weights not in i'th fold
    ), 
    .options = furrr_options(
      packages = c("dplyr", "grf", "future", "furrr"),
      globals = c(
        "num_threads", "num_clusters", "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "temp_path",
        "num.trees", "thiazide_limit", "horizon", "cluster_name"
        ),
      seed = TRUE
    ),
    .f = function(idx, fold, sample_weights, setup_name, data_index, cohort, continuous, discrete, tune_parameters) {
      if (length(list.files(temp_path, pattern = paste0("^Y_hat_", fold, "_", data_index, "_", setup_name, ".rds$"))) > 0) {
        return(NULL)
      }
      data_obs <- filter(cohort[-idx,], !is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon))))
      continuous <- continuous$exp_out
      discrete <- discrete$exp_out
      tune_parameters <- tune_parameters$exp_out
      forest_Y <- GRFAnalysisWrapper(
        data_obs,
        type = "regression",
        continuous = !!continuous,
        discrete = !!discrete,
        Y = paste0("Y_", thiazide_limit, "_", horizon),
        sample.weights = sample_weights,
        clusters = data_obs[[cluster_name]],
        num.trees = num.trees, alpha = 0.0002,
        tune.parameters = tune_parameters,
        tune.num.trees = 200, tune.num.reps = 100, tune.num.draws = 1000,
        seed = sample(-999999999:999999999, 1),
        num.threads = max(floor(num_threads / num_clusters / num_fold_par), 1)
      )
      newdata <- cohort |> 
        slice(idx) |>
        filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
        select(!!continuous) |>
        bind_cols(
          cohort |>
            slice(idx) |>
            filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
            select(!!discrete) |>
            DiscreteCovariatesToOneHot()
        )
      Y_hat <- predict(
        forest_Y, 
        newdata = newdata
      )$predictions
      saveRDS(Y_hat, paste0(temp_path, "Y_hat_", fold, "_", data_index, "_", setup_name, ".rds"))
      return(NULL)
    },
    cohort = cohort,
    continuous = continuous,
    discrete = discrete,
    tune_parameters = tune_parameters,
    setup_name = setup_name,
    data_index = data_index
  )
  W_hat_dyn <- furrr::future_pmap(
    .l = list(
      idx = indices_dyn,
      fold = seq_along(indices_dyn),
      sample_weights = sample_weights_dyn_inv # need all sample weights not in i'th fold
    ),
    .options = furrr_options(
      packages = c("dplyr", "grf", "future", "furrr"),
      globals = c(
        "num_threads", "num_clusters", "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "temp_path",
        "num.trees", "thiazide_limit", "horizon", "cluster_name"
        ),
      seed = TRUE
    ),
    .f = function(idx, fold, sample_weights, setup_name, data_index, cohort, continuous, discrete, tune_parameters) {
      if (length(list.files(temp_path, pattern = paste0("^W_hat_", fold, "_", data_index, "_", setup_name, ".rds$"))) > 0) {
        return(NULL)
      }
      data_obs <- filter(cohort[-idx,], !is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon))))
      continuous <- continuous$exp_out
      discrete <- discrete$exp_out
      tune_parameters <- tune_parameters$exp_out
      forest_W <- GRFAnalysisWrapper(
        data_obs,
        type = "regression",
        continuous = !!continuous,
        discrete = !!discrete,
        Y = "W",
        sample.weights = sample_weights,
        clusters = data_obs[[cluster_name]],
        num.trees = num.trees, alpha = 0.0002,
        tune.parameters = tune_parameters,
        tune.num.trees = 200, tune.num.reps = 100, tune.num.draws = 1000,
        seed = sample(-999999999:999999999, 1),
        num.threads = max(floor(num_threads / num_clusters / num_fold_par), 1)
      )
      newdata <- cohort |> 
        slice(idx) |>
        filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
        select(!!continuous) |>
        bind_cols(
          cohort |>
            slice(idx) |>
            filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
            select(!!discrete) |>
            DiscreteCovariatesToOneHot()
        )
      W_hat <- predict(
        forest_W, 
        newdata = newdata
      )$predictions
      saveRDS(W_hat, paste0(temp_path, "W_hat_", fold, "_", data_index, "_", setup_name, ".rds"))
      return(NULL)
    },
    cohort = cohort,
    continuous = continuous,
    discrete = discrete,
    tune_parameters = tune_parameters,
    setup_name = setup_name,
    data_index = data_index
  )
  plan(sequential)
  for (
    j in seq_along(list.files(temp_path, pattern = paste0("^Y_hat_\\d{1,}_", data_index, "_", setup_name, ".rds$")))
  ) {
    Y_hat_dyn[[j]] <- readRDS(paste0(temp_path, "Y_hat_", j, "_", data_index, "_", setup_name, ".rds"))
  }
  for (
    j in seq_along(list.files(temp_path, pattern = paste0("^W_hat_\\d{1,}_", data_index, "_", setup_name, ".rds$")))
  ) {
    W_hat_dyn[[j]] <- readRDS(paste0(temp_path, "W_hat_", j, "_", data_index, "_", setup_name, ".rds"))
  }
  
  # causal forests for dynamic subgroups -----------------------------------------------------------------------------
  # Compute tau_hat and rankings within each fold using causal forest trained on the remaining folds
  Y_hat_dyn_inv <- vector("list", length(Y_hat_dyn))
  W_hat_dyn_inv <- vector("list", length(W_hat_dyn))
  for(i in seq_along(Y_hat_dyn)) Y_hat_dyn_inv[[i]] <- unlist(Y_hat_dyn[-i])
  for(i in seq_along(W_hat_dyn)) W_hat_dyn_inv[[i]] <- unlist(W_hat_dyn[-i])
  
  if (num_fold_par > 1) {
    plan(multisession, workers = min(num_fold_par, length(indices_dyn)))
  }
  cfw_4mo_dyn <- furrr::future_pmap(
    .l = list(
      idx = indices_dyn,
      fold = seq_along(indices_dyn),
      sample_weights_train = sample_weights_dyn_inv,
      sample_weights_test = sample_weights_dyn,
      Y_hat_train = Y_hat_dyn_inv,
      W_hat_train = W_hat_dyn_inv,
      Y_hat_test = Y_hat_dyn,
      W_hat_test = W_hat_dyn
    ),
    .options = furrr_options(
      packages = c("stringr", "purrr", "dplyr", "grf", "future", "furrr", "DescTools"),
      globals = c("num_threads", "num_clusters", "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "temp_path", 
                  "cf_append", "num.trees", "thiazide_limit", "horizon", "cluster_name"),
      seed = TRUE
    ),
    .f = function(idx, fold, sample_weights_train, sample_weights_test, 
             Y_hat_train, W_hat_train, Y_hat_test, W_hat_test, 
             cohort, continuous, discrete, tune_parameters, 
             n_rankings, setup_name, data_index) {
      if (
        length(list.files(temp_path, pattern = paste0("^dynamic_subgroups_", fold, "_", data_index, "_", setup_name, "_", cf_append, ".rds$"))) > 0
      ) {
        return(NULL)
      }
      data_obs_train <- filter(cohort[-idx,], !is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]]))
      data_obs_test <- cohort |> 
        slice(idx) |>
        filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]]))
      continuous <- continuous$cf
      discrete <- discrete$cf
      tune_parameters <- tune_parameters$cf
      cf <- GRFAnalysisWrapper(
        data_obs_train,
        type = "causal",
        continuous = !!continuous,
        discrete = !!discrete,
        Y = paste0("Y_", thiazide_limit, "_", horizon), W = "W",
        Y.hat = Y_hat_train, W.hat = W_hat_train,
        sample.weights = sample_weights_train,
        clusters = data_obs_train[[cluster_name]],
        num.trees = num.trees, alpha = 0.0002,
        tune.parameters = tune_parameters,
        tune.num.trees = 200, tune.num.reps = 100, tune.num.draws = 1000,
        num.threads = max(floor(num_threads / num_clusters / num_fold_par), 1)
      )
      if("X_03_SjÃ¦lland" %in% names(cf$X.orig)) {
        names(cf$X.orig)[
          which(names(cf$X.orig) == "X_03_SjÃ¦lland")
        ] <- "X_03_Sjælland"
      }
      newdata <- cohort |> 
        slice(idx) |>
        filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
        select(!!continuous) |>
        bind_cols(
          cohort |>
            slice(idx) |>
            filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
            select(!!discrete) |>
            DiscreteCovariatesToOneHot()
        )
      
      pred_cf <- predict(
        object = cf,
        newdata = newdata,
        estimate.variance = TRUE,
        num.threads = max(floor(num_threads / num_clusters / num_fold_par), 1)
      )
      
      tau_hat <- pred_cf$predictions
      tau_hat_var <- pred_cf$variance.estimates
      
      mu_hat_0 <- Y_hat_test - W_hat_test * tau_hat
      mu_hat_1 <- Y_hat_test + (1 - W_hat_test) * tau_hat
      
      aipw_scores <- 
        tau_hat +
        data_obs_test[["W"]] / W_hat_test * (data_obs_test[[paste0("Y_", thiazide_limit, "_", horizon)]] - mu_hat_1) -
        (1 - data_obs_test[["W"]]) / (1 - W_hat_test) * (data_obs_test[[paste0("Y_", thiazide_limit, "_", horizon)]] - mu_hat_0)
      
      # collect output
      out <- dplyr::tibble(
        id = idx[which(idx %in% which(!is.na(cohort[[paste0("Y_", thiazide_limit, "_", horizon)]])))],
        sample_weights = sample_weights_test,
        W_hat = W_hat_test,
        Y_hat = Y_hat_test,
        mu_hat_0 = mu_hat_0,
        mu_hat_1 = mu_hat_1,
        tau_hat = tau_hat,
        tau_hat_var = tau_hat_var,
        aipw_scores = aipw_scores,
        fold = fold
      )
      
      # add rankings to output
      for (i in seq_along(n_rankings)) {
        # rank observations by cate in held-out fold
        tau_hat_quantiles <- DescTools::Quantile(
          x = tau_hat,
          weights = sample_weights_test,
          probs = seq(0, 1, by = 1 / n_rankings[i])
        )
        # if quantiles are not unique, manually sort and cut into appropriate 
        # groups
        if (length(tau_hat_quantiles) == length(unique(tau_hat_quantiles))) {
          ranking <- cut(
            x = tau_hat,
            breaks = tau_hat_quantiles,
            include.lowest = TRUE,
            labels = seq_len(n_rankings[i])
          )
        } else {
          len <- length(tau_hat)
          ranking <- tau_hat |>
            (\(x) {
              dplyr::tibble(
                id = seq_along(x),
                tau_hat = x
              )
            }
            )() |>
            dplyr::arrange(.data$tau_hat) |>
            dplyr::mutate(
              id_2 = seq_along(.data$tau_hat),
              rank = cut(
                x = .data$id_2,
                breaks = c(
                  seq(0, len %% n_rankings[i], by = 1) *
                    (len %/% n_rankings[i] + 1),
                  seq(
                    len %% n_rankings[i] + 1,
                    n_rankings[i],
                    length.out = n_rankings[i] - len %% n_rankings[i]
                  ) *
                    len %/% n_rankings[i] + len %% n_rankings[i]
                ),
                include.lowest = TRUE,
                labels = seq_len(n_rankings[i])
              )
            ) |>
            dplyr::arrange(.data$id) |>
            dplyr::pull(.data$rank)
        }
        out <- out |>
          mutate("ranking_{n_rankings[i]}" := .env$ranking)
      }
      
      saveRDS(out, paste0(temp_path, "dynamic_subgroups_", fold, "_", data_index, "_", setup_name, "_", cf_append, ".rds"))
      return(NULL)
    },
    cohort = cohort,
    continuous = continuous,
    discrete = discrete,
    tune_parameters = tune_parameters,
    n_rankings = n_rankings, 
    setup_name = setup_name,
    data_index = data_index
  )
  plan(sequential)
  
  # read in CATE estimates and ranks from dynamic subgroups in each fold
  for (
    j in seq_along(list.files(temp_path, pattern = paste0("^dynamic_subgroups_\\d{1,}_", data_index, "_", setup_name, "_", cf_append, ".rds$")))
  ) {
    cfw_4mo_dyn[[j]] <- readRDS(paste0(temp_path, "dynamic_subgroups_", j, "_", data_index, "_", setup_name, "_", cf_append,  ".rds"))
  }
  
  # bind together results from folds
  cfw_4mo_dyn <- cfw_4mo_dyn |>
    list_rbind() |>
    arrange(id)
  
  # save CATE estimates and ranks from dynamic subgroups
  saveRDS(cfw_4mo_dyn, paste0(path, model_name, "_dyn_", data_index, ".rds"))
  
  # create data for heatmap of covariate values in each rank
  X_heatmap <- cohort |> 
    filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
    select(!!continuous$cf) |>
    bind_cols(
      cohort |>
        filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
        select(!!discrete$cf) |>
        DiscreteCovariatesToOneHot()
    )
  cfw_4mo_dyn <- list(rank_tbl = cfw_4mo_dyn)
  for (i in seq_along(n_rankings)) {
    heatmap_data <- purrr::map(
      names(X_heatmap),
      \(covariate, cohort, X_heatmap) {
        # calculate the average and standard error of each covariate
        fmla <- formula(glue::glue("`{covariate}` ~ 0 + ranking_{n_rankings[i]}"))
        data_forest <- X_heatmap |>
          dplyr::as_tibble() |>
          dplyr::mutate(
            Y = cohort[[paste0("Y_", thiazide_limit, "_", horizon)]],
            W = cohort[["W"]],
            "ranking_{n_rankings[i]}" := factor(cfw_4mo_dyn[["rank_tbl"]][[glue::glue("ranking_{n_rankings[i]}")]])
          )
        ols <- lm(formula = fmla, data = data_forest)
        
        # results
        avg <- coef(ols)
        stderr <- sqrt(diag(vcovHC(ols)))
        
        # collect results in table
        dplyr::tibble(
          covariate = covariate,
          avg = avg,
          stderr = stderr,
          ranking = paste0("Q", seq_len(n_rankings[i])),
          scaling = pnorm((avg - mean(avg)) / sd(avg)),
          # standard error between groups normalized by
          # standard error between all observations:
          variation = sd(avg) / sd(data_forest[[{{ covariate }}]]),
          labels = paste0(
            sprintf("%.3f", avg),
            "\n",
            "(",
            sprintf("%.3f", stderr),
            ")"
          )
        )
      },
      cohort = cohort |> filter(!is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |> as_tibble(),
      X_heatmap = X_heatmap
    ) |>
      purrr::list_rbind()
    
    cfw_4mo_dyn <- c(
      cfw_4mo_dyn,
      rlang::list2("heatmap_data_{n_rankings[i]}" := heatmap_data)
    )
  }
  
  # save CATE estimates, ranks and heatmap data from dynamic subgroups
  saveRDS(cfw_4mo_dyn, paste0(path, model_name, "_dyn_", data_index, ".rds"))
  
  # remove all temporary files
  for (
    j in seq_along(list.files(temp_path, pattern = paste0("^dynamic_subgroups_\\d{1,}_", data_index, "_", setup_name, "_", cf_append, ".rds$")))
  ) {
    if (remove_temp %in% c("all", "dynamic_subgroups")) {
      file.remove(paste0(temp_path, "dynamic_subgroups_", j, "_", data_index, "_", setup_name, "_", cf_append, ".rds"))
    }
    if (remove_temp == "all") {
      file.remove(paste0(temp_path, "sample_weights_", j, "_", data_index, "_", setup_name, ".rds"))
      file.remove(paste0(temp_path, "Y_hat_", j, "_", data_index, "_", setup_name, ".rds"))
      file.remove(paste0(temp_path, "W_hat_", j, "_", data_index, "_", setup_name, ".rds"))
    }
  }
  
  return(NULL)
}