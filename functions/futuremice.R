futuremice <- function (data, m = 5, parallelseed = NA, n.core = NULL, seed = NA, 
          use.logical = TRUE, future.plan = "multisession", packages = NULL, 
          globals = NULL, available = parallelly::availableCores(logical = use.logical), ...) 
{
  data <- mice:::check.dataform(data)
  m <- mice:::check.m(m)
  if (sum(is.na(data)) == 0) {
    stop("Data has no missing values")
  }
  n.core <- mice:::check.cores(n.core, available, m)
  if (n.core > 1) {
    dist.core <- cut(1:m, n.core, labels = paste0("core", 1:n.core))
  }
  else {
    dist.core <- rep("core1", m)
  }
  n.imp.core <- as.vector(table(dist.core))
  if (!is.na(seed)) {
    if (n.core > 1) {
      if (interactive()) {
        msg <- "Be careful; specifying seed rather than parallelseed results in duplicate imputations.\nDo you want to continue?\n"
        ask <- askYesNo(msg, prompts = getOption("askYesNo", 
                                                 gettext(c("Yes", "No, ignore seed", "Cancel"))))
        if (isTRUE(ask)) {
          seed <- seed
          warning("Be careful; the imputations will be the same over the cores.")
        }
        else if (isFALSE(ask)) {
          seed <- NA
          message("Parallelseed is specified for you, and is accessible in the output object under $parallelseed.")
        }
        else if (is.na(ask)) {
          stop("You stopped futuremice. To obtain unique, but reproducible imputations, specify parallelseed.")
        }
      }
      else {
        warning("Be careful; the imputations will be identical over the cores. Perhaps you want to specify parallelseed, for unique, but reproducible results.")
      }
    }
  }
  if (!is.na(parallelseed)) {
    set.seed(parallelseed)
  }
  else {
    if (!exists(".Random.seed")) {
      set.seed(NULL)
    }
    parallelseed <- get(".Random.seed", envir = globalenv(), 
                        mode = "integer", inherits = FALSE)
  }
  future::plan(future.plan, workers = n.core)
  imps <- furrr::future_map(
    n.imp.core, 
    function(x) {
    mice(data = data, m = x, printFlag = FALSE, seed = seed, ...)
  }, 
  .options = furrr::furrr_options(seed = TRUE, globals = globals, packages = packages)
  )
  future::plan(future::sequential)
  imp <- imps[[1]]
  if (length(imps) > 1) {
    for (i in 2:length(imps)) {
      imp <- ibind(imp, imps[[i]])
    }
  }
  for (i in 1:length(imp$imp)) {
    colnames(imp$imp[[i]]) <- 1:imp$m
  }
  imp$parallelseed <- parallelseed
  return(imp)
}
