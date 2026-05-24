library(parallel)
library(doParallel)
library(foreach)

root_path <- "."
data_path <- file.path(root_path, "realdata")

source(file.path(root_path, "admm.R"))

target_datasets <- c(
  "crabs_train",
  "iris_train",
  "penguins_train",
  "seeds_train",
  "wdbc_train",
  "wine_train"
)

rds_files <- list.files(
  path = data_path,
  pattern = "\\.[Rr][Dd][Ss]$",
  full.names = FALSE
)

data_files <- rds_files[
  grepl("_train\\.[Rr][Dd][Ss]$", rds_files) &
    !grepl("_true\\.[Rr][Dd][Ss]$", rds_files)
]

dataset_names <- sub("\\.[Rr][Dd][Ss]$", "", data_files)
dataset_names <- intersect(dataset_names, target_datasets)

data_list <- list()

for (dataset in dataset_names) {
  X <- readRDS(file.path(data_path, paste0(dataset, ".Rds")))

  if (is.list(X)) {
    X <- do.call(cbind, X)
  } else {
    X <- as.matrix(X)
  }

  data_list[[dataset]] <- X
}

lambda_grid <- round(exp(seq(log(10), log(10000), length.out = 60)), 0)
gamma_grid <- 10

task_grid <- expand.grid(
  dataset = names(data_list),
  lambda = lambda_grid,
  gamma = gamma_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

num_str <- function(x) {
  format(x, scientific = FALSE, trim = TRUE)
}

rho <- 1000
a <- 1e-7

result_dir <- file.path(data_path, "results_ours")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- detectCores() - 1
n_cores <- min(n_cores, nrow(task_grid))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(
  cl,
  varlist = c("root_path"),
  envir = environment()
)

clusterEvalQ(
  cl,
  source(file.path(root_path, "admm.R"))
)

invisible(
  foreach(
    i = seq_len(nrow(task_grid)),
    .packages = character(0),
    .export = c("data_list", "task_grid", "rho", "a", "result_dir", "num_str")
  ) %dopar% {
    dataset <- task_grid$dataset[i]
    lambda <- task_grid$lambda[i]
    gamma <- task_grid$gamma[i]

    X <- data_list[[dataset]]

    n <- nrow(X)
    mu <- n^2 * rho * ((1 - a) / a)

    save_path <- file.path(result_dir, dataset)
    dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

    res <- run_admm(
      X = X,
      rho = rho,
      mu = mu,
      lambda = lambda,
      gamma = gamma,
      max_iter = 10000000,
      tol = 1e-6,
      save_path = save_path,
      scen_nm = dataset,
      warm_start = NULL,
      save_every = 20
    )

    NULL
  }
)

stopCluster(cl)