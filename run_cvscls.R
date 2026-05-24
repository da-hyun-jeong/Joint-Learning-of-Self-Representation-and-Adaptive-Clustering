library(cvxclustr)
library(parallel)
library(doParallel)
library(foreach)

target_datasets <- c(
  "crabs_train",
  "iris_train",
  "penguins_train",
  "seeds_train",
  "wdbc_train",
  "wine_train"
)

root_dir <- "realdata"

rds_files <- list.files(
  path = root_dir,
  pattern = "\\.[Rr][Dd][Ss]$",
  full.names = FALSE
)

x_files <- rds_files[
  grepl("_train\\.[Rr][Dd][Ss]$", rds_files) &
    !grepl("_true\\.[Rr][Dd][Ss]$", rds_files)
]

y_files <- rds_files[
  grepl("_train_true\\.[Rr][Dd][Ss]$", rds_files)
]

x_names <- sub("\\.[Rr][Dd][Ss]$", "", x_files)
y_names <- sub("_true\\.[Rr][Dd][Ss]$", "", y_files)

dataset_names <- sort(intersect(x_names, y_names))
dataset_names <- intersect(dataset_names, target_datasets)

data_list <- list()

for (dataset in dataset_names) {
  x_file <- x_files[x_names == dataset][1]
  X <- readRDS(file.path(root_dir, x_file))
  data_list[[dataset]] <- data.matrix(X)
}

n_size <- sapply(data_list, nrow)
data_list <- data_list[order(n_size)]

knn_grid <- 10
gamma_grid <- unique(round(exp(seq(log(0.01), log(20000), length.out = 1000)), 4))

result_dir <- file.path(root_dir, "results_baseline")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

task_grid <- expand.grid(
  dataset = names(data_list),
  knn = knn_grid,
  stringsAsFactors = FALSE
)

n_cores <- max(1, detectCores() - 1)
n_cores <- min(n_cores, nrow(task_grid))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

invisible(
  foreach(
    i = seq_len(nrow(task_grid)),
    .packages = c("cvxclustr"),
    .export = c("task_grid", "data_list", "gamma_grid", "result_dir")
  ) %dopar% {
    dataset <- task_grid$dataset[i]
    knn <- task_grid$knn[i]

    X <- t(data_list[[dataset]])
    n <- ncol(X)

    phi <- 0.5

    w <- kernel_weights(X, phi)
    w <- knn_weights(w, knn, n)

    save_dir <- file.path(result_dir, dataset)
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

    for (gamma in gamma_grid) {
      res <- cvxclust(
        X = X,
        w = w,
        gamma = gamma,
        method = "ama",
        nu = 0.001,
        tol = 0.0001,
        max_iter = 100000,
        type = 2,
        accelerate = TRUE
      )

      file_name <- paste0(
        dataset,
        "_cvxcls_k",
        knn,
        "_gam",
        gamma,
        ".Rds"
      )

      saveRDS(res, file.path(save_dir, file_name))
    }

    NULL
  }
)

stopCluster(cl)