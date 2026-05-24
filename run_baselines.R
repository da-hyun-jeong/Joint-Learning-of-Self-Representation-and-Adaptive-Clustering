library(cluster)
library(kohonen)
library(cclust)
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
  y_file <- y_files[y_names == dataset][1]

  X <- readRDS(file.path(root_dir, x_file))
  y <- readRDS(file.path(root_dir, y_file))

  data_list[[dataset]] <- list(
    X = data.matrix(X),
    y = as.vector(unlist(y))
  )
}

k_range <- 2:50

seed_vec <- c(
  1001, 1002, 1003, 1004, 1005,
  1006, 1007, 1008, 1009, 1010,
  1011, 1012, 1013, 1014, 1015,
  1016, 1017, 1018, 1019, 1020
)

result_dir <- file.path(root_dir, "results_baseline")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

task_grid <- expand.grid(
  dataset = names(data_list),
  k = k_range,
  seed = seed_vec,
  stringsAsFactors = FALSE
)

n_cores <- detectCores() - 1

cl <- makeCluster(n_cores)
registerDoParallel(cl)

res_log <- foreach(
  i = seq_len(nrow(task_grid)),
  .combine = rbind,
  .packages = c("cluster", "kohonen", "cclust"),
  .export = c("task_grid", "data_list", "result_dir")
) %dopar% {
  dataset <- task_grid$dataset[i]
  k <- task_grid$k[i]
  seed <- task_grid$seed[i]

  X <- data_list[[dataset]]$X

  save_dir <- file.path(result_dir, dataset)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  set.seed(seed)
  km_res <- kmeans(
    X,
    centers = k,
    nstart = 1,
    iter.max = 10000
  )

  saveRDS(
    km_res,
    file.path(save_dir, paste0(dataset, "_kmeans_k", k, "_seed", seed, ".Rds"))
  )

  set.seed(seed)
  pm_res <- pam(
    X,
    k = k,
    cluster.only = FALSE,
    do.swap = TRUE
  )

  saveRDS(
    pm_res,
    file.path(save_dir, paste0(dataset, "_kmedoids_k", k, "_seed", seed, ".Rds"))
  )

  set.seed(seed)
  som_grid <- somgrid(
    xdim = k,
    ydim = 1,
    topo = "rectangular"
  )

  som_res <- som(
    X,
    grid = som_grid,
    rlen = 10000
  )

  saveRDS(
    som_res,
    file.path(save_dir, paste0(dataset, "_som_k", k, "_seed", seed, ".Rds"))
  )

  set.seed(seed)
  ng_res <- cclust(
    X,
    centers = k,
    method = "neuralgas",
    iter.max = 10000
  )

  saveRDS(
    ng_res,
    file.path(save_dir, paste0(dataset, "_neuralgas_k", k, "_seed", seed, ".Rds"))
  )

  data.frame(
    dataset = dataset,
    k = k,
    seed = seed,
    status = "done",
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)

saveRDS(res_log, file.path(result_dir, "baseline_run_log.Rds"))