root_path = "D:/acvx"
dir.create(root_path, showWarnings = FALSE, recursive = TRUE)

setwd(root_path)

result_dir = file.path(root_path, "simulation")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(root_path, "admm.R"))


df_fs = c("moon.R", "spiral.R", "ring.R")

lambda_grid = round(exp(seq(log(10), log(10000), length.out = 60)), 0)
gamma_grid = 10

rho = 1000
a = 1e-7


for (df_f in df_fs) {
  df = sub("\\.R$", "", df_f)
  
  source(file.path(root_path, df_f))

  n = nrow(X)
  mu = n^2 * rho * ((1 - a) / a)
  
  save_path = file.path(result_dir, df)
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  
  for (gamma in gamma_grid) {
    for (lambda in lambda_grid) {
      
      cat( "dataset:", df, "| lambda:", lambda, "| gamma:", gamma, "\n")
      
      res = run_admm(
        X = X,
        rho = rho,
        mu = mu,
        lambda = lambda,
        gamma = gamma,
        max_iter = 10000000,
        tol = 1e-6,
        save_path = save_path,
        scen_nm = df,
        warm_start = NULL,
        save_every = 100
      )
    }
  }
}