root_dir = getwd()

set.seed(123)

n_each = 100

r_mean = c(0.45, 1.15)
r_sd = c(0.055, 0.070)

theta1 = runif(n_each, 0, 2 * pi)
theta2 = runif(n_each, 0, 2 * pi)

r1 = rnorm(n_each, mean = r_mean[1], sd = r_sd[1])
r2 = rnorm(n_each, mean = r_mean[2], sd = r_sd[2])

X_inner = cbind(r1 * cos(theta1), r1 * sin(theta1))
X_outer = cbind(r2 * cos(theta2), r2 * sin(theta2))

X = rbind(X_inner, X_outer)

z = runif(nrow(X), min = -1.5, max = 1.5)
X = cbind(X, z)

colnames(X) = c("x", "y", "z")

y = c(
  rep(1, n_each),
  rep(2, n_each)
)

saveRDS(X, file = file.path(root_dir, "ring.Rds"))
saveRDS(y, file = file.path(root_dir, "ring.true.Rds"))