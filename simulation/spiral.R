rm(list=ls())
root_dir = 'D:/acvx/simulation'

fig_dir  = file.path(root_dir, 'figure')
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)

set.seed(1)

n_per_curve <- 100

r0 <- 1.0
a <- 0.0
b <- 0.56
turns <- 3.05
phi <- 0
noise_sd <- 0.0

height <- 6

# ------------------------------------------------------------
# angle and radius
# ------------------------------------------------------------

theta0 <- (r0 - a) / b
theta1 <- theta0 + turns * pi

th <- seq(theta0, theta1, length.out = n_per_curve)

r <- a + b * th

# ------------------------------------------------------------
# three spiral arms
# ------------------------------------------------------------

ang_1 <- th - theta0 + phi
ang_2 <- ang_1 + 2 * pi / 3
ang_3 <- ang_1 + 4 * pi / 3

x1 <- r * cos(ang_1)
y1 <- r * sin(ang_1)

x2 <- r * cos(ang_2)
y2 <- r * sin(ang_2)

x3 <- r * cos(ang_3)
y3 <- r * sin(ang_3)

# ------------------------------------------------------------
# optional noise
# ------------------------------------------------------------

if (noise_sd > 0) {
  x1 <- x1 + rnorm(n_per_curve, 0, noise_sd)
  y1 <- y1 + rnorm(n_per_curve, 0, noise_sd)
  
  x2 <- x2 + rnorm(n_per_curve, 0, noise_sd)
  y2 <- y2 + rnorm(n_per_curve, 0, noise_sd)
  
  x3 <- x3 + rnorm(n_per_curve, 0, noise_sd)
  y3 <- y3 + rnorm(n_per_curve, 0, noise_sd)
}

# ------------------------------------------------------------
# z-axis path
# ------------------------------------------------------------

z_path <- seq(0, height, length.out = n_per_curve)

# ------------------------------------------------------------
# final data
# ------------------------------------------------------------

M1 <- cbind(x1, y1, z_path)
M2 <- cbind(x2, y2, z_path)
M3 <- cbind(x3, y3, z_path)

X_spiral <- rbind(M1, M2, M3)

colnames(X_spiral) <- c("x", "y", "z")

y_spiral <- c(
  rep(1, n_per_curve),
  rep(2, n_per_curve),
  rep(3, n_per_curve)
)

spiral <- list(
  X = X_spiral,
  y = y_spiral
)

saveRDS(X_spiral, file = file.path(root_dir, 'spiral.Rds'))
saveRDS(y_spiral, file = file.path(root_dir, 'spiral_true.Rds'))
