root_dir = getwd()

set.seed(123)

n = 100
sd_noise = 0.02

theta1 = runif(n, 0, pi)
theta2 = runif(n, 0, pi)

x1 = cos(theta1)
y1 = sin(theta1)

x2 = 1 - cos(theta2)
y2 = -sin(theta2) + 0.5

x1 = x1 + rnorm(n, 0, sd_noise)
y1 = y1 + rnorm(n, 0, sd_noise)
x2 = x2 + rnorm(n, 0, sd_noise)
y2 = y2 + rnorm(n, 0, sd_noise)

X = rbind(cbind(x1, y1), cbind(x2, y2))

colnames(X) = c("x", "y")

y = rep(1:2, each = n)

saveRDS(X, file = file.path(root_dir, "moon.Rds"))
saveRDS(y, file = file.path(root_dir, "moon_true.Rds"))