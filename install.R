packages <- c(
  "cvxclustr",
  "cluster",
  "kohonen",
  "cclust",
  "foreach",
  "doParallel"
)

installed <- rownames(installed.packages())
missing <- setdiff(packages, installed)

if (length(missing) > 0) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org",
    type = "binary"
  )
}