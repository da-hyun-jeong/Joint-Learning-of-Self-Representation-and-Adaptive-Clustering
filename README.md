# Reproducibility Code

This repository contains the code used to fit the proposed adaptive convex clustering model and the baseline prototype-based methods used in the experiments.

## Method

Given a data matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$, the proposed method jointly learns a prototype matrix $U$ and an adaptive proximity matrix $Z$.

<p align="center">
  <img src="method.png" alt="Optimization problem" width="900">
</p>

Here, $\lambda,\gamma > 0$, $U_{i\cdot}$ denotes the $i$-th row of $U$, and $Z_{ij}$ denotes the $(i,j)$-th element of $Z$. The unique rows of $U$ are used as learned prototypes.

## Repository Structure

```text
.
├── admm.R
├── install.R
├── run_ours.R
├── run_cvscls.R
├── run_baselines.R
├── realdata/
└── simulation/
```

## Scripts

### `install.R`

Installs the R packages required for fitting the proposed method and the baseline methods.

### `admm.R`

Implements the ADMM solver for the proposed optimization problem.

### `run_ours.R`

Runs the proposed method on the real datasets stored in `realdata/`. The results are saved under

```text
realdata/results_ours/
```

### `run_cvscls.R`

Runs the standard convex clustering baseline using the `cvxclustr` package. The results are saved under

```text
realdata/results_baseline/
```

### `run_baselines.R`

Runs the following prototype-based baseline methods:

- k-means
- k-medoids
- self-organizing map
- neural gas

The results are saved under

```text
realdata/results_baseline/
```

### `simulation/moon.R`

Generates the two-moon synthetic dataset. The generated files are

```text
simulation/moon.Rds
simulation/moon_true.Rds
```

### `simulation/ring.R`

Generates the 3D ring synthetic dataset. The generated files are

```text
simulation/ring.Rds
simulation/ring_true.Rds
```

### `simulation/spiral.R`

Generates the 3D spiral synthetic dataset. The generated files are

```text
simulation/spiral.Rds
simulation/spiral_true.Rds
```

## Data Format

Each dataset is stored as an `.Rds` file. The corresponding ground-truth labels are stored using the suffix `_true.Rds`.

For example:

```text
iris_train.Rds
iris_train_true.Rds
iris_test.Rds
iris_test_true.Rds
```

The data files contain numeric matrices, and the label files contain class labels.

## Running the Code

All scripts are intended to be run from the repository root directory.

Install required packages:

```bash
Rscript install.R
```

Generate synthetic datasets:

```bash
Rscript simulation/moon.R
Rscript simulation/ring.R
Rscript simulation/spiral.R
```

Run the proposed method:

```bash
Rscript run_ours.R
```

Run the convex clustering baseline:

```bash
Rscript run_cvscls.R
```

Run the remaining baseline methods:

```bash
Rscript run_baselines.R
```

## Main Parameters

The proposed method uses the following default settings:

```text
rho = 1000
a = 1e-7
mu = n^2 * rho * ((1 - a) / a)
lambda grid = 60 logarithmically spaced values from 10 to 10000
gamma = 10
```

The convex clustering baseline uses the following settings:

```text
k-nearest-neighbor graph sizes = 5, 10, 20
gamma grid = 60 logarithmically spaced values from 10 to 10000
```

The other baseline methods are evaluated over the following settings:

```text
number of prototypes k = 2, ..., 50
random seeds = 1001, ..., 1020
```

## Notes

This repository provides fitting scripts for reproducing the experimental results. Plotting and post-processing scripts are not included in this minimal reproducibility package.
