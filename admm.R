
frob_norm = function(A) sqrt(sum(A * A))

init_state = function(X) {
  n = nrow(X)
  p = ncol(X)
  
  x = rep(1:n, each = n)
  y = rep(1:n, times = n)
  edges = cbind(x, y) 
  edges = edges[edges[,2] > edges[,1],]
  L = nrow(edges)
  
  I_n = diag(n)
  S = (matrix(1, n, n) - I_n) / (n - 1)
  
  st = list()
  st$X = X
  st$edges = edges
  st$U = X
  st$U1 = X
  st$U2 = X
  st$U3 = X
  
  st$Z = S
  st$Z1 = S
  st$Z2 = S
  st$Z3 = I_n
  st$Z4 = I_n
  
  st$A = matrix(0, L, p)
  
  st$E1 = matrix(0, n, n)
  st$E2 = matrix(0, n, n)
  st$E3 = matrix(0, n, p)
  st$E4 = I_n - S
  st$E5 = matrix(0, n, p)
  st$E6 = matrix(0, L, p)
  st$E7 = matrix(0, n, p)
  st$E8 = matrix(0, n, p)
  st$E9 = I_n - S

  st$Psi1 = matrix(0, n, n)
  st$Psi2 = matrix(0, n, n)
  st$Psi3 = matrix(0, n, p)
  st$Psi4 = matrix(0, n, n)
  st$Psi5 = matrix(0, n, p)
  st$Psi6 = matrix(0, L, p)
  st$Psi7 = matrix(0, n, p)
  st$Psi8 = matrix(0, n, p)
  st$Psi9 = matrix(0, n, n)
  st
}

project_to_simplex = function(y) {
  u = sort(y, decreasing = TRUE)
  s = cumsum(u)
  idx = seq_along(u)
  k = max(which(u + (1 - s) / idx > 0))
  theta = (1 - s[k]) / k
  pmax(y + theta, 0)
}

update_Z = function(st) {
  st$Z = st$rho / (2 * st$gamma + 3 * st$rho) * (0.5 * (st$Z1 + t(st$Z2)) + st$Z3 + st$Z4 - st$E1 - st$E4 - st$E9 + st$Psi1 + st$Psi4 + st$Psi9)
  st
}

update_Z1 = function(st) {
  n = nrow(st$Z1)
  Z1_tilde = (4 / 5) * (0.5 * st$Z + st$Z2 - 0.25 * t(st$Z2) + 0.5 * st$E1 + st$E2 - 0.5 * st$Psi1 - st$Psi2)
  Z1_new = matrix(0, n, n)
  
  for (i in seq_len(n)) {
    Z1_new[i, -i] = project_to_simplex(Z1_tilde[i, -i])
  }
  
  st$Z1 = Z1_new
  st
}

update_Z2 = function(st) {
  n = nrow(st$Z2)
  Z2_tilde = (4 / 5) * (st$Z1 - 0.25 * t(st$Z1) + 0.5 * t(st$Z) + 0.5 * t(st$E1) - st$E2 - 0.5 * t(st$Psi1) + st$Psi2)
  Z2_new = matrix(0, n, n)
  
  for (j in seq_len(n)) {
    Z2_new[-j, j] = project_to_simplex(Z2_tilde[-j, j])
  }
  
  st$Z2 = Z2_new
  st
}

update_Z3 = function(st) {
  n = nrow(st$Z3)
  st$Z3 = ((st$U + st$E3 - st$Psi3) %*% t(st$U1) + st$Z + st$E4 - st$Psi4) %*% solve(st$U1 %*% t(st$U1) + diag(n))
  st
}

update_Z4 = function(st) {
  L = nrow(st$edges)
  Z4_new = st$Z + st$E9 - st$Psi9
  
  for (l in seq_len(L)) {
    i = st$edges[l, 1]
    j = st$edges[l, 2]
    d = st$U2[i, ] - st$U3[j, ]
    
    Z4_new[i, j] = (
      sum((st$A[l, ] + st$E6[l, ] - st$Psi6[l, ]) * d) +
        st$Z[i, j] + st$E9[i, j] - st$Psi9[i, j]
    ) / (sum(d * d) + 1)
  }
  
  st$Z4 = Z4_new
  st
}

update_U = function(st) {
  st$U = (st$X + st$rho * (st$Z3 %*% st$U1 + st$U1 + st$U2 + st$U3 - st$E3 - st$E5 - st$E7 - st$E8 + st$Psi3 + st$Psi5 + st$Psi7 + st$Psi8)) / (1 + 4 * st$rho)
  st
}

update_U1 = function(st) {
  n = nrow(st$U1)
  st$U1 = solve(t(st$Z3) %*% st$Z3 + diag(n), t(st$Z3) %*% (st$U + st$E3 - st$Psi3) + st$U + st$E5 - st$Psi5)
  st
}

update_U2 = function(st) {
  n = nrow(st$U2)
  U2_new = matrix(0, n, ncol(st$U2))
  
  for (i in seq_len(n)) {
    idx = which(st$edges[, 1] == i)
    num = st$U[i, ] + st$E7[i, ] - st$Psi7[i, ]
    den = 1
    
    if (length(idx) > 0) {
      for (l in idx) {
        j = st$edges[l, 2]
        z = st$Z4[st$edges[l, 1], st$edges[l, 2]]
        
        num = num +
          z^2 * st$U3[j, ] +
          z * st$A[l, ] +
          z * st$E6[l, ] -
          z * st$Psi6[l, ]
        
        den = den + z^2
      }
    }
    
    U2_new[i, ] = num / den
  }
  
  st$U2 = U2_new
  st
}

update_U3 = function(st) {
  n = nrow(st$U3)
  U3_new = matrix(0, n, ncol(st$U3))
  
  for (j in seq_len(n)) {
    idx = which(st$edges[, 2] == j)
    num = st$U[j, ] + st$E8[j, ] - st$Psi8[j, ]
    den = 1
    
    if (length(idx) > 0) {
      for (l in idx) {
        i = st$edges[l, 1]
        z = st$Z4[st$edges[l, 1], st$edges[l, 2]]
        
        num = num +
          z^2 * st$U2[i, ] -
          z * st$A[l, ] -
          z * st$E6[l, ] +
          z * st$Psi6[l, ]
        
        den = den + z^2
      }
    }
    
    U3_new[j, ] = num / den
  }
  
  st$U3 = U3_new
  st
}

update_A = function(st) {
  L = nrow(st$edges)
  A_new = matrix(0, L, ncol(st$A))
  
  for (l in seq_len(L)) {
    i = st$edges[l, 1]
    j = st$edges[l, 2]
    v = st$Z4[i, j] * (st$U2[i, ] - st$U3[j, ]) - st$E6[l, ] + st$Psi6[l, ]
    nv = sqrt(sum(v * v))
    
    if (nv > st$lambda / st$rho) {
      A_new[l, ] = ((st$rho * nv - st$lambda) / ((st$lambda + st$rho) * nv)) * v
    }
  }
  
  st$A = A_new
  st
}

update_E1 = function(st) {
  st$E1 = st$rho / (st$mu + st$rho) * (0.5 * (st$Z1 + t(st$Z2)) - st$Z + st$Psi1)
  st
}

update_E2 = function(st) {
  st$E2 = st$rho / (st$mu + st$rho) * (st$Z1 - st$Z2 + st$Psi2)
  st
}

update_E3 = function(st) {
  st$E3 = st$rho / (st$mu + st$rho) * (st$Z3 %*% st$U1 - st$U + st$Psi3)
  st
}

update_E4 = function(st) {
  st$E4 = st$rho / (st$mu + st$rho) * (st$Z3 - st$Z + st$Psi4)
  st
}

update_E5 = function(st) {
  st$E5 = st$rho / (st$mu + st$rho) * (st$U1 - st$U + st$Psi5)
  st
}

update_E6 = function(st) {
  L = nrow(st$edges)
  E6_new = matrix(0, L, ncol(st$E6))
  
  for (l in seq_len(L)) {
    i = st$edges[l, 1]
    j = st$edges[l, 2]
    
    E6_new[l, ] = st$rho / (st$mu + st$rho) * (
      st$Z4[i, j] * (st$U2[i, ] - st$U3[j, ]) - st$A[l, ] + st$Psi6[l, ]
    )
  }
  
  st$E6 = E6_new
  st
}

update_E7 = function(st) {
  st$E7 = st$rho / (st$mu + st$rho) * (st$U2 - st$U + st$Psi7)
  st
}

update_E8 = function(st) {
  st$E8 = st$rho / (st$mu + st$rho) * (st$U3 - st$U + st$Psi8)
  st
}

update_E9 = function(st) {
  st$E9 = st$rho / (st$mu + st$rho) * (st$Z4 - st$Z + st$Psi9)
  st
}

update_duals = function(st) {
  L = nrow(st$edges)
  Psi6_new = st$Psi6
  
  for (l in seq_len(L)) {
    i = st$edges[l, 1]
    j = st$edges[l, 2]
    Psi6_new[l, ] = Psi6_new[l, ] + st$Z4[i, j] * (st$U2[i, ] - st$U3[j, ]) - st$A[l, ] - st$E6[l, ]
  }
  
  st$Psi1 = st$Psi1 + 0.5 * (st$Z1 + t(st$Z2)) - st$Z - st$E1
  st$Psi2 = st$Psi2 + st$Z1 - st$Z2 - st$E2
  st$Psi3 = st$Psi3 + st$Z3 %*% st$U1 - st$U - st$E3
  st$Psi4 = st$Psi4 + st$Z3 - st$Z - st$E4
  st$Psi5 = st$Psi5 + st$U1 - st$U - st$E5
  st$Psi6 = Psi6_new
  st$Psi7 = st$Psi7 + st$U2 - st$U - st$E7
  st$Psi8 = st$Psi8 + st$U3 - st$U - st$E8
  st$Psi9 = st$Psi9 + st$Z4 - st$Z - st$E9
  st
}

stopping_criteria = function(st, st_prev) {
  n = nrow(st$Z)
  p = ncol(st$U)
  L = nrow(st$edges)
  
  R1 = 0.5 * (st$Z1 + t(st$Z2)) - st$Z - st$E1
  R2 = st$Z1 - st$Z2 - st$E2
  R3 = st$Z3 %*% st$U1 - st$U - st$E3
  R4 = st$Z3 - st$Z - st$E4
  R5 = st$U1 - st$U - st$E5
  R7 = st$U2 - st$U - st$E7
  R8 = st$U3 - st$U - st$E8
  R9 = st$Z4 - st$Z - st$E9
  
  D6 = matrix(0, L, p)
  R6 = matrix(0, L, p)
  
  for (l in seq_len(L)) {
    i = st$edges[l, 1]
    j = st$edges[l, 2]
    
    D6[l, ] = st$Z4[i, j] * (st$U2[i, ] - st$U3[j, ])
    R6[l, ] = D6[l, ] - st$A[l, ] - st$E6[l, ]
  }
  
  r = sqrt(
    sum(R1 * R1) + sum(R2 * R2) + sum(R3 * R3) +
      sum(R4 * R4) + sum(R5 * R5) + sum(R6 * R6) +
      sum(R7 * R7) + sum(R8 * R8) + sum(R9 * R9)
  )
  
  lhs_norm = sqrt(
    sum((0.5 * (st$Z1 + t(st$Z2)))^2) +
      sum((st$Z1 - st$Z2)^2) +
      sum((st$Z3 %*% st$U1)^2) +
      sum(st$Z3^2) +
      sum(st$U1^2) +
      sum(D6^2) +
      sum(st$U2^2) +
      sum(st$U3^2) +
      sum(st$Z4^2)
  )
  
  rhs_norm = sqrt(
    sum((st$Z + st$E1)^2) +
      sum(st$E2^2) +
      sum((st$U + st$E3)^2) +
      sum((st$Z + st$E4)^2) +
      sum((st$U + st$E5)^2) +
      sum((st$A + st$E6)^2) +
      sum((st$U + st$E7)^2) +
      sum((st$U + st$E8)^2) +
      sum((st$Z + st$E9)^2)
  )
  
  m_eq = 4 * n^2 + 4 * n * p + L * p
  eps_pri = st$tol * sqrt(m_eq) + st$tol * max(lhs_norm, rhs_norm)
  
  delta = sqrt(
    sum((st$Z - st_prev$Z)^2) +
      sum((st$Z1 - st_prev$Z1)^2) +
      sum((st$Z2 - st_prev$Z2)^2) +
      sum((st$Z3 - st_prev$Z3)^2) +
      sum((st$Z4 - st_prev$Z4)^2) +
      sum((st$U- st_prev$U)^2) +
      sum((st$U1 - st_prev$U1)^2) +
      sum((st$U2 - st_prev$U2)^2) +
      sum((st$U3 - st_prev$U3)^2) +
      sum((st$A - st_prev$A)^2) +
      sum((st$E1 - st_prev$E1)^2) +
      sum((st$E2 - st_prev$E2)^2) +
      sum((st$E3 - st_prev$E3)^2) +
      sum((st$E4 - st_prev$E4)^2) +
      sum((st$E5 - st_prev$E5)^2) +
      sum((st$E6 - st_prev$E6)^2) +
      sum((st$E7 - st_prev$E7)^2) +
      sum((st$E8 - st_prev$E8)^2) +
      sum((st$E9 - st_prev$E9)^2)
  )
  
  h_now = sqrt(
    sum(st$Z^2) +
      sum(st$Z1^2) + sum(st$Z2^2) +
      sum(st$Z3^2) + sum(st$Z4^2) +
      sum(st$U^2) +
      sum(st$U1^2) + sum(st$U2^2) + sum(st$U3^2) +
      sum(st$A^2) +
      sum(st$E1^2) + sum(st$E2^2) + sum(st$E3^2) +
      sum(st$E4^2) + sum(st$E5^2) + sum(st$E6^2) +
      sum(st$E7^2) + sum(st$E8^2) + sum(st$E9^2)
  )
  
  h_prev = sqrt(
    sum(st_prev$Z^2) +
      sum(st_prev$Z1^2) + sum(st_prev$Z2^2) +
      sum(st_prev$Z3^2) + sum(st_prev$Z4^2) +
      sum(st_prev$U^2) +
      sum(st_prev$U1^2) + sum(st_prev$U2^2) + sum(st_prev$U3^2) +
      sum(st_prev$A^2) +
      sum(st_prev$E1^2) + sum(st_prev$E2^2) + sum(st_prev$E3^2) +
      sum(st_prev$E4^2) + sum(st_prev$E5^2) + sum(st_prev$E6^2) +
      sum(st_prev$E7^2) + sum(st_prev$E8^2) + sum(st_prev$E9^2)
  )
  
  m_var = 9 * n^2 + 8 * n * p + 2 * L * p
  eps_chg = st$tol * sqrt(m_var) + st$tol * max(h_now, h_prev)
  
  list(
    residual = r,
    tolerance = eps_pri,
    delta = delta,
    tolerance_chg = eps_chg,
    converge = (r <= eps_pri) && (delta <= eps_chg)
  )
}

run_admm = function(X, rho, mu, lambda, gamma,
                    max_iter = 2000, tol = 1e-6,
                    save_path = getwd(), scen_nm = 'admm',
                    warm_start = NULL, save_every = 300) {
  
  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  
  lam_str = format(lambda, scientific = FALSE, trim = TRUE)
  gam_str = format(gamma, scientific = FALSE, trim = TRUE)
  fnm = sprintf("%s_lam%s_gam%s", scen_nm, lam_str, gam_str)
  out_file = file.path(save_path, paste0(fnm, ".Rds"))
  
  if (is.null(warm_start)) {
    st = init_state(X)
  } else {
    st = warm_start
  }
  
  st$X = X
  st$rho = rho
  st$mu = mu
  st$lambda = lambda
  st$gamma = gamma
  st$tol = tol
  
  for (k in seq_len(max_iter)) {
    
    st_prev = st
    
    st = update_Z(st)
    st = update_Z1(st)
    st = update_Z2(st)
    st = update_Z3(st)
    st = update_Z4(st)
    
    st = update_U(st)
    st = update_U1(st)
    st = update_U2(st)
    st = update_U3(st)
    
    st = update_A(st)
    
    st = update_E1(st)
    st = update_E2(st)
    st = update_E3(st)
    st = update_E4(st)
    st = update_E5(st)
    st = update_E6(st)
    st = update_E7(st)
    st = update_E8(st)
    st = update_E9(st)
    
    st = update_duals(st)
    
    sc = stopping_criteria(st, st_prev)
    st$converge = sc$converge
    st$iter = k
    st$pri_residual = sc$residual
    st$pri_tolerance = sc$tolerance
    st$chg_residual = sc$delta
    st$chg_tolerance = sc$tolerance_chg
    
    cat(sprintf(
      "Iter %d | pri=%.5f | eps_pri=%.5f | chg=%.5f | eps_chg=%.5f\n",
      k, sc$residual, sc$tolerance, sc$delta, sc$tolerance_chg
    ))
    
    if (sc$converge) {
      cat(sprintf("Converged at iter %d\n", k))
      saveRDS(st, file = out_file)
      return(st)
    }
    
    if (k == max_iter) {
      cat(sprintf("Reached max_iter=%d (not converged)\n", k))
      saveRDS(st, file = out_file)
      return(st)
    }
    
    if (k %% save_every == 0) {
      saveRDS(st, file = out_file)
    }
  }
  
  st
}