library("gtools")
library(MASS)
library(Matrix)

data_generate = function(N, G, C,
                         beta, # G*C
                         V_alpha, # G*G matrix 
                         V_GC, # C * G * G
                         V_delta # C * G * G
                         ) {

  # Beta
  Beta <- array(rep(beta, each = N), dim = c(N, nrow(beta), ncol(beta)))

  # \alpha
  R = chol(V_alpha)
  Z <- matrix(rnorm(N * G), nrow = N, ncol = G)
  alpha <- Z %*% R 
  Alpha <- array(rep(alpha, C), dim = c(nrow(alpha), ncol(alpha), C))
  
  # generate \Gamma
  Gamma <- array(NA, dim = c(N, G, C))
  for (c in 1:C) { 
    res <- try(chol(V_GC[c,,]), silent = TRUE)
    if (inherits(res, "try-error")) {     
      M <- V_GC[c,,]
      min_eig <- min(eigen(M, symmetric = TRUE)$values)
      ridge <- max(0, -min_eig) + 1e-8
      diag(M) <- diag(M) + ridge
      V_GC[c,,] <- M
      R = chol(V_GC[c,,])
    } else { R <- res }
    Z <- matrix(rnorm(N * G), nrow = N, ncol = G)
    gamma_c <- Z %*% R # N*G
    Gamma[,,c] = gamma_c
  }

  # generate delta -- for simplicity, assume delta_ic iid across i
  Delta <- array(NA, dim = c(N, G, C))
  for (c in 1:C) {
    R = chol(V_delta[c,,])
    Z <- matrix(rnorm(N * G), nrow = N, ncol = G)
    delta_c <- Z %*% R # N*G
    Delta[,,c] = delta_c
  }

  # generate Y:N*G*C
  Y = Beta + Alpha + Gamma + Delta
  
  return(Y)
}