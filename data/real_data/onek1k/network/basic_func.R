
# used basic functions
pd_matrix = function(par){  
  n = (sqrt(1+8*length(par))-1) / 2
  L = matrix(0, nrow = n, ncol = n)
  L[lower.tri(L, diag = TRUE)] = par
  V = L %*% t(L)
  return(V)
}

sym_matrix = function(par) {
  n = (sqrt(1+8*length(par))-1) / 2
  V = matrix(0, n, n)
  V[lower.tri(V, diag = TRUE)] = par
  V = V + t(V) - diag(diag(V))
  return(V)
}

vectorized = function(tensor) {
  # tensor 2 vector
  # a new vectorization method
  tensor_permuted = aperm(tensor, perm = c(3, 2, 1))
  vectorized = as.vector(tensor_permuted)
  return(vectorized)
}

tensor_reconstructed = function(vec, N, G, C){
  # vector 2 tensor
  tensor_permuted_reconstructed = array(vec, dim = c(C, G, N))
  tensor_reconstructed = aperm(tensor_permuted_reconstructed, perm = c(3, 2, 1))
  return(tensor_reconstructed)
}

only_one = function(d,i) {
  v = matrix(0, d, d)
  v[i,i] = 1
  return(v)
}
