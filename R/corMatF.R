corMatF <- function(alpha, ni, corstr){
  switch(corstr,
         "exchangeable" = toeplitz(c(1, rep(alpha, ni-1))),
         "ar1" = toeplitz(alpha^(0:(ni-1))),
         "unstructured" = vec2uMat(alpha=alpha[1:(ni*(ni-1)/2)], ni),
         "independence" = diag(ni),
  )
}

vec2uMat <- function(alpha, ni){
  x <- matrix(1, ni, ni)
  x[upper.tri(x)] <- alpha
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  return(x)
}
